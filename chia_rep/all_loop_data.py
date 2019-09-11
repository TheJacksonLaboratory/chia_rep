from collections import OrderedDict

from prettytable import PrettyTable

from .chrom_loop_data import ChromLoopData
import numpy as np
import os
import logging
import math

VERSION = 8
log = logging.getLogger()
log_bin = logging.getLogger('bin')

# Missing in many miseq peak files
CHROMS_TO_IGNORE = ['chrY', 'chrM']

DEFAULT_NUM_PEAKS = 10


class AllLoopData:

    """
    A class used to represent a sample

    Attributes
    ----------
    species_name : str
        Name of the sample's species (hg38, mm10, ...)
    sample_name : str
        Name of the sample (LHH0061, LHH0061_0061H, ...)
    peak_dict : dict(str, list)
        Key: Name of chromosome (chr1, chr2, ...)
        Value: List of peaks in chromosome
        Peak format: [start, end, length, max_value, mean_value (optional)]
    chrom_dict : dict(str, ChromLoopData)
        Key: Name of chromosome
        Value: ChromLoopData object
    """

    def __init__(self, chrom_size_file_path, loop_file_path, bedgraph,
                 peak_dict, chrom_to_load=None, min_loop_value=0):
        """
        Parameters
        ----------
        chrom_size_file_path : str
            File containing the base pair size of each chromosome to use
        loop_file_path : str
            File containing loops in format:
            chr1  start   end chr1  start   end pet_count
        bedgraph : BedGraph
            The bedgraph file for this sample (from pyBedGraph)
        peak_dict : dict(str, list)
            Key: Name of chromosome (chr1, chr2, ...)
            Value: List of peaks in chromosome
            Peak format: [start, end, length]
        chrom_to_load : str, optional
             Name of chromosome to load (default is None)
        min_loop_value : int, optional
            Minimum loop value (PET count) to include (default is 0)
        """

        # Prints peak_dict which is too large to be meaningful
        # log.debug(locals())

        self.is_hiseq = True
        self.species_name = os.path.basename(chrom_size_file_path).split('.')[0]
        self.sample_name = os.path.basename(loop_file_path).split('.')[0]

        for chrom_name, peak_chrom in peak_dict.items():
            if not bedgraph.has_chrom(chrom_name):
                log.error(f'{bedgraph.name} does not have {chrom_name}')
                continue

            bedgraph.load_chrom_data(chrom_name)
            start_list = [x[0] for x in peak_chrom]
            end_list = [x[1] for x in peak_chrom]
            max_list = \
                bedgraph.stats(start_list=start_list, end_list=end_list,
                               chrom_name=chrom_name, stat='max')
            mean_list = \
                bedgraph.stats(start_list=start_list, end_list=end_list,
                               chrom_name=chrom_name, stat='mean')
            for i in range(max_list.size):
                peak_chrom[i].append(max_list[i])
                peak_chrom[i].append(mean_list[i])
            bedgraph.free_chrom_data(chrom_name)
        self.peak_dict = peak_dict

        # Initialize all chromosomes to be loaded
        self.chrom_dict = {}
        with open(chrom_size_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chrom_to_load and chrom_name != chrom_to_load:
                    continue

                if chrom_name in CHROMS_TO_IGNORE:
                    continue

                if chrom_name not in peak_dict:
                    continue

                chrom_size = int(line[1])

                self.chrom_dict[chrom_name] = \
                    ChromLoopData(chrom_name, chrom_size, self.sample_name)

        with open(loop_file_path) as in_file:
            loop_anchor_list = []
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chrom_name not in self.chrom_dict:
                    continue

                loop_value = int(line[6])
                if loop_value < min_loop_value:
                    continue

                loop_end1 = int(line[4])
                loop_end2 = int(line[5])
                loop_start1 = int(line[1])
                loop_start2 = int(line[2])

                self.chrom_dict[chrom_name].add_loop(loop_start1, loop_start2,
                                                     loop_end1, loop_end2,
                                                     loop_value)

                start_interval = loop_start2 - loop_start1
                end_interval = loop_end2 - loop_end1

                loop_anchor_list.append(start_interval)
                loop_anchor_list.append(end_interval)

            log.debug(f'Anchor mean width: {np.mean(loop_anchor_list)}')

        to_remove = []
        for chrom_name in self.chrom_dict:
            if not self.chrom_dict[chrom_name].finish_init(bedgraph):
                to_remove.append(chrom_name)

        # Chromosomes with no loops or other random problems
        for chrom_name in to_remove:
            del self.chrom_dict[chrom_name]

    def preprocess(self, num_peaks=DEFAULT_NUM_PEAKS, both_peak_support=False,
                   kept_dir=None):
        """
        Preprocess all the chromosomes in this object.

        Removes all problematic chromosomes (not enough loops, etc...). Keeps
        only num_peaks peaks in each list in peak_dict

        Parameters
        ----------
        num_peaks : int, optional
            The number of peaks to use when filtering
            (default is DEFAULT_NUM_PEAKS)
        both_peak_support : bool, optional
            Whether to only keep loops that have peak support on both sides
            (default is False)
        kept_dir : str, optional
            Directory to output kept peaks and filters
            (default is None)

        Returns
        ------
        int
            Number of kept loops in all chromosomes
        """

        to_remove = []
        for name, chrom in self.chrom_dict.items():
            if not chrom.preprocess(num_peaks, self.peak_dict[name],
                                    both_peak_support=both_peak_support,
                                    kept_dir=kept_dir):
                to_remove.append(name)

            # Keep only used peaks
            self.peak_dict[name] = self.peak_dict[name][:num_peaks]

        # Remove problematic chromosomes
        for name in to_remove:
            del self.chrom_dict[name]

        total_kept_loops = 0
        for name, chrom_loop in self.chrom_dict.items():
            total_kept_loops += chrom_loop.filtered_numb_values
        return total_kept_loops

    def compare(self, o_loop_data, bin_size, window_size, window_index=None,
                chroms_to_compare=None, is_rep=False):
        """
        Compares this sample against another sample

        Gets the values for each window in each chromosome and combines that
        into a genome-wide value.

        Parameters
        ----------
        o_loop_data : AllLoopData
            The sample to compare against
        bin_size : int
            Determines which loops are the same by putting them into bins
        window_size : int
            Cuts off all loops that are not entirely inside the window
        window_index : int, optional
            Determines a specific window index to look at (Default is None)
        chroms_to_compare : list(str), optional
            A list of chromosomes to compare (Default is All)
        is_rep : bool, optional
            Debugging purposes

        Returns
        ------
        dict
            graph_type : str
                Information about which comparison method was used
            rep : str
                Combined unweighted value of all chromosomes
            w_rep : str
                Combined weighted value of all chromosomes
            main : str
                The reproducibility value of this comparison (w_rep)
        """

        # my_kept_loops = self.preprocess(peak_percent_kept)
        # o_kept_loops = o_loop_data.preprocess(peak_percent_kept)
        #
        # my_loops = (self, my_kept_loops)
        # o_loops = (o_loop_data, o_kept_loops)
        # lesser_numb = min(my_loops, o_loops, key=lambda x: x[1])
        # greater_numb = max(my_loops, o_loops, key=lambda x: x[1])
        # if lesser_numb[1] / greater_numb[1] < 0.75:
        #     lesser_numb[0].preprocess(peak_percent_kept, greater_numb[1])

        # Default: Compare all the chromosomes
        if chroms_to_compare is None:
            chroms_to_compare = list(self.chrom_dict.keys())

        chrom_value_table = \
            PrettyTable(['chrom', 'graph_type', 'rep', 'w_rep'])

        chrom_value_list = []
        log.info(f'Chroms to compare: {chroms_to_compare}')
        for chrom_name in chroms_to_compare:

            if chrom_name not in self.chrom_dict:
                log.warning(f'{chrom_name} is not in {self.sample_name}. '
                            f'Skipping {chrom_name}')
                continue

            if chrom_name not in o_loop_data.chrom_dict:
                log.warning(f'{chrom_name} is in {self.sample_name} but '
                            f'not in {o_loop_data.sample_name}. Skipping '
                            f'{chrom_name}')
                continue

            log.info(f"Comparing {chrom_name} ...")

            # Compare for all windows in chrom
            chrom_size = self.chrom_dict[chrom_name].size
            value_dict_list = []
            numb_windows = math.ceil(chrom_size / window_size)
            for k in range(numb_windows):

                # If there is a specified window, just compare that
                if window_index is not None:
                    k = window_index

                window_start = window_size * k
                window_end = window_size * (k + 1)
                if window_end > chrom_size:
                    window_end = chrom_size

                value_dict_list.append(
                    self.chrom_dict[chrom_name].compare(
                        o_loop_data.chrom_dict[chrom_name], window_start,
                        window_end, bin_size, is_rep=is_rep))

                if window_index is not None:
                    break

            values = [x['rep'] for x in value_dict_list]
            chrom_value = {
                'graph_type': value_dict_list[0]['graph_type'],
                'rep': np.mean(values),
            }
            # for value_dict in value_dict_list:
            #    assert chrom_value['graph_type'] == value_dict['graph_type']
            try:  # Weigh value from each bin according to max loop in graph
                chrom_value['w_rep'] = \
                    np.average(values, weights=[x['w'] for x in
                                                value_dict_list])
            except ZeroDivisionError:  # sum of weights == 0
                log.exception(f"No loops were found in either graphs. Skipping"
                              f"{chrom_name}")
                continue

            # chrom_value = self.chrom_dict[chrom_name].compare(
            #     o_loop_data.chrom_dict[chrom_name],
            #     self.peak_dict[chrom_name] + o_loop_data.peak_dict[chrom_name],
            #     is_rep=is_rep)

            chrom_value_list.append(chrom_value)

            log.debug(chrom_value)
            chrom_value_table.add_row([chrom_name] + list(chrom_value.values()))
            log_bin.info(chrom_value_table)
            chrom_value_table.clear_rows()

        log.debug(chrom_value_list)

        # Weigh value from each chromosome equally
        avg_value = {
            'graph_type': chrom_value_list[0]['graph_type'],
            'rep': np.mean([x['rep'] for x in chrom_value_list]),
            'w_rep': np.mean([x['w_rep'] for x in chrom_value_list])
        }
        avg_value['main'] = avg_value['w_rep']
        log.debug(avg_value)

        return avg_value
