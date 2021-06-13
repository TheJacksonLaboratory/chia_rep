import numpy as np
import os
import logging
import math
from pyBedGraph import BedGraph
from typing import Dict, List

from .chrom_loop_data import ChromLoopData, PEAK_MAX_VALUE_INDEX

log = logging.getLogger()

# Missing in many miseq peak files
CHROMS_TO_IGNORE = ['chrY', 'chrM']


class GenomeLoopData:
    """
    A class used to represent a sample

    Attributes
    ----------
    species_name : str
        Name of the sample's species (hg38, mm10, ...)
    sample_name : str
        Name of the sample (LHH0061, LHH0061_0061H, ...)
    peak_dict : dict[str, list]
        Key: Name of chromosome (chr1, chr2, ...)
        Value: List of peaks in chromosome
        Peak format: [start, end, length, max_value, mean_value (optional)]
    chrom_dict : dict[str, ChromLoopData]
        Key: Name of chromosome
        Value: ChromLoopData object
    """

    def __init__(
        self,
        chrom_size_file: str,
        loop_file: str,
        bedgraph: BedGraph,
        peak_dict: Dict[str, list],
        chroms_to_load: List[str] = None,
        min_loop_value: int = 0
    ):
        """
        Initializes all chromosomes and adds loops to them from given file.

        Finds peak max from bedgraph

        Parameters
        ----------
        chrom_size_file : str
            File containing the base pair size of each chromosome to use
        loop_file : str
            File containing loops in format:
            chrom1  start1   end1 chrom2  start2   end2 pet_count
        bedgraph : BedGraph
            The bedgraph file for this sample (from pyBedGraph)
        peak_dict : dict[str, list]
            Key: Name of chromosome (chr1, chr2, ...)
            Value: List of peaks in chromosome
            Peak format: [start, end, length]
        chroms_to_load : list, optional
             List of names of chromosome to load (default is None)
        min_loop_value : int, optional
            Minimum loop value (PET count) to include (default is 0)
        """

        # Prints peak_dict which is too large to be meaningful
        # log.debug(locals())

        self.species_name = os.path.basename(chrom_size_file).split('.')[0]
        self.sample_name = os.path.basename(loop_file).split('.')[0]

        self.total_samples = 0

        self.peak_dict = {}

        # Find values for each peak since peak caller is not accurate sometimes
        for chrom_name, peak_chrom in peak_dict.items():
            if not bedgraph.has_chrom(chrom_name):
                log.warning(f'{bedgraph.name} does not have {chrom_name}')
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

            self.peak_dict[chrom_name] = peak_dict[chrom_name]

        # Initialize all chromosomes to be loaded
        self.chrom_dict = {}
        with open(chrom_size_file) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chroms_to_load and chrom_name not in chroms_to_load:
                    continue

                if chrom_name in CHROMS_TO_IGNORE:
                    continue

                if chrom_name not in peak_dict:
                    continue

                chrom_size = int(line[1])

                self.chrom_dict[chrom_name] = \
                    ChromLoopData(chrom_name, chrom_size, self.sample_name)

        with open(loop_file) as in_file:
            loop_anchor_list = []
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chrom_name not in self.chrom_dict:
                    continue

                loop_value = int(line[6])
                if loop_value < min_loop_value:
                    continue

                # head interval
                loop_start1 = int(line[1])
                loop_end1 = int(line[2])

                # tail anchor
                loop_start2 = int(line[4])
                loop_end2 = int(line[5])

                self.chrom_dict[chrom_name].add_loop(loop_start1, loop_end1,
                                                     loop_start2, loop_end2,
                                                     loop_value)

                head_interval = loop_end1 - loop_start1
                tail_interval = loop_end2 - loop_start2

                loop_anchor_list.append(head_interval)
                loop_anchor_list.append(tail_interval)

            log.debug(f'Anchor mean width: {np.mean(loop_anchor_list)}')

        # Get rid of chroms that had problems initializing
        to_remove = []
        for chrom_name in self.chrom_dict:
            if self.chrom_dict[chrom_name].finish_init(bedgraph):
                self.total_samples += \
                    np.sum(self.chrom_dict[chrom_name].value_list)
            else:
                to_remove.append(chrom_name)

        # Chromosomes with no loops or other random problems
        for chrom_name in to_remove:
            del self.chrom_dict[chrom_name]

    def filter_peaks(
        self,
        num_peaks: int,
        base_chrom: str
    ) -> None:
        """
        Keeps only the top num_peaks peaks. Takes the ratio of kept peaks in the
        base_chrom and uses it for all other chromosomes to calculate which
        peaks are kept.

        Parameters
        ----------
        num_peaks
        base_chrom

        Returns
        -------
        None
        """
        # Find percentage of num_peaks in base_chrom to use in other chroms
        if base_chrom not in self.chrom_dict or \
                base_chrom not in self.peak_dict or num_peaks < 1:
            log.warning(f'Unable to filter peaks since {base_chrom} is not '
                        f'available or num_peaks is not positive: {num_peaks}')
            return

        base_peak_list = self.peak_dict[base_chrom]

        if num_peaks > len(base_peak_list):
            num_peaks = len(base_peak_list)

        # while num_peaks - 1 > 0 and base_peak_list[num_peaks - 1][PEAK_MAX_VALUE_INDEX] < MIN_PEAK_VALUE:
        #     num_peaks -= 1

        log.debug(f'Number of peaks: {num_peaks}')

        min_peak_value = \
            self.peak_dict[base_chrom][num_peaks - 1][PEAK_MAX_VALUE_INDEX]
        chrom_peak_ratio = num_peaks / len(base_peak_list)
        log.debug(f'Min peak value from {base_chrom}: {min_peak_value}')
        log.debug(f'Chrom peak ratio from {base_chrom}: {chrom_peak_ratio}')

        to_remove = []
        for chrom_name in self.chrom_dict:
            peak_list = self.peak_dict[chrom_name]
            # Filter by minimum peak value found in base_chrom
            # if min_peak_value > 0:
            #     for i in range(len(peak_list)):
            #         if peak_list[i][PEAK_MAX_VALUE_INDEX] < min_peak_value:
            #             log.debug(
            #                 f'Filtered out {len(peak_list) - i} peaks')
            #             self.peak_dict[chrom_name] = self.peak_dict[chrom_name][:i]
            #             break
            self.peak_dict[chrom_name] = peak_list[:int(
                len(peak_list) * chrom_peak_ratio)]

            if len(self.peak_dict[chrom_name]) == 0:
                log.warning(
                    f'Removing {chrom_name} since it has no peaks above {min_peak_value}')
                # log.warning(f'{name} total peaks: {len(peak_list)}, '
                #             f'ratio: {peak_num_ratio}')
                to_remove.append(chrom_name)
                continue

        for name in to_remove:
            del self.chrom_dict[name]

    def preprocess(
        self,
        num_peaks: int = None,
        both_peak_support: bool = False,
        output_dir: str = 'output',
        base_chrom: str = 'chr1'
    ) -> None:
        """
        Preprocess all the chromosomes in chrom_dict attribute of this object.

        Removes all problematic chromosomes (not enough loops, etc...). Keeps
        only num_peaks peaks in each list in peak_dict

        Parameters
        ----------
        num_peaks : int, optional
            The number of peaks to use when filtering. Only to be used with chr1
            since other chromosomes will be dependent on min peak used from chr1
            (default is DEFAULT_NUM_PEAKS)
        both_peak_support : bool, optional
            Whether to only keep loops that have peak support on both sides
            (default is False)
        output_dir : str, optional
            Directory to output found peaks and filters
        base_chrom : str, optional

        Returns
        ------
        None
        """

        for peak_list in self.peak_dict.values():
            peak_list.sort(key=lambda x: x[PEAK_MAX_VALUE_INDEX], reverse=True)

        skip_peak_filter = False
        num_peak_param = num_peaks
        if num_peaks is None:
            num_peak_param = 'all'
            num_peaks = len(self.peak_dict[base_chrom])
            skip_peak_filter = True

        if num_peaks == -1:
            log.error(f'num_peaks is not positive')
            return

        if not skip_peak_filter:
            self.filter_peaks(num_peaks, base_chrom)

        to_remove = []
        for chrom_name, chrom_data in self.chrom_dict.items():
            success = chrom_data.preprocess(self.peak_dict[chrom_name],
                                            both_peak_support=both_peak_support)
            if not success:
                to_remove.append(chrom_name)

        # Remove problematic chromosomes
        for chrom_name in to_remove:
            del self.chrom_dict[chrom_name]

        os.makedirs(f'{output_dir}/peaks', exist_ok=True)
        os.makedirs(f'{output_dir}/loops', exist_ok=True)

        with open(
                f'{output_dir}/peaks/{self.sample_name}.{num_peak_param}.peaks',
                'w') as out_file:
            for chrom_name, peak_list in self.peak_dict.items():
                for peak in peak_list:
                    out_file.write(f'{chrom_name}\t{peak[0]}\t{peak[1]}\t'
                                   f'{peak[PEAK_MAX_VALUE_INDEX]}\n')

        with open(
                f'{output_dir}/loops/{self.sample_name}.{num_peak_param}.loops',
                'w') as out_file:
            for chrom_name, chrom_data in self.chrom_dict.items():
                for i in chrom_data.kept_indexes:
                    out_file.write(
                        f'{chrom_name}\t'
                        f'{chrom_data.start_anchor_list[0][i]}\t'
                        f'{chrom_data.start_anchor_list[1][i]}\t'
                        f'{chrom_name}\t'
                        f'{chrom_data.end_anchor_list[0][i]}\t'
                        f'{chrom_data.end_anchor_list[1][i]}\t'
                        f'{chrom_data.pet_count_list[i]}\t'
                        f'{chrom_data.value_list[i]}\n')

    def compare(
        self,
        o_loop_data: 'GenomeLoopData',
        window_size: int,
        bin_size: int,
        num_peaks: any,
        chroms_to_compare: List[str] = None,
        output_dir: str = 'output',
        do_output_graph: bool = False
    ) -> Dict[str, float]:
        """
        Compares this sample against another sample

        Gets the comparison values for each window in each chromosome and
        combines that into a genome-wide comparison value. Each window is given
        a weight based on the highest loop in it. Each chromosome is weighted
        equally.

        Parameters
        ----------
        o_loop_data : GenomeLoopData
            The sample to compare this sample against
        window_size : int
            Splits the chromosome into ceil(chrom_size / window_size) windows.
            Loops that start in one window and end in another are filter out.
        bin_size : int
            Splits each window into window_size / bin_size bins which determines
            the start and end of each loop
        num_peaks : any
        chroms_to_compare : list, optional
            A list of chromosomes to compare (Default is All)
        output_dir : str, optional
            Directory to output data
        do_output_graph : bool, optional
            Whether to output graph used for comparison (Default is False)

        Returns
        ------
        dict[str, float]
            Contains comparison values based on EMD and Jensen-Shannon formulas
        """
        # Default: Compare all the chromosomes
        if chroms_to_compare is None:
            chroms_to_compare = list(self.chrom_dict.keys())

        comparison_name = f'{self.sample_name}_{o_loop_data.sample_name}'
        param_str = f'{window_size}.{bin_size}.{num_peaks}'
        os.makedirs(f'{output_dir}/{param_str}/comparisons/{comparison_name}',
                    exist_ok=True)
        os.makedirs(f'{output_dir}/{param_str}/scores/windows', exist_ok=True)
        os.makedirs(f'{output_dir}/{param_str}/scores/chromosomes',
                    exist_ok=True)

        if do_output_graph:
            with open(f'{output_dir}/{param_str}/comparisons/{comparison_name}/'
                      f'chr1_graph.txt', 'w') as out_file:
                out_file.write(
                    f'sample_name\tchrom_name\twindow_start\twindow_end\tgraph\n')

        chrom_score_dict = {}
        log.info(f'Chromosomes to compare: {chroms_to_compare}')
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
            chrom_size = self.chrom_dict[chrom_name].size
            value_dict_list = []  # Contains comparison values for each window
            numb_windows = math.ceil(chrom_size / window_size)
            window_starts = []
            window_ends = []

            with open(f'{output_dir}/{param_str}/comparisons/{comparison_name}/'
                      f'{chrom_name}.txt', 'w') as out_file:
                out_file.write(
                    f'sample_name\tchrom_name\twindow_start\twindow_end\tnumb_nonzeros\tnumb_zeros\n')

            # Compare all windows in chromosome
            for k in range(numb_windows):
                window_start = window_size * k
                window_end = window_size * (k + 1)
                if window_end > chrom_size:
                    window_end = chrom_size
                window_starts.append(window_start)
                window_ends.append(window_end)

                value_dict_list.append(
                    self.chrom_dict[chrom_name].compare(
                        o_loop_data.chrom_dict[chrom_name], window_start,
                        window_end, window_size, bin_size, num_peaks,
                        output_dir=output_dir, do_output_graph=do_output_graph))

            emd_values = [x['emd_value'] for x in value_dict_list]
            j_values = [x['j_value'] for x in value_dict_list]
            weights = [x['w'] for x in value_dict_list]
            chrom_comp_values = []
            try:  # Weigh value from each window according to max loop in graph
                chrom_comp_values.append(
                    np.average(emd_values, weights=weights))
                chrom_comp_values.append(np.average(j_values, weights=weights))
            except ZeroDivisionError:  # sum of weights == 0
                log.exception(f"No loops were found in either graphs. Skipping"
                              f"{chrom_name}")
                continue

            chrom_score_dict[chrom_name] = chrom_comp_values
            log.debug(f'{chrom_name} comp values: {chrom_comp_values}')

            with open(f'{output_dir}/{param_str}/scores/windows/'
                      f'{comparison_name}_{chrom_name}.txt', 'w') as out_file:
                out_file.write(
                    f'sample_name\tchrom_name\twindow_start\twindow_end\t'
                    f'emd_value\tj_value\tweight\temd_dist\tj_divergence\n')
                for i, value_dict in enumerate(value_dict_list):
                    window_start = window_starts[i]
                    window_end = window_ends[i]
                    emd_value = value_dict['emd_value']
                    j_value = value_dict['j_value']
                    emd_dist = value_dict['emd_dist']
                    j_divergence = value_dict['j_divergence']
                    weight = value_dict['w']
                    out_file.write(f'{chrom_name}\t{window_start}\t{window_end}'
                                   f'\t{emd_value}\t{j_value}\t{weight}\t'
                                   f'{emd_dist}\t{j_divergence}\n')

        with open(f'{output_dir}/{param_str}/scores/chromosomes/'
                  f'{comparison_name}.txt', 'w') as out_file:
            out_file.write(f'chrom_name\temd_value\tj_value\n')
            for chrom_name, score_dict in chrom_score_dict.items():
                emd_value = score_dict[0]
                j_value = score_dict[1]
                out_file.write(f'{chrom_name}\t{emd_value}\t{j_value}\n')

        # Weigh value from each chromosome equally
        avg_value = {
            'emd_value': np.mean([x[0] for x in chrom_score_dict.values()]),
            'j_value': np.mean([x[1] for x in chrom_score_dict.values()])
        }
        log.debug(avg_value)

        return avg_value
