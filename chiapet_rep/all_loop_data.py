from collections import OrderedDict

from prettytable import PrettyTable

from .chrom_loop_data import ChromLoopData
import numpy as np
import os
import logging

VERSION = 2
log = logging.getLogger()
log_bin = logging.getLogger('bin')

# Missing in many miseq peak files
CHROM_TO_IGNORE = 'chrY'

DEFAULT_PEAK_PERCENT = 0.2


class AllLoopData:

    def __init__(self, chrom_size_file_path, in_file_path, peak_file_path,
                 bedgraph, is_hiseq, peak_percent_kept=DEFAULT_PEAK_PERCENT,
                 wanted_chroms=None, min_loop_value=0):

        log.debug(locals())

        self.is_hiseq = is_hiseq

        self.sample_name = os.path.basename(in_file_path).split('.')[0]

        self.chrom_dict = {}
        with open(chrom_size_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if wanted_chroms and chrom_name not in wanted_chroms:
                    continue

                if chrom_name == CHROM_TO_IGNORE:
                    continue

                self.chrom_dict[chrom_name] = \
                    ChromLoopData(chrom_name, line[1], self.sample_name)

        # Read in loops
        with open(in_file_path) as in_file:
            loop_anchor_list = []
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chrom_name not in self.chrom_dict:
                    continue

                loop_end1 = int(line[4])
                loop_end2 = int(line[5])
                loop_start1 = int(line[1])
                loop_start2 = int(line[2])
                loop_value = int(line[6])

                if loop_value < min_loop_value:
                    continue

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
            if not self.chrom_dict[chrom_name].finish_init(bedgraph,
                                                           peak_file_path,
                                                           peak_percent_kept):
                to_remove.append(chrom_name)

        # Chromosomes with no loops or other random problems
        for chrom_name in to_remove:
            del self.chrom_dict[chrom_name]

    # To speed up testing process by avoiding loading bedgraphs every time
    # Untested
    def preprocess(self, peak_file_path,
                   peak_percent_kept=DEFAULT_PEAK_PERCENT):
        for chrom in self.chrom_dict.values():
            chrom.preprocess(peak_file_path, peak_percent_kept)

    def compare(self, o_loop_data, bin_size, window_size, window_index=None,
                wanted_chroms=None):

        # Default: Compare all the chromosomes
        if wanted_chroms is None:
            wanted_chroms = list(self.chrom_dict.keys())

        chrom_value_table = \
            PrettyTable(['chrom', 'graph_type', 'rep', 'w_rep'])

        chrom_value_list = []
        for chrom_name in wanted_chroms:

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
            numb_windows = int(chrom_size / window_size)
            if numb_windows == 0:
                numb_windows = 1
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
                        window_end, bin_size,
                        self.is_hiseq == o_loop_data.is_hiseq))

                if window_index is not None:
                    break

            values = [x['rep'] for x in value_dict_list]
            chrom_value = {
                'graph_type': value_dict_list[0]['graph_type'],
                'rep': np.mean(values),
            }
            try:  # Weigh value from each bin according to max loop in graph
                chrom_value['w_rep'] = \
                    np.average(values, weights=[x['w'] for x in
                                                value_dict_list])
            except ZeroDivisionError:  # sum of weights == 0
                log.exception(f"No loops were found in either graphs. Skipping"
                              f"{chrom_name}")
                continue

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
        log.debug(avg_value)

        return avg_value
