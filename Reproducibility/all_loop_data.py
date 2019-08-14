from .chrom_loop_data import ChromLoopData
import numpy as np
import os
import logging

VERSION = 1
log = logging.getLogger(__name__.split('.')[-1])

# Missing in many miseq peak files
CHROM_TO_IGNORE = 'chrY'


class AllLoopData:

    def __init__(self, chrom_size_file_path, in_file_path, peak_file_path,
                 bedgraph, is_hiseq, peak_percent_kept=0.2,
                 wanted_chroms=None, min_loop_value=0):

        log.info(locals())

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

            log.info(f'Anchor mean width: {np.mean(loop_anchor_list)}')

        to_remove = []
        for chrom_name in self.chrom_dict:
            if not self.chrom_dict[chrom_name].finish_init(bedgraph,
                                                           peak_file_path,
                                                           peak_percent_kept):
                to_remove.append(chrom_name)

        for chrom_name in to_remove:
            del self.chrom_dict[chrom_name]

    def compare(self, o_loop_data, window_start, window_end,
                bin_size, wanted_chroms=None):

        # Compare all the chromosomes
        if wanted_chroms is None:
            wanted_chroms = list(self.chrom_dict.keys())

        chrom_rep = {}

        for chrom_name in wanted_chroms:

            if chrom_name not in self.chrom_dict:
                log.warning(f'{chrom_name} is not in {self.sample_name}')
                continue

            # Check that chrom name is also in other
            if chrom_name not in o_loop_data.chrom_dict:
                log.warning(f'{chrom_name} is in {self.sample_name} but'
                            f'not in {o_loop_data.sample_name}')
                continue

            chrom_rep[chrom_name] = \
                self.chrom_dict[chrom_name].compare(
                    o_loop_data.chrom_dict[chrom_name], window_start,
                    window_end, bin_size, self.is_hiseq == o_loop_data.is_hiseq)

        value_dict_list = []
        for chrom in chrom_rep:
            log.debug(chrom)
            for i in range(0, len(chrom_rep[chrom])):
                if len(value_dict_list) == i:
                    value_dict_list.append([])
                value_dict_list[i].append(chrom_rep[chrom][i])

        log.debug(value_dict_list)
        genome_rep_values = []
        for value_type_list in value_dict_list:
            log.debug(value_type_list)
            avg_value = dict.fromkeys(value_type_list[0].keys(), 0)
            for value_dict in value_type_list:
                for k in value_dict:
                    if isinstance(value_dict[k], str):
                        avg_value[k] = value_dict[k]
                        continue
                    avg_value[k] += value_dict[k]

            for k in avg_value:
                if isinstance(avg_value[k], str):
                    continue

                avg_value[k] /= len(value_type_list)

            log.debug(avg_value)
            genome_rep_values.append(avg_value)

        log.debug(genome_rep_values)
        return genome_rep_values
