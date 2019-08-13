from chrom_loop_data import ChromLoopData
import numpy as np
import os
import logging as log

WINDOW_SIZE = 3000000
BIN_SIZE = 1000
DEFAULT_CHROM_NAME = 'chr1'
HISEQ_MIN_LOOP_VALUE = 10
MISEQ_MIN_LOOP_VALUE = 10

PEAK_FILE = '/media/hirwo/extra/jax/data/chia_pet/peaks/{0}.for.BROWSER.spp.z6.broadPeak'

VERSION = 3


class GenomeLoopData:

    def __init__(self, chrom_size_file_path, in_file_path, bin_size=BIN_SIZE,
                 window_size=WINDOW_SIZE, wanted_chrom_name=DEFAULT_CHROM_NAME,
                 min_loop_value=None, hiseq=False, is_CTCF=False):

        log.info(f"Reading in {in_file_path} ...")
        log.info(f"Min Loop Value: {min_loop_value}")
        log.info(f"Is CTCF: {is_CTCF}")

        self.name = os.path.basename(in_file_path).split('.')[0]
        self.sample_name = self.name

        self.chrom_dict = {}
        with open(chrom_size_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if wanted_chrom_name != chrom_name:
                    continue
                self.chrom_dict[chrom_name] = ChromLoopData(chrom_name, line[1],
                                                            self.name)

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

                if min_loop_value is not None and loop_value < min_loop_value:
                    continue

                self.chrom_dict[chrom_name].add_loop(loop_start1, loop_start2,
                                                     loop_end1, loop_end2,
                                                     loop_value)

                start_interval = loop_start2 - loop_start1
                end_interval = loop_end2 - loop_end1

                loop_anchor_list.append(start_interval)
                loop_anchor_list.append(end_interval)

            log.info(f'Anchor mean width: {np.mean(loop_anchor_list)}')

        for chrom_name in self.chrom_dict:
            self.chrom_dict[chrom_name].finish_init()

    def find_loop_anchor_points(self, bedGraph, chrom_name=DEFAULT_CHROM_NAME):
        self.chrom_dict[chrom_name].find_loop_anchor_points(bedGraph)

    def adjust_with_bedGraph(self, bedGraph):
        return
        log.info(f'Adjusting {self.name} with bedgraph file')
        for name in self.chrom_dict:
            self.chrom_dict[name].adjust_with_bedGraph(bedGraph.get_chrom(name))

    def filter_with_bedGraph(self, chrom_name):
        log.info(f'Filtering {chrom_name} of {self.name} with bedgraph file')
        self.chrom_dict[chrom_name].filter_with_bedGraph(
            PEAK_FILE.format(self.name))

    def compare(self, o_loop_data, window_start, window_end,
                bin_size, chrom_name=DEFAULT_CHROM_NAME):
        if chrom_name is None:
            for chrom_name in self.chrom_dict:
                if chrom_name not in o_loop_data.chrom_dict:
                    continue
                self.chrom_dict[chrom_name].compare(
                    o_loop_data.chrom_dict[chrom_name])

        return self.chrom_dict[chrom_name].compare(
            o_loop_data.chrom_dict[chrom_name], window_start, window_end,
            bin_size)
