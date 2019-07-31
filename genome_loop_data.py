from chrom_loop_data import ChromLoopData
import numpy as np
import os

WINDOW_SIZE = 5000000
BIN_SIZE = 1000
DEFAULT_CHROM_NAME = 'chr1'


class GenomeLoopData:

    def __init__(self, chrom_size_file_path, in_file_path, bin_size=BIN_SIZE,
                 window_size=WINDOW_SIZE, wanted_chrom_name=DEFAULT_CHROM_NAME):
        print(f"Bin size: {bin_size}")

        self.name = os.path.basename(in_file_path).split('.')[0]

        self.chrom_dict = {}
        with open(chrom_size_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if wanted_chrom_name != chrom_name:
                    continue
                self.chrom_dict[chrom_name] = ChromLoopData(chrom_name,
                                                            int(line[1]),
                                                            bin_size,
                                                            window_size)

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

                if loop_value < 10:
                    continue
                    pass

                self.chrom_dict[chrom_name].add_loop(loop_start1, loop_start2,
                                                     loop_end1, loop_end2,
                                                     loop_value)

                start_interval = loop_start2 - loop_start1
                end_interval = loop_end2 - loop_end1
                if start_interval > bin_size:
                    # print("Over: ", start_interval, loop_value)
                    pass

                if end_interval > bin_size:
                    # print("Over: ", end_interval, loop_value)
                    pass

                loop_anchor_list.append(start_interval)
                loop_anchor_list.append(end_interval)

            print(f'Anchor mean width: {np.mean(loop_anchor_list)}')

    def adjust_with_bedGraph(self, bedGraph):
        print(f'Adjusting {self.name} with bedgraph file')
        for name in self.chrom_dict:
            self.chrom_dict[name].adjust_with_bedGraph(bedGraph.get_chrom(name))
