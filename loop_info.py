from chrom_loops import ChromLoopData
import os


class LoopInfo:

    def __init__(self, chrom_size_file_path, in_file_path, bin_size=50000,
                 wanted_chrom_name='chr1'):
        self.name = os.path.basename(in_file_path)

        self.chrom_dict = {}
        with open(chrom_size_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if wanted_chrom_name != chrom_name:
                    continue
                self.chrom_dict[chrom_name] = ChromLoopData(chrom_name,
                                                            int(line[1]),
                                                            bin_size)

        with open(in_file_path) as in_file:
            for line in in_file:
                line = line.strip().split()
                chrom_name = line[0]
                if chrom_name not in self.chrom_dict:
                    continue

                loop_start = (int(line[1]) + int(line[2])) / 2
                loop_end = (int(line[4]) + int(line[5])) / 2
                loop_value = line[6]

                self.chrom_dict[chrom_name].add_loop(loop_start, loop_end,
                                                     loop_value)

        for chrom_name in self.chrom_dict:
            self.chrom_dict[chrom_name].create_graph()
