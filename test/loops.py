import sys
sys.path.append('..')
from genome_loop_data import GenomeLoopData, WINDOW_SIZE
import os
import numpy as np

DATA_DIR = '../data'

chr1_size = 248956422
chr2_size = 242193529
chr_size = chr1_size
CHROM_NAME = 'chr1'


def compare_loops(loop_info_dict):
    len_loop = len(loop_info_dict)
    keys = list(loop_info_dict.keys())
    for i in range(len_loop):
        for j in range(i + 1, len_loop, 1):
            key1 = keys[i]
            key2 = keys[j]
            loop1 = loop_info_dict[key1]
            loop2 = loop_info_dict[key2]
            chrom_loop1 = loop1.chrom_dict[CHROM_NAME]
            chrom_loop2 = loop2.chrom_dict[CHROM_NAME]
            print(f'Comparing {loop1.name} vs. {loop2.name}:')
            out_file = open(f'{loop1.name}_V_{loop2.name}.txt', 'w+')
            out_file.write(f'Comparing {loop1.name} vs. {loop2.name}:\n')

            j_values = []
            w_values = []
            for k in range(int(chr_size / WINDOW_SIZE)):
                values = chrom_loop1.compare(chrom_loop2, k)
                # print(f'Window Index: {k}\n'
                #       f'Jensen-Shannon: {values[0]}')
                #       f'Wasserstein: {values[1]}')
                out_file.write(f'Window Index: {k}\n'
                               f'Jensen-Shannon: {values[0]}\n')
                if values[0] != -1:
                    j_values.append(values[0])

                if len(values) > 1:
                    w_values.append(values[1])

            print(f'Jensen-Shannon Average: {np.mean([x for x in j_values])}')
            out_file.write(f'Jensen-Shannon Average: {np.mean(j_values)}\n')
            if len(w_values) > 0:
                print(f'Wasserstein Average: {np.mean(w_values)}')

            out_file.close()

        print()


def read_loops(data_directory, bin_size):
    loop_info_dict = {}
    for filename in os.listdir(data_directory):
        if filename.endswith('.BE3'):
            file_path = os.path.join(data_directory, filename)
            loopInfo = GenomeLoopData(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                      file_path, bin_size, wanted_chrom_name=CHROM_NAME)
            loop_info_dict[loopInfo.name] = loopInfo

    return loop_info_dict


def main():
    for datatype in ['CTCF', 'RNAPII']:
        loop_info_list = read_loops(DATA_DIR + f'/{datatype}/small')
        compare_loops(loop_info_list)

        loop_info_list = read_loops(DATA_DIR + f'/{datatype}/big')
        compare_loops(loop_info_list)
        break


if __name__ == '__main__':
    main()
