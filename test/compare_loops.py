import os
import numpy as np
import csv
import sys

sys.path.append('..')
from genome_loop_data import GenomeLoopData, WINDOW_SIZE

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'
CHROM_SIZES_FILE = f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes'

chr1_size = 248956422
chr2_size = 242193529
chr_size = chr1_size
CHROM_NAME = 'chr1'

BIN_SIZE = 10000


def compare_loops(loop_1, loop_2, bin_size, window_size, window_index=None):
    j_values = []
    w_values = []
    for k in range(int(chr_size / window_size)):

        if window_index is not None:
            k = window_index

        window_start = window_size * k
        window_end = window_size * (k + 1)
        if window_end > chr1_size:
            window_end = chr1_size

        values = loop_1.compare(loop_2, window_start, window_end, bin_size)
        j_values.append(values[0])

        if len(values) > 1:
            w_values.append(values[1])

        if window_index is not None:
            break

    # Make j_value range from -1 to 1
    j_value = 2 * (1 - np.mean([x for x in j_values])) - 1

    return j_value


def read_loops(datatype, bin_size=BIN_SIZE, hiseq=False, min_loop_value=None):
    loop_info_dict = {}
    if hiseq:
        data_directory = DATA_DIR + f'/{datatype}/hiseq'
    else:
        data_directory = DATA_DIR + f'/{datatype}/miseq'

    for filename in os.listdir(data_directory):
        if filename.endswith('.BE3'):
            file_path = os.path.join(data_directory, filename)
            loopInfo = GenomeLoopData(CHROM_SIZES_FILE, file_path, bin_size,
                                      wanted_chrom_name=CHROM_NAME, hiseq=hiseq,
                                      min_loop_value=min_loop_value)
            loop_info_dict[loopInfo.name] = loopInfo

    return loop_info_dict


def main():
    for datatype in ['CTCF', 'RNAPII']:
        loop_info_list = read_loops(datatype)
        # compare_loops(loop_info_list)

        loop_info_list = read_loops(datatype, hiseq=True)
        # compare_loops(loop_info_list)
        break


if __name__ == '__main__':
    main()
