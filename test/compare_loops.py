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


def output_csv(scores, output_file):
    with open(output_file, 'w+') as out_file:
        header = list(scores.keys())
        header.insert(0, 'Sample Name')
        writer = csv.DictWriter(out_file, fieldnames=header)

        writer.writeheader()
        for sample_name, sample_scores in scores.items():
            writer.writerow(sample_scores)


def compare_loops(loop_info_dict, interval_size, window_size, window_index=None):
    len_loop = len(loop_info_dict)
    keys = list(loop_info_dict.keys())
    scores = {}
    for key in keys:
        scores[key] = {}
        scores[key][key] = 1
        scores[key]['Sample Name'] = key

    for i in range(len_loop):
        for j in range(i + 1, len_loop, 1):
            key1 = keys[i]
            key2 = keys[j]
            print(f'Comparing {key1} vs. {key2}:')
            loop1 = loop_info_dict[key1]
            loop2 = loop_info_dict[key2]

            j_values = []
            w_values = []
            if window_index is None:
                # misses comparing part of the chromosome
                for k in range(int(chr_size / window_size)):
                    values = loop1.compare(loop2, window_size * k,
                                           window_size * (k + 1), interval_size)
                    # print(f'Window Index: {k}\n'
                    #       f'Jensen-Shannon: {values[0]}')
                    #       f'Wasserstein: {values[1]}')
                    # if values[0] != -1:
                    #     pass
                    j_values.append(values[0])

                    if len(values) > 1:
                        w_values.append(values[1])
            else:
                values = loop1.compare(loop2, window_size * window_index,
                                       window_size * (window_index + 1),
                                       interval_size)
                j_values.append(values[0])

            j_value = 2 * (1 - np.mean([x for x in j_values])) - 1
            if len(w_values) > 0:
                w_mean = np.mean(w_values)
                # print(f'Wasserstein Average: {w_mean}')

            # print(f'Jensen-Shannon value: {j_value}')

            scores[key1][key2] = j_value
            scores[key2][key1] = j_value

        # print()

    return scores


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
        compare_loops(loop_info_list)

        loop_info_list = read_loops(datatype, hiseq=True)
        compare_loops(loop_info_list)
        break


if __name__ == '__main__':
    main()
