import sys

sys.path.append('../../pyBedGraph')
sys.path.append('..')
from pyBedGraph import BedGraph
import numpy as np
import os
import remove_noise
import matplotlib.pyplot as plt

RANDOM_SEED = 1
NUM_TESTS = 1000000
INTERVAL_SIZE = 1000
chr1_size = 248956422

CHROM_START = 74443235
CHROM_END = 75673623

CHROM_START = 0
CHROM_END = chr1_size

MAX_NOISE_VALUE = 100

DATA_DIR = '../data'
BIG_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'


def create_test_cases(interval_size=INTERVAL_SIZE):

    test_cases = np.empty(int(chr1_size / interval_size), dtype=np.int32)
    for i in range(test_cases.size):
        test_cases[i] = interval_size * i
    test_cases = np.vstack((test_cases, test_cases + interval_size))

    return test_cases


def get_stats(bedGraph, test_cases):
    chrom = bedGraph.chromosome_map['chr1']
    '''print("Searching for noise threshold ...")
    freq_array = np.zeros(MAX_NOISE_VALUE, dtype=np.int32)
    for i in range(chrom.num_intervals):
        index = chrom.value_map[i]
        if index < MAX_NOISE_VALUE:
            index = int(index)
            freq = chrom.intervals[1][i] - chrom.intervals[0][i]
            freq_array[index] += freq

    noise_level = 0
    median_value = np.median(freq_array)
    max_value = np.max(freq_array)
    max_index = np.where(freq_array == max_value)[0][0]
    for i in range(max_index, MAX_NOISE_VALUE, 1):

        if freq_array[i] / median_value > 10:
            continue

        noise_level = i
        break
    if median_value > 10000:
        noise_level = 25
    elif median_value < 2000:
        noise_level = 5
    else:
        print(f"unknown median value: {median_value}")

    if noise_level == 0:
        print("Noise value was not found")
        plt.plot([x for x in range(MAX_NOISE_VALUE)], freq_array)
        plt.show()
        exit(-1)

    print(f"Signals at or below {noise_level} are considered noise")
    print(f'Median value was {median_value}')
    plt.plot([x for x in range(MAX_NOISE_VALUE)], freq_array)
    plt.show()

    # zero out noise
    for i in range(chrom.num_intervals):
        if chrom.value_map[i] <= noise_level:
            chrom.value_map[i] = 0'''

    bedGraph_stats = {
        'num_samples': chrom.num_samples,
        'mean_list': bedGraph.stats(start_list=test_cases[0],
                                    end_list=test_cases[1],
                                    chrom_name='chr1'),
        'name': bedGraph.name,
        'max_value': np.max(chrom.value_map)
    }
    return bedGraph_stats


def compare_bedGraphs(bedGraph1_stat, bedGraph2_stat, out_file, test_cases):
    info_str = f"Comparing {bedGraph1_stat['name']} vs. {bedGraph2_stat['name']}"
    print(info_str)
    out_file.write(info_str + '\n')

    mean1_list = bedGraph1_stat['mean_list']
    mean2_list = bedGraph2_stat['mean_list']

    assert len(mean1_list) == len(mean2_list)
    mean_length = len(mean1_list)

    ratio = bedGraph1_stat['num_samples'] / bedGraph2_stat['num_samples']

    similarity_list = []
    ms_error_list = []
    num_counted = 0
    for i in range(mean_length):
        m1 = mean1_list[i]
        m2 = mean2_list[i]

        # account for bedGraphs that were sampled more
        m2 *= ratio
        max_m = max(m1, m2)

        if max_m != 0:
            similarity_error = abs(m1 - m2) * max_m
            mse = (m1 - m2) * (m1 - m2) * max_m * max_m

            # test_case_str = f'{test_cases[0][i]} {test_cases[1][i]}\n{m1} ' \
            #                f'{m2}\n{percentage_error} {mse}\n'
            # print(test_case_str)
            # out_file.write(test_case_str + '\n')

            similarity_list.append(similarity_error)
            ms_error_list.append(mse)

            num_counted += 1

    num_counted_str = f'Num Counted: {num_counted}'
    print(num_counted_str)
    out_file.write(num_counted_str + '\n')

    stat_str = f'Ratio: {ratio}\nSimilarity (inverted): {np.mean(similarity_list)}' \
               f'\nMean Squared Error (ish): {np.mean(ms_error_list)}\n'
    print(stat_str)
    out_file.write(stat_str + '\n')


def compare_dataset(bedGraph_list, datatype, test_cases=None):

    if test_cases is None:
        test_cases = create_test_cases()

    bedGraph_stats_list = []
    for bedGraph in bedGraph_list:
        remove_noise.remove_noise([bedGraph])

        bedGraph_stats = get_stats(bedGraph, test_cases)
        bedGraph_stats_list.append(bedGraph_stats)

        bedGraph.free_chrom_data('chr1')

    out_file = open(f'{datatype}_results.txt', 'w+')
    for i in range(len(bedGraph_list)):
        for j in range(i + 1, len(bedGraph_list), 1):
            compare_bedGraphs(bedGraph_stats_list[i], bedGraph_stats_list[j],
                              out_file, test_cases)
    out_file.close()


def read_bedGraphs(data_directory, min_value=-1):
    bedGraph_dict = {}
    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph'):
            file_path = data_directory + filename
            bedGraph = BedGraph(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path, 'chr1', ignore_missing_bp=False,
                                min_value=min_value)
            bedGraph_dict[bedGraph.name] = bedGraph

    return bedGraph_dict


def read_small_bedGraphs(datatype):
    data_directory = f'{DATA_DIR}/{datatype}/small/fixed/'

    #  return read_bedGraphs(data_directory, min_value=6)
    return read_bedGraphs(data_directory)


def read_big_bedGraphs(datatype):
    data_directory = f'{BIG_DATA_DIR}/{datatype}/big/'

    #  return read_bedGraphs(data_directory, min_value=26)
    return read_bedGraphs(data_directory)


def main():
    test_cases = create_test_cases()

    for datatype in ['CTCF', 'RNAPII']:
        bedGraph_list = read_small_bedGraphs(datatype)
        # bedGraph_list += read_big_bedGraphs(datatype)

        compare_dataset(bedGraph_list, datatype, test_cases)
        break


if __name__ == "__main__":
    main()
