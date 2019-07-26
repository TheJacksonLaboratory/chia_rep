import sys

sys.path.append('../pyBedGraph')
from pyBedGraph import BedGraph
import numpy as np
import os
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

# DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'
DATA_DIR = '/home/hirwo/Documents/chia_pet'


def create_test_cases(num_tests, interval_size):
    np.random.seed(RANDOM_SEED)

    test_cases = np.random.randint(CHROM_START, CHROM_END - interval_size,
                                   num_tests, dtype=np.int32)
    test_cases = np.vstack((test_cases, test_cases + interval_size))

    return test_cases


def get_stats(bedGraph, test_cases):
    print("Searching for noise threshold ...")
    chrom = bedGraph.chromosome_map['chr1']
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
            chrom.value_map[i] = 0

    num_samples = 0
    for i in range(chrom.num_intervals):
        num_samples += chrom.intervals[1][i] - chrom.intervals[0][i]

    mean_list = bedGraph.stats(start_list=test_cases[0], end_list=test_cases[1],
                               chrom_name='chr1')

    return chrom.num_samples, mean_list


def compare_bedGraphs(bedGraph1_stat, bedGraph2_stat, out_file, test_cases):
    info_str = f"Comparing {bedGraph1_stat['name']} vs. {bedGraph2_stat['name']}"
    print(info_str)
    similarity_list = []
    ms_error_list = []
    num_counted = 0
    for i in range(mean_length):
        m1 = mean1_list[i]
        m2 = mean2_list[i]

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


def clean_bedGraphs(bedGraph_list, no_modify=False):

    for bedGraph in bedGraph_list:
        bedGraph.load_chrom_data('chr1')

        remove_noisy_peaks(bedGraph, no_modify)

        if not no_modify:
            with open(f'{DATA_DIR}/hg38.blacklist.bed') as in_file:
                for line in in_file:
                    line = line.strip().split()
                    if line[0] not in bedGraph.chromosome_map:
                        continue

                    chrom = bedGraph.chromosome_map[line[0]]
                    start = int(line[1])
                    end = int(line[2])

                    index = chrom.index_list[start]
                    chrom.value_map[index] = 0
                    chrom.index_list[start:end] = -1

        if len(bedGraph_list) > 1:
            bedGraph.free_chrom_data('chr1')


def compare_dataset(bedGraph_list, test_cases, datatype):

    bedGraph_stats_list = []
    for bedGraph in bedGraph_list:
        bedGraph.load_chrom_data('chr1')
        clean_bedGraphs([bedGraph])

        num_samples, mean_list = get_stats(bedGraph, test_cases)
        bedGraph_stats = {
            'num_samples': num_samples,
            'mean_list': mean_list,
            'name': bedGraph.name
        }
        bedGraph_stats_list.append(bedGraph_stats)

        bedGraph.free_chrom_data('chr1')

    out_file = open(f'{datatype}_results.txt', 'w+')
    for i in range(len(bedGraph_list)):
        for j in range(i + 1, len(bedGraph_list), 1):
            compare_bedGraphs(bedGraph_stats_list[i], bedGraph_stats_list[j],
                              out_file, test_cases)
    out_file.close()


def read_small_bedGraphs(datatype):
    data_directory = f'{DATA_DIR}/{datatype}/small/'
    bedGraph_list = []
    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph'):
            file_path = data_directory + filename
            bedGraph = BedGraph(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path, 'chr1', ignore_missing_bp=False,
                                min_value=6)
            bedGraph_list.append(bedGraph)

    return bedGraph_list


def read_big_bedGraphs(datatype):
    data_directory = f'{DATA_DIR}/{datatype}/big/'
    bedGraph_list = []
    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph'):
            file_path = data_directory + filename
            bedGraph = BedGraph(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path, 'chr1', ignore_missing_bp=False,
                                min_value=26)
            bedGraph_list.append(bedGraph)

    return bedGraph_list


# noisy peaks have small groups of intervals
def remove_noisy_peaks(bedGraph, no_modify=False):
    print(bedGraph.name)
    chrom = bedGraph.chromosome_map['chr1']
    average_value = chrom.total_value / chrom.num_intervals
    min_threshold = average_value * 1.5
    num_noisy_peaks = 0
    interval_group = [0]
    interval_group_values = [chrom.value_map[0]]
    for i in range(1, chrom.num_intervals, 1):

        # new group of intervals
        if chrom.intervals[0][i] != chrom.intervals[1][interval_group[-1]]:
            group_size = len(interval_group)
            group_max_value = np.max(interval_group_values)

            if group_size < 20 and group_max_value > min_threshold:
                max_index = -1
                for j in range(group_size):
                    if interval_group_values[j] == group_max_value:
                        max_index = j
                        break

                interval_index = interval_group[max_index]
                interval_start = chrom.intervals[0][interval_index]
                interval_end = chrom.intervals[1][interval_index]

                try:
                    before_index = chrom.index_list[interval_start - 25]
                except IndexError:
                    before_index = -1
                try:
                    after_index = chrom.index_list[interval_end + 25]
                except IndexError:
                    after_index = -1

                if (before_index == -1 or chrom.value_map[before_index] == 0)\
                    and (after_index == -1 or chrom.value_map[after_index] == 0):
                    print([chrom.intervals[0][x] for x in interval_group], group_max_value)

                    if not no_modify:
                        for x in interval_group:
                            chrom.value_map[x] = 0
                    num_noisy_peaks += 1

            interval_group.clear()
            interval_group_values.clear()

        interval_group.append(i)
        interval_group_values.append(chrom.value_map[i])

    print(f'Noisy Peaks: {num_noisy_peaks}\n')


def main():
    test_cases = create_test_cases(NUM_TESTS, INTERVAL_SIZE)

    for datatype in ['CTCF', 'RNAPII']:
        bedGraph_list = read_small_bedGraphs(datatype)
        # bedGraph_list += read_big_bedGraphs(datatype)

        compare_dataset(bedGraph_list, test_cases, datatype)
        break


if __name__ == "__main__":
    main()
