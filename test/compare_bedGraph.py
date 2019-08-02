import sys

sys.path.append('../../pyBedGraph')
sys.path.append('..')
from pyBedGraph import BedGraph
import numpy as np
import os
import matplotlib.pyplot as plt
from math import sqrt

WINDOW_SIZE = 100000
INTERVAL_SIZE = 1000
chr1_size = 248956422

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'


def create_test_cases(interval_size=INTERVAL_SIZE, start=0, end=chr1_size,
                      step=INTERVAL_SIZE):
    test_cases = np.empty(int((end - start) / step) - interval_size,
                          dtype=np.int32)
    for i in range(test_cases.size):
        test_cases[i] = i * step + start
    test_cases = np.vstack((test_cases, test_cases + interval_size))

    return test_cases


def get_stats(bedGraph_dict, test_cases):

    bedGraph_stats_dict = {}
    for name in bedGraph_dict:
        print(f'Getting stats for {name}')
        bedGraph = bedGraph_dict[name]
        bedGraph.load_chrom_data('chr1')
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

        bedGraph_stats_dict[name] = {
            'num_samples': chrom.num_samples,
            'mean_list': bedGraph.stats(start_list=test_cases[0],
                                        end_list=test_cases[1],
                                        chrom_name='chr1',
                                        stat='max'),
            'name': bedGraph.name,
            'max_value': np.max(chrom.value_map)
        }
        bedGraph.free_chrom_data('chr1')

    return bedGraph_stats_dict


def compare_bedGraphs_with_window(bedGraph1_stat, bedGraph2_stat, out_file, test_cases,
                      window_index):
    info_str = f"Comparing {bedGraph1_stat['name']} vs. {bedGraph2_stat['name']}"
    # print(info_str)
    out_file.write(info_str + '\n')

    mean1_list = bedGraph1_stat['mean_list']
    mean2_list = bedGraph2_stat['mean_list']

    window_start = window_index * WINDOW_SIZE
    window_end = (window_index + 1) * WINDOW_SIZE

    mean_start = int(window_start / INTERVAL_SIZE)
    mean_end = int(window_end / INTERVAL_SIZE)
    mean_length = mean_end - mean_start

    assert len(mean1_list) == len(mean2_list)

    # ratio = bedGraph1_stat['num_samples'] / bedGraph2_stat['num_samples']
    mean1_sample_size = np.mean(mean1_list[mean_start:mean_end])
    mean2_sample_size = np.mean(mean2_list[mean_start:mean_end])
    if mean2_sample_size == 0 and mean1_sample_size == 0:
        return 1
    elif mean2_sample_size == 0:
        return -1

    mean1_window_max = np.max(mean1_list[mean_start:mean_end])
    mean2_window_max = np.max(mean2_list[mean_start:mean_end])
    mean_window_max = max(mean1_window_max, mean2_window_max)
    mean_window_max = max(bedGraph1_stat['max_value'], bedGraph2_stat['max_value'])

    ratio = mean1_sample_size / mean2_sample_size

    similarity_list = []
    ms_error_list = []
    num_counted = 0
    mean1_bigger = 0
    for i in range(mean_start, mean_end, 1):
        m1 = mean1_list[i]
        m2 = mean2_list[i]

        # account for bedGraphs that were sampled more
        m2 *= ratio
        max_m = max(m1, m2)

        if m1 > m2:
            mean1_bigger += 1

        if max_m != 0:
            # higher value == more error
            similarity_error = sqrt(abs(m1 - m2)) * \
                               (mean_window_max - max_m) / \
                               (sqrt(mean_window_max) * mean_window_max)
            mse = (m1 - m2) * (m1 - m2) * max_m

            if bedGraph1_stat['name'] == 'LHH0083H' and similarity_error > 0.5:
                out_file.write(
                    f'{test_cases[0][i]} - {test_cases[1][i]}, {similarity_error}\n')

            # test_case_str = f'{test_cases[0][i]} {test_cases[1][i]}\n{m1} ' \
            #                f'{m2}\n{percentage_error} {mse}\n'
            # print(test_case_str)
            # out_file.write(test_case_str + '\n')
        else:
            similarity_error = 0
            mse = 0

        num_counted += 1
        similarity_list.append(similarity_error)
        ms_error_list.append(mse)

    stat_str = f'{bedGraph1_stat["name"]} bigger: {mean1_bigger / mean_length}\n' \
               f'Num Counted: {num_counted}\n' \
               f'Ratio: {ratio}\n' \
               f'Window: {window_start} - {window_end}\n' \
               f'Max in window: {mean_window_max}\n' \
               f'Similarity: {2 * (1 - np.mean(similarity_list)) - 1}\n' \
               f'Mean Squared Error (ish): {np.mean(ms_error_list)}\n'
    # print(stat_str)
    out_file.write(stat_str + '\n')

    return 2 * (1 - np.median(similarity_list)) - 1


def compare_bedGraphs(bedGraph1_stat, bedGraph2_stat, out_file, test_cases):
    info_str = f"Comparing {bedGraph1_stat['name']} vs. {bedGraph2_stat['name']}"
    # print(info_str)
    out_file.write(info_str + '\n')

    mean1_list = bedGraph1_stat['mean_list']
    mean2_list = bedGraph2_stat['mean_list']

    assert len(mean1_list) == len(mean2_list)
    mean_len = len(mean1_list)

    ratio = bedGraph1_stat['num_samples'] / bedGraph2_stat['num_samples']

    similarity_list = []
    ms_error_list = []
    num_counted = 0
    mean1_bigger = 0
    for i in range(mean_len):
        m1 = mean1_list[i]
        m2 = mean2_list[i]

        # account for bedGraphs that were sampled more
        m2 *= ratio
        max_m = max(m1, m2)

        if m1 > m2:
            mean1_bigger += 1

        if max_m != 0:
            # higher value == more error
            similarity_error = abs(m1 - m2) / max_m
            mse = (m1 - m2) * (m1 - m2) * max_m

            if bedGraph1_stat['name'] == 'LHH0083H' and similarity_error > 0.5:
                out_file.write(
                    f'{test_cases[0][i]} - {test_cases[1][i]}, {similarity_error}\n')

            # test_case_str = f'{test_cases[0][i]} {test_cases[1][i]}\n{m1} ' \
            #                f'{m2}\n{percentage_error} {mse}\n'
            # print(test_case_str)
            # out_file.write(test_case_str + '\n')
        else:
            similarity_error = 0
            mse = 0

        num_counted += 1
        similarity_list.append(similarity_error)
        ms_error_list.append(mse)

    stat_str = f'{bedGraph1_stat["name"]} bigger: {mean1_bigger / mean_len}\n' \
               f'Num Counted: {num_counted}\n' \
               f'Ratio: {ratio}\n' \
               f'Similarity: {2 * (1 - np.mean(similarity_list)) - 1}\n' \
               f'Mean Squared Error (ish): {np.mean(ms_error_list)}\n'
    # print(stat_str)
    out_file.write(stat_str + '\n')

    return 2 * (1 - np.median(similarity_list)) - 1


def compare_dataset(bedGraph_dict, datatype, test_cases=None):
    if test_cases is None:
        test_cases = create_test_cases()

    bedGraph_stats_dict = get_stats(bedGraph_dict, test_cases)

    out_file = open(f'{datatype}_results.txt', 'w+')

    keys = list(bedGraph_dict.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys), 1):
            name1 = keys[i]
            name2 = keys[j]

            '''similarity_values = []
            for k in range(int(chr1_size / WINDOW_SIZE)):
                value = compare_bedGraphs(bedGraph_stats_dict[name1],
                                          bedGraph_stats_dict[name2], out_file,
                                          test_cases, k)
                similarity_values.append(value)
            print(f"{name1} - {name2} : {np.mean(similarity_values)}")'''
            value = compare_bedGraphs(bedGraph_stats_dict[name1],
                                      bedGraph_stats_dict[name2], out_file,
                                      test_cases)
            print(f"{name1} - {name2} : {value}")

        print()

    out_file.close()


def read_bedGraphs(datatype, hiseq=False, min_value=-1):
    bedGraph_dict = {}
    if hiseq:
        data_directory = f'{DATA_DIR}/{datatype}/hiseq/'
    else:
        data_directory = f'{DATA_DIR}/{datatype}/miseq/'

    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph'):
            file_path = data_directory + filename
            bedGraph = BedGraph(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path, 'chr1', ignore_missing_bp=False,
                                min_value=min_value)
            bedGraph_dict[bedGraph.name] = bedGraph

    return bedGraph_dict


def main():
    test_cases = create_test_cases()

    for datatype in ['CTCF', 'RNAPII']:
        bedGraph_list = read_bedGraphs(datatype)
        # bedGraph_list += read_bedGraphs(datatype, big=True)

        compare_dataset(bedGraph_list, datatype, test_cases)
        break


if __name__ == "__main__":
    main()
