import sys

sys.path.append('../../pyBedGraph')
sys.path.append('..')
from pyBedGraph import BedGraph
from scipy.stats import pearsonr, spearmanr
import numpy as np
import os
import matplotlib.pyplot as plt
from math import sqrt
from copy import deepcopy
import logging as log
from prettytable import PrettyTable

WINDOW_SIZE = 100000
INTERVAL_SIZE = 1000
chr1_size = 248956422

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'

VERSION = 1


def create_test_cases(interval_size=INTERVAL_SIZE, start=0, end=chr1_size,
                      step=INTERVAL_SIZE):
    num_test_cases = int((end - start) / step) - interval_size
    test_cases = np.arange(start, start + step * num_test_cases, step, dtype=np.int32)
    test_cases = np.vstack((test_cases, test_cases + interval_size))

    return test_cases


def get_stats(bedGraph_dict, test_cases, stat='mean'):

    bedGraph_stats_dict = {}
    for name in bedGraph_dict:
        print(f'Getting stats for {name}')
        bedGraph = bedGraph_dict[name]
        bedGraph.load_chrom_data('chr1')
        chrom = bedGraph.chromosome_map['chr1']

        bedGraph_stats_dict[name] = {
            'num_samples': chrom.num_samples,
            'mean_list': bedGraph.stats(start_list=test_cases[0],
                                        end_list=test_cases[1],
                                        chrom_name='chr1',
                                        stat=stat),
            'name': bedGraph.name,
            'max_value': np.max(chrom.value_map)
        }
        bedGraph.free_chrom_data('chr1')

    return bedGraph_stats_dict


def linear_regression(stat_lists):
    n = stat_lists[0].size

    squared_sum = np.sum(stat_lists * stat_lists, axis=1)
    x2_sum = squared_sum[0]
    xy_sum = np.sum(stat_lists[0] * stat_lists[1])
    sums = np.sum(stat_lists, axis=1)
    x_sum = sums[0]
    y_sum = sums[1]

    # y = a + b * x
    b = (n * xy_sum - x_sum * y_sum) / (n * x2_sum - x_sum * x_sum)
    a = (y_sum - b * x_sum) / n

    return a, b


def get_weighted_linear_error(a, b, stats_lists):
    error_arr = (stats_lists[1] - (a + b * stats_lists[0]))
    error_arr *= error_arr

    weights = np.empty(stats_lists[0].size)
    for i in range(stats_lists[1].size):
        if stats_lists[1][i] > stats_lists[0][i]:
            if stats_lists[0][i] == 0:
                stats_lists[0][i] = 0.01
            weights[i] = stats_lists[1][i] / stats_lists[0][i]
        else:
            if stats_lists[1][i] == 0:
                stats_lists[1][i] = 0.01
            weights[i] = stats_lists[0][i] / stats_lists[1][i]

    error_arr *= weights

    return np.mean(error_arr)


def get_adjusted_r2(r, num_points, num_variables=2):
    r2 = r * r
    return 1 - ((1 - r2) * (num_points - 1) / (num_points - num_variables - 1))


def my_pearson(stat_list1, stat_list2):
    mean1 = np.mean(stat_list1)
    mean2 = np.mean(stat_list2)

    xy_sum = 0
    x_sum = 0
    y_sum = 0
    for i in range(stat_list1.size):
        xy_sum += (stat_list1[i] - mean1) * (stat_list2[i] - mean2)
        x_sum += (mean1 - stat_list1[i]) * (mean1 - stat_list1[i])
        y_sum += (mean2 - stat_list2[i]) * (mean2 - stat_list2[i])

    x_sum = sqrt(x_sum)
    y_sum = sqrt(y_sum)

    try:
        assert x_sum * y_sum != 0
    except AssertionError:
        log.exception(f'{x_sum} {y_sum}')
        return 0

    return xy_sum / (x_sum * y_sum)


def get_coefficients(bedGraph1_stat, bedGraph2_stat, min_value):

    stat_list1 = bedGraph1_stat['mean_list']
    stat_list2 = bedGraph2_stat['mean_list']
    assert stat_list1.size == stat_list2.size

    num_samples1 = bedGraph1_stat['num_samples']
    num_samples2 = bedGraph2_stat['num_samples']
    ratio = num_samples1 / num_samples2
    modified_min_value = min_value * ratio

    to_remove = []
    for i in range(stat_list1.size):
        if stat_list2[i] < min_value and stat_list1[i] < modified_min_value:
            to_remove.append(i)
    stat_list1 = np.delete(stat_list1, to_remove)
    stat_list2 = np.delete(stat_list2, to_remove)

    orig_pearson_value, pearson_p = pearsonr(stat_list1, stat_list2)
    pearson_value = my_pearson(stat_list1, stat_list2)
    a, b = linear_regression(np.array([stat_list1, stat_list2]))
    weighted_error = get_weighted_linear_error(a, b, np.array([stat_list1, stat_list2]))
    adjusted_r2 = get_adjusted_r2(pearson_value, stat_list1.size)
    table = PrettyTable(['orig_r', 'my_r', 'adj_r2', 'error'])
    table.add_row([round(orig_pearson_value, 4), round(pearson_value, 4),
                  round(adjusted_r2, 4), round(weighted_error/1000000000, 4)])
    print(table)

    scatter_plot_dir = f'scatter_plots_{min_value}'
    if not os.path.isdir(scatter_plot_dir):
        os.mkdir(scatter_plot_dir)

    name1 = bedGraph1_stat['name']
    name2 = bedGraph2_stat['name']
    plt.scatter(stat_list1, stat_list2)
    plt.title(f'{name1} vs. {name2}, r={round(pearson_value, 3)}')
    plt.savefig(f'{scatter_plot_dir}/{name1}__v__{name2}')
    plt.close()

    return pearson_value


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


def compare_bedGraph_stats(bedGraph1_stat, bedGraph2_stat, min_value=0):
    pearson_value = get_coefficients(bedGraph1_stat, bedGraph2_stat, min_value)
    results = {'pearson_value': pearson_value}

    results['main_value'] = results['pearson_value']
    return results


def read_bedGraphs(datatype=None, hiseq=False, data_directory=None, min_value=-1):
    bedGraph_dict = {}
    if data_directory is None:
        try:
            assert datatype is not None
        except AssertionError:
            log.error("Missing data_directory and datatype")
            return bedGraph_dict

        if hiseq:
            data_directory = f'{DATA_DIR}/{datatype}/hiseq/'
        else:
            data_directory = f'{DATA_DIR}/{datatype}/miseq/'

    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph'):
            file_path = os.path.join(data_directory, filename)
            bedGraph = BedGraph(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path, 'chr1', ignore_missing_bp=False,
                                min_value=min_value)
            bedGraph_dict[bedGraph.name] = bedGraph

    return bedGraph_dict
