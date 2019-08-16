from scipy.stats import pearsonr
import numpy as np
import os
import matplotlib.pyplot as plt
from math import sqrt
from copy import deepcopy
import logging
from prettytable import PrettyTable
from pyBedGraph import BedGraph

TEST_CASE_INTERVAL_SIZE = 1000

VERSION = 1

log = logging.getLogger()


def create_test_cases(end, interval_size=TEST_CASE_INTERVAL_SIZE, start=0,
                      step=TEST_CASE_INTERVAL_SIZE):
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


def compare(bedGraph1_stat, bedGraph2_stat, min_value=0):
    pearson_value = get_coefficients(bedGraph1_stat, bedGraph2_stat, min_value)
    results = {'pearson_value': pearson_value}

    results['main_value'] = results['pearson_value']
    return results

