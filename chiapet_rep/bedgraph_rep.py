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
DEFAULT_CHROM = 'chr1'

VERSION = 1

log = logging.getLogger()


def create_test_cases(end, interval_size=TEST_CASE_INTERVAL_SIZE, start=0,
                      step=TEST_CASE_INTERVAL_SIZE):
    log.info(f"Made test cases from {start} - {end} with interval sizes: "
             f"{interval_size} and step size: {step}")
    num_test_cases = int((end - start) / step) - interval_size
    test_cases = np.arange(start, start + step * num_test_cases, step, dtype=np.int32)
    test_cases = np.vstack((test_cases, test_cases + interval_size))

    return test_cases


def get_stats(bedgraph_dict, test_cases, stat='mean'):

    bedgraph_stats_dict = {}
    for name in bedgraph_dict:
        log.info(f'Getting {stat} for {name}')
        bedgraph = bedgraph_dict[name]
        bedgraph.load_chrom_data(DEFAULT_CHROM)
        chrom = bedgraph.chromosome_map[DEFAULT_CHROM]

        bedgraph_stats_dict[name] = {
            'num_samples': chrom.num_samples,
            'stat_list': bedgraph.stats(start_list=test_cases[0],
                                        end_list=test_cases[1],
                                        chrom_name=DEFAULT_CHROM,
                                        stat=stat),
            'name': bedgraph.name,
            'max_value': np.max(chrom.value_map)
        }
        bedgraph.free_chrom_data(DEFAULT_CHROM)

    return bedgraph_stats_dict


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


# min_value > 0 cuts a square in the bottom left corner of scatter plot
# Rationale: It's too dense there
def get_coefficients(bedgraph1_stat, bedgraph2_stat, min_value):

    stat_list1 = bedgraph1_stat['stat_list']
    stat_list2 = bedgraph2_stat['stat_list']
    assert stat_list1.size == stat_list2.size

    # Cut a rectangle instead of a square?
    num_samples1 = bedgraph1_stat['num_samples']
    num_samples2 = bedgraph2_stat['num_samples']
    ratio = num_samples1 / num_samples2
    modified_min_value = min_value * ratio

    to_remove = []
    for i in range(stat_list1.size):
        if stat_list2[i] < min_value and stat_list1[i] < modified_min_value:
            to_remove.append(i)
    stat_list1 = np.delete(stat_list1, to_remove)
    stat_list2 = np.delete(stat_list2, to_remove)

    # scipy pearson
    orig_pearson_value, pearson_p = pearsonr(stat_list1, stat_list2)

    # An attempt at a weighted value
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

    name1 = bedgraph1_stat['name']
    name2 = bedgraph2_stat['name']
    plt.scatter(stat_list1, stat_list2)
    plt.title(f'{name1} vs. {name2}, r={round(pearson_value, 3)}')
    plt.savefig(f'{scatter_plot_dir}/{name1}__v__{name2}')
    plt.close()

    return pearson_value


def compare(bedgraph1_stat, bedgraph2_stat, min_value=0):
    pearson_value = get_coefficients(bedgraph1_stat, bedgraph2_stat, min_value)
    results = {
        'graph_type': 'pearson',
        'rep': pearson_value,
        'w_rep': pearson_value}

    return results

