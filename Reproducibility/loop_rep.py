import os
from collections import OrderedDict
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import logging

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'
CHROM_SIZES_FILE = f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes'

chr1_size = 248956422
chr2_size = 242193529
chr_size = chr1_size
CHROM_NAME = 'chr1'

VERSION = 5

log = logging.getLogger(__name__.split('.')[-1])


def compare(loop_1, loop_2, bin_size, window_size, window_index=None):
    total_values = []
    for k in range(int(chr_size / window_size)):

        if window_index is not None:
            k = window_index

        window_start = window_size * k
        window_end = window_size * (k + 1)
        if window_end > chr1_size:
            window_end = chr1_size

        avg_result_list = loop_1.compare(loop_2, window_start, window_end,
                                         bin_size)
        for i in range(len(avg_result_list)):
            if i == len(total_values):
                total_values.append([])
            total_values[i].append(avg_result_list[i])

        if window_index is not None:
            break

    avg_result_list = []
    for result_type_list in total_values:
        if not result_type_list:
            continue
        avg_result = OrderedDict()
        avg_result['graph_type'] = result_type_list[0]['graph_type']
        weights = [x['w'] for x in result_type_list]
        for key in result_type_list[0]:
            if key == 'w' or key == 'graph_type':
                continue
            stat_list = [x[key] for x in result_type_list]
            avg_result['uw_' + key] = np.mean(stat_list)
            avg_result['w_' + key] = np.average(stat_list, weights=weights)
        avg_result['main_value'] = avg_result['w_j-s_value']
        avg_result_list.append(avg_result)

    # avg_result_list.sort(key=lambda x: x['w_j-s_value'], reverse=True)

    final_results = PrettyTable(list(avg_result_list[0].keys()))
    for avg_result in avg_result_list:
        final_results.add_row(list(avg_result.values()))

    log.info(f"{loop_1.sample_name}_{loop_2.sample_name}" +
             final_results.get_string())

    '''plt.close()
    plt.hist(j_values, bins=20)
    plt.xlabel('Jensen-Shannon window values')
    plt.ylabel('Frequency')
    plt.xlim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - Reproducibility: '
              f'{round(j_value, 2)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_hist.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(j_values))],
             j_values)
    plt.xlabel('Jensen-Shannon Window')
    plt.ylabel('j-value')
    plt.ylim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - Reproducibility: '
              f'{round(j_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_j.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(emd_values))],
             emd_values)
    plt.xlabel('Window')
    plt.ylabel('Earth Mover\'s Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - Reproducibility: '
              f'{round(emd_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_emd.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(w_values))],
             w_values)
    plt.xlabel('Window')
    plt.ylabel('Wasserstein Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - Reproducibility: '
              f'{round(w_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_w.png')
    plt.close()'''

    return avg_result_list
