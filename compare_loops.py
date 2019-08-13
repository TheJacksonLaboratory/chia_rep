import os
from collections import OrderedDict
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import logging as log

sys.path.append('..')
from genome_loop_data import GenomeLoopData, WINDOW_SIZE

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'
CHROM_SIZES_FILE = f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes'

chr1_size = 248956422
chr2_size = 242193529
chr_size = chr1_size
CHROM_NAME = 'chr1'

VERSION = 4


def compare_loops(loop_1, loop_2, bin_size, window_size, window_index=None):
    total_values = [{}, {}, {}]
    for k in range(int(chr_size / window_size)):

        if window_index is not None:
            k = window_index

        window_start = window_size * k
        window_end = window_size * (k + 1)
        if window_end > chr1_size:
            window_end = chr1_size

        values_dict_list = loop_1.compare(loop_2, window_start, window_end,
                                          bin_size)
        for i in range(len(values_dict_list)):
            values_dict = values_dict_list[i]
            for key in values_dict:
                if key not in total_values[i]:
                    total_values[i][key] = []
                total_values[i][key].append(values_dict[key])

        if window_index is not None:
            break

    results_list = []
    for value_dict in total_values:
        if not value_dict:
            break
        results = OrderedDict()
        results['type'] = value_dict['type'][0]
        for key in value_dict:
            if key == 'w' or key == 'type':
                continue
            results['uw_' + key] = np.mean(value_dict[key])
            results['w_' + key] = np.average(value_dict[key],
                                             weights=value_dict['w'])
        results_list.append(results)

    results_list.sort(key=lambda x: x['w_j-s_value'], reverse=True)

    final_results = PrettyTable(list(results_list[0].keys()))
    for results in results_list:
        final_results.add_row(list(results.values()))

    table_str = f"{loop_1.name}_{loop_2.name}" + \
                final_results.get_string()
    log.info(table_str)

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

    best_result = results_list[0]
    best_result['main_value'] = best_result['w_j-s_value']
    return best_result


def read_loops(datatype, hiseq=False, min_loop_value=None):
    loop_info_dict = {}
    if hiseq:
        data_directory = DATA_DIR + f'/{datatype}/hiseq'
    else:
        data_directory = DATA_DIR + f'/{datatype}/miseq'

    for filename in os.listdir(data_directory):
        if filename.endswith('.BE3') or filename.endswith('.cis'):
            file_path = os.path.join(data_directory, filename)
            loopInfo = GenomeLoopData(CHROM_SIZES_FILE, file_path,
                                      wanted_chrom_name=CHROM_NAME, hiseq=hiseq,
                                      min_loop_value=min_loop_value)
            loop_info_dict[loopInfo.name] = loopInfo

    return loop_info_dict
