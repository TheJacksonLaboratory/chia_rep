# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from copy import deepcopy
import os

sys.path.append('..')
from chia_rep import reproducibility
from chia_rep import GenomeLoopData
from chia_rep import ChromLoopData

# MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/prob_mouse'
MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/mouse'
PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/peaks'
MY_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_peaks'
MY_FULL_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_full_peaks'
# HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/prob_human'
HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/human'
BEDGRAPH_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/bedgraphs'
CHROM_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/chrom_sizes'

LOG_MAIN_FORMAT = '%(levelname)s - %(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(message)s'
LOG_BIN_FORMAT = '%(asctime)s - %(filename)s:%(lineno)d\n%(message)s'
LOG_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

# To see log info statements (optional)
from logging.config import fileConfig

fileConfig('chia_rep.conf')
main_formatter = logging.Formatter(LOG_MAIN_FORMAT, datefmt=LOG_TIME_FORMAT)

loop_dict = reproducibility.read_data(loop_data_dir=HUMAN_DATA_DIR,
                                      chrom_size_file=f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                      peak_data_dir=PEAK_DATA_DIR)

loop_dict.update(reproducibility.read_data(loop_data_dir=MOUSE_DATA_DIR,
                                           chrom_size_file=f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                                           bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                           peak_data_dir=PEAK_DATA_DIR))

# for name in loop_dict:
#     for chrom in loop_dict[name].peak_dict:
#         peaks = []
#         chrom_peak_dict = loop_dict[name].peak_dict[chrom]
#         for x in range(chrom_peak_dict['num']):
#             peaks.append([chrom_peak_dict['start_list'][x],
#                           chrom_peak_dict['end_list'][x],
#                           chrom_peak_dict['peak_len'][x],
#                           chrom_peak_dict['value_list'][x]])
#         loop_dict[name].peak_dict[chrom] = peaks


def comparison():
    parent_dir = 'all_complete'
    # parent_dir = 'small_complete'
    for i in [10, 20, 30, 40, 50, 100, 500]:  # bin size
    # for i in [1]:
        i *= 1000
        for j in [10, 20]:  # window size:
        # for j in [3]:
            j *= 1000000
            for k in [20, 30, 40, 50, 60]:  # Peaks kept
            # for k in [80, 100, 200]:  # Peaks kept
                temp_str = f'half.{k}peaks.{i}.{j}'

                # if os.path.isfile(f'{parent_dir}/results/{temp_str}.emd_value.csv'):
                #     print(f'Skipping {temp_str}')
                #     continue

                main_handler = logging.FileHandler(f'{parent_dir}/tmp/main.{temp_str}.log', mode='w')
                main_handler.setFormatter(main_formatter)
                main_handler.setLevel(logging.DEBUG)

                stream_handler = logging.StreamHandler(sys.stdout)
                stream_handler.setFormatter(main_formatter)
                stream_handler.setLevel(logging.INFO)

                log = logging.getLogger()
                log_all = logging.getLogger('all')

                log.handlers = [main_handler, stream_handler]
                log_all.handlers = [main_handler]

                l = deepcopy(loop_dict)
                reproducibility.preprocess(l, num_peaks=k, kept_dir=f'kept')

                rep, non_rep, emd_scores, j_scores = \
                    reproducibility.compare(l, bin_size=i, window_size=j)
                # reproducibility.output_results(rep, non_rep)
                reproducibility.output_results(rep, non_rep,
                                               f'{parent_dir}/results/info',
                                               temp_str)
                reproducibility.output_to_csv(emd_scores,
                                              f'{parent_dir}/results/{temp_str}.emd_value.csv')
                reproducibility.output_to_csv(j_scores,
                                              f'{parent_dir}/results/{temp_str}.j_value.csv')

comparison()
