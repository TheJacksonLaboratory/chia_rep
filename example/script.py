"""
Use this file for testing purposes in an interactive python session. This is to
avoid the time-consuming process of reloading bedgraphs each time a comparison
is done.
"""

# coding: utf-8
# Automatically updates modules that have been imported
# Still needs to reload/reset if extra methods/fields are added to class objects
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from copy import deepcopy
import os

# Include the following line if chia_rep is not installed
sys.path.append('..')
import chia_rep

PARENT_DIR = '/media/hyrial/KESU/jax/data'
DATA_DIR = f'{PARENT_DIR}/to_use'
CHROM_SIZE_PATH = f'{PARENT_DIR}/chrom_sizes/hg38.chrom.sizes'

LOG_MAIN_FORMAT = '%(levelname)s - %(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(message)s'
LOG_BIN_FORMAT = '%(asctime)s - %(filename)s:%(lineno)d\n%(message)s'
LOG_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

# To see log info statements (optional)
from logging.config import fileConfig

fileConfig('log.conf')
main_formatter = logging.Formatter(LOG_MAIN_FORMAT, datefmt=LOG_TIME_FORMAT)

if len(sys.argv) == 2:
    DATA_DIR = sys.argv[1]
    CHROM_SIZE_PATH = f'{DATA_DIR}/hg38.chrom.sizes'

# Make the loop_dict a global variable to be untouched so it doesn't have to be reloaded
loop_dict = chia_rep.read_data(loop_data_dir=DATA_DIR,
                               chrom_size_file=CHROM_SIZE_PATH,
                               bedgraph_data_dir=DATA_DIR,
                               peak_data_dir=DATA_DIR,
                               chroms_to_load=['chr1'])


# Simply copy/paste following method into interactive session and run it if changed
def comparison():
    parent_dir = 'test_results'
    # parent_dir = 'small_complete'

    directories = [
        parent_dir,
        f'{parent_dir}/results',
        f'{parent_dir}/results/info',
        f'{parent_dir}/tmp'
    ]

    for directory in directories:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    for bin_size in [5, 10]:  # bin size kb
        # for i in [1]:
        bin_size *= 1000
        for window_size in [3, 5, 10]:  # window size mb
            # for j in [3]:
            window_size *= 1000000
            temp_str = f'half.{bin_size}_bin.{window_size}_window'

            if os.path.isfile(f'{parent_dir}/results/{temp_str}.emd_value.csv'):
                print(f'Skipping {temp_str}')
                continue

            main_handler = logging.FileHandler(
                f'{parent_dir}/tmp/main.{temp_str}.log', mode='w')
            main_handler.setFormatter(main_formatter)
            main_handler.setLevel(logging.DEBUG)

            stream_handler = logging.StreamHandler(sys.stdout)
            stream_handler.setFormatter(main_formatter)
            stream_handler.setLevel(logging.INFO)

            log = logging.getLogger()
            log_all = logging.getLogger('all')

            log.handlers = [main_handler, stream_handler]
            log_all.handlers = [main_handler]

            # Make sure not to modify the original dict
            l = deepcopy(loop_dict)
            chia_rep.preprocess(l, saved_loops_dir=f'kept')

            rep, non_rep, emd_scores, j_scores = \
                chia_rep.compare(l, bin_size=bin_size, window_size=window_size)
            # chia_rep.output_results(rep, non_rep)
            chia_rep.output_results(rep, non_rep, f'{parent_dir}/results/info',
                                    temp_str)
            chia_rep.output_to_csv(emd_scores,
                                   f'{parent_dir}/results/{temp_str}.emd_value.csv')
            chia_rep.output_to_csv(j_scores,
                                   f'{parent_dir}/results/{temp_str}.j_value.csv')
