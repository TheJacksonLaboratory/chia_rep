# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from logging.config import fileConfig
from copy import deepcopy
import time

sys.path.append('..')
from chia_rep import all_loop_data
from chia_rep import chrom_loop_data
from chia_rep import reproducibility
from chia_rep import loop_rep

MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/mouse'
PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/peaks'
MY_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_peaks'
HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/human'
BEDGRAPH_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/bedgraphs'
CHROM_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/chrom_sizes'
DATA_DIR = 'data'
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 3000000  # 3mb

LOG_MAIN_FORMAT = '%(levelname)s - %(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(message)s'
LOG_BIN_FORMAT = '%(asctime)s - %(filename)s:%(lineno)d\n%(message)s'
LOG_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

# To see log info statements (optional)
from logging.config import fileConfig

fileConfig('chia_rep.conf')
main_formatter = logging.Formatter(LOG_MAIN_FORMAT, datefmt=LOG_TIME_FORMAT)
bin_formatter = logging.Formatter(LOG_BIN_FORMAT, datefmt=LOG_TIME_FORMAT)

main_handler = logging.FileHandler(f'tmp/input.log', mode='w')
main_handler.setFormatter(main_formatter)
main_handler.setLevel(logging.DEBUG)

log = logging.getLogger()
log.handlers = []

# Since reading in bedgraph file can take a long time, load them first if in an interactive session
# bedgraph_dict = reproducibility.read_bedgraphs(DATA_DIR, f'{DATA_DIR}/hg38.chrom.sizes')

loop_dict = reproducibility.read_data(loop_data_dir=HUMAN_DATA_DIR,
                                      chrom_size_file=f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      is_hiseq=True,
                                      bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                      chrom_to_read='chr1')

loop_dict.update(reproducibility.read_data(loop_data_dir=MOUSE_DATA_DIR,
                                           chrom_size_file=f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                                           is_hiseq=True,
                                           bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                           chrom_to_read='chr1'))

'''loop_dict = reproducibility.read_data(loop_data_dir=DATA_DIR,
                                      chrom_size_file=f'{DATA_DIR}/hg38.chrom.sizes',
                                      is_hiseq=False,  # Determines if match comparison is made (TODO)
                                      bedgraph_data_dir=DATA_DIR)'''

for i in [1, 5, 10, 50]:
    i *= 1000  # bin size
    for j in [1, 3, 5, 10]:
        j *= 1000000  # window size
        for k in [0.005, 0.01, 0.025, 0.05, 0.1]:

            temp_str = f'both_peak_weight.{k * 100}.{i}.{j}'
            main_handler = logging.FileHandler(f'tmp/main.{temp_str}.log', mode='w')
            main_handler.setFormatter(main_formatter)
            main_handler.setLevel(logging.DEBUG)

            bin_handler = logging.FileHandler(f'tmp/bin.{temp_str}.log', mode='w')
            bin_handler.setFormatter(bin_formatter)
            bin_handler.setLevel(logging.DEBUG)

            stream_handler = logging.StreamHandler(sys.stdout)
            stream_handler.setFormatter(main_formatter)
            stream_handler.setLevel(logging.INFO)

            log = logging.getLogger()
            log_bin = logging.getLogger('bin')
            log_all = logging.getLogger('all')

            log.handlers = [main_handler, stream_handler]
            log_all.handlers = [main_handler, bin_handler]
            log_bin.handlers = [bin_handler]

            l = deepcopy(loop_dict)
            reproducibility.preprocess(l, MY_PEAK_DATA_DIR,
                                       peak_percentage_kept=k,
                                       window_size=None)

            rep, non_rep, scores = reproducibility.compare(l, loop_rep.compare,
                                                           bin_size=i,
                                                           window_size=j)
            # reproducibility.output_results(rep, non_rep)
            reproducibility.output_results(rep, non_rep, f'results/info/{temp_str}.txt')
            reproducibility.output_to_csv(scores, f'results/{temp_str}.csv')
