# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from logging.config import fileConfig
sys.path.append('..')
from chiapet_rep import all_loop_data
from chiapet_rep import chrom_loop_data
from chiapet_rep import reproducibility
from chiapet_rep import loop_rep

TEST_DATA_DIR = 'data'
DATA_PATH = '/media/hirwo/extra/jax/data/chia_pet'
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 3000000  # 3mb

# To see log info statements (optional)
from logging.config import fileConfig
fileConfig('test.conf')


# Since reading in bedgraph file can take a long time, load them first if in an interactive session
bedgraph_dict = reproducibility.read_bedGraphs(TEST_DATA_DIR, f'{TEST_DATA_DIR}/hg38.chrom.sizes')

loop_dict = reproducibility.read_data(loop_data_dir=TEST_DATA_DIR,
                                      peak_data_dir=TEST_DATA_DIR,
                                      chrom_size_file=f'{TEST_DATA_DIR}/hg38.chrom.sizes',
                                      is_hiseq=False,  # Determines if match comparison is made (TODO)
                                      bedgraph_dict=bedgraph_dict)

rep, non_rep, scores = reproducibility.compare(loop_dict, loop_rep.compare,
                                               bin_size=BIN_SIZE,
                                               window_size=WINDOW_SIZE)

reproducibility.output_results(rep, non_rep)
