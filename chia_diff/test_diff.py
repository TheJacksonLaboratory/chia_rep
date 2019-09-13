# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from copy import deepcopy

sys.path.append('..')
from chia_rep import reproducibility
from chia_rep import GenomeLoopData
from chia_rep import ChromLoopData

MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/prob_mouse'
MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/mouse'
PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/peaks'
MY_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_peaks'
MY_FULL_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_full_peaks'
HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/prob_human'
HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/human'
BEDGRAPH_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/bedgraphs'
CHROM_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/chrom_sizes'
CHIA_DIFF_DIR = '/media/hirwo/extra/jax/data/chia_pet/chia_diff'

LOG_MAIN_FORMAT = '%(levelname)s - %(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(message)s'
LOG_BIN_FORMAT = '%(asctime)s - %(filename)s:%(lineno)d\n%(message)s'
LOG_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

# To see log info statements (optional)
from logging.config import fileConfig

fileConfig('chia_rep.conf')
main_formatter = logging.Formatter(LOG_MAIN_FORMAT, datefmt=LOG_TIME_FORMAT)

loop_dict = reproducibility.read_data(loop_data_dir=CHIA_DIFF_DIR,
                                      chrom_size_file=f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                      peak_data_dir=PEAK_DATA_DIR)


def find_diff():
    l = deepcopy(loop_dict)
    keys = list(l.keys())
    for k in [60]:
        reproducibility.preprocess(l, num_peaks=k, both_peak_support=True)
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                sample1 = l[keys[i]]
                sample2 = l[keys[j]]
                sample1.find_diff_loops(sample2)
