# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from logging.config import fileConfig
from copy import deepcopy
import time

from create_graphs import *

sys.path.append('..')
from chia_rep import all_loop_data
from chia_rep import chrom_loop_data
from chia_rep import reproducibility
from chia_rep import loop_rep

loop_dict = reproducibility.read_data(HUMAN_DATA_DIR,
                                      f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      BEDGRAPH_DATA_DIR, chrom_to_load='chr1')
loop_dict.update(reproducibility.read_data(MOUSE_DATA_DIR,
                                           f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                                           BEDGRAPH_DATA_DIR,
                                           chrom_to_load='chr1'))
main()
