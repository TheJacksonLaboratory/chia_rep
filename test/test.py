# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
import logging
from logging.config import fileConfig
sys.path.append('..')
from Reproducibility import reproducibility
from Reproducibility import loop_rep
from Reproducibility import all_loop_data
from Reproducibility import chrom_loop_data

DATA_PATH = '/media/hirwo/extra/jax/data/chia_pet'

fileConfig('test.conf')

'''bedgraph_dict = reproducibility.read_bedGraphs(f'{DATA_PATH}/bedgraphs',
                                               'data/hg38.chrom.sizes',
                                               chrom_to_load=None)

loop_dict = reproducibility.read_data(loop_data_dir=f'{DATA_PATH}/loops',
                                      peak_data_dir=f'{DATA_PATH}/peaks',
                                      chrom_size_file='data/hg38.chrom.sizes',
                                      bedgraph_dict=bedgraph_dict)'''

bedgraph_dict = reproducibility.read_bedGraphs(f'data/small',
                                               'data/hg38.chrom.sizes',
                                               chrom_to_load=None)

loop_dict = reproducibility.read_data(loop_data_dir=f'data/small',
                                      peak_data_dir=f'data/small',
                                      chrom_size_file='data/hg38.chrom.sizes',
                                      bedgraph_dict=bedgraph_dict,
                                      is_hiseq=False)

rep, non_rep, scores = reproducibility.compare(loop_dict, loop_rep.compare,
                                               bin_size=10000,
                                               window_size=3000000)

reproducibility.output_results(rep, non_rep)
