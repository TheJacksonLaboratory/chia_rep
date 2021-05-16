# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

import sys
from copy import deepcopy

sys.path.append('..')
from chia_rep import chia_rep
from misc import bedgraph_rep

MOUSE_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/mouse'
PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/peaks'
MY_PEAK_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/my_peaks'
HUMAN_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/human'
BEDGRAPH_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/bedgraphs'
PEAK_BEDGRAPH_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/bedgraph_from_peak'
CHROM_DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet/chrom_sizes'
DATA_DIR = 'data'
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 3000000  # 3mb

# To see log info statements (optional)
from logging.config import fileConfig

fileConfig('chia_rep.conf')

# Since reading in bedgraph file can take a long time, load them first if in an interactive session
bedgraph_dict = chia_rep.read_bedgraphs(PEAK_BEDGRAPH_DATA_DIR + '/human',
                                               f'{DATA_DIR}/hg38.chrom.sizes',
                                        chrom_to_read='chr1')
bedgraph_dict.update(
    chia_rep.read_bedgraphs(PEAK_BEDGRAPH_DATA_DIR + '/mouse',
                                   f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                            chrom_to_read='chr1'))

'''sample_dict = reproducibility.read_data(loop_data_dir=DATA_DIR,
                                      chrom_size_file=f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      is_hiseq=True,
                                      bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                      chrom_to_read='chr1')

sample_dict.update(reproducibility.read_data(loop_data_dir=MOUSE_DATA_DIR,
                                           chrom_size_file=f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                                           is_hiseq=True,
                                           bedgraph_data_dir=BEDGRAPH_DATA_DIR,
                                           chrom_to_read='chr1'))'''

'''sample_dict = reproducibility.read_data(loop_data_dir=PARENT_DIR,
                                      chrom_size_file=f'{PARENT_DIR}/hg38.chrom.sizes',
                                      is_hiseq=False,  # Determines if match comparison is made (TODO)
                                      bedgraph_data_dir=PARENT_DIR)'''

test_cases = bedgraph_rep.create_test_cases(195471971)
bedgraph_stats = bedgraph_rep.get_stats(bedgraph_dict, test_cases, 'max')

for i in [1]:
    b_stats = deepcopy(bedgraph_stats)
    rep, non_rep, scores = chia_rep.compare(bedgraph_stats,
                                            bedgraph_rep.compare,
                                            min_value=i)
    temp_str = f'pearson_r_{i}'
    chia_rep.check_results(rep, non_rep, f'results/info/{temp_str}.txt')
    chia_rep.output_to_csv(scores, f'results/{temp_str}.csv')
