# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')

from create_graphs import *

sys.path.append('..')
from chia_rep import reproducibility

loop_dict = reproducibility.read_data(HUMAN_DATA_DIR,
                                      f'{CHROM_DATA_DIR}/hg38.chrom.sizes',
                                      BEDGRAPH_DATA_DIR, chrom_to_load='chr1')
loop_dict.update(reproducibility.read_data(MOUSE_DATA_DIR,
                                           f'{CHROM_DATA_DIR}/mm10.chrom.sizes',
                                           BEDGRAPH_DATA_DIR,
                                           chrom_to_load='chr1'))
main()
