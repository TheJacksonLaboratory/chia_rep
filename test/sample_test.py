import sys
import os
sys.path.append('..')

from chia_rep import reproducibility

TEST_DATA_DIR = 'data'
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 10000000  # 10mb
NUM_PEAKS = 20

# To see log info statements (optional)
from logging.config import fileConfig
fileConfig('chia_rep.conf')

if not os.path.isdir('sample_test'):
    os.mkdir('sample_test')

loop_dict = reproducibility.read_data(loop_data_dir=TEST_DATA_DIR,
                                      chrom_size_file=f'{TEST_DATA_DIR}/hg38.chrom.sizes',
                                      bedgraph_data_dir=TEST_DATA_DIR,
                                      peak_data_dir=TEST_DATA_DIR)

reproducibility.preprocess(loop_dict, num_peaks=NUM_PEAKS, kept_dir='sample_test')

rep, non_rep, emd_scores, j_scores = \
    reproducibility.compare(loop_dict, bin_size=BIN_SIZE,
                            window_size=WINDOW_SIZE)

reproducibility.output_results(rep, non_rep, 'sample_test')
reproducibility.output_to_csv(emd_scores, f'sample_test/emd_value.csv')
reproducibility.output_to_csv(j_scores, f'sample_test/j_value.csv')

