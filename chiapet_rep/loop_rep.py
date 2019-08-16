import os
from collections import OrderedDict
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import logging

DATA_DIR = '/media/hirwo/extra/jax/data/chia_pet'
CHROM_SIZES_FILE = f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes'

chr1_size = 248956422
chr2_size = 242193529
chr_size = chr1_size
VERSION = 5

log = logging.getLogger()


def compare(loop_1, loop_2, bin_size, window_size, window_index=None,
            wanted_chroms=None):

    result = loop_1.compare(loop_2, bin_size, window_size, window_index,
                            wanted_chroms)

    '''plt.close()
    plt.hist(j_values, bins=20)
    plt.xlabel('Jensen-Shannon window values')
    plt.ylabel('Frequency')
    plt.xlim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chiapet_rep: '
              f'{round(j_value, 2)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_hist.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(j_values))],
             j_values)
    plt.xlabel('Jensen-Shannon Window')
    plt.ylabel('j-value')
    plt.ylim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chiapet_rep: '
              f'{round(j_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_j.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(emd_values))],
             emd_values)
    plt.xlabel('Window')
    plt.ylabel('Earth Mover\'s Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chiapet_rep: '
              f'{round(emd_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_emd.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(w_values))],
             w_values)
    plt.xlabel('Window')
    plt.ylabel('Wasserstein Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chiapet_rep: '
              f'{round(w_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_w.png')
    plt.close()'''

    return result
