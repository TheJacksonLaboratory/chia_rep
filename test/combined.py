import os
import sys
sys.path.append('..')
import remove_noise
from bedGraph import *
from loops import *
import matplotlib.pyplot as plt
import numpy as np


DATA_DIR = '../data'
BIN_SIZE = 10000


def read_data(datatype):
    genome_loop_data_dict = read_loops(DATA_DIR + f'/{datatype}/small', BIN_SIZE)
    bedGraph_dict = read_small_bedGraphs(datatype)
    return genome_loop_data_dict, bedGraph_dict


def main():
    for datatype in ['CTCF', 'RNAPII']:
        genome_loop_data_dict, bedGraph_dict = read_data(datatype)
        remove_noise.remove_noise(list(bedGraph_dict.values()))

        for name in bedGraph_dict:
            genome_loop_data_dict[name].adjust_with_bedGraph(bedGraph_dict[name])
            pass

        compare_loops(genome_loop_data_dict)
        break


if __name__ == '__main__':
    main()
