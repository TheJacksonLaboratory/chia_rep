import os
import sys
sys.path.append('..')
import compare_bedGraph
import compare_loops
import matplotlib.pyplot as plt
import numpy as np


CHROM_START = 36000000
CHROM_END = 39000000
WINDOW_SIZE = CHROM_END - CHROM_START
INTERVAL_SIZE = 5000


def compare_with_diff_min_value(bedGraph_dict):
    for min_value in [x for x in range(2, 20, 2)]:
        genome_loop_data_dict = read_loops('CTCF', INTERVAL_SIZE, True, min_value)
        genome_loop_data_dict.update(read_loops('RNAPII', INTERVAL_SIZE, True, min_value))

        for name in bedGraph_dict:
            genome_loop_data_dict[name].adjust_with_bedGraph(bedGraph_dict[name])

        results = compare_loops(genome_loop_data_dict)
        compare_loops.output_csv(results, f'scores/avg/min_value_{min_value}.csv')


def read_loop_data(hiseq=False):
    genome_loop_dict = compare_loops.read_loops('CTCF', INTERVAL_SIZE, hiseq=hiseq)
    genome_loop_dict.update(compare_loops.read_loops('RNAPII', INTERVAL_SIZE, hiseq=hiseq))

    return genome_loop_dict


def read_bedGraph_data(hiseq=False):
    bedGraph_dict = compare_bedGraph.read_bedGraphs('CTCF', hiseq=hiseq)
    bedGraph_dict.update(compare_bedGraph.read_bedGraphs('RNAPII', hiseq=hiseq))

    return bedGraph_dict


def main():
    genome_loop_dict = read_loop_data()
    bedGraph_dict = read_bedGraph_data()

    test_cases = compare_bedGraph.create_test_cases(interval_size=INTERVAL_SIZE,
                                                    start=CHROM_START,
                                                    end=CHROM_END,
                                                    step=1)
    bedGraph_stats_dict = compare_bedGraph.get_stats(bedGraph_dict, test_cases)

    for name in bedGraph_dict:
        genome_loop_dict[name].filter_with_bedGraph(test_cases, 'chr1',
                                                    bedGraph_stats_dict[name]['mean_list'])

    results = compare_loops.compare_loops(genome_loop_dict, INTERVAL_SIZE,
                                          WINDOW_SIZE, int(CHROM_START / WINDOW_SIZE))
    compare_loops.output_csv(results, 'results.csv')


if __name__ == '__main__':
    main()
