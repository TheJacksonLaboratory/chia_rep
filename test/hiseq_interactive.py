# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
import sys
sys.path.append('..')
sys.path.append('../../pyBedGraph')
import combined
import compare_loops
import compare_bedGraph
import genome_loop_data
import chrom_loop_data
import copy
import pyBedGraph
from pyBedGraph import BedGraph, Chrom_Data, include_missing_bp
g = combined.read_loop_data(True)
b = combined.read_bedGraph_data(True)
test_cases = compare_bedGraph.create_test_cases(combined.BIN_SIZE, combined.CHROM_START, combined.CHROM_END, 1)
stats = compare_bedGraph.get_stats(b, test_cases)
g1 = copy.deepcopy(g)
for name in b:
    g1[name].find_loop_anchor_points(b[name])
    g1[name].filter_with_bedGraph(test_cases, 'chr1', stats[name]['mean_list'])
scores = combined.compare(g1, compare_loops.compare_loops,
                          bin_size=combined.BIN_SIZE,
                          window_size=combined.WINDOW_SIZE,
                          window_index=int(combined.CHROM_START / combined.WINDOW_SIZE))
combined.output_to_csv(scores, 'csv_results/results.csv')

