# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
import sys

sys.path.append('..')
sys.path.append('../../pyBedGraph')
import combined
import compare_bedGraph
import copy
import pyBedGraph
from pyBedGraph import BedGraph, Chrom_Data, include_missing_bp
import logging as log

log_path = 'logs/pearson'
parameters = f'pearson'
log.basicConfig(
    level=log.INFO,
    handlers=[
        log.FileHandler("{0}/{1}".format(log_path, parameters), mode='a+'),
        log.StreamHandler()
    ]
)

#b = combined.read_bedGraph_data(True)
#b.update(combined.read_bedGraph_data())
b = combined.read_loop_bedgraph()
complete_test_cases = compare_bedGraph.create_test_cases(1000)
n_overlap_stats = compare_bedGraph.get_stats(b, complete_test_cases)
for x in [-1]:
    scores = combined.compare(n_overlap_stats,
                              compare_bedGraph.compare_bedGraph_stats,
                              min_value=x)
    combined.output_to_csv(scores, f'pearson_{x}.csv')
