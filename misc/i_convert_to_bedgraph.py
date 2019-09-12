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
import convert_loop_to_bedGraph
import copy
import pyBedGraph
from pyBedGraph import BedGraph, Chrom_Data, include_missing_bp
import logging as log

log_path = 'logs'
parameters = f'convert'
file_handler = log.StreamHandler()
file_handler.terminator = ''
log.basicConfig(
    level=log.INFO,
    handlers=[
        log.FileHandler("{0}/{1}".format(log_path, parameters), mode='a+'),
        file_handler
    ]
)

g = combined.read_loop_data(True)
b = combined.read_bedGraph_data(True)
g1 = copy.deepcopy(g)
for name in b:
    g1[name].find_loop_anchor_points(b[name])
    g1[name].filter_with_peaks('chr1')
    arr = convert_loop_to_bedGraph.get_bedgraph_array(g1[name].chrom_dict['chr1'])
    intervals = convert_loop_to_bedGraph.get_intervals(arr)
    convert_loop_to_bedGraph.output_as_bedgraph(intervals,
                                                f'loop_bedgraphs/{name}.bedgraph')


