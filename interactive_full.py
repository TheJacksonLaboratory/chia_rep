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
import logging as log

g = combined.read_loop_data(True)
g.update(combined.read_loop_data())
b = combined.read_bedGraph_data(True)
b.update(combined.read_bedGraph_data())

