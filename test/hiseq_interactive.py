# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
import combined
import compare_loops
import compare_bedGraph
import copy
g = combined.read_loop_data(True)
b = combined.read_bedGraph_data(True)
test_cases = compare_bedGraph.create_test_cases(combined.INTERVAL_SIZE, combined.CHROM_START, combined.CHROM_END, 1)
stats = compare_bedGraph.get_stats(b, test_cases)
g1 = copy.deepcopy(g)
for name in b:
    g1[name].filter_with_bedGraph(test_cases, 'chr1', stats[name]['mean_list'])
compare_loops.compare_loops(g1, combined.INTERVAL_SIZE, combined.WINDOW_SIZE, int(combined.CHROM_START / combined.WINDOW_SIZE))

