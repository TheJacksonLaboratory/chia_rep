# coding: utf-8
g1 = copy.deepcopy(g)
for name in b:
    g1[name].filter_with_bedGraph(test_cases, 'chr1', stats[name]['mean_list'])
    g1[name].adjust_with_bedGraph(b[name])
compare_loops.compare_loops(g1, combined.INTERVAL_SIZE, combined.WINDOW_SIZE, int(combined.CHROM_START / combined.WINDOW_SIZE))

