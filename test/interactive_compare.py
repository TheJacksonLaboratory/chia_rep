# coding: utf-8
g1 = copy.deepcopy(g)
for name in b:
    g1[name].find_loop_anchor_points(b[name])
    g1[name].filter_with_bedGraph(test_cases, 'chr1', stats[name]['mean_list'])
scores = combined.compare(g1, compare_loops.compare_loops,
                          bin_size=combined.BIN_SIZE,
                          window_size=combined.WINDOW_SIZE,
                          window_index=int(combined.CHROM_START / combined.WINDOW_SIZE))
