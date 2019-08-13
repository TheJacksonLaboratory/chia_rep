import logging as log
import copy
import combined

for bin_size in [10000]:
    combined.BIN_SIZE = bin_size
    log_path = 'logs'
    parameters = f'{combined.BIN_SIZE}_{chrom_loop_data.PERCENT_PEAK_KEPT}_' \
                 f'{combined.WINDOW_SIZE}'
    log.basicConfig(
        level=log.INFO,
        handlers=[
            log.FileHandler("{0}/{1}.log".format(log_path, parameters),
                            mode='w+'),
            log.StreamHandler()
        ]
    )

    # test_cases = compare_bedGraph.create_test_cases(combined.BIN_SIZE, combined.CHROM_START, combined.CHROM_END, 10)
    # stats = compare_bedGraph.get_stats(b, test_cases)
    log.info(f'BIN_SIZE: {combined.BIN_SIZE}')
    log.info(f'WINDOW_SIZE: {combined.WINDOW_SIZE}')
    log.info(f'PEAK PERCENTAGE KEPT: {chrom_loop_data.PERCENT_PEAK_KEPT}')
    log.info(f"combined: {combined.VERSION}")
    log.info(f"compare_loops: {compare_loops.VERSION}")
    log.info(f"ChromLoop: {chrom_loop_data.VERSION}")
    log.info(f"GenomeLoop: {genome_loop_data.VERSION}")
    g1 = copy.deepcopy(g)
    for name in b:
        g1[name].find_loop_anchor_points(b[name])
        g1[name].filter_with_bedGraph('chr1')
    rep, non_rep, scores = combined.compare(g1, compare_loops.compare_loops,
                                            bin_size=combined.BIN_SIZE,
                                            window_size=combined.WINDOW_SIZE,
                                            window_index=None)
                                            # window_index=int(combined.CHROM_START / combined.WINDOW_SIZE))
    combined.output_results(rep, non_rep)
    combined.output_to_csv(scores, f'csv_results/{parameters}.csv')
    log.info('--------------------------------------------------------------------'
             '--------------------------------------------------------------------'
             '--------------------------------------------------------------------')
