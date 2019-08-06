
b = combined.read_bedGraph_data(True)
complete_test_cases = compare_bedGraph.create_test_cases(1000)
n_overlap_stats = compare_bedGraph.get_stats(b, complete_test_cases)
for x in [0, 100, 200]:
    scores = combined.compare(n_overlap_stats, compare_bedGraph.compare_bedGraph_stats, min_value=x)
    combined.output_to_csv(scores, f'pearson_{x}.csv')