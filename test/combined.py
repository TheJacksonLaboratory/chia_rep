import os
import sys

sys.path.append('..')
import compare_bedGraph
import compare_loops
import matplotlib.pyplot as plt
import numpy as np
import csv

CHROM_START = 36000000
CHROM_END = 39000000
WINDOW_SIZE = CHROM_END - CHROM_START
BIN_SIZE = 3000

REPLICATES = [
    ['84', '86'],
    ['83', '85'],
    ['48', '54'],
    ['58', '61']
]


def compare_with_diff_min_value(bedGraph_dict):
    for min_value in [x for x in range(2, 20, 2)]:
        genome_loop_data_dict = read_loops('CTCF', BIN_SIZE, True, min_value)
        genome_loop_data_dict.update(
            read_loops('RNAPII', BIN_SIZE, True, min_value))

        for name in bedGraph_dict:
            genome_loop_data_dict[name].adjust_with_bedGraph(
                bedGraph_dict[name])

        results = compare_loops(genome_loop_data_dict)
        compare_loops.output_csv(results,
                                 f'scores/avg/min_value_{min_value}.csv')


def read_loop_data(hiseq=False):
    genome_loop_dict = compare_loops.read_loops('CTCF', BIN_SIZE, hiseq=hiseq)
    genome_loop_dict.update(
        compare_loops.read_loops('RNAPII', BIN_SIZE, hiseq=hiseq))

    return genome_loop_dict


def read_bedGraph_data(hiseq=False):
    bedGraph_dict = compare_bedGraph.read_bedGraphs('CTCF', hiseq=hiseq)
    bedGraph_dict.update(compare_bedGraph.read_bedGraphs('RNAPII', hiseq=hiseq))

    return bedGraph_dict


def output_to_csv(scores, output_file):
    with open(output_file, 'w+') as out_file:
        header = list(scores.keys())
        header.insert(0, 'Sample Name')
        writer = csv.DictWriter(out_file, fieldnames=header)

        writer.writeheader()
        for sample_name, sample_scores in scores.items():
            writer.writerow(sample_scores)


def compare(data_dict, compare_func, **kwargs):
    print(f"Comparing {data_dict}\n"
          f"Function: {compare_func}\n"
          f"Arguments: {kwargs}")

    non_replicates = {}
    replicates = {}

    len_dict = len(data_dict)
    keys = list(data_dict.keys())
    scores = {}
    for key in keys:
        scores[key] = {}
        scores[key][key] = 1
        scores[key]['Sample Name'] = key

    for i in range(len_dict):
        for j in range(i + 1, len_dict):
            key1 = keys[i]
            key2 = keys[j]
            print(f'{key1} vs. {key2}:')
            data1 = data_dict[key1]
            data2 = data_dict[key2]
            value = compare_func(data1, data2, **kwargs)

            scores[key1][key2] = value
            scores[key2][key1] = value

            comb_keys = key1 + '_' + key2
            is_rep = False
            for rep in REPLICATES:
                if rep[0] in comb_keys and rep[1] in comb_keys:
                    replicates[comb_keys] = value
                    is_rep = True
                    break
            if not is_rep:
                non_replicates[comb_keys] = value

            print(f'Reproducibility: {value}')
            print()

    replicate_values = list(replicates.values())
    non_replicate_values = list(non_replicates.values())

    min_diff = np.min(replicate_values) - np.max(non_replicate_values)
    avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
    min_rep = np.min(replicate_values)
    max_non_rep = np.max(non_replicate_values)
    print(f"Number of replicates found: {len(replicates)}")
    print(f"Number of non-replicates found: {len(non_replicates)}")
    print(f"Max non-replicate value: "
          f"{max(non_replicates, key=non_replicates.get)} - {max_non_rep}")
    print(f"Min replicate value: "
          f"{min(replicates, key=replicates.get)} - {min_rep}")
    print(f"Min diff between replicates and non-replicates: {min_diff}")
    print(f"Diff between replicate and non-replicate average: {avg_diff}")

    return scores


def main():
    bedGraph_dict = read_bedGraph_data()
    # genome_loop_dict = read_loop_data()

    test_cases = compare_bedGraph.create_test_cases(interval_size=BIN_SIZE,
                                                    start=CHROM_START,
                                                    end=CHROM_END,
                                                    step=1)
    bedGraph_stats_dict = compare_bedGraph.get_stats(bedGraph_dict, test_cases)

    bedGraph_scores = compare(bedGraph_stats_dict,
                              compare_bedGraph.compare_bedGraph_stats)
    output_to_csv(bedGraph_scores, 'pearson.csv')

    '''for name in bedGraph_dict:
        genome_loop_dict[name].find_loop_anchor_points(bedGraph_dict[name])
        genome_loop_dict[name].filter_with_bedGraph(test_cases, 'chr1',
                                                    bedGraph_stats_dict[name]['mean_list'])

    results = compare_loops.compare_loops(genome_loop_dict, INTERVAL_SIZE,
                                          WINDOW_SIZE, int(CHROM_START / WINDOW_SIZE))
    compare_loops.output_csv(results, 'results.csv')'''


if __name__ == '__main__':
    main()
