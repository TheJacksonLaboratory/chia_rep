import os
import sys

sys.path.append('..')
import compare_bedGraph
import compare_loops
import matplotlib.pyplot as plt
import numpy as np
import csv
from prettytable import PrettyTable
import logging as log

CHROM_START = 36000000
CHROM_END = 39000000
CHROM_START = 0
CHROM_END = compare_bedGraph.chr1_size
WINDOW_SIZE = CHROM_END - CHROM_START
WINDOW_SIZE = 3000000
BIN_SIZE = 10000

VERSION = 2

REPLICATES = [
    ['LHH0054H', 'LHH0048H'],
    ['LHH0054H', 'LHH0048'],
    ['LHH0054H', 'LHH0054L'],
    ['LHH0048H', 'LHH0048'],
    ['LHH0048H', 'LHH0054L'],
    ['LHH0048', 'LHH0054L'],
    ['LHH0084H', 'LHH0086V'],
    ['LHH0084H', 'LHH0086'],
    ['LHH0084H', 'LHH0084'],
    ['LHH0086V', 'LHH0084'],
    ['LHH0086V', 'LHH0086'],
    ['LHH0086', 'LHH0084'],
    ['LHH0083H', 'LHH0085V'],
    ['LHH0083H', 'LHH0085'],
    ['LHH0083H', 'LHH0083'],
    ['LHH0085V', 'LHH0083'],
    ['LHH0085V', 'LHH0085'],
    ['LHH0085', 'LHH0083'],
    ['LHH0058H', 'LHH0061H'],
    ['LHH0058H', 'LHH0061'],
    ['LHH0058H', 'LHH0058'],
    ['LHH0061H', 'LHH0058'],
    ['LHH0061H', 'LHH0061'],
    ['LHH0061', 'LHH0058']
]

TO_CHECK = [
               ['LHH0054H', 'LHH0084H'],
               ['LHH0054H', 'LHH0086V'],
               ['LHH0054H', 'LHH0086'],
               ['LHH0054H', 'LHH0084'],
               ['LHH0048H', 'LHH0084H'],
               ['LHH0048H', 'LHH0086V'],
               ['LHH0048H', 'LHH0086'],
               ['LHH0048H', 'LHH0084'],
               ['LHH0083H', 'LHH0058H'],
               ['LHH0083H', 'LHH0061H'],
               ['LHH0085V', 'LHH0058H'],
               ['LHH0085V', 'LHH0061H'],
               ['LHH0058H', 'LHH0083'],
               ['LHH0058H', 'LHH0085'],
               ['LHH0061H', 'LHH0085'],
               ['LHH0086', 'LHH0054L'],
               ['LHH0084', 'LHH0054L']
           ] + REPLICATES


def read_loop_data(hiseq=False):
    genome_loop_dict = compare_loops.read_loops('CTCF', hiseq=hiseq)
    genome_loop_dict.update(
        compare_loops.read_loops('RNAPII', hiseq=hiseq))

    return genome_loop_dict


def read_loop_bedgraph():
    bedGraph_dict = \
        compare_bedGraph.read_bedGraphs(data_directory='loop_bedgraphs')
    return bedGraph_dict


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
    log.info(f"Comparing {[x for x in data_dict]}\n"
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
            combined_keys = key1 + '_' + key2

            want_to_check = False
            for to_check in TO_CHECK:
                if key1 in to_check and key2 in to_check:
                    want_to_check = True
                    break
            if not want_to_check:
                continue

            log.info(f'{key1} vs. {key2}:')
            data1 = data_dict[key1]
            data2 = data_dict[key2]
            reproducibility_values = compare_func(data1, data2, **kwargs)
            main_value = reproducibility_values['main_value']

            scores[key1][key2] = main_value
            scores[key2][key1] = main_value

            is_rep = False
            for rep in REPLICATES:
                if key1 in rep and key2 in rep:
                    replicates[combined_keys] = reproducibility_values
                    is_rep = True
                    break
            if not is_rep:
                non_replicates[combined_keys] = reproducibility_values

            log.info(
                f'{combined_keys} Reproducibility: {main_value}\n'
                '-----------------------------------------------------------'
                '-----------------------------------------------------------')

    comp_str = f"Compared {[x for x in data_dict]}\n" \
               f"Function: {compare_func}\n" \
               f"Arguments: {kwargs}"
    log.info(comp_str)

    return replicates, non_replicates, scores


def output_results(replicates, non_replicates):
    value_keys = list(list(replicates.values())[0].keys())
    value_table = PrettyTable(['output'] + value_keys)
    value_table.sortby = 'main_value'
    value_table.reversesort = True
    for key, reproducibility_value in replicates.items():
        value_table.add_row([key] + [x if isinstance(x, str) else round(x, 6)
                                     for x in reproducibility_value.values()])
    log.info(value_table)
    value_table.clear_rows()
    for key, reproducibility_value in non_replicates.items():
        value_table.add_row([key] + [x if isinstance(x, str) else round(x, 6)
                                     for x in reproducibility_value.values()])
    log.info(value_table)

    replicate_values = [x['main_value'] for x in replicates.values()]
    non_replicate_values = [x['main_value'] for x in non_replicates.values()]

    replicate_main_values = {key: x['main_value'] for key, x in
                             replicates.items()}
    non_replicate_main_values = {key: x['main_value'] for key, x in
                                 non_replicates.items()}

    min_diff = np.min(replicate_values) - np.max(non_replicate_values)
    avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
    min_rep = np.min(replicate_values)
    max_non_rep = np.max(non_replicate_values)
    log.info(f"Number of replicates found: {len(replicates)}")
    log.info(f"Number of non-replicates found: {len(non_replicates)}")
    log.info(f"Max non-replicate value: "
             f"{max(non_replicate_main_values, key=non_replicate_main_values.get)}"
             f" -> {max_non_rep}")
    log.info(f"Min replicate value: "
             f"{min(replicate_main_values, key=replicate_main_values.get)}"
             f" -> {min_rep}")
    log.info(f"Min diff between replicates and non-replicates: {min_diff}")
    log.info(f"Diff between replicate and non-replicate average: {avg_diff}")
