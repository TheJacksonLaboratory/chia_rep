import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from prettytable import PrettyTable
import logging
from pyBedGraph import BedGraph
from .all_loop_data import AllLoopData

log = logging.getLogger()

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
    ['LHH0061', 'LHH0058'],
    ['sampleA1', 'sampleA2']
]

# Non-replicates that had high reproducibility values (> -0.5)
NON_REP_PROB = [
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
]

TO_CHECK = REPLICATES + NON_REP_PROB

TO_CHECK_SET = set([sample for pair in TO_CHECK for sample in pair])


def test_log_func():
    log.info('info log')
    log.warning('warning log')
    print("normal print statement")


def output_to_csv(scores, output_file):
    with open(output_file, 'w+') as out_file:
        header = list(scores.keys())
        header.insert(0, 'Sample Name')
        writer = csv.DictWriter(out_file, fieldnames=header)

        writer.writeheader()
        for sample_name, sample_scores in scores.items():
            writer.writerow(sample_scores)


# Make compare_all == False to save time by not comparing known reproducibility
def compare(data_dict, compare_func, compare_all=True, **kwargs):
    log.info(f"Comparing {[x for x in data_dict]}\n"
             f"Function: {compare_func}\n"
             f"Arguments: {kwargs}")

    non_replicates = {}
    replicates = {}

    len_dict = len(data_dict)
    keys = list(data_dict.keys())
    scores = {}  # To easily output in .csv format
    for key in keys:
        scores[key] = {}
        scores[key][key] = 1
        scores[key]['Sample Name'] = key

    for i in range(len_dict):
        for j in range(i + 1, len_dict):
            key1 = keys[i]
            key2 = keys[j]
            combined_keys = key1 + '_' + key2

            if not compare_all:
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
            rep_dict = compare_func(data1, data2, **kwargs)
            main_value = rep_dict['w_rep']

            scores[key1][key2] = main_value
            scores[key2][key1] = main_value

            # Separate reproducibility values into replicates and non-replicates
            is_rep = False
            for rep in REPLICATES:
                if key1 in rep and key2 in rep:
                    replicates[combined_keys] = rep_dict
                    is_rep = True
                    break
            if not is_rep:
                non_replicates[combined_keys] = rep_dict

            log.info(f'{combined_keys} Reproducibility: {main_value}')

    comp_str = f"Compared {[x for x in data_dict]}\n" \
               f"Function: {compare_func}\n" \
               f"Arguments: {kwargs}"
    log.info(comp_str)

    return replicates, non_replicates, scores


def output_results(rep, non_rep, out_file_path=None):
    print(f"Number of replicates found: {len(rep)}")
    print(f"Number of non-replicates found: {len(non_rep)}")

    out_file = None
    if out_file_path:
        out_file = open(out_file_path, 'w')

    # Dictionaries of dictionaries
    value_keys = ['graph_type', 'rep', 'w_rep']
    value_table = PrettyTable(['output'] + value_keys)
    value_table.sortby = 'w_rep'
    value_table.reversesort = True

    for k, value_dict in rep.items():
        value_table.add_row([k] + [x if isinstance(x, str) else round(x, 6)
                                   for x in value_dict.values()])
    temp_str = 'Replicates\n' + value_table.get_string()
    if out_file:
        out_file.write(temp_str + '\n')
    print(temp_str)

    value_table.clear_rows()
    for k, value_dict in non_rep.items():
        value_table.add_row([k] + [x if isinstance(x, str) else round(x, 6)
                                   for x in value_dict.values()])
    temp_str = 'Non-Replicates\n' + value_table.get_string()
    if out_file:
        out_file.write(temp_str + '\n')
    print(temp_str)

    replicate_values = [x['w_rep'] for x in rep.values()]
    non_replicate_values = [x['w_rep'] for x in non_rep.values()]

    if len(non_replicate_values) == 0 or len(replicate_values) == 0:
        return

    min_diff = np.min(replicate_values) - np.max(non_replicate_values)
    avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
    min_rep = np.min(replicate_values)
    max_non_rep = np.max(non_replicate_values)
    temp_str = f"Min replicate value: " \
               f"{min(rep, key=lambda x: rep[x]['w_rep'])} -> {min_rep}\n" \
               f"Max non-replicate value: " \
               f"{max(non_rep, key=lambda x: non_rep[x]['w_rep'])} -> {max_non_rep}\n" \
               f"Min diff between replicates and non-replicates: {min_diff}\n" \
               f"Diff between replicate and non-replicate average: {avg_diff}"
    print(temp_str)
    if out_file:
        out_file.write(temp_str + '\n')
        log.info(f"Results have been written to {out_file_path}")


# Only uses loop files with .cis and .BE3 endings
def read_data(loop_data_dir, peak_data_dir, chrom_size_file, is_hiseq,
              bedgraph_data_dir=None, bedgraph_dict=None,
              peak_percentage_kept=0.2, min_loop_value=0,
              min_bedgraph_value=-1):
    if not os.path.isfile(chrom_size_file):
        log.error(f"{chrom_size_file} is not a valid file")
        return

    loop_info_dict = {}

    if bedgraph_dict is None:
        if bedgraph_data_dir is None:
            log.error(f"Missing directory for bedgraph_data_dir")
            return
        bedgraph_dict = read_bedgraphs(bedgraph_data_dir, chrom_size_file,
                                       min_bedgraph_value)

    for loop_file_name in os.listdir(loop_data_dir):
        if loop_file_name.endswith('.BE3') or loop_file_name.endswith('.cis'):
            sample_name = loop_file_name.split('.')[0]
            log.info(f'Reading {sample_name}')

            loop_file_path = os.path.join(loop_data_dir, loop_file_name)

            peak_file_path = None
            for peak_file_name in os.listdir(peak_data_dir):
                if sample_name in peak_file_name and \
                        'peak' in peak_file_name.lower():
                    peak_file_path = os.path.join(peak_data_dir, peak_file_name)
                    break

            if peak_file_path is None:
                log.error(f"Sample name from loop file |{sample_name}| not "
                          f"found in peak directory: \n"
                          f"{os.listdir(peak_data_dir)}")
                continue

            if sample_name not in bedgraph_dict:
                log.error(f"Sample name from loop file |{sample_name}| not "
                          f"found in bedgraph dictionary: \n"
                          f"{list(bedgraph_dict.keys())}")
                continue

            loop_data = AllLoopData(chrom_size_file, loop_file_path,
                                    peak_file_path, bedgraph_dict[sample_name],
                                    is_hiseq,
                                    peak_percent_kept=peak_percentage_kept,
                                    min_loop_value=min_loop_value)
            loop_info_dict[loop_data.sample_name] = loop_data

    return loop_info_dict


def read_bedgraphs(data_directory, chrom_size_file, min_value=-1,
                   chrom_to_load=None):
    bedgraph_dict = {}

    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph') or filename.endswith('.bedGraph'):
            file_path = os.path.join(data_directory, filename)
            bedgraph = BedGraph(chrom_size_file, file_path, chrom_to_load,
                                ignore_missing_bp=False, min_value=min_value)
            bedgraph_dict[bedgraph.name] = bedgraph

    return bedgraph_dict
