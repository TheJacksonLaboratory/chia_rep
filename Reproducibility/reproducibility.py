import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from prettytable import PrettyTable
import logging
from pyBedGraph import BedGraph
from .all_loop_data import AllLoopData

log = logging.getLogger(__name__.split('.')[-1])

CHROM_START = 36000000
CHROM_END = 39000000
CHROM_START = 0
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

TO_CHECK = REPLICATES + [
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

TO_CHECK += [['LHH0048', 'LHH0058']]

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
            main_value = reproducibility_values[0]['main_value']

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


def output_results(rep, non_rep):
    log.info(f"Number of replicates found: {len(rep)}")
    log.info(f"Number of non-replicates found: {len(non_rep)}")
    num_stats = len(list(rep.values())[0])
    # Dictionaries with value as list of dictionaries
    for i in range(num_stats):
        log.info(f'Stat_{i}')
        value_keys = list(list(rep.values())[0][i].keys())
        value_table = PrettyTable(['output'] + value_keys)
        value_table.sortby = 'main_value'
        value_table.reversesort = True

        rep_value_dict = {x: rep[x][i] for x in rep}
        for k, rep_value in rep_value_dict.items():
            value_table.add_row([k] + [x if isinstance(x, str) else round(x, 6)
                                       for x in rep_value.values()])
        log.info('Replicates\n' + value_table.get_string())
        value_table.clear_rows()
        non_rep_value_dict = {x: non_rep[x][i] for x in non_rep}
        for k, rep_value in non_rep_value_dict.items():
            value_table.add_row([k] + [x if isinstance(x, str) else round(x, 6)
                                       for x in rep_value.values()])
        log.info('Non-Replicates\n' + value_table.get_string())

        replicate_values = [x['main_value'] for x in
                            list(rep_value_dict.values())]
        non_replicate_values = [x['main_value'] for x in
                                list(non_rep_value_dict.values())]

        replicate_main_values = {key: x['main_value'] for key, x in
                                 rep_value_dict.items()}
        non_replicate_main_values = {key: x['main_value'] for key, x in
                                     non_rep_value_dict.items()}

        if len(non_replicate_values) == 0 or len(replicate_values) == 0:
            return

        min_diff = np.min(replicate_values) - np.max(non_replicate_values)
        avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
        min_rep = np.min(replicate_values)
        max_non_rep = np.max(non_replicate_values)
        log.info(f"Min replicate value: "
                 f"{min(replicate_main_values, key=replicate_main_values.get)}"
                 f" -> {min_rep}")
        log.info(f"Max non-replicate value: "
                 f"{max(non_replicate_main_values, key=non_replicate_main_values.get)}"
                 f" -> {max_non_rep}")
        log.info(f"Min diff between replicates and non-replicates: {min_diff}")
        log.info(
            f"Diff between replicate and non-replicate average: {avg_diff}\n")


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
        bedgraph_dict = read_bedGraphs(bedgraph_data_dir, chrom_size_file,
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


def read_bedGraphs(data_directory, chrom_size_file, min_value=-1,
                   chrom_to_load=None):
    bedGraph_dict = {}

    for filename in os.listdir(data_directory):
        if filename.endswith('.bedgraph') or filename.endswith('.bedGraph'):
            file_path = os.path.join(data_directory, filename)
            bedGraph = BedGraph(chrom_size_file, file_path, chrom_to_load,
                                ignore_missing_bp=False, min_value=min_value)
            bedGraph_dict[bedGraph.name] = bedGraph

    return bedGraph_dict
