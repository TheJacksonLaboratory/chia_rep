import os
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import csv
from prettytable import PrettyTable
import logging
from pyBedGraph import BedGraph
from .all_loop_data import AllLoopData, DEFAULT_PEAK_PERCENT

log = logging.getLogger()

VERSION = 8

REPLICATES = [
    ['LHH0054', 'LHH0048'],
    ['LHH0084', 'LHH0086'],
    ['LHH0083', 'LHH0085'],
    ['LHH0058', 'LHH0061'],
    ['LME0028', 'LME0034'],
    # ['LME0031', 'LME0033'],
    ['LHM0008', 'LHM0012'],
    ['LHM0010', 'LHM0007'],
    ['LHM0011', 'LHM0014'],
    ['LHM0009', 'LHM0013'],
    ['LMP0001', 'LMP0002'],
    ['LMP0003', 'LMP0004'],
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

    log.info(f"Results have been written to {output_file}")


# Make compare_all == False to save time by not comparing known reproducibility
def compare(data_dict, compare_func, compare_all=True, **kwargs):
    log.info(f"Comparing {[x for x in data_dict]}\n"
             f"Function: {compare_func}\n"
             f"Arguments: {kwargs}")

    non_replicates = {}
    replicates = {}

    len_dict = len(data_dict)
    keys = list(data_dict.keys())
    scores = OrderedDict()  # To easily output in .csv format
    for key in keys:
        scores[key] = OrderedDict()
        scores[key][key] = 1
        scores[key]['Sample Name'] = key

    for i in range(len_dict):
        for j in range(i + 1, len_dict):
            key1 = keys[i]
            key2 = keys[j]
            combined_keys = key1 + '_' + key2
            data1 = data_dict[key1]
            data2 = data_dict[key2]

            if data2.species_name != data1.species_name:
                continue

            if not compare_all:
                want_to_check = False
                for to_check in TO_CHECK:
                    if key1 in to_check and key2 in to_check:
                        want_to_check = True
                        break
                if not want_to_check:
                    continue

            is_rep = False
            for rep in REPLICATES:
                if (rep[0] in key1 or rep[0] in key2) and (rep[1] in key1
                                                           or rep[1] in key2):
                    is_rep = True
                    break

            log.info(f'{key1} vs. {key2}:')

            rep_dict = compare_func(data1, data2, **kwargs, is_rep=is_rep)
            main_value = rep_dict['main']

            scores[key1][key2] = main_value
            scores[key2][key1] = main_value

            # Separate reproducibility values into replicates and non-replicates
            is_rep = False
            for rep in REPLICATES:
                if (rep[0] in key1 or rep[0] in key2) and (rep[1] in key1
                                                           or rep[1] in key2):
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
    log.info(f"Number of replicates found: {len(rep)}")
    log.info(f"Number of non-replicates found: {len(non_rep)}")

    out_file = None
    if out_file_path:
        out_file = open(out_file_path, 'w')

    # Dictionaries of dictionaries
    value_table = PrettyTable(['output', 'graph_type', 'main'])
    value_table.sortby = 'main'
    value_table.reversesort = True

    for rep_value in [rep, non_rep]:
        for k, value_dict in rep_value.items():
            value_table.add_row([k, value_dict['graph_type'],
                                 round(value_dict['main'], 5)])
        temp_str = 'Non-Replicates\n' + value_table.get_string()
        if out_file:
            out_file.write(temp_str + '\n')
        log.info(temp_str)
        value_table.clear_rows()

    replicate_values = [x['main'] for x in rep.values()]
    non_replicate_values = [x['main'] for x in non_rep.values()]

    if len(non_replicate_values) == 0 or len(replicate_values) == 0:
        return

    min_diff = np.min(replicate_values) - np.max(non_replicate_values)
    avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
    min_rep = np.min(replicate_values)
    max_non_rep = np.max(non_replicate_values)
    temp_str = f"Min replicate value: " \
               f"{min(rep, key=lambda x: rep[x]['main'])} -> {min_rep}\n" \
               f"Max non-replicate value: " \
               f"{max(non_rep, key=lambda x: non_rep[x]['main'])} -> {max_non_rep}\n" \
               f"Min diff between replicates and non-replicates: {min_diff}\n" \
               f"Diff between replicate and non-replicate average: {avg_diff}"
    log.info(temp_str)
    if out_file:
        out_file.write(temp_str + '\n')
        log.info(f"Results have been written to {out_file_path}")


# Only uses loop files with .cis and .BE3 endings
def read_data(loop_data_dir, chrom_size_file, bedgraph_data_dir, is_hiseq=True,
              min_loop_value=0, min_bedgraph_value=-1, chrom_to_read=None):
    if not os.path.isfile(chrom_size_file):
        log.error(f"{chrom_size_file} is not a valid file")
        return

    loop_info_dict = OrderedDict()

    log.info(os.listdir(loop_data_dir))
    log.info(len(os.listdir(loop_data_dir)))
    for loop_file_name in os.listdir(loop_data_dir):
        if loop_file_name.endswith('.BE3') or loop_file_name.endswith('.cis'):
            sample_name = loop_file_name.split('.')[0]
            log.info(f'Reading {sample_name}')

            loop_file_path = os.path.join(loop_data_dir, loop_file_name)

            bedgraph = None
            for bedgraph_file_name in os.listdir(bedgraph_data_dir):
                if f'{sample_name}.' in bedgraph_file_name and \
                        'bedgraph' in bedgraph_file_name.lower():
                    bedgraph_file_path = os.path.join(bedgraph_data_dir,
                                                      bedgraph_file_name)

                    bedgraph = BedGraph(chrom_size_file, bedgraph_file_path,
                                        chrom_wanted=chrom_to_read,
                                        ignore_missing_bp=False,
                                        min_value=min_bedgraph_value)
                    break

            if bedgraph is None:
                log.error(f"{sample_name} not in {bedgraph_data_dir}. Skipping "
                          f"...")
                continue

            loop_data = AllLoopData(chrom_size_file, loop_file_path, bedgraph,
                                    is_hiseq, min_loop_value=min_loop_value,
                                    wanted_chroms=[chrom_to_read])
            loop_info_dict[loop_data.sample_name] = loop_data

    return loop_info_dict


def preprocess(loop_dict, peak_data_dir,
               peak_percentage_kept=DEFAULT_PEAK_PERCENT, window_size=None):
    for sample_name in loop_dict:
        peak_file_path = None
        for peak_file_name in os.listdir(peak_data_dir):
            if f'{sample_name}.' in peak_file_name and \
                    'peak' in peak_file_name.lower():
                peak_file_path = os.path.join(peak_data_dir, peak_file_name)
                break

        if peak_file_path is None:
            log.error(f"Sample name from loop dict |{sample_name}| not "
                      f"found in peak directory: \n"
                      f"{os.listdir(peak_data_dir)}")
            continue

        loop_dict[sample_name].preprocess(peak_file_path,
                                          peak_percentage_kept, window_size)


def read_bedgraphs(data_directory, chrom_size_file, min_value=-1,
                   chrom_to_read=None):
    bedgraph_dict = {}

    for filename in os.listdir(data_directory):
        if filename.lower().endswith('.bedgraph') or filename.endswith('peaks'):
            file_path = os.path.join(data_directory, filename)
            bedgraph = BedGraph(chrom_size_file, file_path, chrom_to_read,
                                ignore_missing_bp=False, min_value=min_value)
            bedgraph_dict[bedgraph.name] = bedgraph

    return bedgraph_dict
