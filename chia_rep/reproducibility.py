import os
from collections import OrderedDict
import multiprocessing

import matplotlib.pyplot as plt
import numpy as np
import csv
from prettytable import PrettyTable
import logging
from pyBedGraph import BedGraph
from .all_loop_data import AllLoopData, DEFAULT_NUM_PEAKS

log = logging.getLogger()

VERSION = 9

REPLICATES = [
    ['LHM0011_0011H', 'LHM0014_0014H'],
    ['LHM0009_0009H', 'LHM0013_0013H'],
    ['LMP0002_0002V', 'LMP0001_0001V'],
    ['LHH0083_0083H', 'LHH0085_0085V'],
    ['LMP0004_0004V', 'LMP0003_0003V'],
    ['LHM0007_0007H', 'LHM0010H'],
    ['LHM0008_0008H', 'LHM0012_0012H'],
    ['LME0034_0034V', 'LME0028_0028N'],
    ['LHH0084_0084H', 'LHH0086_0086V'],
    ['LHH0048_0048H', 'LHH0054L_0054H'],
    ['LHH0061_0061H', 'LHH0058_0058H']
]

# Non-replicates that had high reproducibility values (> -0.5)
NON_REP_PROB = [
    'LMP0002_0002V', 'LME0033_0033V',
    'LMP0001_0001V', 'LME0033_0033V',
    'LMP0003_0003V', 'LME0028_0028N'
]

TO_CHECK = REPLICATES + NON_REP_PROB

TO_CHECK_SET = set([sample for pair in TO_CHECK for sample in pair])


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
def compare(data_dict, compare_all=True, **kwargs):
    log.info(f"Comparing {[x for x in data_dict]}\n"
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
                if key1 in rep and key2 in rep:
                    is_rep = True
                    break

            log.info(f'{key1} vs. {key2}:')
            rep_dict = data1.compare(data2, **kwargs, is_rep=is_rep)

            # rep_dict = compare_func(data1, data2, **kwargs, is_rep=is_rep)
            main_value = rep_dict['main']

            scores[key1][key2] = main_value
            scores[key2][key1] = main_value

            # Separate reproducibility values into replicates and non-replicates
            if is_rep:
                replicates[combined_keys] = rep_dict
            else:
                non_replicates[combined_keys] = rep_dict

            log.info(f'{combined_keys} Reproducibility: {main_value}')

    comp_str = f"Compared {[x for x in data_dict]}\n" \
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
def read_data(loop_data_dir, chrom_size_file, bedgraph_data_dir, peak_data_dir,
              min_loop_value=0, min_bedgraph_value=0, chrom_to_load=None):
    if not os.path.isfile(chrom_size_file):
        log.error(f"Chrom size file: {chrom_size_file} is not a valid file")
        return

    if not os.path.isdir(loop_data_dir):
        log.error(f"Loop dir: {loop_data_dir} is not a valid directory")
        return

    if not os.path.isdir(bedgraph_data_dir):
        log.error(f"bedgraph dir: {bedgraph_data_dir} is not a valid directory")
        return

    if not os.path.isdir(peak_data_dir):
        log.error(f"Peak dir: {peak_data_dir} is not a valid directory")
        return

    loop_info_dict = OrderedDict()

    log.info(os.listdir(loop_data_dir))
    log.info(f'Number of samples: {len(os.listdir(loop_data_dir))}')
    for loop_file_name in os.listdir(loop_data_dir):
        if loop_file_name.endswith('.BE3') or loop_file_name.endswith('.cis'):
            sample_name = loop_file_name.split('.')[0]
            log.info(f'Loading {sample_name} ...')

            loop_file_path = os.path.join(loop_data_dir, loop_file_name)

            peak_dict = None
            for peak_file_name in os.listdir(peak_data_dir):
                if 'peak' not in peak_file_name.lower() or \
                        sample_name not in peak_file_name:
                    continue

                is_narrowPeak = False
                if 'narrowpeak' in peak_file_name.lower():
                    is_narrowPeak = True
                elif 'broadpeak' in peak_file_name.lower():
                    is_narrowPeak = False
                else:
                    log.error(f"{peak_file_name} is an unknown peak file")
                peak_file_path = os.path.join(peak_data_dir, peak_file_name)
                peak_dict = read_peak_file(peak_file_path, is_narrowPeak)
                break

            if peak_dict is None:
                log.error(f"{sample_name} not in {peak_data_dir}. Skipping ...")
                continue

            bedgraph = None
            for bedgraph_file_name in os.listdir(bedgraph_data_dir):
                if f'{sample_name}.' in bedgraph_file_name and \
                        'bedgraph' in bedgraph_file_name.lower():
                    bedgraph_file_path = os.path.join(bedgraph_data_dir,
                                                      bedgraph_file_name)

                    bedgraph = BedGraph(chrom_size_file, bedgraph_file_path,
                                        chrom_wanted=chrom_to_load,
                                        ignore_missing_bp=False,
                                        min_value=min_bedgraph_value)
                    break

            if bedgraph is None:
                log.error(f"{sample_name}. not in {bedgraph_data_dir}. "
                          f"Skipping")
                continue

            loop_data = AllLoopData(chrom_size_file, loop_file_path, bedgraph,
                                    peak_dict, min_loop_value=min_loop_value,
                                    chrom_to_load=chrom_to_load)
            loop_info_dict[loop_data.sample_name] = loop_data

    return loop_info_dict


def preprocess(loop_dict, num_peaks=DEFAULT_NUM_PEAKS):
    for sample_name in loop_dict:
        loop_dict[sample_name].preprocess(num_peaks)


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


def read_peak_file(peak_file_path, is_narrowPeak):
    peak_dict = {}

    if is_narrowPeak:
        value_index = 8
    else:  # .broadPeak files
        value_index = 6

    with open(peak_file_path) as peak_file:
        for line in peak_file:
            data = line.split()
            chrom_name = data[0]

            try:
                peak_start = int(data[1])
                peak_end = int(data[2])
                peak_value = float(data[value_index])
            except ValueError:
                print("Error:", line)
                continue

            if chrom_name not in peak_dict:
                peak_dict[chrom_name] = []

            peak_dict[chrom_name].append([peak_start, peak_end,
                                          peak_end - peak_start])

    return peak_dict
