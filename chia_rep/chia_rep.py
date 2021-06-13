import os
from collections import OrderedDict
import time
import csv
import logging
from pyBedGraph import BedGraph
from typing import Dict, List

from .genome_loop_data import GenomeLoopData
import numpy as np

log = logging.getLogger()
# score_dict = OrderedDict[str, OrderedDict[str, float]]
score_dict = Dict[str, Dict[str, float]]


def output_score(
    scores: score_dict,
    file_path: str
) -> None:
    """
    Output scores to a specified file path.

    Parameters
    ----------
    scores
    file_path

    Returns
    -------
    None
    """
    with open(f'{file_path}', 'w') as out_file:
        header = ['Sample Name'] + list(scores.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)
        writer.writeheader()
        for sample_name, sample_scores in scores.items():
            writer.writerow(sample_scores)


def output_to_csv(
    emd_scores: score_dict,
    j_scores: score_dict,
    window_size: int,
    bin_size: int,
    num_peaks: any,
    output_dir: str = 'output'
) -> None:
    """
    Outputs emd_scores and j_scores to a specified file path.

    Parameters
    ----------
    emd_scores
    j_scores
    window_size
    bin_size
    num_peaks
    output_dir
        Directory to output scores

    Returns
    ------
    None
    """
    param_str = f'{window_size}.{bin_size}.{num_peaks}'
    score_dir = f'{output_dir}/{param_str}/scores'
    os.makedirs(score_dir, exist_ok=True)
    output_score(emd_scores, f'{score_dir}/emd_complete.csv')
    output_score(j_scores, f'{score_dir}/j_complete.csv')

    log.info(f"Results have been written to {score_dir}")


def compare(
    sample_dict: OrderedDict,
    num_peaks: any,
    compare_list: list = None,
    compare_list_file: str = None,
    output_dir: str = 'output',
    window_size: int = 3000000,
    bin_size: int = 5000,
    do_output_graph: bool = False
) -> (score_dict, score_dict):
    """
    Compares specified samples against each other. Specify comparisons in either
    compare_list or compare_list_file.

    Parameters
    ----------
    sample_dict : OrderedDict
        (Key: sample name, value: sample data)
    num_peaks : any
        Number of peaks to use when creating folder name for results.
    compare_list : list, optional
        List of comparisons to make. Shape is n x 2 where n is the number of
        comparisons
    compare_list_file : str, optional
        File that contains a list of comparisons to make.
        Format:
        sample1_name    sample2_name
        sample3_name    sample4_name
        ...
    output_dir : str
        Directory to output data
    window_size : int
        Split each chromosome into windows with specified size.
    bin_size : int
        Split each window into bins with with specified size.
    do_output_graph : bool
        Output generated matrix to file. Outputs two arrays for each window in
        each comparison of length (window_size / bin_size)^2.

    Returns
    -------
    tuple
        Two OrderedDicts containing emd_scores and j_scores
    """
    total_start_time = time.time()
    os.makedirs(f'{output_dir}/timings', exist_ok=True)
    sample_list = list(sample_dict.keys())

    if compare_list is None:
        if compare_list_file is None or not os.path.isfile(compare_list_file):
            log.error(f"{compare_list_file} is not a valid file")
            return OrderedDict(), OrderedDict()

        to_compare_list = []
        with open(compare_list_file) as in_file:
            for line in in_file:
                comparison = line.split()
                if len(comparison) != 2:
                    log.error(f'Invalid number of columns in {compare_list_file}')
                    return OrderedDict(), OrderedDict()

                to_compare_list.append(comparison)
    else:
        for comparison in compare_list:
            if len(comparison) != 2:
                log.error(f'Invalid list length in {comparison}')
                return OrderedDict(), OrderedDict()

        to_compare_list = compare_list

    # To easily output in .csv format
    emd_scores = OrderedDict()
    j_scores = OrderedDict()
    for key in sample_list:
        emd_scores[key] = OrderedDict()
        emd_scores[key][key] = 1
        emd_scores[key]['Sample Name'] = key

        j_scores[key] = OrderedDict()
        j_scores[key][key] = 1
        j_scores[key]['Sample Name'] = key

    comparison_timings = OrderedDict()
    for comparison in to_compare_list:
        comparison_start_time = time.time()

        sample1_name, sample2_name = comparison
        sample1 = sample_dict[sample1_name]
        sample2 = sample_dict[sample2_name]
        comparison_name = f'{sample1_name}_{sample2_name}'
        log.info(f'Compare {comparison_name}')

        if sample1.species_name != sample2.species_name:
            log.error('Tried to compare two different species. Skipping')

        value_dict = sample1.compare(sample2, window_size, bin_size, num_peaks,
                                     output_dir=output_dir,
                                     do_output_graph=do_output_graph)

        # Save values in OrderedDict
        emd_value = value_dict['emd_value']
        j_value = value_dict['j_value']
        emd_scores[sample1_name][sample2_name] = emd_value
        emd_scores[sample2_name][sample1_name] = emd_value
        j_scores[sample1_name][sample2_name] = j_value
        j_scores[sample2_name][sample1_name] = j_value
        log.info(f'{comparison_name} emd_value: {emd_value}')
        log.info(f'{comparison_name} j_value: {j_value}')

        comparison_timings[comparison_name] = time.time() - comparison_start_time

    param_str = f'{window_size}.{bin_size}.{num_peaks}'
    with open(f'{output_dir}/timings/comparison.{param_str}.txt',
              'w') as out_file:
        out_file.write(f'comparison\ttime_taken\n')
        for comparison_name, compare_timing in comparison_timings.items():
            out_file.write(f'{comparison_name}\t{compare_timing}\n')
        out_file.write(f'total\t{time.time() - total_start_time}\n')

    return emd_scores, j_scores


def check_results(rep, non_rep, out_file_dir=None, desc_str=None):
    """
    Deprecated:
    Outputs results in readable text file table format

    Parameters
    ----------
    rep : dict(str, dict)
        Key: Combined comparison name (LHH0048_LHH0054L)
        Value: dict containing keys: 'emd_value' and/or 'j_value'
        Dictionary containing information for replicate comparisons
    non_rep : dict(str, dict)
        Key: Combined comparison name (LHH0048_LHH0054L)
        Value: dict containing keys: 'emd_value' and/or 'j_value'
        Dictionary containing information for non-replicate comparisons
    out_file_dir : str, optional
        Default is None
    desc_str : str, optional
        Used as file name to describe settings/parameters for this comparison
        Not optional if out_file_dir is not None
        Default is None
    """

    log.info(f"Number of known replicates: {len(rep)}")
    log.info(f"Number of non-replicates or unknown: {len(non_rep)}")

    for value_type in ['emd_value', 'j_value']:
        out_file = None
        out_file_path = None
        if out_file_dir:
            if not os.path.isdir(out_file_dir):
                os.mkdir(out_file_dir)

            if desc_str is None:
                out_file_path = os.path.join(out_file_dir, f'{value_type}.txt')
            else:
                out_file_path = os.path.join(out_file_dir,
                                             f'{desc_str}.{value_type}.txt')

            out_file = open(out_file_path, 'w')

        rep_table = PrettyTable(['Comparison', 'emd_value', 'j_value'])
        rep_table.sortby = value_type
        rep_table.reversesort = True

        for comparison_value in [rep, non_rep]:
            if len(comparison_value) == 0:
                continue

            for k, value_dict in comparison_value.items():
                rep_table.add_row([k, round(value_dict['emd_value'], 5),
                                   round(value_dict['j_value'], 5)])

            if comparison_value == rep:
                temp_str = f'Replicates sorted by {value_type}'
            else:
                temp_str = f'Non-Replicates sorted by {value_type}'
            temp_str += f'\n{rep_table.get_string()}'
            if out_file:
                out_file.write(temp_str + '\n')
            log.info(temp_str)
            rep_table.clear_rows()

        replicate_values = [x[value_type] for x in rep.values()]
        non_replicate_values = [x[value_type] for x in non_rep.values()]

        # No more statistics can be made without knowing replicates
        if len(non_replicate_values) == 0 or len(replicate_values) == 0:
            return

        min_diff = np.min(replicate_values) - np.max(non_replicate_values)
        avg_diff = np.mean(replicate_values) - np.mean(non_replicate_values)
        min_rep = np.min(replicate_values)
        max_non_rep = np.max(non_replicate_values)
        temp_str = f"Min replicate value: " \
                   f"{min(rep, key=lambda x: rep[x][value_type])} -> {min_rep}\n" \
                   f"Max non-replicate value: " \
                   f"{max(non_rep, key=lambda x: non_rep[x][value_type])} -> {max_non_rep}\n" \
                   f"Min diff between replicates and non-replicates: {min_diff}\n" \
                   f"Diff between replicate and non-replicate average: {avg_diff}"
        log.info(temp_str)

        if out_file_path:
            out_file.write(temp_str + '\n')
            log.info(f"Results have been written to {out_file_path}")
            out_file.close()


def read_data(
    input_data_file: str,
    chrom_size_file: str,
    min_loop_value: int = 1,
    min_bedgraph_value: int = 1,
    chroms_to_load: List[str] = None,
    use_bigwig: bool = False,
    output_dir: str = 'output'
) -> Dict[str, GenomeLoopData]:
    """
    Reads all samples that are found in loop_data_dir.

    loop_data_dir/peak_data_dir/bedgraph_data_dir do not have to be separate
    directories.

    Parameters
    ----------
    input_data_file : str
        File with file paths to all necessary input files.
        Format:
        sample1_name bedgraph1_file   peak1_file   loop2_file
        sample2_name bedgraph2_file   peak2_file   loop1_file
        ...
    chrom_size_file : str
        Path to chromosome size file
    min_loop_value : int, optional
        Minimum loop value accepted by GenomeLoopData/ChromLoopData
    min_bedgraph_value : int, optional
        Minimum value accepted by BedGraph obj from pyBedGraph
    chroms_to_load : list, optional
        Specify specific chromosomes to load instead of the entire genome
    use_bigwig : bool, optional
        Specify if input_file is bigwig or not. Not implemented yet.
    output_dir : str
        Directory to output data

    Returns
    -------
    OrderedDict[str, GenomeLoopData]
    """
    total_start_time = time.time()
    os.makedirs(f'{output_dir}/timings', exist_ok=True)
    sample_data_dict = OrderedDict()

    if not os.path.isfile(chrom_size_file):
        log.error(f"Chrom size file: {chrom_size_file} is not a valid file")
        return sample_data_dict

    if not os.path.isfile(input_data_file):
        log.error(f"Data file: {input_data_file} is not a valid file")
        return sample_data_dict

    # Get input file names
    input_sample_files = []
    with open(input_data_file) as in_file:
        for line in in_file:
            sample_files = line.split()
            if len(sample_files) != 4:
                log.error(f"Invalid number of columns in {input_data_file}")
                return sample_data_dict
            input_sample_files.append(sample_files)

    sample_timings = OrderedDict()
    for sample_files in input_sample_files:
        sample_start_time = time.time()

        sample_name = sample_files[0]
        bedgraph_file = sample_files[1]
        peak_file = sample_files[2]
        loop_file = sample_files[3]

        # Check for file validity
        invalid_file = False
        for i in range(1, 4):
            if not os.path.isfile(sample_files[i]):
                log.error(f"Data file: {sample_files[i]} is not a valid file")
                invalid_file = True
                break
        if invalid_file:
            continue

        log.info(f'Loading {sample_name} ...')

        peak_dict = read_peak_file(peak_file)
        bedgraph = BedGraph(chrom_size_file, bedgraph_file,
                            chroms_to_load=chroms_to_load,
                            ignore_missing_bp=False,
                            min_value=min_bedgraph_value)

        gld = GenomeLoopData(chrom_size_file, loop_file, bedgraph,
                             peak_dict, min_loop_value=min_loop_value,
                             chroms_to_load=chroms_to_load)
        sample_data_dict[sample_name] = gld
        sample_timings[sample_name] = time.time() - sample_start_time

    with open(f'{output_dir}/timings/read_data.txt', 'w') as out_file:
        out_file.write(f'sample_name\ttime_taken\n')
        for sample_name, sample_timing in sample_timings.items():
            out_file.write(f'{sample_name}\t{sample_timing}\n')
        out_file.write(f'total\t{time.time() - total_start_time}\n')

    return sample_data_dict


def preprocess(
    sample_dict: OrderedDict,
    num_peaks: int = None,
    both_peak_support: bool = False,
    output_dir: str = 'output',
    base_chrom: str = 'chr1'
) -> None:
    """
    Takes num_peaks top peaks from the peak file for base_chrom and uses the
    same ratio for each other chromosome.
    Preprocess all the chromosomes in this object.

    Removes all problematic chromosomes (not enough loops, etc...).

    Parameters
    ----------
    sample_dict : OrderedDict[str, GenomeLoopData]
        Samples to compare
    num_peaks : int, optional
        The number of peaks to use when filtering. Only to be used with chr1
        since other chromosomes will be dependent on min peak used from chr1
    both_peak_support : bool, optional
        Whether to only keep loops that have peak support on both sides
    output_dir : str, optional
        Directory to output found peaks and filters
    base_chrom : str, optional
        Chromosome to use with num_peaks to equalize number of peaks used for
        each chromosome.

    Returns
    ------
    None
    """
    for sample_data in sample_dict.values():
        sample_data.preprocess(num_peaks, both_peak_support=both_peak_support,
                               output_dir=output_dir, base_chrom=base_chrom)


def read_peak_file(
    peak_file_path: str
) -> Dict[str, list]:
    """
    Finds the start and ends of every peak in chromosome for one sample. Find
    the max value within the interval for each peak using the bedgraph file.

    File format must have at least 3 columns with the first 3 being:
    chrom_name   start   end

    Parameters
    ----------
    peak_file_path : str
        File path of peak file

    Returns
    -------
    dict[str, list]
        Dictionary containing list of peaks start and ends for every chromosome
    """
    peak_dict = {}

    with open(peak_file_path) as peak_file:
        for line in peak_file:
            data = line.split()
            chrom_name = data[0]

            try:
                peak_start = int(data[1])
                peak_end = int(data[2])
            except ValueError:  # Sometimes have peaks with 1+E08 as a value
                log.error(f'Invalid peak: {line}')
                continue

            if chrom_name not in peak_dict:
                peak_dict[chrom_name] = []

            peak_dict[chrom_name].append([peak_start, peak_end,
                                          peak_end - peak_start])
    return peak_dict
