import sys
import os
import pytest

sys.path.append('.')
import chia_rep


@pytest.fixture
def loop_dict():
    data_dir = 'test/test_files'
    chrom_size_path = f'test/test_files/hg38.chrom.sizes'

    print(data_dir)
    print(chrom_size_path)

    # Make the loop_dict a global variable to be untouched so it doesn't have to be reloaded
    return chia_rep.read_data(loop_data_dir=data_dir,
                              chrom_size_file=chrom_size_path,
                              bedgraph_data_dir=data_dir,
                              peak_data_dir=data_dir,
                              chroms_to_load=['chr1'])


def test_filter_peaks(loop_dict):
    for sample in loop_dict:
        loop_dict[sample].filter_peaks(60, 'chr1')
        assert len(loop_dict[sample].peak_dict['chr1']) == 60

        # if sample == 'sampleA1':
        #     print(len(loop_dict[sample].peak_dict['chr2']))
        #     assert len(loop_dict[sample].peak_dict['chr2']) == 30


def test_package(loop_dict):

    parent_dir = 'test_results'
    directories = [
        parent_dir,
        f'{parent_dir}/results',
        f'{parent_dir}/results/info',
        f'{parent_dir}/tmp'
    ]

    for directory in directories:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    for bin_size in [5]:  # bin size kb
        # for i in [1]:
        bin_size *= 1000
        for window_size in [3]:  # window size mb
            # for j in [3]:
            window_size *= 1000000
            temp_str = f'{bin_size}_bin.{window_size}_window'

            if os.path.isfile(f'{parent_dir}/results/{temp_str}.emd_value.csv'):
                print(f'Skipping {temp_str}')
                continue

            # Make sure not to modify the original dict
            l = loop_dict
            chia_rep.preprocess(l, extra_data_dir=f'extra_data')

            rep, non_rep, emd_scores, j_scores = \
                chia_rep.compare(l, bin_size=bin_size, window_size=window_size)
            # chia_rep.output_results(rep, non_rep)
            chia_rep.output_results(rep, non_rep, f'{parent_dir}/results/info',
                                    temp_str)
            chia_rep.output_to_csv(emd_scores,
                                   f'{parent_dir}/results/{temp_str}.emd_value.csv')
            chia_rep.output_to_csv(j_scores,
                                   f'{parent_dir}/results/{temp_str}.j_value.csv')
