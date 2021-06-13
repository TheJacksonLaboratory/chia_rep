import sys
import os
import shutil

sys.path.append('.')
import chia_rep


def test_filter_peaks():
    sample_dict = chia_rep.read_data('test/sample_input_file.txt',
                                     'test/test_files/hg38.chrom.sizes',
                                     output_dir='test/output')
    for sample in sample_dict:
        sample_dict[sample].filter_peaks(60, 'chr1')
        assert len(sample_dict[sample].peak_dict['chr1']) == 60

        # if sample == 'sampleA1':
        #     print(len(sample_dict[sample].peak_dict['chr2']))
        #     assert len(sample_dict[sample].peak_dict['chr2']) == 30


def test_package():
    bin_size = 5000
    window_size = 3000000

    shutil.rmtree('test/output')

    sample_dict = chia_rep.read_data('test/sample_input_file.txt',
                                     'test/test_files/hg38.chrom.sizes',
                                     output_dir='test/output')

    l = sample_dict
    chia_rep.preprocess(l, output_dir='test/output')

    emd_scores, j_scores = chia_rep.compare(l, 'all',
                                            compare_list_file='test/pairs.txt',
                                            bin_size=bin_size,
                                            window_size=window_size,
                                            output_dir='test/output')
    chia_rep.output_to_csv(emd_scores, j_scores, window_size, bin_size, 'all',
                           output_dir='test/output')

    assert os.path.isfile('test/output/loops/sampleA1.all.loops')
    assert os.path.isfile('test/output/loops/sampleA2.all.loops')
    assert os.path.isfile('test/output/loops/sampleB1.all.loops')

    assert os.path.isfile('test/output/peaks/sampleA1.all.peaks')
    assert os.path.isfile('test/output/peaks/sampleA2.all.peaks')
    assert os.path.isfile('test/output/peaks/sampleB1.all.peaks')

    param = f'{window_size}.{bin_size}.all'
    assert os.path.isfile(f'test/output/{param}/scores/emd_complete.csv')
    assert os.path.isfile(f'test/output/{param}/scores/j_complete.csv')

    assert os.path.isfile(f'test/output/timings/comparison.{param}.txt')
    assert os.path.isfile(f'test/output/timings/read_data.txt')

    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA1_sampleA2_chr1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA1_sampleB1_chr1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA2_sampleB1_chr1.txt')

    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA1_sampleA2.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA1_sampleB1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA2_sampleB1.txt')


def test_package2():
    bin_size = 5000
    window_size = 3000000

    shutil.rmtree('test/output')

    sample_dict = chia_rep.read_data('test/sample_input_file.txt',
                                     'test/test_files/hg38.chrom.sizes',
                                     output_dir='test/output')
    l = sample_dict
    chia_rep.preprocess(l, output_dir='test/output')

    comparison_list = [
        ['sampleA1', 'sampleA2'],
        ['sampleA1', 'sampleB1'],
        ['sampleA2', 'sampleB1']
    ]

    emd_scores, j_scores = chia_rep.compare(l, 'all',
                                            compare_list=comparison_list,
                                            bin_size=bin_size,
                                            window_size=window_size,
                                            output_dir='test/output')

    chia_rep.output_to_csv(emd_scores, j_scores, window_size, bin_size, 'all',
                           output_dir='test/output')

    assert os.path.isfile('test/output/loops/sampleA1.all.loops')
    assert os.path.isfile('test/output/loops/sampleA2.all.loops')
    assert os.path.isfile('test/output/loops/sampleB1.all.loops')

    assert os.path.isfile('test/output/peaks/sampleA1.all.peaks')
    assert os.path.isfile('test/output/peaks/sampleA2.all.peaks')
    assert os.path.isfile('test/output/peaks/sampleB1.all.peaks')

    param = f'{window_size}.{bin_size}.all'
    assert os.path.isfile(f'test/output/{param}/scores/emd_complete.csv')
    assert os.path.isfile(f'test/output/{param}/scores/j_complete.csv')

    assert os.path.isfile(f'test/output/timings/comparison.{param}.txt')
    assert os.path.isfile(f'test/output/timings/read_data.txt')

    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA1_sampleA2_chr1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA1_sampleB1_chr1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/windows/sampleA2_sampleB1_chr1.txt')

    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA1_sampleA2.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA1_sampleB1.txt')
    assert os.path.isfile(
        f'test/output/{param}/scores/chromosomes/sampleA2_sampleB1.txt')
