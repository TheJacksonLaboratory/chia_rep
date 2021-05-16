import sys
sys.path.append('.')
import chia_rep


def test_read_peak_file():
    sampleA1_peaks = [
        [30801330, 30802398, 0],
        [226407716, 226408581, 1],
        [150281286, 150281708, 1252]
    ]
    peak_dict = chia_rep.read_peak_file('test/test_files/sampleA1.broadPeak')
    assert 'chr1' in peak_dict
    for peak in sampleA1_peaks:
        peak_start = peak[0]
        peak_end = peak[1]
        peak_len = peak_end - peak_start
        assert peak_dict['chr1'][peak[2]] == [peak_start, peak_end, peak_len]


def test_read_data():
    sample_dict = chia_rep.read_data('test/sample_input_file.txt',
                                     'test/test_files/hg38.chrom.sizes')

    assert 'sampleA1' in sample_dict
    assert 'sampleA2' in sample_dict
    assert 'sampleB1' in sample_dict

    assert 'chr1' in sample_dict['sampleA1'].chrom_dict
    assert 'chr1' in sample_dict['sampleA2'].chrom_dict
    assert 'chr1' in sample_dict['sampleB1'].chrom_dict
    assert sample_dict['sampleA1'].chrom_dict['chr1'].numb_loops == 232525

    pet_count_values = [[39051, 6],
                        [46456, 6],
                        [57845, 7],
                        [69292, 6],
                        [139492, 7],
                        [140486, 6],
                        [140488, 7],
                        [169264, 7],
                        [171005, 8],
                        [189477, 7],
                        [213365, 7]]
    for value in pet_count_values:
        i = value[0]
        pet_count = value[1]
        assert sample_dict['sampleA1'].chrom_dict['chr1'].pet_count_list[i] == pet_count
