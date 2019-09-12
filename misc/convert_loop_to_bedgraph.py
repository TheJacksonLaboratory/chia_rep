import sys
import numpy as np
sys.path.append('..')
from genome_loop_data import GenomeLoopData
import combined
from compare_bedGraph import chr1_size
import logging as log

VERSION = 1


def get_intervals(bedgraph_array):
    intervals = []
    prev_value = -1
    prev_start = -1
    count = 0
    for i in range(bedgraph_array.size):
        if bedgraph_array[i] == 0:
            continue

        if prev_start == -1:
            prev_start = i
            prev_value = bedgraph_array[i]
            continue

        if prev_value != bedgraph_array[i]:
            intervals.append([prev_start, i, prev_value])
            count += 1
            prev_start = -1
            prev_value = -1

            if count % 500 == 0:
                print(count)

    log.info(len(intervals))
    return intervals


def get_bedgraph_array(chrom_loop, min_value=0):
    log.info(f"Getting bedgraph for {chrom_loop.sample_name}")

    loop_bedgraph_array = np.zeros(chr1_size, dtype=np.uint32)

    for i in range(chrom_loop.numb_values):
        if min_value > chrom_loop.value_list[i]:
            continue

        start = chrom_loop.start_list[i]
        end = chrom_loop.end_list[i]
        loop_bedgraph_array[start:end] += chrom_loop.value_list[i]

    return loop_bedgraph_array


def output_as_bedgraph(intervals, out_file, chrom_name='chr1'):
    with open(out_file, 'w+') as out_file:
        for interval in intervals:
            out_file.write(f'{chrom_name}\t{interval[0]}\t{interval[1]}\t'
                           f'{interval[2]}\n')


def main():
    loop_dict = combined.read_loop_data(True)

    for name in loop_dict:
        loop_bedgraph_array = get_bedgraph_array(loop_dict[name].chrom_dict['chr1'])

        output_as_bedgraph(loop_bedgraph_array)


if __name__ == '__main__':
    main()
