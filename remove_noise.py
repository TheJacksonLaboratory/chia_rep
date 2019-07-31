import numpy as np
import os


# Removes noise from:
# Blacklist file
# Noisy Peaks
# arg1 - list of BedGraph objects from pyBedGraph
# arg2 - Set to True to see potentially removed peaks
def remove_noise(bedGraph_list, no_modify=False):
    for bedGraph in bedGraph_list:
        for chrom_name in bedGraph.chromosome_map:
            bedGraph.load_chrom_data(chrom_name)

            chrom = bedGraph.chromosome_map[chrom_name]

            remove_noisy_peaks(bedGraph.name, chrom, no_modify=no_modify)

            # Don't free yet since it is used later (prob)
            if len(bedGraph_list) > 1:
                bedGraph.free_chrom_data(chrom_name)
                print()


def remove_noisy_peaks(bedGraph_name, chrom, no_modify=False):
    if not no_modify:
        print(f"Removing noisy peaks from {bedGraph_name}")
    else:
        print(f"Potentially noisy peaks in {bedGraph_name}")

    # average_value == sum of values / # of peaks
    # TODO: could potentially weight each value based on interval width
    average_value = np.sum(chrom.value_map) / chrom.num_intervals
    min_threshold = average_value * 3  # change to whatever

    # Personal tests found most peaks to have groups of less than 15 for miseq
    # Made 20 to be safe
    # TODO: Could be different for hiseq
    min_group_size = 20

    # Initialize starting group
    interval_group = [0]  # holds the index of the intervals
    interval_group_values = [chrom.value_map[0]]  # holds the value of the intervals

    num_noisy_peaks = 0
    for i in range(1, chrom.num_intervals, 1):

        # There is a new group of intervals
        if chrom.intervals[0][i] != chrom.intervals[1][interval_group[-1]]:
            group_size = len(interval_group)
            group_max_value = np.max(interval_group_values)

            # Conditions for a potential noisy peak
            if group_size < min_group_size and group_max_value > min_threshold:

                max_index = -1
                for j in range(group_size):
                    if interval_group_values[j] == group_max_value:
                        max_index = j
                        break

                # Get interval bounds for the largest in the group
                interval_index = interval_group[max_index]
                interval_start = chrom.intervals[0][interval_index]
                interval_end = chrom.intervals[1][interval_index]

                # Noisy peaks should be alone within +-25 base pairs
                # TODO: Will need to change if BedGraph files don't have minimum values
                try:
                    before_index = chrom.index_list[interval_start - 25]
                except IndexError:
                    before_index = -1
                try:
                    after_index = chrom.index_list[interval_end + 25]
                except IndexError:
                    after_index = -1

                # Rationale: value in index_list should be -1 or pointing to an
                # interval that has value of ~0
                if (before_index == -1 or chrom.value_map[before_index] == 0) \
                        and (after_index == -1 or chrom.value_map[after_index] == 0):

                    if no_modify:
                        print([f"{chrom.intervals[0][x]}, {chrom.intervals[1][x]}"
                               for x in interval_group], group_max_value)
                    else:
                        # nullify entry from BedGraph object for each interval
                        # in the group
                        # print([f"{chrom.intervals[0][x]}, {chrom.intervals[1][x]}"
                        #        for x in interval_group], group_max_value)
                        for x in interval_group:
                            chrom.value_map[x] = 0

                    num_noisy_peaks += 1

            # start a new group
            interval_group.clear()
            interval_group_values.clear()

        interval_group.append(i)
        interval_group_values.append(chrom.value_map[i])

    print(f'Noisy Peaks: {num_noisy_peaks}')


def output_bedGraph(bedGraph, output_folder_path):
    output_file_path = os.path.join(output_folder_path, bedGraph.name)

    with open(output_file_path, 'w+') as out_file:
        for chrom_name in bedGraph.chromosome_map:
            chrom = bedGraph.chromosome_map[chrom_name]
            for i in range(chrom.num_intervals):
                if chrom.value_map[i] > 0:
                    out_file.write(f'{chrom.name}\t{chrom.intervals[0][i]}\t'
                                   f'{chrom.intervals[1][i]}\t{chrom.value_map[i]}\n')
