import math
from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt
from copy import deepcopy
import logging
from .chiapet_rep_util import *

from prettytable import PrettyTable

log = logging.getLogger()
log_bin = logging.getLogger('bin')
log_all = logging.getLogger('all')

EMD_WEIGHT = 0.5
J_WEIGHT = 0.5
MAX_LOOP_LEN = 1000000

VERSION = 0


# Second value returned is the weight of the specific Earth Mover's Distance
def emd(p, q):
    assert p.size == q.size
    p_sum = p.sum()
    q_sum = q.sum()

    if p_sum == 0 and q_sum == 0:
        return 0, 0

    if p_sum == 0:
        return q.size, q_sum

    if q_sum == 0:
        return p.size, p_sum

    p = p / p_sum
    q = q / q_sum
    return c_emd(p, q, p.size), q_sum + p_sum


def jensen_shannon_divergence(p, q, base=2):
    assert p.size == q.size
    p_sum = p.sum()
    q_sum = q.sum()

    if p_sum == 0 and q_sum == 0:
        return 0

    if p_sum == 0:
        return 1

    if q_sum == 0:
        return 1

    p, q = p / p_sum, q / q_sum
    # log.info(np.info(p))
    # log.info(np.info(q))
    # log.info(p[p > 0.00000001])

    m = (p + q) / 2
    return sp.entropy(p, m, base=base) / 2 + sp.entropy(q, m, base=base) / 2


def get_matching_graphs(p, q):
    assert p.size == q.size
    q1 = deepcopy(q)
    p1 = deepcopy(p)
    match_graphs(p, q, p1, q1, p[0].size)

    # Could probably be made simpler: return p1, q1
    return [[p1, p], [q1, q]]


# Finds the places where loops overlapped in one of the samples
# Uses dtype unsigned char because bool doesn't work as well in Cython
# Same memory usage in a numpy array
def get_removed_area(chrom_size, to_remove):
    removed_area = np.zeros(chrom_size, dtype=np.uint8)
    for i in range(len(to_remove[0])):
        start = to_remove[0][i]
        end = to_remove[1][i]
        removed_area[start:end] = 1

    return removed_area


class ChromLoopData:

    def __init__(self, chrom_name, chrom_size, sample_name):
        self.name = chrom_name
        self.size = int(chrom_size)
        self.sample_name = sample_name

        self.start_anchor_list = [[], []]
        self.end_anchor_list = [[], []]
        self.start_list = None  # Later initialized with bedgraph
        self.end_list = None
        self.value_list = []
        self.numb_values = 0
        self.removed_intervals = [[], []]  # Start/end anchors are on same peak

        self.max_loop_value = 0

    def add_loop(self, loop_start1, loop_start2, loop_end1, loop_end2,
                 loop_value):

        self.start_anchor_list[0].append(loop_start1)
        self.start_anchor_list[1].append(loop_start2)
        self.end_anchor_list[0].append(loop_end1)
        self.end_anchor_list[1].append(loop_end2)

        self.value_list.append(int(loop_value))

        self.numb_values += 1

    # If return False -> remove this chromosome from containing object
    def finish_init(self, bedgraph, peak_file_path, peak_percent_kept):
        self.value_list = np.asarray(self.value_list, dtype=np.uint32)
        self.start_anchor_list = np.asarray(self.start_anchor_list,
                                            dtype=np.int32)
        self.end_anchor_list = np.asarray(self.end_anchor_list, dtype=np.int32)

        if self.numb_values == 0:
            log.debug(f"No loops for {self.name}")
            return False

        # if not bedgraph.has_chrom(self.name):
        if self.name not in bedgraph.chromosome_map:
            log.warning(f"{self.name} was not found in corresponding bedgraph: "
                        f"{bedgraph.name}")
            return False

        self.find_loop_anchor_points(bedgraph)
        return self.filter_with_peaks(peak_file_path, peak_percent_kept)

    def find_loop_anchor_points(self, bedgraph):

        log.info(f'Finding anchor points for {self.sample_name}\'s {self.name}'
                 f' from {bedgraph.name}')

        bedgraph.load_chrom_data(self.name)

        # Get index of peaks in every anchor interval
        self.start_list = bedgraph.stats(start_list=self.start_anchor_list[0],
                                         end_list=self.start_anchor_list[1],
                                         chrom_name=self.name, stat='max_index')
        self.end_list = bedgraph.stats(start_list=self.end_anchor_list[0],
                                       end_list=self.end_anchor_list[1],
                                       chrom_name=self.name, stat='max_index')

        # Get peak value for every anchor interval
        start_list_peaks = bedgraph.stats(start_list=self.start_anchor_list[0],
                                          end_list=self.start_anchor_list[1],
                                          chrom_name=self.name, stat='max')
        end_list_peaks = bedgraph.stats(start_list=self.end_anchor_list[0],
                                        end_list=self.end_anchor_list[1],
                                        chrom_name=self.name, stat='max')
        bedgraph.free_chrom_data(self.name)

        for i in range(self.numb_values):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]

            # Remove anchors that have the same* peak
            if not loop_start < loop_end:
                self.value_list[i] = 0

                # Removed interval goes from
                # (start of start anchor, end of end anchor)
                self.removed_intervals[0].append(self.start_anchor_list[0][i])
                self.removed_intervals[1].append(self.end_anchor_list[1][i])
                continue

            # Weigh each loop based on its corresponding bedgraph peak
            peak_value = max(start_list_peaks[i], end_list_peaks[i])
            self.value_list[i] *= peak_value

            # Remove loops over a given threshold
            loop_length = int(loop_end) - int(loop_start)
            if loop_length > MAX_LOOP_LEN:
                self.value_list[i] = 0

        self.max_loop_value = np.max(self.value_list)

    # To speed up testing process by avoiding loading bedgraphs every time
    def preprocess(self, peak_file_path, peak_percent_kept):
        return self.filter_with_peaks(peak_file_path, peak_percent_kept)

    def filter_with_peaks(self, peak_file, peak_percent_kept, start=0,
                          end=None):
        if end is None:
            end = self.size

        log.info(f"Filtering {self.sample_name} {self.name}:{start}-{end} ...")

        start_time = time.time()
        peaks = []
        peak_length = []
        with open(peak_file) as in_file:
            log.debug(f"Reading {peak_file} ...")
            for line in in_file:
                line = line.split()
                chrom_name = line[0]
                if chrom_name != self.name:
                    continue

                peak_start = int(line[1])
                peak_end = int(line[2])
                peak_value = float(line[6])
                peaks.append([peak_start, peak_end, peak_value])
                peak_length.append(peak_end - peak_start)

        median_peak_len = np.median(peak_length)
        log.debug(f"Median length of all peaks: {median_peak_len}")
        log.debug(f'Number of peaks: {len(peaks)}')

        num_wanted_peaks = int(len(peaks) * peak_percent_kept)

        if num_wanted_peaks == 0:
            log.error(f"Not enough peaks were found for {self.name} in "
                      f"{peak_file}. Only {peak_percent_kept * 100}% "
                      f"of {len(peaks)} were kept")
            return False

        wanted_peaks = peaks[:num_wanted_peaks]
        peak_length = peak_length[:num_wanted_peaks]
        min_peak_value = wanted_peaks[-1][2]  # Since peak file comes sorted

        log.debug(f"Median length of wanted peaks: {np.median(peak_length)}")
        log.debug(f"Num wanted peaks: {num_wanted_peaks}")
        log.debug(f'Min peak value: {min_peak_value}')

        # Get the coverage of each wanted peak
        index_array = np.full(self.size, False, dtype=bool)
        for i in range(num_wanted_peaks):
            peak_start = wanted_peaks[i][0]
            peak_end = wanted_peaks[i][1]
            index_array[peak_start:peak_end] = True

        log.debug(f'Time: {time.time() - start_time}')

        '''plt.close()
        plt.hist(peak_length, bins=20)
        plt.xlabel('Peak Length')
        plt.ylabel('Frequency')
        plt.title(f'{self.sample_name} Peak Length Frequency: '
                  f'{median_peak_len}')
        plt.savefig(f'{self.sample_name}_peak_length.png')
        plt.close()

        dist_between_peaks = []
        i = peak_start
        while i < index_array.size:
            if index_array[i]:
                peak_dist = i - peak_start
                if peak_dist < 0:
                    log.info(i, peak_start)
                    raise ValueError
                if peak_dist != 0:
                    dist_between_peaks.append(peak_dist)
                while index_array[i]:
                    i += 1
                peak_start = i
            i += 1

        median_peak_dist = np.median(dist_between_peaks)
        log.info(f"Median length of distance between peaks: "
              f"{median_peak_dist}")

        plt.close()
        plt.hist(np.log10(dist_between_peaks), bins=20)
        plt.xlabel('log10(Distance between peaks)')
        plt.ylabel('Frequency')
        plt.title(f'{self.sample_name} Distance between peaks: '
                  f'{median_peak_dist}')
        plt.savefig(f'{self.sample_name}_peak_dist.png')
        plt.close()'''

        # Should move below to Cython -> debug information will be lost
        numb_deleted = 0
        total_loops = 0
        removed_loop_lengths = []
        removed_pet_count = []
        kept_loop_lengths = []
        kept_pet_count = []
        to_remove = []
        for i in range(self.numb_values):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]
            loop_value = self.value_list[i]

            total_loops += 1

            # Both of loop anchor points are not in a wanted peak
            if loop_value == 0 or loop_start < start or loop_end > end or \
                    (not index_array[loop_start] and not index_array[loop_end]):
                removed_pet_count.append(self.value_list[i])
                to_remove.append(i)
                self.value_list[i] = 0
                numb_deleted += 1
                removed_loop_lengths.append(
                    self.end_list[i] - self.start_list[i])
                continue

            kept_pet_count.append(self.value_list[i])
            kept_loop_lengths.append(self.end_list[i] - self.start_list[i])

        self.value_list = np.delete(self.value_list, to_remove)
        self.start_list = np.delete(self.start_list, to_remove)
        self.end_list = np.delete(self.end_list, to_remove)
        self.start_anchor_list = \
            np.delete(self.start_anchor_list, to_remove, -1)
        self.end_anchor_list = np.delete(self.end_anchor_list, to_remove, -1)
        self.numb_values -= len(to_remove)

        log.debug(f'Total loops: {total_loops}')
        log.debug(f"Number of loops removed: {numb_deleted}")
        log.debug(f"Number of loops kept: {total_loops - numb_deleted}")
        log.debug(f'Avg loop length removed: {np.mean(removed_loop_lengths)}')
        log.debug(f'Avg loop length kept: {np.mean(kept_loop_lengths)}')
        log.debug(f'Avg PET count removed: {np.mean(removed_pet_count)}')
        log.debug(f'Avg PET count kept: {np.mean(kept_pet_count)}')
        log.debug(f'Time taken: {time.time() - start_time}\n')

        # Success
        return True

    def temp_func(self):
        log.info("HI")
        pass

    def create_graph(self, loops, num_loops, bin_size, window_size,
                     random=False):
        graph_len = ceil(window_size / bin_size)
        graph = np.zeros((graph_len, graph_len), dtype=np.float64)

        # So comparing a sample with no loops vs. a sample with loops doesn't
        # result in high reproducibility
        if num_loops == 0:
            num_loops = len(loops)

        indexes = None
        if random:
            indexes = np.random.choice(len(loops), num_loops, replace=False)
        log.debug(f"Number of loops to use: {num_loops}")

        num_loops_used = 0
        for i in range(num_loops):

            if indexes:
                loop_index = loops[indexes[i]]
            else:
                loop_index = loops[i]

            value = self.value_list[loop_index]
            start = self.start_list[loop_index]
            end = self.end_list[loop_index]

            start = start % window_size
            end = end % window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            # Get the other side of the graph as well for emd calculation
            graph[bin_start][bin_end] += value
            graph[bin_end][bin_start] += value

            # Also get areas surrounding this loop
            # May not be needed with emd calculation
            # Helps with finding jensen-shannon
            for j in range(bin_start - 1, bin_start + 2):
                if j < 0 or j == graph_len:
                    continue
                for k in range(bin_end - 1, bin_end + 2):
                    if k < 0 or k == graph_len:
                        continue
                    graph[j][k] += value
                    graph[k][j] += value

            num_loops_used += 1

        # log.info(f"Number of loops in {self.sample_name} graph: {num_loops_used}")

        # plt.plot([x for x in range(len(loop_len))], [np.log(x) for x in loop_len])
        # plt.show()

        # log.info(f'Max value in graph: {np.max(graph)}')
        # return graph, total_PET_count / self.total_loop_value
        return graph

    def get_stats(self, graph, o_graph, o_chromLoopData, graph_type, table):

        graph_flat = graph.flatten()
        o_graph_flat = o_graph.flatten()

        max_emd_dist = graph[0].size - 1
        log.debug(f'Graph size: {max_emd_dist + 1}')

        result = {'graph_type': graph_type}

        start_time = time.time()
        j_divergence = jensen_shannon_divergence(graph_flat, o_graph_flat)
        log.debug(f'Shannon time: {time.time() - start_time}')

        # Make j_value range from -1 to 1
        j_value = 2 * (1 - j_divergence) - 1

        # Calculate emd for all rows and columns -> Take weighted average
        emd_distance_list = []
        emd_weight_list = []
        start_time = time.time()
        for k in range(graph[0].size):
            emd_value, emd_weight = emd(graph[k], o_graph[k])
            emd_distance_list.append(emd_value)
            emd_weight_list.append(emd_weight)

            emd_value, emd_weight = emd(graph[:, k], o_graph[:, k])
            emd_distance_list.append(emd_value)
            emd_weight_list.append(emd_weight)
        log.debug(f'EMD time: {time.time() - start_time}')

        '''if np.sum(emd_weight_list) == 0:
            emd_distance = 0
        else:
            emd_distance = np.average(emd_distance_list,
                                      weights=emd_weight_list)'''
        emd_dist = np.mean(emd_distance_list)

        # Scale from -1 to 1
        # Since most values are in bottom half, don't scale linearly
        emd_value = 2 * math.pow(emd_dist - max_emd_dist, 2) / \
                    math.pow(max_emd_dist, 2) - 1

        # linear scale
        # emd_value = 1 - 2 / max_emd * emd_dist

        total_weight = np.max(graph) / self.max_loop_value + \
                       np.max(o_graph) / o_chromLoopData.max_loop_value

        result['rep'] = emd_value * EMD_WEIGHT + j_value * J_WEIGHT
        result['w'] = total_weight
        table.add_row([graph_type, j_value, emd_value, total_weight])

        return result

    def compare_randomly_choose(self, loops, o_loops, o_chrom, bin_size,
                                window_size, inner_table, num_iterations=10):
        start_time = time.time()

        num_loops_to_choose = min(len(loops), len(o_loops))

        value_table = PrettyTable(['graph_type', 'rep', 'weight'])
        results = []
        for _ in range(num_iterations):
            # Create a different graph for each iteration
            graph = self.create_graph(loops, num_loops_to_choose, bin_size,
                                      window_size, True)
            o_graph = o_chrom.create_graph(o_loops, num_loops_to_choose,
                                           bin_size, window_size, True)
            result = self.get_stats(graph, o_graph, o_chrom, 'random',
                                    inner_table)
            results.append(result)

            value_table.add_row([x for x in list(result.values())])

        log_bin.info(value_table)

        # Take average of all iterations
        final_result = {
            'graph_type': 'random',
            'w': sum(x['w'] for x in results) / num_iterations,
            'rep': sum(x['rep'] for x in results) / num_iterations
        }

        log.debug(f'Random time: {time.time() - start_time} ------------------')
        return final_result

    def compare_subsample(self, loops, o_loops, o_chrom, bin_size,
                          window_size, inner_table):
        start_time = time.time()

        num_loops_to_choose = min(len(loops), len(o_loops))

        # Sort the index list that has more to take the top loops
        if len(loops) > num_loops_to_choose:
            loops = loops.tolist()
            loops.sort(key=lambda x: self.value_list[x], reverse=True)
            log.debug(f"Taking the top {num_loops_to_choose} of "
                      f"{self.sample_name}")
        if len(o_loops) > num_loops_to_choose:
            o_loops = o_loops.tolist()
            o_loops.sort(key=lambda x: o_chrom.value_list[x], reverse=True)
            log.debug(f"Taking the top {num_loops_to_choose} of "
                      f"{o_chrom.sample_name}")

        value_table = PrettyTable(['graph_type', 'rep', 'weight'])
        graph = self.create_graph(loops, num_loops_to_choose, bin_size,
                                  window_size)
        o_graph = o_chrom.create_graph(o_loops, num_loops_to_choose,
                                       bin_size, window_size)
        result = self.get_stats(graph, o_graph, o_chrom, 'subsample',
                                inner_table)

        value_table.add_row([x for x in list(result.values())])
        log_bin.info(value_table)
        log.debug(f'Sample time: {time.time() - start_time} ------------------')

        return result

    def compare_only_match(self, graph, o_graph, o_chrom, inner_table):
        start_time = time.time()

        graph_type = ['mod', 'orig']

        # Returns [(mod1, orig1), (mod2, orig2)]
        new_graph_list = get_matching_graphs(graph, o_graph)
        value_table = PrettyTable(['graph_type', 'rep', 'weight'])
        results = []

        # Compare only mod vs. orig
        for i in range(2):
            j = 1 - i

            new_graph = new_graph_list[0][i]
            new_o_graph = new_graph_list[1][j]

            result = self.get_stats(new_graph, new_o_graph, o_chrom,
                                    f'{graph_type[i]}_{graph_type[j]}',
                                    inner_table)
            results.append(result)

            value_table.add_row([x for x in list(result.values())])

        log_bin.info(value_table)
        log.debug(f'Match time: {time.time() - start_time} ------------------')

        # Not sure what to return: max, avg?
        return results[-1]

    def compare(self, o_chrom, window_start, window_end, bin_size, same_seq):
        start_time = time.time()

        log_all.info(f'{self.sample_name} vs. {o_chrom.sample_name} '
                     f'{self.name}:{window_start} - {window_end}')
        window_size = window_end - window_start

        # Get areas removed due to overlapping start/end anchors
        remove_time = time.time()
        to_remove = np.array(
            [self.removed_intervals[0] + o_chrom.removed_intervals[0],
             self.removed_intervals[1] + o_chrom.removed_intervals[1]],
            dtype=np.int32)
        combined_removed = get_removed_area(self.size, to_remove)
        log.debug(f'Removed Area time: {time.time() - remove_time}')

        # Get loop indexes in the window
        loop_time = time.time()
        loops = get_loops(window_start, window_end, combined_removed,
                          self.start_list, self.end_list, self.value_list)
        o_loops = get_loops(window_start, window_end, combined_removed,
                            o_chrom.start_list, o_chrom.end_list,
                            o_chrom.value_list)
        log.debug(f'Get Loop time: {time.time() - loop_time}')

        log.debug(f"Numb of loops in {self.sample_name}: {len(loops)}")
        log.debug(f"Numb of loops in {o_chrom.sample_name}: "
                  f"{len(o_loops)}")

        value_table = PrettyTable(['graph_type', 'rep', 'weight'])
        inner_value_table = \
            PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])

        all_results = []

        # Randomly subsample from larger sample to better match miseq vs. hiseq
        '''result = self.compare_randomly_choose(loops, o_loops, o_chromLoopData,
                                              bin_size, window_size, 
                                              inner_value_table)
        log_bin.info(inner_value_table)
        inner_value_table.clear_rows()
        all_results.append(result)
        value_table.add_row([x for x in list(result.values())])'''

        # Subsample top loops from larger sample to better match miseq vs. hiseq
        sub_result = self.compare_subsample(loops, o_loops, o_chrom, bin_size,
                                            window_size, inner_value_table)
        log_bin.info(inner_value_table)
        inner_value_table.clear_rows()
        all_results.append(sub_result)
        value_table.add_row([x for x in list(sub_result.values())])

        # Make graphs using all loops in the window
        graph = self.create_graph(loops, 0, bin_size, window_size)
        o_graph = \
            o_chrom.create_graph(o_loops, 0, bin_size, window_size)

        # Normal comparison
        norm_time = time.time()
        norm_result = self.get_stats(graph, o_graph, o_chrom, 'norm',
                                     inner_value_table)
        log.debug(f'Norm time: {time.time() - norm_time} ------------------')
        log_bin.info(inner_value_table)
        inner_value_table.clear_rows()
        all_results.append(norm_result)
        value_table.add_row([x for x in list(norm_result.values())])

        # Match graphs only where loops overlap
        '''if not same_seq:  # Don't do this for hiseq-hiseq or miseq-miseq?
            result = self.compare_only_match(graph, o_graph, o_chromLoopData,
                                             inner_value_table)
            log_bin.info(inner_value_table)
            inner_value_table.clear_rows()
            all_results.append(result)
            value_table.add_row([x for x in list(result.values())])'''

        log_bin.info(value_table)
        log.debug(f'Total time: {time.time() - start_time}')
        log_all.info(
            '-----------------------------------------------------------------')

        return sub_result
