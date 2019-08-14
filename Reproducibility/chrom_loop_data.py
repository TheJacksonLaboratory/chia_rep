from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt
from copy import deepcopy
import logging
from .reproducibility_util import *

from prettytable import PrettyTable

log = logging.getLogger(__name__.split('.')[-1])

MAX_LOOP_LEN = 1000000

VERSION = 0


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
    # normalize p, q to probabilities
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


def normalize(p):
    graph_sum = p.sum()
    if graph_sum == 0:
        return p
    else:
        return p / graph_sum


def get_graphs(p, q):
    assert p.size == q.size
    q1 = deepcopy(q)
    p1 = deepcopy(p)
    match_graphs(p, q, p1, q1, p[0].size)

    return [[p1, p], [q1, q]]


class ChromLoopData:

    def __init__(self, chrom_name, chrom_size, sample_name):
        self.name = chrom_name
        self.size = int(chrom_size)
        self.sample_name = sample_name

        self.start_anchor_list = [[], []]
        self.end_anchor_list = [[], []]
        self.start_list = None
        self.end_list = None
        self.value_list = []
        self.numb_values = 0
        self.removed_intervals = [[], []]

        self.total_loop_value = 0
        self.max_loop_value = 0

        self.largest_loop_len = -1

    def add_loop(self, loop_start1, loop_start2, loop_end1, loop_end2,
                 loop_value):

        self.start_anchor_list[0].append(loop_start1)
        self.start_anchor_list[1].append(loop_start2)
        self.end_anchor_list[0].append(loop_end1)
        self.end_anchor_list[1].append(loop_end2)

        self.value_list.append(int(loop_value))

        self.numb_values += 1

    def finish_init(self, bedgraph, peak_file_path, peak_percent_kept):
        self.value_list = np.asarray(self.value_list, dtype=np.int32)
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

    def find_loop_anchor_points(self, bedGraph):

        log.info(f'Finding anchor points for {self.sample_name}\'s {self.name}')

        '''self.start_list = np.array(
            (self.start_anchor_list[0] + self.start_anchor_list[1]) / 2,
            dtype=np.int32)
        self.end_list = np.array(
            (self.end_anchor_list[0] + self.end_anchor_list[1]) / 2,
            dtype=np.int32)

        return'''

        bedGraph.load_chrom_data(self.name)
        self.start_list = bedGraph.stats(start_list=self.start_anchor_list[0],
                                         end_list=self.start_anchor_list[1],
                                         chrom_name=self.name, stat='max_index')

        self.end_list = bedGraph.stats(start_list=self.end_anchor_list[0],
                                       end_list=self.end_anchor_list[1],
                                       chrom_name=self.name, stat='max_index')

        start_list_peaks = bedGraph.stats(start_list=self.start_anchor_list[0],
                                          end_list=self.start_anchor_list[1],
                                          chrom_name=self.name, stat='max')

        end_list_peaks = bedGraph.stats(start_list=self.end_anchor_list[0],
                                        end_list=self.end_anchor_list[1],
                                        chrom_name=self.name, stat='max')

        for i in range(self.numb_values):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]

            # Weight each loop based on a corresponding bedgraph peak
            peak_value = max(start_list_peaks[i], end_list_peaks[i])
            self.value_list[i] *= peak_value

            # Remove overlapping anchors
            if not loop_start < loop_end:
                self.value_list[i] = 0
                self.start_list[i] = -1
                self.end_list[i] = -1
                self.removed_intervals[0].append(self.start_anchor_list[0][i])
                self.removed_intervals[1].append(self.end_anchor_list[1][i])
                continue

            self.total_loop_value += self.value_list[i]

            if self.max_loop_value < self.value_list[i]:
                self.max_loop_value = self.value_list[i]

            loop_length = int(loop_end) - int(loop_start)

            # Remove loops over a given threshold
            if loop_length > MAX_LOOP_LEN:
                self.value_list[i] = 0

        bedGraph.free_chrom_data(self.name)

    def filter_with_peaks(self, peak_file, peak_percent_kept, start=0,
                          end=None):
        if end is None:
            end = self.size

        start_time = time.time()
        peaks = []
        peak_length = []
        with open(peak_file) as in_file:
            log.info(f"Reading {peak_file} ...")
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
        log.info(f"Median length of all peaks: {median_peak_len}")
        log.info(f'Number of peaks: {len(peaks)}')

        num_wanted_peaks = int(len(peaks) * peak_percent_kept)

        if num_wanted_peaks == 0:
            log.error(f"Not enough peaks were found for {self.name} in "
                      f"{peak_file}. Only {peak_percent_kept * 100}% "
                      f"of {len(peaks)} were kept")
            return False

        wanted_peaks = peaks[:num_wanted_peaks]
        peak_length = peak_length[:num_wanted_peaks]
        min_peak_value = wanted_peaks[-1][2]

        log.info(f"Median length of wanted peaks: {np.median(peak_length)}")
        log.info(f"Num wanted peaks: {num_wanted_peaks}")
        log.info(f'Min peak value: {min_peak_value}')

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

            if loop_start == -1 or loop_end == -1:
                continue

            if loop_start < start or loop_end > end:
                continue

            total_loops += 1

            if index_array[loop_start - start] or index_array[loop_end - start]:
                kept_pet_count.append(self.value_list[i])
                kept_loop_lengths.append(self.end_list[i] - self.start_list[i])
                continue

            removed_pet_count.append(self.value_list[i])
            to_remove.append(i)
            self.value_list[i] = 0
            numb_deleted += 1
            removed_loop_lengths.append(self.end_list[i] - self.start_list[i])

        self.value_list = np.delete(self.value_list, to_remove)
        self.start_list = np.delete(self.start_list, to_remove)
        self.end_list = np.delete(self.end_list, to_remove)
        self.start_anchor_list = np.delete(self.start_anchor_list, to_remove,
                                           -1)
        self.end_anchor_list = np.delete(self.end_anchor_list, to_remove, -1)
        self.numb_values -= len(to_remove)

        log.info(f'Total loops: {total_loops}')
        log.info(f"Number of loops removed: {numb_deleted}")
        log.info(f"Number of loops kept: {total_loops - numb_deleted}")
        log.info(f'Avg loop length removed: {np.mean(removed_loop_lengths)}')
        log.info(f'Avg loop length kept: {np.mean(kept_loop_lengths)}')
        log.info(f'Avg PET count removed: {np.mean(removed_pet_count)}')
        log.info(f'Avg PET count kept: {np.mean(kept_pet_count)}')
        log.debug(f'Time taken: {time.time() - start_time}\n')

        # Success
        return True

    def temp_func(self):
        log.info("HI")
        pass

    def create_graph(self, loops, num_loops, bin_size, window_size):
        graph_len = ceil(window_size / bin_size)
        graph = np.zeros((graph_len, graph_len), dtype=np.float64)

        if num_loops == 0:
            num_loops = len(loops)

        num_loops_used = 0
        indexes = np.random.choice(len(loops), num_loops, replace=False)
        for i in indexes:
            loop_index = loops[i]

            value = self.value_list[loop_index]
            start = self.start_list[loop_index]
            end = self.end_list[loop_index]

            num_loops_used += 1

            start = start % window_size
            end = end % window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            graph[bin_start][bin_end] += value
            graph[bin_end][bin_start] += value
            for j in range(bin_start - 1, bin_start + 2):
                if j < 0 or j == graph_len:
                    continue
                for k in range(bin_end - 1, bin_end + 2):
                    if k < 0 or k == graph_len:
                        continue
                    graph[j][k] += value
                    graph[k][j] += value

        # log.info(f"Number of loops in {self.sample_name} graph: {num_loops_used}")

        # plt.plot([x for x in range(len(loop_len))], [np.log(x) for x in loop_len])
        # plt.show()

        # log.info(f'Max value in graph: {np.max(graph)}')
        # return graph, total_PET_count / self.total_loop_value
        return graph

    def get_loops(self, window_start, window_end, removed_area):
        loops = []
        for i in range(self.numb_values):
            value = self.value_list[i]
            start = self.start_list[i]
            end = self.end_list[i]

            if start < window_start or end > window_end or value == 0:
                continue

            if removed_area[start] or removed_area[end]:
                continue

            if start > window_end:
                break

            loops.append(i)

        return loops

    def get_stats(self, graph, o_graph, o_chromLoopData, graph_type):
        start_time = time.time()

        graph_flat = graph.flatten()
        o_graph_flat = o_graph.flatten()

        result = {'graph_type': graph_type}

        j_divergence = jensen_shannon_divergence(graph_flat, o_graph_flat)
        log.debug(f'Shannon time: {time.time() - start_time}')

        # Make j_value range from -1 to 1
        result['j-s_value'] = 2 * (1 - j_divergence) - 1

        emd_distance_list = []
        emd_weight_list = []
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
        emd_distance = np.mean(emd_distance_list)
        result['earth_mover'] = emd_distance

        total_weight = np.max(graph) / self.max_loop_value + \
                       np.max(o_graph) / o_chromLoopData.max_loop_value

        result['w'] = total_weight

        return result

    def compare_randomly_choose(self, loops, o_loops, o_chromLoopData, bin_size,
                                window_size, num_iterations=10):
        start_time = time.time()

        num_loops_to_choose = min(len(loops), len(o_loops))

        value_table = PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])
        results = []
        for _ in range(num_iterations):
            graph = self.create_graph(loops, num_loops_to_choose, bin_size,
                                      window_size)
            o_graph = o_chromLoopData.create_graph(o_loops, num_loops_to_choose,
                                                   bin_size, window_size)
            result = self.get_stats(graph, o_graph, o_chromLoopData, 'random')
            results.append(result)

            value_table.add_row([x for x in list(result.values())])

        log.info(value_table)

        final_result = {'graph_type': 'random'}
        for result in results:
            for key in result:
                if key == 'graph_type':
                    continue

                if key not in final_result:
                    final_result[key] = 0
                final_result[key] += result[key]

        for key in final_result:
            if key == 'graph_type':
                continue
            final_result[key] /= len(results)

        log.debug(f'Random time: {time.time() - start_time}')
        return final_result

    def compare_only_match(self, graph, o_graph, o_chromLoopData):
        start_time = time.time()

        graph_type = ['mod', 'orig']
        new_graph_list = get_graphs(graph, o_graph)
        value_table = PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])
        results = []

        for i in range(2):
            # Only compare originals for now
            # if i == 0:
            #    continue

            for j in range(2):
                # Don't compare mod vs. mod
                if i == 0 and j == 0:
                    continue

                # Only compare originals for now
                # if j == 0:
                #   continue

                new_graph = new_graph_list[0][i]
                new_o_graph = new_graph_list[1][j]

                result = self.get_stats(new_graph, new_o_graph, o_chromLoopData,
                                        f'{graph_type[i]}_{graph_type[j]}')
                results.append(result)

                value_table.add_row([x for x in list(result.values())])

        log.info(value_table)
        log.debug(f'Match time: {time.time() - start_time}')
        return results[-1]

    @staticmethod
    def get_removed_area(chrom_size, to_remove):
        removed_area = np.empty(chrom_size, dtype=bool)
        for i in range(len(to_remove[0])):
            start = to_remove[0][i]
            end = to_remove[1][i]
            removed_area[start:end] = True

        return removed_area

    def compare(self, o_chromLoopData, window_start, window_end, bin_size,
                same_seq):
        start_time = time.time()

        log.info(f'{self.sample_name} Window: {window_start} - {window_end}')
        window_size = window_end - window_start

        # Get areas removed due to overlapping start/end anchors
        remove_time = time.time()
        to_remove = np.array(
            [self.removed_intervals[0] + o_chromLoopData.removed_intervals[0],
             self.removed_intervals[1] + o_chromLoopData.removed_intervals[1]],
            dtype=np.int32)
        combined_removed = ChromLoopData.get_removed_area(self.size, to_remove)
        log.debug(f'Removed Area time: {time.time() - remove_time}')

        # Get loops in the window
        loop_time = time.time()
        loops = self.get_loops(window_start, window_end, combined_removed)
        o_loops = o_chromLoopData.get_loops(window_start, window_end,
                                            combined_removed)
        log.debug(f'Get Loop time: {time.time() - loop_time}')

        log.debug(f"Numb of loops in {self.sample_name}: {len(loops)}")
        log.debug(f"Numb of loops in {o_chromLoopData.sample_name}: "
                  f"{len(o_loops)}")

        value_table = PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])
        all_results = []

        # Randomly subsample from larger sample to better match miseq vs. hiseq
        '''result = self.compare_randomly_choose(loops, o_loops, o_chromLoopData,
                                              bin_size, window_size)
        all_results.append(result)
        value_table.add_row([x for x in list(result.values())])'''

        # Normal comparison
        norm_time = time.time()
        graph = self.create_graph(loops, 0, bin_size, window_size)
        o_graph = \
            o_chromLoopData.create_graph(o_loops, 0, bin_size, window_size)
        result = self.get_stats(graph, o_graph, o_chromLoopData, 'norm')
        all_results.append(result)
        value_table.add_row([x for x in list(result.values())])
        log.debug(f'Norm time: {time.time() - norm_time}')

        # Match graphs only where loops overlap
        if not same_seq:
            pass
        result = self.compare_only_match(graph, o_graph, o_chromLoopData)
        all_results.append(result)
        value_table.add_row([x for x in list(result.values())])

        log.info(value_table)
        log.debug(f'Total time: {time.time() - start_time}')
        log.info(
            '---------------------------------------------------------------\n')

        return all_results
