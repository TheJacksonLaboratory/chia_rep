from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt
from copy import deepcopy
import logging as log

from prettytable import PrettyTable

MAX_LOOP_LEN = 70000

VERSION = 0

PERCENT_PEAK_KEPT = 0.20


def emd(p, q):
    assert p.size == q.size
    p_sum = p.sum()
    q_sum = q.sum()

    if p_sum == 0 and q_sum == 0:
        return 0

    if p_sum == 0:
        return p.size

    if q_sum == 0:
        return q.size

    p = p / p_sum
    q = q / q_sum
    to_move = 0
    cost = 0
    for i in range(p.size):
        to_move += p[i] - q[i]
        cost += abs(to_move)
    return cost


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
    for i in range(p[0].size):
        for j in range(p[0].size):
            if p[i][j] == 0:
                q1[i][j] = 0

    p1 = deepcopy(p)
    for i in range(p[0].size):
        for j in range(p[0].size):
            if q[i][j] == 0:
                p1[i][j] = 0

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

    def finish_init(self):
        self.value_list = np.asarray(self.value_list, dtype=np.int32)
        self.start_anchor_list = np.asarray(self.start_anchor_list,
                                            dtype=np.int32)
        self.end_anchor_list = np.asarray(self.end_anchor_list, dtype=np.int32)

    def find_loop_anchor_points(self, bedGraph):

        log.info(f'Finding anchor points for {self.sample_name}')

        '''self.start_list = np.array(
            (self.start_anchor_list[0] + self.start_anchor_list[1]) / 2,
            dtype=np.int32)
        self.end_list = np.array(
            (self.end_anchor_list[0] + self.end_anchor_list[1]) / 2,
            dtype=np.int32)

        return'''

        bedGraph.load_chrom_data('chr1')
        self.start_list = bedGraph.stats(start_list=self.start_anchor_list[0],
                                         end_list=self.start_anchor_list[1],
                                         chrom_name='chr1', stat='max_index')

        self.end_list = bedGraph.stats(start_list=self.end_anchor_list[0],
                                       end_list=self.end_anchor_list[1],
                                       chrom_name='chr1', stat='max_index')

        start_list_peaks = bedGraph.stats(start_list=self.start_anchor_list[0],
                                          end_list=self.start_anchor_list[1],
                                          chrom_name='chr1', stat='max')

        end_list_peaks = bedGraph.stats(start_list=self.end_anchor_list[0],
                                        end_list=self.end_anchor_list[1],
                                        chrom_name='chr1', stat='max')

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

        bedGraph.free_chrom_data('chr1')

    def filter_with_bedGraph(self, peak_file, start=0, end=None):
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

        num_wanted_peaks = int(len(peaks) * PERCENT_PEAK_KEPT)
        assert len(peaks) >= num_wanted_peaks

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
            if not index_array[loop_start - start] and \
                    not index_array[loop_end - start]:
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
        self.start_anchor_list = np.delete(self.start_anchor_list, to_remove, -1)
        self.end_anchor_list = np.delete(self.end_anchor_list, to_remove, -1)
        self.numb_values -= len(to_remove)

        log.info(f'Total loops: {total_loops}')
        log.info(f"Number of loops removed: {numb_deleted}")
        log.info(f"Number of loops kept: {total_loops - numb_deleted}")
        log.info(f'Avg loop length removed: {np.mean(removed_loop_lengths)}')
        log.info(f'Avg loop length kept: {np.mean(kept_loop_lengths)}')
        log.info(f'Avg PET count removed: {np.mean(removed_pet_count)}')
        log.info(f'Avg PET count kept: {np.mean(kept_pet_count)}')
        log.info(f'Time taken: {time.time() - start_time}\n')

    def temp_func(self):
        log.info("HI")
        pass

    @staticmethod
    def get_removed_area(chrom_size, to_remove):
        removed_area = np.full(chrom_size, False, dtype=bool)
        for i in range(len(to_remove[0])):
            start = to_remove[0][i]
            end = to_remove[1][i]
            removed_area[start:end] = True

        return removed_area

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

    def compare(self, o_chromLoopData, window_start, window_end, bin_size):

        log.info(f'Window: {window_start} - {window_end}')

        window_size = window_end - window_start

        start_time = time.time()
        to_remove = \
            [self.removed_intervals[0] + o_chromLoopData.removed_intervals[0],
             self.removed_intervals[1] + o_chromLoopData.removed_intervals[1]]
        combined_removed = ChromLoopData.get_removed_area(self.size, to_remove)

        loops = self.get_loops(window_start, window_end, combined_removed)
        o_loops = o_chromLoopData.get_loops(window_start, window_end,
                                            combined_removed)

        num_loops = min(len(loops), len(o_loops))
        log.info(f"Numb of loops in {self.sample_name}: {len(loops)}")
        log.info(f"Numb of loops in {o_chromLoopData.sample_name}: "
                 f"{len(o_loops)}")
        log.info(f"Numb of loops to use: {num_loops}")
        log.info(f'Prep time: {time.time() - start_time}')

        value_table = PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])
        all_results = []
        for _ in range(10):
            orig_graph = self.create_graph(loops, num_loops, bin_size, window_size)
            orig_o_graph = o_chromLoopData.create_graph(o_loops, num_loops,
                                                        bin_size, window_size)

            graph_type = ['mod', 'orig']
            new_graph_list = get_graphs(orig_graph, orig_o_graph)
            for i in range(2):
                for j in range(2):
                    # Don't compare mod vs. mod
                    if i == 0 and j == 0:
                        continue

                    # Only compare originals for now
                    if i == 0 or j == 0:
                        continue

                    new_graph = new_graph_list[0][i]
                    new_o_graph = new_graph_list[1][j]
                    graph_flat = new_graph.flatten()
                    o_graph_flat = new_o_graph.flatten()

                    results = {}

                    j_divergence = jensen_shannon_divergence(graph_flat,
                                                             o_graph_flat)
                    # Make j_value range from -1 to 1
                    j_value = 2 * (1 - j_divergence) - 1
                    results['j-s_value'] = j_value

                    emd_distance_list = []
                    for k in range(new_graph[0].size):
                        norm_graph = normalize(new_graph[k])
                        o_norm_graph = normalize(new_o_graph[k])
                        emd_distance_list.append(emd(norm_graph, o_norm_graph))

                        norm_graph = normalize(new_graph[:, k])
                        o_norm_graph = normalize(new_o_graph[:, k])
                        emd_distance_list.append(emd(norm_graph, o_norm_graph))

                    emd_distance = np.mean(emd_distance_list)
                    results['earth_mover'] = emd_distance

                    total_weight = np.max(new_graph) / self.max_loop_value + \
                                np.max(new_o_graph) / o_chromLoopData.max_loop_value
                    # total_weight = max(weight, o_weight)
                    if total_weight == 0:
                        pass
                        # total_weight = 100000
                    results['w'] = total_weight
                    results['type'] = graph_type[i] + '_' + graph_type[j]

                    value_table.add_row([graph_type[i] + '_' + graph_type[j],
                                         round(j_value, 6), round(emd_distance, 6),
                                         round(total_weight, 6)])
                    all_results.append(results)

        combined_result = {'type': 'orig_orig'}
        for result in all_results:
            for key in result:
                if key == 'type':
                    continue

                if key not in combined_result:
                    combined_result[key] = 0
                combined_result[key] += result[key]

        for key in combined_result:
            if key == 'type':
                continue
            combined_result[key] /= len(all_results)
        value_table.add_row([combined_result[key] for key in combined_result])

        log.info(value_table)
        log.info(f'Total time: {time.time() - start_time}')
        log.info(
            '-----------------------------------------------------------------')
        return [combined_result]
