from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt
import heapq

MAX_LOOP_LEN = 70000


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
    # print(np.info(p))
    # print(np.info(q))
    # print(p[p > 0.00000001])

    m = (p + q) / 2
    return sp.entropy(p, m, base=base) / 2 + sp.entropy(q, m, base=base) / 2


def wasserstein(p, q):
    return sp.wasserstein_distance(p, q)


def dot_product(p, q):
    return np.sum(np.dot(p, q)) / (np.sum(p) * np.sum(q))


def graph_similarity(p, q, base=2):
    pass


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
        self.start_anchor_list = np.asarray(self.start_anchor_list, dtype=np.int32)
        self.end_anchor_list = np.asarray(self.end_anchor_list, dtype=np.int32)

    def find_loop_anchor_points(self, bedGraph):

        print(f'Finding anchor points for {self.sample_name}')

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

            # peak_value = max(start_list_peaks[i], end_list_peaks)

            if not loop_start < loop_end:
                self.value_list[i] = 0
                self.start_list[i] = -1
                self.end_list[i] = -1
                continue

            loop_length = int(loop_end) - int(loop_start)

            # loops over some threshold are considered noise
            if loop_length > MAX_LOOP_LEN:
                self.value_list[i] = 0

        bedGraph.free_chrom_data('chr1')

    def temp_func(self):
        print("HI")
        pass

    def get_removed_area(self):
        removed_area = np.full(self.size, False, dtype=bool)
        for i in range(self.numb_values):
            if self.start_list[i] == -1 or self.end_list[i] == -1:
                start = self.start_anchor_list[0][i]
                end = self.end_anchor_list[1][i]
                removed_area[start:end] = True

        return removed_area

    def create_graph(self, window_start, window_end, bin_size, removed_area):
        window_size = window_end - window_start
        graph_len = ceil(window_size / bin_size)

        graph = np.zeros((graph_len, graph_len), dtype=np.float64)

        # loop_len_len = ceil((self.largest_loop_len + 1) / self.bin_size)
        # loop_len = np.zeros(loop_len_len, dtype=np.uint16)
        # print(loop_len_len)

        num_loops_used = 0
        for i in range(self.numb_values):
            start = self.start_list[i]
            end = self.end_list[i]
            value = self.value_list[i]

            if removed_area[start] or removed_area[end]:
                continue

            if start < window_start or end > window_end or value == 0:
                continue

            num_loops_used += 1

            start = start % window_size
            end = end % window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            if self.sample_name == 'LHH0058H' and value == 36:
                print(self.start_list[i], self.end_list[i])
                print(graph[bin_start - 1][bin_end])
                print(graph[bin_start][bin_end])
                print()

            graph[bin_start][bin_end] += value
            for j in range(bin_start - 1, bin_start + 2, 1):
                if j < 0 or j == graph_len:
                    continue
                for k in range(bin_end - 1, bin_end + 2, 1):
                    if k < 0 or k == graph_len:
                        continue
                    graph[j][k] += value

            if self.sample_name == 'LHH0058H' and value == 36:
                print(self.start_list[i], self.end_list[i])
                print(graph[bin_start - 1][bin_end])
                print(graph[bin_start][bin_end])
                print()

        print(f"Number of loops in {self.sample_name} graph: {num_loops_used}")

        # plt.plot([x for x in range(len(loop_len))], [np.log(x) for x in loop_len])
        # plt.show()

        # print(f'Max value in graph: {np.max(graph)}')
        return graph

    def adjust_with_bedGraph(self, bedGraph_chrom):
        bedGraph_chrom.load_index_array()

        for i in range(2):
            self.start_anchor_list[i] = \
                np.asarray(self.start_anchor_list[i], dtype=np.int32)
            self.end_anchor_list[i] = \
                np.asarray(self.end_anchor_list[i], dtype=np.int32)

        start_mean_values = bedGraph_chrom.get_max(
            self.start_anchor_list[0], self.start_anchor_list[1])
        end_mean_values = bedGraph_chrom.get_max(
            self.end_anchor_list[0], self.end_anchor_list[1])

        for i in range(self.numb_values):
            # peak_value = (start_mean_values[i] + end_mean_values[i]) / 2
            peak_value = max(start_mean_values[i], end_mean_values[i])
            change = peak_value
            # print(f'Change: {change}')

            self.value_list[i] *= change

        bedGraph_chrom.free_index_list()

    def filter_with_bedGraph(self, test_cases, bedGraph_mean_list):

        start_time = time.time()
        assert len(test_cases[0]) == len(bedGraph_mean_list)
        num_test_cases = len(test_cases[0])

        '''plt.close()
        plt.hist([np.log10(x) for x in bedGraph_mean_list], bins=20)
        plt.xlabel('log10(Max value in interval)')
        plt.ylabel('Frequency')
        plt.title(self.sample_name)
        plt.savefig(f'{self.sample_name}.png')'''

        peaks = []
        for i in range(num_test_cases):
            peaks.append({
                'start': test_cases[0][i],
                'end': test_cases[1][i],
                'value': bedGraph_mean_list[i]
            })

        num_wanted_peaks = int(len(peaks) / 200)
        wanted_peaks = heapq.nlargest(num_wanted_peaks, peaks, key=lambda k: k['value'])
        min_peak_value = min(wanted_peaks, key=lambda k: k['value'])['value']

        start = test_cases[0][0]
        end = test_cases[1][-1]

        # eventually want to weight each loop with its corresponding peak
        # index_array = np.full(end - start, False, dtype=bool)
        index_array = np.full(end - start, -1, dtype=np.float64)

        print(num_wanted_peaks, min_peak_value)
        for i in range(num_wanted_peaks):
            peak_start = wanted_peaks[i]['start'] - start
            peak_end = wanted_peaks[i]['end'] - start
            index_array[peak_start:peak_end] = wanted_peaks[i]['value']

        numb_deleted = 0
        total_loops_inside_window = 0
        removed_loop_lengths = []
        removed_pet_count = []
        kept_loop_lengths = []
        kept_pet_count = []
        for i in range(self.numb_values):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]

            if loop_start == -1 or loop_end == -1:
                continue

            if loop_start < start or loop_end > end:
                continue

            total_loops_inside_window += 1
            if index_array[loop_start - start] == -1 and \
                    index_array[loop_end - start] == -1:
                removed_pet_count.append(self.value_list[i])
                self.value_list[i] = 0
                numb_deleted += 1
                removed_loop_lengths.append(
                    self.end_list[i] - self.start_list[i])
                continue

            kept_pet_count.append(self.value_list[i])
            kept_loop_lengths.append(self.end_list[i] - self.start_list[i])
            '''else:
                loop_strength = max(index_array[loop_start - start],
                                    index_array[loop_end - start])
                self.value_list[i] *= loop_strength'''

        print(f'Time taken: {time.time() - start_time}')
        print(f"Number of loops removed: {numb_deleted}")
        print(f'Total loops: {total_loops_inside_window}')
        print(f'Average loop length removed: {np.mean(removed_loop_lengths)}')
        print(f'Average PET count removed: {np.mean(removed_pet_count)}')
        print(f'Average loop length kept: {np.mean(kept_loop_lengths)}')
        print(f'Average PET count kept: {np.mean(kept_pet_count)}')
        print()

    def compare(self, o_chromLoopData, window_start, window_end, bin_size):
        removed1 = self.get_removed_area()
        removed2 = o_chromLoopData.get_removed_area()
        combined_removed = removed1 | removed2

        graph = self.create_graph(window_start, window_end, bin_size,
                                  combined_removed)
        o_graph = o_chromLoopData.create_graph(window_start, window_end,
                                               bin_size, combined_removed)

        graph = graph
        o_graph = o_graph
        graph_len = graph.size
        if self.sample_name == 'LHH0058H' and o_chromLoopData.sample_name == 'LHH0061H':
            average_error = 0
            for i in range(graph[0].size):
                for j in range(graph[0].size):
                    average_error += abs(graph[i][j] - o_graph[i][j])
                    if abs(graph[i][j] - o_graph[i][j]) > 50:
                        print(i * 3000 + 36000000, j * 3000 + 36000000, graph[i][j],
                              o_graph[i][j])

            print(f'Average error: {average_error / (graph_len * graph_len)}')
        graph_flat = graph.flatten()
        o_graph_flat = o_graph.flatten()
        j_divergence = jensen_shannon_divergence(graph_flat, o_graph_flat)
        # w_value = wasserstein(p_array, q_array)
        # print(f'Dot product: {dot_product(p, q)}')

        return [j_divergence]  # , w_value
