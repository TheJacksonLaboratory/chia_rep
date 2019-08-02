from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt

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
    '''for i in range(p[0].size):
        for j in range(p[0].size):
            if ((p[i][j] == 0 and q[i][j] != 0) or
                    (q[i][j] == 0 and p[i][j] != 0)) and \
                    abs(p[i][j] - q[i][j]) > 100:
                print(i * 1000 + 36000000, j * 1000 + 36000000, p[i][j],
                      q[i][j])'''
    p_array = p.flatten()
    q_array = q.flatten()
    j_divergence = jensen_shannon_divergence(p_array, q_array, base)
    # w_value = wasserstein(p_array, q_array)
    # print(f'Dot product: {dot_product(p, q)}')

    return [j_divergence]  # , w_value


class ChromLoopData:

    def __init__(self, chrom_name, chrom_size, sample_name, bedGraph):
        self.name = chrom_name
        self.size = int(chrom_size)
        self.sample_name = sample_name
        self.bedGraph = bedGraph

        self.start_interval_list = [[], []]
        self.end_interval_list = [[], []]
        self.start_list = []
        self.end_list = []
        self.value_list = []
        self.numb_values = 0

        self.largest_loop_len = -1

    def add_loop(self, loop_start1, loop_start2, loop_end1, loop_end2,
                 loop_value):

        loop_start = (loop_start2 + loop_start1) / 2
        loop_end = (loop_end1 + loop_end2) / 2

        self.start_interval_list[0].append(loop_start1)
        self.start_interval_list[1].append(loop_start2)
        self.end_interval_list[0].append(loop_end1)
        self.end_interval_list[1].append(loop_end2)

        loop_length = int(loop_end) - int(loop_start)

        # loops over 1Mb are considered noise
        if loop_length > MAX_LOOP_LEN:
            return

        if int(loop_value) < 4:
            pass
            # return

        if int(loop_end) - int(loop_start) > self.largest_loop_len:
            self.largest_loop_len = int(loop_end) - int(loop_start)

        self.start_list.append(int(loop_start))
        self.end_list.append(int(loop_end))
        self.value_list.append(int(loop_value))

        self.numb_values += 1

    def create_graph(self, window_start, window_end, bin_size):
        window_size = window_end - window_start
        graph_len = ceil(window_size / bin_size)

        graph = np.zeros((graph_len, graph_len), dtype=np.float64)

        # loop_len_len = ceil((self.largest_loop_len + 1) / self.bin_size)
        # loop_len = np.zeros(loop_len_len, dtype=np.uint16)
        # print(loop_len_len)

        for i in range(self.numb_values):
            start = self.start_list[i]
            end = self.end_list[i]
            value = self.value_list[i]

            if start < window_start or end > window_end or value == 0:
                continue

            start = start % window_size
            end = end % window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            value = value * value
            graph[bin_start][bin_end] += value
            for j in range(bin_start - 1, bin_start + 2, 1):
                if j < 0 or j == graph_len:
                    continue
                for k in range(bin_end - 1, bin_end + 2, 1):
                    if k < 0 or k == graph_len:
                        continue
                    graph[j][k] += value

        # plt.plot([x for x in range(len(loop_len))], [np.log(x) for x in loop_len])
        # plt.show()

        # print(f'Max value in graph: {np.max(graph)}')
        return graph

    def adjust_with_bedGraph(self, bedGraph_chrom):
        bedGraph_chrom.load_index_array()

        for i in range(2):
            self.start_interval_list[i] = \
                np.asarray(self.start_interval_list[i], dtype=np.int32)
            self.end_interval_list[i] = \
                np.asarray(self.end_interval_list[i], dtype=np.int32)

        start_mean_values = bedGraph_chrom.get_max(
            self.start_interval_list[0], self.start_interval_list[1])
        end_mean_values = bedGraph_chrom.get_max(
            self.end_interval_list[0], self.end_interval_list[1])

        for i in range(self.numb_values):

            # peak_value = (start_mean_values[i] + end_mean_values[i]) / 2
            peak_value = max(start_mean_values[i], end_mean_values[i])
            change = peak_value
            # print(f'Change: {change}')

            self.value_list[i] *= change

        bedGraph_chrom.free_index_list()

    def filter_with_bedGraph(self, test_cases, bedGraph_mean_list):
        assert len(test_cases[0]) == len(bedGraph_mean_list)
        num_test_cases = len(test_cases[0])

        plt.close()
        plt.hist([np.log10(x) for x in bedGraph_mean_list], bins=20)
        plt.xlabel('log10(Max value in interval)')
        plt.ylabel('Frequency')
        plt.title(self.sample_name)
        plt.savefig(f'{self.sample_name}.png')

        peaks = []
        for i in range(num_test_cases):
            peaks.append({
                'start': test_cases[0][i],
                'end': test_cases[1][i],
                'value': bedGraph_mean_list[i]
            })

        peaks.sort(key=lambda k: k['value'], reverse=True)

        start = test_cases[0][0]
        end = test_cases[1][-1]

        # eventually want to weight each loop with its corresponding peak
        # index_array = np.full(end - start, False, dtype=bool)
        index_array = np.full(end - start, -1, dtype=np.float64)

        top_percentile = int(len(peaks) / 50)
        print(top_percentile, peaks[top_percentile])
        for i in range(top_percentile):
            if peaks[i]['value'] < 50:
                break
            peak_start = peaks[i]['start'] - start
            peak_end = peaks[i]['end'] - start
            index_array[peak_start:peak_end] = peaks[i]['value']

        numb_deleted = 0
        total_loops_inside_window = 0
        removed_loop_lengths = []
        removed_pet_count = []
        kept_loop_lengths = []
        kept_pet_count = []
        for i in range(self.numb_values):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]

            if loop_start < start or loop_end > end:
                continue

            total_loops_inside_window += 1
            if index_array[loop_start - start] == -1 and \
                    index_array[loop_end - start] == -1:
                removed_pet_count.append(self.value_list[i])
                self.value_list[i] = 0
                numb_deleted += 1
                removed_loop_lengths.append(self.end_list[i] - self.start_list[i])
                continue
            kept_pet_count.append(self.value_list[i])
            kept_loop_lengths.append(self.end_list[i] - self.start_list[i])
            '''else:
                loop_strength = max(index_array[loop_start - start],
                                    index_array[loop_end - start])
                self.value_list[i] *= loop_strength'''

        print(f"Number of loops removed: {numb_deleted}")
        print(f'Total loops: {total_loops_inside_window}')
        print(f'Average loop length removed: {np.mean(removed_loop_lengths)}')
        print(f'Average PET count removed: {np.mean(removed_pet_count)}')
        print(f'Average loop length kept: {np.mean(kept_loop_lengths)}')
        print(f'Average PET count kept: {np.mean(kept_pet_count)}')

    def compare(self, o_chromLoopData, window_start, window_end, bin_size):
        graph = self.create_graph(window_start, window_end, bin_size)
        o_graph = o_chromLoopData.create_graph(window_start, window_end,
                                               bin_size)
        return graph_similarity(graph, o_graph)
