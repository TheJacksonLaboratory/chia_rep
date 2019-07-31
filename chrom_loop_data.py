from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt


def jensen_shannon_divergence(p, q, base=2):
    # normalize p, q to probabilities
    p_sum = p.sum()
    q_sum = q.sum()

    if p_sum == 0 and q_sum == 0:
        return -1

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
    # print(p.shape)
    p_array = p.flatten()
    q_array = q.flatten()
    j_divergence = jensen_shannon_divergence(p_array, q_array, base)
    # w_value = wasserstein(p_array, q_array)
    # print(f'Dot product: {dot_product(p, q)}')

    return [j_divergence]# , w_value


class ChromLoopData:

    def __init__(self, chrom_name, chrom_size, bin_size, window_size):
        self.name = chrom_name
        self.size = chrom_size
        self.bin_size = bin_size
        self.window_size = window_size

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
        if loop_length > 1000000:
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

    def create_graph(self, window_index):
        graph_len = ceil(self.window_size / self.bin_size)

        # graph = np.zeros((graph_len, graph_len), dtype=np.uint32)
        graph = np.zeros((graph_len, graph_len), dtype=np.float64)
        # self.graph.fill(1)

        # loop_len_len = ceil((self.largest_loop_len + 1) / self.bin_size)
        # loop_len = np.zeros(loop_len_len, dtype=np.uint16)
        # print(loop_len_len)

        window_start = window_index * self.window_size
        window_end = (window_index + 1) * self.window_size

        for i in range(self.numb_values):
            start = self.start_list[i]
            end = self.end_list[i]
            value = self.value_list[i]

            if start < window_start or end > window_end:
                continue

            start = start % self.window_size
            end = end % self.window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / self.bin_size)
            bin_end = int(end / self.bin_size)

            graph[bin_start][bin_end] += value * value
            graph[bin_end][bin_start] += value * value

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

        start_mean_values = bedGraph_chrom.get_exact_mean(
            self.start_interval_list[0], self.start_interval_list[1])
        end_mean_values = bedGraph_chrom.get_exact_mean(
            self.end_interval_list[0], self.end_interval_list[1])

        for i in range(self.numb_values):

            mean_value = np.mean([start_mean_values[i], end_mean_values[i]])
            change = mean_value / bedGraph_chrom.avg_chrom_value
            # print(f'Change: {change}')

            self.value_list[i] *= change

        bedGraph_chrom.free_index_list()

    def compare(self, o_chromLoopData, window_index):
        o_graph = o_chromLoopData.create_graph(window_index)
        graph = self.create_graph(window_index)
        return graph_similarity(graph, o_graph)
