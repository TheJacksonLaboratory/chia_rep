from math import ceil
import numpy as np
import scipy.stats as sp
import time
import matplotlib.pyplot as plt


def jensen_shannon_divergence(p, q, base=2):
    # normalize p, q to probabilities
    p, q = p / p.sum(), q / q.sum()
    # print(np.info(p))
    # print(np.info(q))

    m = (p + q) / 2
    return sp.entropy(p, m, base=base) / 2 + sp.entropy(q, m, base=base) / 2


class ChromLoopData:

    def __init__(self, chrom_name, chrom_size, bin_size):
        self.name = chrom_name
        self.size = chrom_size
        self.bin_size = bin_size

        self.start_list = []
        self.end_list = []
        self.value_list = []
        self.numb_values = 0

        self.largest_loop_len = -1

        self.graph = None

    def add_loop(self, loop_start, loop_end, loop_value):

        loop_length = int(loop_end) - int(loop_start)

        # loops over 10Mb are considered noise
        if loop_length > 10000000:
            return

        if int(loop_end) - int(loop_start) > self.largest_loop_len:
            self.largest_loop_len = int(loop_end) - int(loop_start)

        self.start_list.append(int(loop_start))
        self.end_list.append(int(loop_end))
        self.value_list.append(int(loop_value))

        self.numb_values += 1

    def create_graph(self):
        graph_len = ceil(self.size / self.bin_size)

        self.graph = np.zeros((graph_len, graph_len), dtype=np.uint16)
        self.graph.fill(1)

        # loop_len_len = ceil((self.largest_loop_len + 1) / 50000)
        # loop_len = np.zeros(loop_len_len, dtype=np.uint16)
        # print(loop_len_len)

        for i in range(self.numb_values):
            start = self.start_list[i]
            end = self.end_list[i]
            value = self.value_list[i]

            # loop_len[int((end - start) / 50000)] += 1

            bin_start = int(start / self.bin_size)
            bin_end = int(end / self.bin_size)

            self.graph[bin_start][bin_end] += value * 10

        # plt.plot([x for x in range(len(loop_len))], loop_len)
        # plt.show()

        print(np.max(self.graph))

    def compare(self, o_chromLoopData):
        o_graph = o_chromLoopData.graph
        return jensen_shannon_divergence(o_graph.flatten(), self.graph.flatten())
