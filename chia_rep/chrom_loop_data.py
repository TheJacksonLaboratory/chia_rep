import math
from math import ceil
import numpy as nl
import numpy as npp
import scipy.stats as sp
from scipy.sparse import coo_matrix, csr_matrix
import time
import matplotlib.pyplot as plt
from copy import deepcopy
import logging
from .chia_rep_util import *

from prettytable import PrettyTable

log = logging.getLogger()
log_bin = logging.getLogger('bin')

EMD_WEIGHT = 1
J_WEIGHT = 1
MIN_NUMB_LOOPS = 5
MAX_LOOP_LEN = 1000000  # 1mb

VERSION = 53

MAX_USHRT = 65535
MIN_RATIO_INCREASE = 1.1

PEAK_START = 0
PEAK_END = 1
PEAK_LEN = 2
PEAK_MAX_VALUE = 3
PEAK_MEAN_VALUE = 4


def emd(p, q):
    """
    Finds the Earth Mover's Distance of two 1D arrays

    Parameters
    ----------
    p : 1D Numpy array
    q : 1D Numpy array

    Returns
    -------
    float/int
        The Earth Mover's Distance
    int
        The weight of this calculation
    """
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
    """
    Finds the jensen-shannon divergence

    Parameters
    ----------
    p : 2D Numpy array
        Graph of sample1
    q : 2D Numpy array
        Graph of sample2
    base : int, optional
        Determines base to be used in calculating scipy.entropy (Default is 2)

    Returns
    -------
    float/int
        Jensen-Shannon divergence
    """
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


def get_removed_area(chrom_size, to_remove):
    """
    Gets removed area from both samples so no loops from that area are compared

    Parameters
    ----------
    chrom_size : int
        Size of chromosome
    to_remove : list[2][]
        2D List containing start and end intervals of removed areas

    Returns
    -------
    1D Numpy array
        Array marked with 1 in removed intervals
    """

    # Uses dtype unsigned char because bool doesn't work as well in Cython
    removed_area = np.zeros(chrom_size, dtype=np.uint8)
    for i in range(len(to_remove[0])):
        start = to_remove[0][i]
        end = to_remove[1][i]
        removed_area[start:end] = 1

    return removed_area


class ChromLoopData:
    """
    A class used to represent a chromosome in a sample

    Attributes
    ----------
    name : str
        Name of chromosome
    size : int
        Size of chromosome
    sample_name : str
        Name of the sample (LHH0061, LHH0061_0061H, ...)
    """

    def __init__(self, chrom_name, chrom_size, sample_name):
        """
        Parameters
        ----------
        chrom_name : str
            Name of chromosome
        chrom_size : int
            Size of chromosome
        sample_name : str
            Name of the sample (LHH0061, LHH0061_0061H, ...)
        """

        self.name = chrom_name
        self.size = chrom_size
        self.sample_name = sample_name

        self.start_anchor_list = [[], []]
        self.end_anchor_list = [[], []]
        self.start_list = None  # Later initialized with bedgraph
        self.end_list = None
        self.value_list = []
        self.pet_count_list = []
        self.numb_loops = 0
        self.removed_intervals = [[], []]  # Start/end anchors are on same peak
        self.start_list_peaks = None
        self.end_list_peaks = None

        self.filtered_start = []
        self.filtered_end = []
        self.filtered_values = []
        self.filtered_numb_values = 0

        # Used in filter_with_peaks, keeps track of peaks for each loop
        self.peak_indexes = []
        self.peaks_used = None

        self.max_loop_value = 0

    def add_loop(self, loop_start1, loop_start2, loop_end1, loop_end2,
                 loop_value):

        self.start_anchor_list[0].append(loop_start1)
        self.start_anchor_list[1].append(loop_start2)
        self.end_anchor_list[0].append(loop_end1)
        self.end_anchor_list[1].append(loop_end2)

        self.value_list.append(loop_value)

        self.numb_loops += 1

    def finish_init(self, bedgraph):
        """
        Finishes the construction of this chromosome. Converts lists to numpy
        arrays and calls find_loop_anchor_points

        Parameters
        ----------
        bedgraph : BedGraph
            Used to find the anchor points of each loop

        Returns
        -------
        bool
            Whether the chromosome was successfully made
        """

        if self.numb_loops == 0:
            return False

        self.pet_count_list = np.asarray(self.value_list, dtype=np.uint32)
        self.value_list = np.asarray(self.value_list, dtype=np.float64)
        self.start_anchor_list = np.asarray(self.start_anchor_list,
                                            dtype=np.int32)
        self.end_anchor_list = np.asarray(self.end_anchor_list, dtype=np.int32)

        log.debug(f"Max PET count: {np.max(self.value_list)}")

        # if not bedgraph.has_chrom(self.name):
        if self.name not in bedgraph.chromosome_map:
            log.warning(f"{self.name} was not found in corresponding bedgraph: "
                        f"{bedgraph.name}")
            return False

        self.find_loop_anchor_points(bedgraph)
        return True

    def find_loop_anchor_points(self, bedgraph):
        """
        Finds the exact loop anchor points. Finds peak values for each anchor
        and weighs the loop. Also finds loops that have overlapping start/end
        indexes due to close and long start/end anchors.

        Parameters
        ----------
        bedgraph : BedGraph
            Used to find the anchor points of each loop
        """

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
        self.start_list_peaks = start_list_peaks
        self.end_list_peaks = end_list_peaks
        bedgraph.free_chrom_data(self.name)

        start_list_peaks = start_list_peaks / start_list_peaks.sum()
        end_list_peaks = end_list_peaks / end_list_peaks.sum()

        # Merge peaks that are close together
        # for i in range(self.numb_values):
        #     for j in range(i, self.numb_values):
        #         pass

        for i in range(self.numb_loops):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]

            # Remove anchors that have the same* peak
            # Keep indexes of loop length to avoid comparisons in interval
            if not loop_start < loop_end:
                self.value_list[i] = 0

                # Removed interval goes from
                # (start of start anchor, end of end anchor)
                self.removed_intervals[0].append(self.start_anchor_list[0][i])
                self.removed_intervals[1].append(self.end_anchor_list[1][i])
                continue

            # Weigh each loop based on its corresponding bedgraph peak
            # peak_value = max(start_list_peaks[i], end_list_peaks[i])
            peak_value = start_list_peaks[i] + end_list_peaks[i]
            self.value_list[i] *= peak_value

            # Remove loops over a given threshold
            # loop_length = int(loop_end) - int(loop_start)
            # if loop_length > MAX_LOOP_LEN:
            #     self.value_list[i] = 0

        self.max_loop_value = np.max(self.value_list)
        log.debug(f"Max loop weighted value: {self.max_loop_value}")

    # May have more pre-processing to do?
    def preprocess(self, num_peaks, peak_dict, both_peak_support=False,
                   kept_dir=None):
        return self.filter_with_peaks(num_peaks, peak_dict, both_peak_support,
                                      kept_dir)

    def filter_with_peaks(self, num_peaks, peaks, both_peak_support=False,
                          kept_dir=None):
        """
        Filters out loops without peak support.

        Parameters
        ----------
        num_peaks : int
            Number of peaks to use when filtering
        peaks : list(list)
            List of peaks to use
        both_peak_support : bool, optional
            Whether to only keep loops that have peak support on both sides
            (default is False)
        kept_dir : str, optional
            Directory to output kept peaks and filters
            (default is None)

        Returns
        -------
        bool
            Whether the chromosome had any problems when filtering
        """

        start_time = time.time()

        log.info(f"Filtering {self.sample_name} {self.name} with "
                 f"{num_peaks} peaks...")

        peaks.sort(key=lambda x: x[PEAK_MAX_VALUE], reverse=True)
        # num_wanted_peaks = math.ceil(len(peaks) * num_peaks)
        num_wanted_peaks = min(len(peaks), num_peaks)

        log.debug(f"Num wanted peaks: {num_wanted_peaks}")
        log.debug(f"Top peaks: {peaks[:5]}")

        wanted_peaks = peaks[:num_wanted_peaks]
        min_peak_value = wanted_peaks[-1][2]
        log.debug(f'Min peak value: {min_peak_value}')

        if kept_dir:
            with open(
                    f'{kept_dir}/{self.sample_name}.{num_peaks}.peaks',
                    'a+') as out_file:
                for peak in wanted_peaks:
                    out_file.write(
                        f'{self.name}\t{peak[0]}\t{peak[1]}\t{peak[2]}\n')

        # Get the coverage of each wanted peak
        # Could be used to find the specific peaks for every loop
        index_array = np.zeros(self.size, dtype=np.uint16)
        assert num_wanted_peaks < MAX_USHRT
        for i in range(num_wanted_peaks):
            peak_start = wanted_peaks[i][0]
            peak_end = wanted_peaks[i][1]
            index_array[peak_start:peak_end] = i + 1

        log.debug(f'Time: {time.time() - start_time}')

        numb_deleted = 0
        removed_loop_lengths = []
        removed_loop_values = []
        kept_loop_lengths = []
        kept_indexes = []
        self.peak_indexes = [[], []]
        self.filtered_start = []
        self.filtered_end = []
        self.filtered_values = []

        # plt.title(f'{self.sample_name} Log10(loop span)')
        # plt.xlabel('log10(Loop Span)')
        # plt.ylabel('Density')
        # loop_spans = [self.end_list[i] - self.start_list[i] for i in range(self.numb_loops)]
        # plt.hist(loop_spans, bins=int(1 + np.log2(len(loop_spans))), density=True)
        # plt.savefig(f'{self.sample_name}_loop_span_all')
        # plt.close()

        for i in range(self.numb_loops):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]
            loop_value = self.value_list[i]

            # Loops that are too long can be considered noise
            if loop_end - loop_start > MAX_LOOP_LEN:
                continue

            # From overlapping anchors
            if loop_value == 0:
                continue

            if both_peak_support:
                to_keep = index_array[loop_start] and index_array[loop_end]
            else:
                to_keep = index_array[loop_start] or index_array[loop_end]

            if not to_keep:
                removed_loop_values.append(loop_value)
                numb_deleted += 1
                removed_loop_lengths.append(loop_end - loop_start)
                continue

            self.filtered_start.append(loop_start)
            self.filtered_end.append(loop_end)
            self.filtered_values.append(loop_value)
            kept_indexes.append(i)

            kept_loop_lengths.append(loop_end - loop_start)

            # Unused for now
            self.peak_indexes[0].append((
                wanted_peaks[index_array[loop_start] - 1][0],
                wanted_peaks[index_array[loop_start] - 1][1]))
            self.peak_indexes[1].append((
                wanted_peaks[index_array[loop_end] - 1][0],
                wanted_peaks[index_array[loop_end] - 1][1]))

        self.filtered_start = np.array(self.filtered_start, dtype=np.int32)
        self.filtered_end = np.array(self.filtered_end, dtype=np.int32)
        self.filtered_values = np.array(self.filtered_values)
        self.filtered_numb_values = self.filtered_start.size

        if kept_dir is not None:
            with open(
                    f'{kept_dir}/{self.sample_name}.{num_peaks}.loops',
                    'a+') as out_file:
                for i in kept_indexes:
                    out_file.write(
                        f'{self.name}\t{self.start_anchor_list[0][i]}\t'
                        f'{self.start_anchor_list[1][i]}\t{self.name}\t'
                        f'{self.end_anchor_list[0][i]}\t'
                        f'{self.end_anchor_list[1][i]}\t'
                        f'{self.pet_count_list[i]}\t'
                        f'{self.value_list[i]}\n')

        log.debug(f'Total loops: {self.numb_loops}')
        log.debug(f"Number of loops removed: {numb_deleted}")
        log.info(f"Number of loops kept: {self.filtered_numb_values}")

        if self.filtered_numb_values == 0:
            log.warning(f"No loops left. Skipping")
            return False

        if numb_deleted > 0:
            log.debug(
                f'Avg loop length removed: {np.mean(removed_loop_lengths)}')
            log.debug(f'Avg loop value removed: {np.mean(removed_loop_values)}')
        else:
            log.debug(f'Avg loop length removed: N/A')
            log.debug(f'Avg loop value removed: N/A')
        log.debug(f'Avg loop length kept: {np.mean(kept_loop_lengths)}')
        log.debug(f'Avg loop value kept: {np.mean(self.filtered_values)}')
        log.debug(f'Largest loop value kept: {np.max(self.filtered_values)}')
        log.debug(f'Time taken: {time.time() - start_time}\n')

        return True

    def create_graph(self, loops, num_loops, bin_size, window_size,
                     random=False, to_debug=False):
        """
        Creates a bin-based graph to easily compare loops

        Parameters
        ----------
        loops : list
            List of loop indexes from self.filtered_*
        num_loops : int
            Number of loops to use when making the graph
        bin_size : int
        window_size : int
        random : bool, optional
            Randomly pick which loops to use (Default is False)
        to_debug : bool, optional
            Log loops used in the graph (Default is False)

        Returns
        -------
        Numpy 2D array
        """

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

            if indexes is not None:
                loop_index = loops[indexes[i]]
            else:
                loop_index = loops[i]

            value = self.filtered_values[loop_index]
            start = self.filtered_start[loop_index]
            end = self.filtered_end[loop_index]
            orig_start = start
            orig_end = end

            start = start % window_size
            end = end % window_size

            # loop_len[int((end - start) / self.bin_size)] += 1

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            graph[bin_start][bin_end] += value

            # Get the other side of the graph as well for emd calculation
            # graph[bin_end][bin_start] += value

            # Avoid double counting the middle
            # if bin_end != bin_start:
            #    graph[bin_end][bin_start] += value

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
                    # graph[k][j] += value
                    # if j != k:
                    #    graph[k][j] += value

            # if to_debug:
            #     log.debug(
            #         f'{self.sample_name}\t{orig_start}\t{orig_end}\t{value}')

            num_loops_used += 1

        # log.info(f"Number of loops in {self.sample_name} graph: {num_loops_used}")

        # plt.plot([x for x in range(len(loop_len))], [np.log(x) for x in loop_len])
        # plt.show()

        # log.info(f'Max value in graph: {np.max(graph)}')
        # return graph, total_PET_count / self.total_loop_value
        return graph

    def get_stats(self, graph, o_graph, o_chrom, graph_type):
        """
        Get the reproducibility stat.

        Jensen Shannon or EMD. Though EMD seems to be better.

        Parameters
        ----------
        graph : 2D Numpy array
            Graph created from window from sample1
        o_graph : 2D Numpy array
            Graph created from window from sample2
        o_chrom : ChromLoopData
            sample2
        graph_type : str
            Description of type of comparison

        Returns
        -------
        dict
            graph_type : str
                Information about which comparison method was used
            rep : str
                Reproducibility statistic for this window comparison
            w : str
                Weight of this window
        """

        max_emd_dist = graph[0].size - 1
        log.debug(f'Graph size: {max_emd_dist + 1}')

        result = {'graph_type': graph_type}

        max_graph = np.max(graph)
        max_o_graph = np.max(o_graph)

        result['w'] = max_graph / self.max_loop_value + \
                      max_o_graph / o_chrom.max_loop_value

        if max_graph == 0 or max_o_graph == 0:
            result['rep'] = -1
            return result

        # with open(f'graphs/{self.sample_name}_{window_start}.csv', 'w') as out:
        #     sparse_mat = coo_matrix(graph)
        #     for row, col, value in zip(sparse_mat.row, sparse_mat.col, sparse_mat.data):
        #         out.write("({0}, {1})\t{2}\n".format(row, col, round(value, 6)))
        # with open(f'graphs/{o_chrom.sample_name}_{window_start}.csv', 'w') as out:
        #     sparse_mat = coo_matrix(graph)
        #     for row, col, value in zip(sparse_mat.row, sparse_mat.col, sparse_mat.data):
        #         out.write("({0}, {1})\t{2}\n".format(row, col, round(value, 6)))

        # total_weight = num_loops + num_o_loops
        # total_weight = graph.sum() + o_graph.sum()

        j_value = 0
        # start_time = time.time()
        # graph_flat = graph.flatten()
        # o_graph_flat = o_graph.flatten()
        # j_divergence = jensen_shannon_divergence(graph_flat, o_graph_flat)
        # log.debug(f'Jensen-Shannon time: {time.time() - start_time}')
        #
        # # Make j_value range from -1 to 1
        # j_value = 2 * (1 - j_divergence) - 1

        # complete_graph(graph, max_emd_dist)
        # complete_graph(o_graph, max_emd_dist)

        # Calculate emd for all rows and columns -> Take weighted average
        emd_value = 0
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
        max_emd_weight = np.max(emd_weight_list)

        if max_emd_weight == 0:
            emd_dist = 0
        else:
            emd_dist = np.average(emd_distance_list, weights=emd_weight_list)
        # emd_dist = np.mean(emd_distance_list)

        # Scale from -1 to 1
        # Since most values are in bottom half, don't scale linearly
        emd_value = 2 * (emd_dist - max_emd_dist) * (emd_dist - max_emd_dist) / \
                    (max_emd_dist * max_emd_dist) - 1

        # linear scale
        # emd_value = 1 - 2 / max_emd_dist * emd_dist

        if max_emd_weight == 0 and result['w'] != 0:
            log.error(f'Total Weight: {result["w"]} with 0 emd dist')

        result['rep'] = emd_value * EMD_WEIGHT + j_value * J_WEIGHT
        # result['rep'] = j_value
        # result['rep'] = emd_value
        # table.add_row([graph_type, j_value, emd_value, result['w']])

        return result

    def compare_randomly_choose(self, loops, o_loops, o_chrom, bin_size,
                                window_size, inner_table, num_iterations=10):
        """
        Randomly choose which loops to compare to better match samples with
        differing sequencing depth

        (Deprecated) (Doesn't work?)
        """

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
        """
        Only compare top loops to better match samples with differing sequencing
         depth

        (Deprecated) (Doesn't work?)
        """

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
        log.debug(f'Sub Sample time: {time.time() - start_time} --------------')

        return result

    def compare_only_match(self, graph, o_graph, o_chrom, inner_table,
                           graph1_smaller):
        """
        Only compare loops that match with other sample to better match
        differing sequencing depth

        (Deprecated) (Doesn't work?)
        """

        start_time = time.time()

        graph_type = ['mod', 'orig']

        # Returns [(mod1, orig1), (mod2, orig2)]
        new_graph_list = get_matching_graphs(graph, o_graph)
        value_table = PrettyTable(['graph_type', 'rep', 'weight'])

        # Compare only mod vs. orig
        for i in range(2):
            j = 1 - i

            if graph1_smaller and i == 0:
                continue

            new_graph = new_graph_list[0][i]
            new_o_graph = new_graph_list[1][j]

            result = self.get_stats(new_graph, new_o_graph, o_chrom,
                                    f'{graph_type[i]}_{graph_type[j]}',
                                    inner_table)

            value_table.add_row([x for x in list(result.values())])

        log_bin.info(value_table)
        log.debug(f'Match time: {time.time() - start_time} ------------------')

        # Not sure what to return: max, avg?
        return result

    def find_diff_loops(self, o_chrom):
        with open(f'diff_loops/{self.sample_name}_{o_chrom.sample_name}.loops',
                  'a+'):
            pass
        pass

    def compare(self, o_chrom, window_start, window_end, bin_size,
                is_rep=False):
        """
        Compare a window of this chromosome of another chromosome from another
        sample

        Parameters
        ----------
        o_chrom : ChromLoopData
            The other chromosome
        window_start : int
            The start of the window
        window_end : int
            The end of the window
        bin_size : int
            Determines which loops are the same by putting them into bins
        is_rep : bool, optional
            Debugging purposes

        Returns
        -------
        dict
            graph_type : str
                Information about which comparison method was used
            rep : str
                Reproducibility statistic for this window comparison
            w : str
                Weight of this window
        """

        # def compare(self, o_chrom, combined_peak_list, is_rep=False):
        start_time = time.time()

        if window_end > self.size:
            window_end = self.size

        if window_start >= self.size:
            log.error(f"Start of window ({window_start}) is larger than "
                      f"{self.name} size: {self.size}")
            return {'graph_type': 'error',
                    'rep': 0,
                    'w': 0
                    }

        # log_all.info(f'{self.sample_name} vs. {o_chrom.sample_name} '
        #             f'{self.name}:{window_start} - {window_end}')

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
                          self.filtered_start, self.filtered_end,
                          self.filtered_values)
        o_loops = get_loops(window_start, window_end, combined_removed,
                            o_chrom.filtered_start, o_chrom.filtered_end,
                            o_chrom.filtered_values)
        # loops = get_loops(0, self.size, combined_removed,
        #                   self.filtered_start, self.filtered_end,
        #                   self.filtered_values)
        # o_loops = get_loops(0, self.size, combined_removed,
        #                     o_chrom.filtered_start, o_chrom.filtered_end,
        #                     o_chrom.filtered_values)
        log.debug(f'Get Loop time: {time.time() - loop_time}')

        num_loops = len(loops)
        num_o_loops = len(o_loops)
        log.debug(f"Numb of loops in {self.sample_name}: {num_loops}")
        log.debug(f"Numb of loops in {o_chrom.sample_name}: {num_o_loops}")

        # value_table = PrettyTable(['graph_type', 'rep', 'weight'])
        # inner_value_table = \
        #     PrettyTable(['graph_type', 'j_value', 'emd', 'weight'])

        all_results = []
        result = None

        if num_loops == 0 and num_o_loops == 0:
            result = {
                'graph_type': 'pass_none',
                'rep': 1,
                'w': 0
            }

        # Randomly subsample from larger sample to better match miseq vs. hiseq
        # result = self.compare_randomly_choose(loops, o_loops, o_chrom,
        #                                       bin_size, window_size,
        #                                       inner_value_table)
        # log_bin.info(inner_value_table)
        # inner_value_table.clear_rows()
        # all_results.append(result)
        # value_table.add_row([x for x in list(result.values())])

        # Subsample top loops from larger sample to better match miseq vs. hiseq
        # result = self.compare_subsample(loops, o_loops, o_chrom, bin_size,
        #                                 window_size, inner_value_table)
        # log_bin.info(inner_value_table)
        # inner_value_table.clear_rows()
        # all_results.append(result)
        # value_table.add_row([x for x in list(result.values())])

        # if greater_numb / lesser_numb < 8:
        if result is None:
            # Make graphs using all loops in the window
            graph = self.create_graph(loops, num_loops, bin_size, window_size,
                                      to_debug=is_rep)
            o_graph = o_chrom.create_graph(o_loops, num_o_loops, bin_size,
                                           window_size, to_debug=is_rep)
            # merged_list = merge_peaks(combined_peak_list)
            # graph = create_peak_graph(merged_list, loops, self)
            # o_graph = create_peak_graph(merged_list, o_loops, o_chrom)

            # if self.numb_values / o_chrom.numb_values > 5 or \
            #        o_chrom.numb_values / self.numb_values > 5:
            # Match graphs only where loops overlap
            # graph1_smaller = self.numb_values < o_chrom.numb_values
            #
            # match_time = time.time()
            # result = self.compare_only_match(graph, o_graph, o_chrom,
            #                                  inner_value_table,
            #                                  graph1_smaller)
            # log.debug(f'Match time: {time.time() - match_time} -----------')
            # log_bin.info(inner_value_table)
            # inner_value_table.clear_rows()
            # all_results.append(result)
            # value_table.add_row([x for x in list(result.values())])

            # Normal comparison
            norm_time = time.time()
            result = self.get_stats(graph, o_graph, o_chrom, 'norm')
            log.debug(f'Norm time: {time.time() - norm_time} -------------')
            # log_bin.info(inner_value_table)
            # inner_value_table.clear_rows()
            all_results.append(result)
            # value_table.add_row([x for x in list(result.values())])

        # log_bin.info(value_table)
        log.debug(f'Total time: {time.time() - start_time}')
        log.debug(f'Reproducibility: {result["rep"]}')
        if is_rep:
            for i in [0.5, 0, -0.5]:
                if result['rep'] < i:
                    log.debug(f'Less than {i}')

        # log_all.info(
        #    '-----------------------------------------------------------------')

        return result


def merge_peaks(combined_peak_list):
    combined_peak_list.sort(key=lambda x: x[PEAK_START])
    merged_peak_list = []
    peak_diffs = []
    for higher_peak in combined_peak_list:
        if not merged_peak_list:
            merged_peak_list.append(higher_peak)
            continue

        lower_peak = merged_peak_list[-1]
        if higher_peak[PEAK_START] <= lower_peak[PEAK_END]:
            if lower_peak[PEAK_END] > higher_peak[PEAK_END]:
                lower_peak[PEAK_MAX_VALUE] += higher_peak[PEAK_MAX_VALUE]
                continue

            dist_diff = higher_peak[PEAK_START] - lower_peak[PEAK_START]
            if dist_diff / lower_peak[PEAK_LEN] >= 0.5 and \
                    dist_diff / higher_peak[PEAK_LEN] >= 0.5:
                peak_diffs.append(dist_diff)
                lower_peak[PEAK_END] = higher_peak[PEAK_END]
                lower_peak[PEAK_MAX_VALUE] += higher_peak[PEAK_MAX_VALUE]
                lower_peak[PEAK_LEN] = lower_peak[PEAK_END] - lower_peak[
                    PEAK_START]
                continue

        merged_peak_list.append(higher_peak)

    log.info(f"Merged {len(combined_peak_list) - len(merged_peak_list)} peaks")
    log.info(f"Avg space between merged peaks: {np.mean(peak_diffs)}")

    return merged_peak_list


def create_peak_graph(merged_list, loop_index_list, chrom_loop_data):
    self = chrom_loop_data

    graph_len = len(merged_list)
    graph = np.zeros((graph_len, graph_len), dtype=np.float64)

    peak_array = np.full(self.size, -1, dtype=np.int16)
    for i, peak in enumerate(merged_list):
        start = peak[PEAK_START]
        end = peak[PEAK_END]
        peak_array[start:end] = i

    for loop_index in loop_index_list:
        value = self.filtered_values[loop_index]
        start = self.filtered_start[loop_index]
        end = self.filtered_end[loop_index]

        start_index = peak_array[start]
        end_index = peak_array[end]

        if start_index == -1 or end_index == -1:
            log.error("ERROR")
            continue

        graph[start_index][end_index] += value

    return graph