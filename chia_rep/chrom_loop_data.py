import sys
from math import ceil
import numpy as np
import scipy.stats as sp
import time
import logging
from pyBedGraph import BedGraph
from typing import Dict, Tuple

from .util import *

log = logging.getLogger()
log_bin = logging.getLogger('bin')

MAX_USHRT = 65535

PEAK_START_INDEX = 0
PEAK_END_INDEX = 1
PEAK_LEN_INDEX = 2
PEAK_MAX_VALUE_INDEX = 3
PEAK_MEAN_VALUE_INDEX = 4

# The length of each normalization call
NORM_LEN = 100


def emd(
    p: np.ndarray,
    q: np.ndarray
) -> Tuple[float, float]:
    """
    Finds the Earth Mover's Distance of two 1D arrays

    Parameters
    ----------
    p : ndarray
    q : ndarray

    Returns
    -------
    tuple[float, float]
        The Earth Mover's Distance; the weight of this calculation
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


def jensen_shannon_divergence(
    p: np.ndarray,
    q: np.ndarray,
    base=2
) -> float:
    """
    Finds the jensen-shannon divergence

    Parameters
    ----------
    p : np.ndarray[np.float64]
        Graph of sample1
    q : np.ndarray[np.float64]
        Graph of sample2
    base : int, optional
        Determines base to be used in calculating scipy.entropy (Default is 2)

    Returns
    -------
    float
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


def output_graph(
    output_dir: str,
    chrom_name: str,
    window_start: int,
    window_end: int,
    graph: np.ndarray,
    sample_name: str,
    do_output_graph: bool = False
) -> None:
    """
    Outputs graph to file

    Parameters
    ----------
    output_dir
    chrom_name
    window_start
    window_end
    graph
    sample_name
    do_output_graph

    Returns
    -------
    None
    """
    with open(f'{output_dir}/{chrom_name}.txt', 'a') as out_file:
        n_non_zeros = np.count_nonzero(graph)
        n_zeros = graph.size - n_non_zeros
        out_file.write(
            f'{sample_name}\t{chrom_name}\t{window_start}\t{window_end}\t'
            f'{n_non_zeros}\t{n_zeros}\n')

    start_time = time.time()
    if chrom_name == 'chr1' and do_output_graph:
        with open(f'{output_dir}/chr1_graph.txt', 'a') as out_file, \
                np.printoptions(threshold=sys.maxsize, precision=6,
                                linewidth=sys.maxsize, suppress=True) as _:
            generated_graph = '\t'.join([str(x) for x in graph.flatten()])

            out_file.write(f'{sample_name}\t{chrom_name}\t{window_start}\t'
                           f'{window_end}\t{generated_graph}\n')
    log.debug(f'Printing graph took {time.time() - start_time}s')


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

    def __init__(
        self,
        chrom_name: str,
        chrom_size: int,
        sample_name: str
    ):
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

        self.norm_list = None

        self.filtered_start = []
        self.filtered_end = []
        self.filtered_values = []
        self.filtered_numb_values = 0
        self.filtered_anchors = []
        self.kept_indexes = []

        # Used in filter_with_peaks, keeps track of peaks for each loop
        self.peak_indexes = []
        self.peaks_used = None

        self.max_loop_value = 0

    def add_loop(
        self,
        loop_start1: int,
        loop_end1: int,
        loop_start2: int,
        loop_end2: int,
        loop_value: int
    ) -> None:

        self.start_anchor_list[0].append(loop_start1)
        self.start_anchor_list[1].append(loop_end1)
        self.end_anchor_list[0].append(loop_start2)
        self.end_anchor_list[1].append(loop_end2)

        self.value_list.append(loop_value)

        self.numb_loops += 1

    def finish_init(
        self,
        bedgraph: BedGraph
    ) -> bool:
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

        log.debug(f"Max PET count: {np.max(self.pet_count_list)}")

        # if not bedgraph.has_chrom(self.name):
        if self.name not in bedgraph.chromosome_map:
            log.warning(f"{self.name} was not found in corresponding bedgraph: "
                        f"{bedgraph.name}")
            return False

        self.find_loop_anchor_points(bedgraph)
        return True

    def find_loop_anchor_points(
        self,
        bedgraph: BedGraph
    ):
        """
        Finds the exact loop anchor points.

        Finds peak values for each anchor and weighs the loop. Also finds loops
        that have overlapping start/end indexes due to close and long start/end
        anchors.

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

        for i in range(self.numb_loops):
            # loop_start = self.start_list[i]
            # loop_end = self.end_list[i]

            # Remove anchors that have the same* peak
            # Keep indexes of loop length to avoid comparisons in interval
            # if not loop_start < loop_end:
            #     self.value_list[i] = 0
            #
            #     # Removed interval goes from
            #     # (start of start anchor, end of end anchor)
            #     self.removed_intervals[0].append(self.start_anchor_list[0][i])
            #     self.removed_intervals[1].append(self.end_anchor_list[1][i])
            #     continue

            # Weigh each loop based on its corresponding bedgraph peak
            # peak_value = max(start_list_peaks[i], end_list_peaks[i])
            peak_value = start_list_peaks[i] + end_list_peaks[i]
            self.value_list[i] *= peak_value

        self.max_loop_value = np.max(self.value_list)

        # Should be very small due to peaks being weighted earlier
        log.debug(f"Max loop weighted value: {self.max_loop_value}")

    # May have more pre-processing to do besides filtering later?
    # Useless extra function otherwise
    def preprocess(
        self,
        peak_list: list,
        both_peak_support: bool = False
    ) -> bool:
        return self.filter_with_peaks(peak_list, both_peak_support)

    def filter_with_peaks(
        self,
        peak_list: list,
        both_peak_support: bool = False
    ) -> bool:
        """
        Filters out loops without peak support.

        Get coverage of peaks that have been chosen to be used. Find loops that
        are not within that coverage and filter them out.

        Parameters
        ----------
        peak_list : list(list)
            List of peaks to use
        both_peak_support : bool, optional
            Whether to only keep loops that have peak support on both sides
            (default is False)

        Returns
        -------
        bool
            Whether the chromosome had any problems when filtering
        """

        start_time = time.time()

        num_peaks = len(peak_list)
        min_peak_value = peak_list[-1][PEAK_MAX_VALUE_INDEX]

        log.info(f"Filtering {self.sample_name} {self.name} with "
                 f"{num_peaks} peaks...")
        log.debug(f"Top peaks: {peak_list[:3]}")
        log.debug(f"Bottom peaks: {peak_list[-3:]}")
        log.debug(f'Min peak value: {min_peak_value}')

        # Get the coverage of each wanted peak
        # Could be used to find the specific peaks for every loop
        index_array = np.zeros(self.size, dtype=np.uint16)

        if num_peaks >= MAX_USHRT:
            log.warning(f'Number of peaks: {num_peaks} is greater than max_unsigned_short: {MAX_USHRT}')
        for i in range(num_peaks):
            peak_start = peak_list[i][0]
            peak_end = peak_list[i][1]
            index_array[peak_start:peak_end] = i + 1

        log.debug(f'Time: {time.time() - start_time}')

        numb_deleted = 0
        removed_loop_lengths = []
        removed_loop_values = []
        kept_loop_lengths = []
        self.kept_indexes = []
        self.peak_indexes = [[], []]
        self.filtered_start = []
        self.filtered_end = []
        self.filtered_values = []
        self.filtered_anchors = []

        for i in range(self.numb_loops):
            loop_start = self.start_list[i]
            loop_end = self.end_list[i]
            loop_value = self.value_list[i]

            if loop_start > loop_end:
                temp_val = loop_start
                loop_start = loop_end
                loop_end = temp_val

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
            self.filtered_anchors.append([self.start_anchor_list[0][i],
                                          self.start_anchor_list[1][i],
                                          self.start_list_peaks[i]])
            self.filtered_anchors.append([self.end_anchor_list[0][i],
                                          self.end_anchor_list[1][i],
                                          self.start_list_peaks[i]])
            self.kept_indexes.append(i)

            kept_loop_lengths.append(loop_end - loop_start)

            # Unused for now
            self.peak_indexes[0].append((
                peak_list[index_array[loop_start] - 1][0],
                peak_list[index_array[loop_start] - 1][1]))
            self.peak_indexes[1].append((
                peak_list[index_array[loop_end] - 1][0],
                peak_list[index_array[loop_end] - 1][1]))

        self.filtered_start = np.array(self.filtered_start, dtype=np.int32)
        self.filtered_end = np.array(self.filtered_end, dtype=np.int32)
        self.filtered_values = np.array(self.filtered_values)
        self.filtered_numb_values = self.filtered_start.size
        self.kept_indexes = np.array(self.kept_indexes, dtype=np.int32)

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

    def create_graph(
        self,
        loops: list,
        bin_size: int,
        window_size: int,
        random: bool = False,
        num_loops: int = 0
    ) -> np.ndarray:
        """
        Creates a bin-based graph to easily compare loops

        Parameters
        ----------
        loops : list
            List of loop indexes from self.filtered_*
        bin_size : int
        window_size : int
        random : bool, optional
            Randomly pick which loops to use (Default is False)
            Useless for now since different sequencing depth is ignored
        num_loops : int, optional
            Number of loops to use when making the graph

        Returns
        -------
        np.ndarray
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

            bin_start = int(start / bin_size)
            bin_end = int(end / bin_size)

            if bin_end < bin_start:
                log.error(
                    f'{orig_start}\t{orig_end}\t{start}\t{end}\t{bin_start}\t{bin_end}')
                temp_val = bin_start
                bin_start = bin_end
                bin_end = temp_val

            graph[bin_start][bin_end] += value

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

            num_loops_used += 1

        return graph

    def get_stats(
        self,
        graph: np.ndarray,
        o_graph: np.ndarray,
        o_chrom: 'ChromLoopData'
    ):
        """
        Get the reproducibility stat.

        Jensen Shannon and EMD. Though EMD seems to be better.

        Parameters
        ----------
        graph : np.ndarray
            2D Numpy array, float64
            Graph created from window from sample1
        o_graph : np.ndarray
            2D Numpy array, float64
            Graph created from window from sample2
        o_chrom : ChromLoopData
            sample2

        Returns
        -------
        dict
            emd_value : float
                Comparison value based on the Earth mover's Distance formula
            j_value : float
                Comparison value based on the Jensen-Shannon formula
            emd_dist : float
                Result of the Earth mover's Distance formula
            j_divergence : float
                Result of the Jensen-Shannon formula
            w : str
                Weight of this window
        """
        result = {}
        max_graph = np.max(graph)
        max_o_graph = np.max(o_graph)
        result['w'] = max_graph / self.max_loop_value + max_o_graph / o_chrom.max_loop_value

        if max_graph == 0 or max_o_graph == 0:
            if max_graph == 0:
                log.debug('No loops in sample A')
            else:
                log.debug('No loops in sample B')

            result['j_divergence'] = 1
            result['j_value'] = -1
            result['emd_value'] = -1
            result['linear_emd_value'] = -1
            result['emd_dist'] = 1
            return result

        log.debug('Loops are found in both samples')

        graph_flat = graph.flatten()
        o_graph_flat = o_graph.flatten()

        j_divergence = jensen_shannon_divergence(graph_flat, o_graph_flat)

        # Make j_value range from -1 to 1
        j_value = 2 * (1 - j_divergence) - 1

        # Calculate emd for all rows and columns -> Take weighted average
        emd_distance_list = []
        emd_weight_list = []
        for k in range(graph[0].size):
            emd_dist, emd_weight = emd(graph[k], o_graph[k])
            emd_distance_list.append(emd_dist)
            emd_weight_list.append(emd_weight)

            emd_dist, emd_weight = emd(graph[:, k], o_graph[:, k])
            emd_distance_list.append(emd_dist)
            emd_weight_list.append(emd_weight)

        max_emd_weight = np.max(emd_weight_list)

        if max_emd_weight == 0:
            overall_emd_dist = 0
        else:
            overall_emd_dist = np.average(emd_distance_list,
                                          weights=emd_weight_list)
        # overall_emd_dist = np.mean(emd_distance_list)

        # Higher emd_dist == samples are more different
        # Lower emd_dist == samples are more similar
        max_emd_dist = graph.shape[0] - 1
        numerator = overall_emd_dist - max_emd_dist
        emd_value = 2 * numerator * numerator / (
                max_emd_dist * max_emd_dist) - 1

        # Linear scale
        # emd_value = 1 - 2 / max_emd_dist * overall_emd_dist
        # linear_emd_value = 1 - 2 * overall_emd_dist

        if max_emd_weight == 0 and result['w'] != 0:
            log.error(f'Total Weight: {result["w"]} with 0 emd dist')

        result['j_divergence'] = j_divergence
        result['j_value'] = j_value
        result['emd_value'] = emd_value
        result['emd_dist'] = overall_emd_dist
        # result['linear_emd_value'] = linear_emd_value

        return result

    def compare(
        self,
        o_chrom: 'ChromLoopData',
        window_start: int,
        window_end: int,
        window_size: int,
        bin_size: int,
        num_peaks: any,
        output_dir: str = 'output',
        do_output_graph: bool = False
    ) -> Dict[str, float]:
        """
        Compare a window of this chromosome to another chromosome from another
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
        window_size : int
            Needed to find the exact loop start/end index in this graph
        num_peaks : any
        output_dir : str, optional
            Directory to output data
        do_output_graph : bool, optional
            Whether to output graph used for comparison (Default is False)

        Returns
        -------
        dict
            emd_value : float
                Comparison value based on the Earth mover's Distance formula
            j_value : float
                Comparison value based on the Jensen-Shannon formula
            emd_dist : float
                Result of the Earth mover's Distance formula
            j_divergence : float
                Result of the Jensen-Shannon formula
            w : float
                Weight of this window
        """

        return_skeleton = {
            'emd_value': 0,
            'j_value': 0,
            'emd_dist': 0,
            'j_divergence': 0,
            'w': 0
        }

        if window_end > self.size:
            window_end = self.size

        if window_start >= self.size:
            log.error(f"Start of window ({window_start}) is larger than "
                      f"{self.name} size: {self.size}")
            return return_skeleton

        log.debug(f'{self.sample_name} vs. {o_chrom.sample_name} '
                  f'{self.name}:{window_start} - {window_end}')

        # Get loop indexes in the window
        loops = get_loops(window_start, window_end, self.filtered_start,
                          self.filtered_end, self.filtered_values)
        o_loops = get_loops(window_start, window_end, o_chrom.filtered_start,
                            o_chrom.filtered_end, o_chrom.filtered_values)
        num_loops = len(loops)
        num_o_loops = len(o_loops)
        log.debug(f"Numb of loops in {self.sample_name}: {num_loops}")
        log.debug(f"Numb of loops in {o_chrom.sample_name}: {num_o_loops}")

        if num_loops == 0 and num_o_loops == 0:
            result = return_skeleton
            log.debug('No loops in either sample')
        else:
            # Make graphs using all loops in the window
            graph = self.create_graph(loops, bin_size, window_size)
            o_graph = o_chrom.create_graph(o_loops, bin_size, window_size)

            comparison_name = f'{self.sample_name}_{o_chrom.sample_name}'
            param_str = f'{window_size}.{bin_size}.{num_peaks}'
            parent_dir = f'{output_dir}/{param_str}/comparisons/{comparison_name}'
            output_graph(parent_dir, self.name, window_start, window_end, graph,
                         self.sample_name, do_output_graph=do_output_graph)
            output_graph(parent_dir, o_chrom.name, window_start, window_end,
                         o_graph, o_chrom.sample_name,
                         do_output_graph=do_output_graph)

            result = self.get_stats(graph, o_graph, o_chrom)

        log.debug(f'emd_dist: {result["emd_dist"]}')
        log.debug(f'emd_value: {result["emd_value"]}')
        log.debug(f'j_divergence: {result["emd_value"]}')
        log.debug(f'j_value: {result["j_value"]}')

        log.debug(
            '-----------------------------------------------------------------')

        return result
