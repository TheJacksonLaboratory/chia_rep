import sys
import click
from logging.config import fileConfig

# Add parent folder to PATH if chia_rep is not installed
sys.path.append('..')
import chia_rep

# To allow logging INFO level statements to stdout
fileConfig('log.conf')


@click.command()
@click.argument('input_data_file', type=click.Path(exists=True))
@click.argument('chrom_size_file', type=click.Path(exists=True))
@click.argument('compare_list_file', type=click.Path(exists=True))
@click.argument('window_size', type=int)
@click.argument('bin_size', type=int)
@click.argument('chroms_to_load', nargs=-1)
@click.option('-l', '--min-loop-value', default=1, type=int)
@click.option('-b', '--min-bedgraph-value', default=1, type=int)
@click.option('-p', '--num-peaks', type=int)
@click.option('-g', '--do-output-graph', type=bool)
def main(input_data_file, chrom_size_file, comparison_list_file, window_size,
         bin_size, chroms_to_load, min_loop_value, min_bedgraph_value,
         num_peaks, do_output_graph):
    if 'all' in chroms_to_load:
        chroms_to_load = None
    sample_data_dict = chia_rep.read_data(input_data_file, chrom_size_file,
                                          min_loop_value, min_bedgraph_value,
                                          chroms_to_load)

    chia_rep.preprocess(sample_data_dict, num_peaks)
    emd_scores, j_scores = \
        chia_rep.compare(sample_data_dict, num_peaks,
                         compare_list_file=comparison_list_file,
                         window_size=window_size, bin_size=bin_size,
                         do_output_graph=do_output_graph)

    # compare_list = [
    #     ['sampleA1', 'sampleA2'],
    #     ['sampleA1', 'sampleB1'],
    #     ['sampleA2', 'sampleB1']
    # ]
    # emd_scores, j_scores = \
    #     chia_rep.compare(sample_data_dict, compare_list=compare_list,
    #                      window_size=window_size, bin_size=bin_size)

    chia_rep.output_to_csv(emd_scores, j_scores, window_size, bin_size, num_peaks)


if __name__ == '__main__':
    main()
