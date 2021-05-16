import itertools
import os
import click


@click.group()
def cli():
    pass


@cli.command()
@click.argument('sample_list_file', type=click.File('r'))
@click.argument('output_pair_file', type=click.File('w'))
def make_pairs(sample_list_file, output_pair_file):
    sample_list = [sample.strip() for sample in sample_list_file]
    for pair in itertools.combinations(sample_list, 2):
        line = "\t".join(pair) + '\n'
        output_pair_file.write(line)


@cli.command()
@click.argument('sample_list_file', type=click.File('r'))
@click.argument('sample_input_file', type=click.File('w'))
@click.argument('sample_data_dir')
def make_sample_input_file(sample_list_file, sample_input_file, sample_data_dir):
    bg_ext = '.bedgraph'
    loop_ext = '.cis.be3'
    peak_ext = '.broadpeak'

    for sample_name in sample_list_file:
        sample_name = sample_name.strip()

        # Get file paths for each needed file
        bg_file_path = None
        loop_file_path = None
        peak_file_path = None
        for file in os.scandir(sample_data_dir):
            if file.name.lower().endswith(bg_ext) and sample_name in file.name:
                bg_file_path = file.path
            elif file.name.lower().endswith(loop_ext) and sample_name in file.name:
                loop_file_path = file.path
            elif file.name.lower().endswith(peak_ext) and sample_name in file.name:
                peak_file_path = file.path

        if not bg_file_path or not loop_file_path or not peak_file_path:
            print(f'Missing files for {sample_name}. Skipping.')
            continue

        sample_input_file.write(f'{sample_name}\t{bg_file_path}\t{peak_file_path}\t{loop_file_path}\n')


if __name__ == '__main__':
    cli()
