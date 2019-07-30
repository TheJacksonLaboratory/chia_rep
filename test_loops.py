from loop_info import LoopInfo
import os

DATA_DIR = 'data'


def compare_loops(loop_info_list):
    len_loop = len(loop_info_list)
    for i in range(len_loop):
        for j in range(i + 1, len_loop, 1):
            loop1 = loop_info_list[i]
            loop2 = loop_info_list[j]
            chrom_loop1 = loop_info_list[i].chrom_dict['chr1']
            chrom_loop2 = loop_info_list[j].chrom_dict['chr1']
            value = chrom_loop1.compare(chrom_loop2)
            print(f'Comparing {loop1.name} vs. {loop2.name}: '
                  f'{value}')

        print()


def read_loops(data_directory):
    loop_info_list = []
    for filename in os.listdir(data_directory):
        if filename.endswith('.BE3'):
            file_path = os.path.join(data_directory, filename)
            loopInfo = LoopInfo(f'{DATA_DIR}/chrom_sizes/hg38.chrom.sizes',
                                file_path)
            loop_info_list.append(loopInfo)

    return loop_info_list


def main():
    for datatype in ['CTCF', 'RNAPII']:
        loop_info_list = read_loops(DATA_DIR + f'/{datatype}')

        compare_loops(loop_info_list)
        break


if __name__ == '__main__':
    main()
