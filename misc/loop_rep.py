import logging

VERSION = 6

log = logging.getLogger()


def compare(loop_1, loop_2, bin_size, window_size, window_index=None,
            wanted_chroms=None, is_rep=False, num_peaks=None):

    result = loop_1.compare(loop_2, bin_size, window_size, window_index=window_index,
                            chroms_to_compare=wanted_chroms, is_rep=is_rep,
                            num_peaks=num_peaks)

    if 'w_rep' in result:
        result['main'] = result['w_rep']
    else:
        result['main'] = result['rep']

    '''plt.close()
    plt.hist(j_values, bins=20)
    plt.xlabel('Jensen-Shannon window values')
    plt.ylabel('Frequency')
    plt.xlim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chia_rep: '
              f'{round(j_value, 2)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_hist.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(j_values))],
             j_values)
    plt.xlabel('Jensen-Shannon Window')
    plt.ylabel('j-value')
    plt.ylim(-1, 1)
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chia_rep: '
              f'{round(j_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_j.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(emd_values))],
             emd_values)
    plt.xlabel('Window')
    plt.ylabel('Earth Mover\'s Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chia_rep: '
              f'{round(emd_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_emd.png')
    plt.close()

    plt.plot([int(x * window_size / 1000000) for x in range(len(w_values))],
             w_values)
    plt.xlabel('Window')
    plt.ylabel('Wasserstein Distance')
    plt.title(f'{loop_1.name} vs. {loop_2.name} - chia_rep: '
              f'{round(w_value, 5)}')
    plt.savefig(f'{loop_1.name}_v_{loop_2.name}_window_w.png')
    plt.close()'''

    return result
