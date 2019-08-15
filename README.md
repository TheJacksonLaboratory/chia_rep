# chiapet_rep
TODO

# Features:
- TODO

# Usage:
### Installation:
```bash
git clone https://github.com/c0ver/chiapet_rep.git
pip install chiapet_rep/
```

### Input data and compare:
```python
from chiapet_rep import reproducibility
from chiapet_rep import loop_rep

# Test data provided in github repository
TEST_DATA_DIR = 'chiapet_rep/test/data'
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 3000000  # 3mb

# To see log info statements (optional)
from logging.config import fileConfig
fileConfig('chiapet_rep/test/test.conf')


# Since reading in bedgraph file can take a long time, load them first if in an interactive session
bedgraph_dict = reproducibility.read_bedGraphs(TEST_DATA_DIR, f'{TEST_DATA_DIR}/hg38.chrom.sizes')

# Takes ~30 seconds for included test data
loop_dict = reproducibility.read_data(loop_data_dir=TEST_DATA_DIR,
                                      peak_data_dir=TEST_DATA_DIR,
                                      chrom_size_file=f'{TEST_DATA_DIR}/hg38.chrom.sizes',
                                      is_hiseq=False,  # Determines if match comparison is made (TODO)
                                      bedgraph_dict=bedgraph_dict)

rep, non_rep, scores = reproducibility.compare(loop_dict, loop_rep.compare,
                                               bin_size=BIN_SIZE,
                                               window_size=WINDOW_SIZE)

reproducibility.output_results(rep, non_rep)
```