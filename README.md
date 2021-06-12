# ChIA-Rep 
A Python package for assessing the reproducibility of ChIA-PET datasets.    
    
## Methods Overview 
### Reading in Data
**Reading peaks**\
Read peak interval from peak file and get the peak value from max value within
the interval from the bedgraph file. 

**Reading loops**\
Since loops come with a start interval and end interval, we assign an anchor
within each interval for each loop based on the largest value within the 
interval. Additionally, each loop is weighted by the anchor intensity of both
the start and end interval.

**Deadzones**\
Remove loops that start and end on the same peak. Or if the loop start is 
somehow past the loop end.

example:\
loop_start_interval = (0, 10)\
loop_end_interval = (5, 8)\
peaks are at index 5 and 9\
peak value at index 9 > peak value at index 5\
Therefore, the loop start anchor is at index 9 and the loop end anchor is at index 5.

We create a "deadzone" from 0 to 10 in the chromosome. When creating graphs to 
compare two chromosomes, we combine the "deadzones" from each chromosome and 
ignore loops from either chromosome in the combined deadzones. Therefore, the 
created graph for each chromosome can be different for each comparison.

### Preprocessing  
- Filter out all peaks that are smaller than a certain value: `num_peaks`
- Find the kept peak ratio from the `base_chrom` and use it for other 
  chromosomes
- Filter out loops that are too long (> 1M)
- Filter out a loop if neither anchor overlaps with a kept peak
    
### Graph representation and comparison (Not exactly correct)
- For a non-overlapping window, bin the loops into bins of fixed size
- Create an adjacency matrix for each window, where index (bin1, bin2) contains 
  the value from the loops going from bin1 to bin2
- Convert each adjacency matrix into a probability vector by reading row-by-row
- Compute the Jensen-Shannon divergence and the Earth Mover's Distance (EMD) \
  between two probability vectors
- Transform each value to be between -1 (dissimilar) and 1 (similar)
- Take the weighted average of values from windows in a chromosome
- Take the average of values from chromosomes to produce a genome-wide 
  reproducibility value
    
### Example
Given two ChIA-PET datasets, create adjacency matrices A1 and A2

**Adjacency matrix A1**   

|         | bin1   | bin2   | bin3   | bin4   |
|-------- |------  |------  |------  |------  | 
| bin1    | 3      | 2      | 0      | 1      |
| bin2    |        | 1      | 5      | 3      |
| bin3    |        |        | 10     | 9      |
| bin4    |        |        |        | 20     |
    
**Adjacency matrix A2**   

|         | bin1   | bin2   | bin3   | bin4   |
|-------- |------  |------  |------  |------  | 
| bin1    | 4      | 5      | 1      | 4      |
| bin2    |        | 3      | 2      | 3      |
| bin3    |        |        | 7      | 9      |
| bin4    |        |        |        | 27     |
    
**Probability vectors p_A1 and p_A2**   
- p_A1 = (0.06, 0.05, 0, 0.02, 0.02, 0.009, 0.06, 0.19, 0.17, 0.37)
- p_A2 = (0.06, 0.08, 0.02, 0.06, 0.05, 0.03, 0.05, 0.11, 0.14, 0.42)

## Results
- ChIA-Rep can clearly distinguish between replicates and non-replicates
- Generally, replicates have positive values and non-replicates have negative values
- Can take 0 as a threshold to determine the similarity
    
## Usage 
### Dependencies:
```
numpy>=1.17.0
scipy>=1.3.1
pybedgraph>=0.5.40
click>=7.0
```
### Installation: 
```bash    
# Install from github
git clone https://github.com/c0ver/chia_rep.git    
pip3 install chia_rep/

# Install from pypi
pip3 install chia_rep
```

### Create Input files
With `example/sample_list.txt` containing the following:
```
LHH0048H
LHH0054H
LHH0084H
LHH0086V
...
```
and `data/` containing bedgraph, peak, and loop files
```bash
cd example
python commands.py --help
python commands.py make-pairs --help
python commands.py make-sample-input-file --help

python commands.py make-pairs sample_list.txt

# Assumes (letter case doesn't matter)
# bedgraph file extension: .bedgraph
# peak files extension: .broadpeak
# loop files extension: .cis.be3
# Creates sample_input_file.txt
python commands.py make-sample-input-file sample_list.txt sample_input_file.txt data/
```

    
### Run script
Example script is included in `example/script.py`.
```bash
cd example
python script.py --help

# Example usages
python script.py sample_input_file.txt hg38.chrom.sizes pairs.txt 3000000 5000 chr1
python script.py sample_input_file.txt hg38.chrom.sizes pairs.txt 3000000 5000 all
python script.py sample_input_file.txt hg38.chrom.sizes pairs.txt 3000000 5000 chr1 chr2
```
  
## Testing  
```bash
pytest  # Runs the tests in test/
```

### Documentation
Included in docs/build/html

## Contact
Contact Minji (minji.kim@jax.org) for general questions, and report software issues in the [Issues](https://github.com/TheJacksonLaboratory/chia_rep/issues) page. 
