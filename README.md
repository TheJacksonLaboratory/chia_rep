# chia_rep 
A Python package for assessing the reproducibility of generated ChIA-PET data.    
    
## Method 
### Preprocessing 
Find the specific start and end anchor within given interval.    
- Get the exact index of largest value in interval from corresponding bedgraph file    
- Weigh each loop value by *(avg_peak_value * PET_count)* - Remove peaks that have overlapping start/end peaks (few)  
   - Do not compare area covered by loop with same location from a chromosome from another sample  
    
Filter out loops without peak support on either start or end anchors.    
 - Take *x* percentage of peaks from peak calling algorithm (narrowPeak, broadPeak)  
 - Better method: Take *x* number of peaks for each chromosome in each sample  
 - Filter out loops that don't have peak support on either side  
    
### Comparison
 **Window Size**  
Using a window size filters out extremely long (supposedly noisy) loops. For example, with window size = 3mb, there will be ~80 windows in chr1 (~250,000,000 bp). Each window produces a reproducibility value and a weight dependendent on the largest loop value inside it. Theses values combine into a single reproducibility value for the chromosome which can be further combined for the genome. The final report will be a value between -1 (dissimilar) and 1 (similar).  
  
Windows without loops for either sample will have a weight of 0 and are not included in the final calculation.   
  
**Bin Size**  
A 2D array of non-overlapping bins within each window is created to make a graph representation of the window where each bin is a node and loops are edges. For example. with window size = 3mb, bin size = 10kb, 300x300 bins are created. Bin (0, 0) will hold loops that start/end in the interval 0-10kb. Bin (1, 1) will hold loops that start/end in the interval 10kb-20kb. Bin (1, 0) will hold loops that start in interval 0-10kb and end in 10kb-20kb.  
  
A smaller bin size may not always be better due to replicate samples not having peaks and loops in the exact same position in the chromosome. As a result, replicate loops may be in slightly different bins which may affect the result.  
    
### Example
Each sample generates a graph similar to the following for each window.    
```    
name  start  end  value     
chr1 1  11  5    
chr1  4  21  5    
chr1  14  26  5    
```    
Assuming bin_size = 10:    
    
**Original window**   


|          |  0 - 9   | 10 - 19  | 20 - 29  |    
|--------- |--------  |--------- |--------- |    
| 0 - 9    | 0        | 5        | 5        |    
| 10 - 19  | 0        | 0        | 5        |    
| 20 - 29  | 0        | 0        | 0        |    
    
**Normalized window**   


|          |     |     |    
|--------- |--------  |--------- |    
| 0        | 0.333    | 0.333    |    
| 0        | 0        | 0.333    |    
| 0        | 0        | 0        |    
    
### Graph Creation 
Normally only use top-right of graph, since 2D array representation is symmetric since loops have a defined start and end.  
    
Since loops from replicate samples may be in slightly different bins due to noise, i.e. 9999 vs. 10001, also add loops to nearby bins so they may still partly overlap. This will also make the graph unsymmetrical if only done for one side.   
    
    
### Comparing Graph 
Jensen-Shannon Divergence    
- Flatten 2D array into 1D array    
- Compares the difference between the distributions    
    
Earth Mover's Distance 1D:    
- Similar to Jensen-Shannon but takes into account distributions that are similar in shape    
- Loops/peaks could be in slightly different bins and be heavily penalized in Jensen-Shannon    
- Take EMD of each row and column in the graph since 2D EMD is much more complicated  
  
Convert EMD to a value to be between (-1, 1)  
Since values tend to be closer to 0 than `max_emd_dist`, use a non-linear scaling.  
```  
emd_value = 2 * (emd_dist - max_emd_dist) * (emd_dist - max_emd_dist) / (max_emd_dist * max_emd_dist) - 1  
```  
    
## Problems 
Many replicate data sets have significantly different number of loops/peaks (sequencing depth).  
  
Potential Fix:  
Taking a defined number of peaks should help with creating a similar number of loops.  
    
## Usage: 
### Dependencies:
```python
# Automatically included when pip installing
install_requires=[
	'numpy>=1.17.0',  
	'scipy>=1.3.1',  
	'prettytable>=0.7.2',  
	'pybedgraph>=0.5.31',  
	'matplotlib>=3.1.1']
```
### Installation: 
```bash    
git clone https://github.com/c0ver/chia_rep.git    
pip3 install chia_rep/
```
    
### Input data and compare:  
```python    
from chia_rep import reproducibility  
  
TEST_DATA_DIR = 'data'  # Folder containing all data
BIN_SIZE = 10000  # 10kb
WINDOW_SIZE = 10000000  # 10mb
NUM_PEAKS = 60

# To see log info statements (optional)  
from logging.config import fileConfig
fileConfig('chia_rep.conf')
  
  
loop_dict = reproducibility.read_data(loop_data_dir=TEST_DATA_DIR,
                                      chrom_size_file=f'{TEST_DATA_DIR}/hg38.chrom.sizes',
                                      bedgraph_data_dir=TEST_DATA_DIR,
                                      peak_data_dir=TEST_DATA_DIR)
  
reproducibility.preprocess(loop_dict, num_peaks=NUM_PEAKS)
  
rep, non_rep, emd_scores, j_scores = reproducibility.compare(loop_dict, bin_size=BIN_SIZE, window_size=WINDOW_SIZE)  

# To compare only certain samples (further documentation provided in reproducibility.py)
rep, non_rep, emd_scores, j_scores = reproducibility.compare(loop_dict, bin_size=BIN_SIZE, window_size=WINDOW_SIZE, 
specified_comparisons=[['sampleA1', 'sampleA2'], ['sampleA1', sampleB1']])

reproducibility.output_results(rep, non_rep, 'sample_test')  
reproducibility.output_to_csv(emd_scores, 'sample_test/sample_test.emd_value.csv')  
reproducibility.output_to_csv(j_scores, 'sample_test/sample_test.j_value.csv')
```  
  
## Testing:  
Use `chia_rep/test/test.py` in an interactive session. Further reasoning provided  
inside file.

