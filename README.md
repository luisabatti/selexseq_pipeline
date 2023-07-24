# SELEX-seq pipeline
A pipeline to process SELEX-seq data and generate a list of enriched motifs.


## Requirements
This script requires the following packages:

### Python
1. argparse
2. numpy
3. pandas
4. eme_selex
5. tqdm

### R
1. tidyverse
2. Biostrings
3. msa
4. ggseqlogo
5. ips
6. argparse

### System-wide
1. [MAFFT](https://mafft.cbrc.jp/alignment/software/)


## Usage
This pipeline requires the following file types:

1) .fasta.gz: These should be the SELEX-seq library files to be analyzed. These should have been previously trimmed according to their insert size. If using paired-end sequencing, merge the two ends before analysis. Files should be placed within a folder with the same name as the library. 

2) selexseq_files_library.csv: a .csv file with information about the SELEX-seq fastq files to be analyzed. It should include the following columns:
	- sample_name: Name of the library. **It needs to match the folder and name of the fastq files to be analyzed.**
	- protein: The name of the protein to be analyzed. For the R0 control, it should be named **none**.
	- replicate: The replicate number for each protein entry.
	- cycle: The cycle number of each library.

## Settings
The selex_python.py script takes the following arguments:
- --input_file: Path to .csv file containing library information. For example: selexseq_files_library.csv;
- --input_extension: Extension of the library files. For example: .fasta.gz
- --input_folder: Path to input folder where library files are stored.
- --output_folder: Path to output folder where processed files will be saved.
- --kmer_min: Minimum kmer size for the analysis. Range: 1-9. Default: 7.
- --kmer_max: Maximum kmer size for the analysis. Range: 1-9. Default: 9. 
- --align_mode: Algorithm to align motifs. Either "MAFFT", or "Clustal".

## Output
This pipeline will generate the following file types:

1) eme_selex_py_results_kmer_k.bz2: This file will store the fragment calculations from all libraries within each kmer.

2) kmer_k_selex_py.csv: This table will contain the results for each protein and kmer size. It will include information such as the motif sequence and fold change above R0.

3) kmer_k_selex_py_top10_motif.txt: This file will contain the consensus motif for each protein and kmer size.

4) kmer_k_selex_py_top10_motifLogo.pdf: This plot will represent the motif logo for each protein and kmer size.


## Credits
Written by Luis Abatti