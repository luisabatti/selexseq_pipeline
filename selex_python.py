import argparse
import os
import shutil
import sys
import numpy as np
import pandas as pd
import pickle
import bz2
import subprocess
from pathlib import Path

from collections import defaultdict
from eme_selex.sequence import canonical
from eme_selex.eme import kmer_fraction_from_file as kf
from tqdm.auto import tqdm

# Create an ArgumentParser object
parser = argparse.ArgumentParser(prog='selex_python', description='This script will process selex seq datasets.')

# Add arguments
parser.add_argument('--input_file', default = False, help='path to .csv file with information')
parser.add_argument('--input_extension', default = '.fasta.gz', help='Extension of the input files')
parser.add_argument('--input_folder', default = False, help='Path to input folder')
parser.add_argument('--output_folder', default = False, help='Path to output folder')
parser.add_argument('--kmer_min', default = 7, help = 'Minimum kmer size')
parser.add_argument('--kmer_max', default = 9, help = 'Maximum kmer size')
parser.add_argument('--align_mode', default = "MAFFT", choices = ['MAFFT', 'Clustal'], help = 'Algorithm to align motifs. Either "MAFFT", or "Clustal"')

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the arguments
input_file = args.input_file
input_extension = args.input_extension
input_folder = args.input_folder
output_folder = args.output_folder
kmer_min = args.kmer_min
kmer_max = args.kmer_max
align_mode = args.align_mode

#Check if the inputted arguments are valid

if os.path.isfile(input_file) == False:
    sys.exit('Incorrect input file')
elif os.path.isdir(input_folder) == False:
    sys.exit('Incorrect input folder')
elif os.path.isdir(output_folder) == False:
    sys.exit('Incorrect output folder')
else:
    print(f'Processing {input_file}. Reading files from {input_folder}. Outputting files to {output_folder}')

counts, fractions, models = defaultdict(dict), defaultdict(dict), defaultdict(dict)

#List of kmers to work with
kmer_list = tuple(range(int(kmer_min), int(kmer_max) + 1))

#Reads in input_file and transform to table
df = pd.read_csv(input_file, usecols = ['sample_name', 'protein', 'replicate', 'cycle']).drop_duplicates().fillna('None').reset_index(drop = True)

#Gets all the sample names
samples = df["sample_name"].values
print(f"Working with the following samples: {samples}")

#Gets all the unique protein names (except to None which will be the control library)
selex_proteins = df["protein"].drop_duplicates().tolist()
selex_proteins.remove('None')
#List of unique proteins in the metadata file

print(f"Working with the following proteins: {selex_proteins}")


#For loop to calculate and dump a .bz2 file with each kmer size
for k in kmer_list:

    pickle_file_path = Path(f"{output_folder}/eme_selex_py_results_kmer_{k}.bz2")

    print("Processing k number",k)

    if os.path.isfile(pickle_file_path):
        print(f'{pickle_file_path} already exists, skipping pickle dump')
        pass
    else:
        print("Calculating kmer",k)
        for sample in tqdm(samples, leave=False):
            c, f, m = kf(f"{input_folder}/{sample}/{sample}{input_extension}", k=k)
            counts[sample] = c
            fractions[sample] = f
            models[sample] = m

        with bz2.open(f"{pickle_file_path}", "wb") as wH:
            pickle.dump(counts, wH)
            pickle.dump(fractions, wH)
            pickle.dump(models, wH)

        #Saves files with kmer data as a pickle file

#Creates dictionary with each replicate and control libraries
_cycle0_libs = df[df["cycle"]==0][["sample_name", "replicate"]].set_index("replicate")["sample_name"].to_dict()
cycle0_libs = {}
for i, j in df[["replicate", "sample_name"]].iterrows():
    l, s = j.values
    cycle0_libs[s] = _cycle0_libs[l]

#For loop to process each kmer size
for k in kmer_list:

    print(f'Running script for kmer: {k}')

    pickle_file_path=(f"{output_folder}/eme_selex_py_results_kmer_{k}.bz2")

    if not os.path.isfile(pickle_file_path):
        raise FileNotFoundError(f"No such file: '{pickle_file_path}'")
    else:
        print(f"Opening {pickle_file_path}...")

    with bz2.open(f"{pickle_file_path}", "rb") as iH:
        counts = pickle.load(iH)
        _fractions = pickle.load(iH)
        models = pickle.load(iH)

    # keep only canonical k-mers in the fractions container
    fractions = defaultdict(dict)
    fold_change = defaultdict(dict)
    for _i, _j in _fractions.items():
        for _k, _l in _j.items():
            fractions[_i][canonical(_k)[0]] = _l

    for _i, _j in fractions.items():
        for _k, _l in _j.items():
            if _i in cycle0_libs.keys():
                fold_change[_i][_k] = _l/fractions[cycle0_libs[_i]][_k]

    # create dict for easy lookup of sample metadata
    protein_dict = df.set_index('sample_name')["protein"].to_dict()
    cycle_dict = df.set_index('sample_name')["cycle"].to_dict()

    def get_at(x):
        return x.count("A") + x.count("T")

    # create a pandas data frame for the foldchange container
    fold_change_df = pd.DataFrame.from_dict(fold_change)
    fold_change_df.index.name = "kmer"

    def melt_df(x):
        melt_df = pd.melt(x.reset_index(), id_vars="kmer")
        # add the information about the number of As and Ts in the k-mer
        melt_df["AT"] = melt_df["kmer"].apply(get_at)

        # add metadata
        melt_df["protein"] = melt_df["variable"].apply(lambda x: protein_dict.get(x, None))    
        melt_df["cycle"] = melt_df["variable"].apply(lambda x: cycle_dict.get(x, None))
        melt_mean = melt_df.dropna().groupby(
            ["cycle", "protein", "kmer", "AT"]).mean(numeric_only=True).reset_index()
        melt_std = melt_df.dropna().groupby(
            ["cycle", "protein", "kmer", "AT"]).std(numeric_only=True).reset_index()
        melt_mean["value_std"] = melt_std["value"]

        return melt_mean

    #melt_fractions_mean = melt_df(fractions_df)

    melt_fold_change_mean = melt_df(fold_change_df)

    # Creates a folder and export a .csv table with the results for each protein

    for protein in selex_proteins:

        selex_run = ['None']

        selex_run.append(protein)
        #Select one protein at a time + None to export

        fold_change_mean_df = melt_fold_change_mean[melt_fold_change_mean['protein'].isin(selex_run)].sort_values(by=["value"], ascending=False)

        filepath = Path(f"{output_folder}/{selex_run[1]}/{selex_run[1]}_kmer_{k}_selex_py.csv")  

        filepath.parent.mkdir(parents=True, exist_ok=True)  

        fold_change_mean_df.to_csv(filepath, index=False)

        print(f"Exported dataframe with {selex_run[0]} and {selex_run[1]}")

print('Running R script to generate motif logos')

subprocess.call(["Rscript", "./selex_motifLogo.R", "--input_file", input_file, "--output_folder", output_folder, "--kmer_min", kmer_min, "--kmer_max", kmer_max, "--align_mode", align_mode])

print('Done!')