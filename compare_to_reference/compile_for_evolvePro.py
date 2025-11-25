import os
import glob
import pandas as pd
from Bio import SeqIO
from statistics import median
import numpy as np
import re

def apply_single_mutation(seq, mutation):
    """Applies a single mutation in the format T233A to the sequence string."""
    try:
        original, rest = mutation[0], mutation[1:]
        nums = ''.join([c for c in rest if c.isdigit()])
        position = int(nums) - 1  # zero-based indexing
        new = rest[len(nums):]

        seq_list = list(seq)
        if seq_list[position] != original:
            print(f"Warning: position {position+1} in WT is {seq_list[position]}, not {original} for mutation {mutation}")
        seq_list[position] = new
        return ''.join(seq_list)
    except Exception as e:
        print(f"Error processing mutation {mutation}: {e}")
        return seq

def generate_mutated_sequence(wt_seq, gene_name):
    """Generates a mutated sequence for single or multiple mutations, or returns WT unchanged."""
    # Case 1: 'WT' in name uses unaltered sequence
    if 'WT' in gene_name.upper():
        return wt_seq

    # Case 2: multiple mutations separated by underscores
    mutations = gene_name.split('_')

    # Check every mutation follows valid pattern (e.g., T233A)
    valid_pattern = re.compile(r'^[A-Z]\d+[A-Z]$')
    mutated_seq = wt_seq

    for mut in mutations:
        mut = mut.strip()
        if valid_pattern.match(mut):
            mutated_seq = apply_single_mutation(mutated_seq, mut)
        else:
            print(f"Skipped invalid mutation format: {mut}")
    return mutated_seq

def compile_csv_data(root_folder, fasta_file, input_file, output_file):
    """Compiles all subfolder CSVs into a combined summary with mutation-based sequences."""
    wt_record = next(SeqIO.parse(fasta_file, "fasta"))
    wt_seq = str(wt_record.seq)

    csv_files = glob.glob(os.path.join(root_folder, "**", input_file), recursive=True)
    combined_data = []

    for file in csv_files:
        df = pd.read_csv(file)
        for gene_name, group in df.groupby('gene_name'):
            med_rmsd = median(group['rmsd'])
            std_rmsd = np.std(group['rmsd'], ddof=1)
            mutated_seq = generate_mutated_sequence(wt_seq, gene_name)
            combined_data.append({
                'gene_name': gene_name,
                'sequence': mutated_seq,
                'median_rmsd': med_rmsd,
                'std_rmsd': std_rmsd
            })

    result_df = pd.DataFrame(combined_data)
    result_df.to_csv(output_file, index=False)
    print(f"Compiled CSV saved to {output_file}")

# Example usage:
folder_path = "/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/previous_csv/"
fasta_path = "/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/prep_work/1LDN_WT.fasta"
input_csv = "rmsd_to_reference_PYR_16NAD_entire_ligand.csv"
file_name = "compiled_summary_for_EvolvePro.csv"
compile_csv_data(folder_path, fasta_path,input_csv, os.path.join(folder_path,file_name))
