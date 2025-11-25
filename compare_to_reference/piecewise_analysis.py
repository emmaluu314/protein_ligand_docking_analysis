import os
import csv
from tqdm import tqdm

# User defined variables
project_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/'
main_path = os.path.join(project_path,'2025_10_30_16NAD_mutants/')
folder_list_file = 'jobs.txt'  # text file with one folder name per line
# User specifies which analyses to run (could come from config or CLI args)
selected_analyses = ["rmsd_to_reference_14NAD_entire_ligand", 
                     "rmsd_to_reference_16NAD_entire_ligand"]

# Function definitinos 
# Define all analyses in a dictionary with metadata
analysis_registry = {
    "rmsd_to_reference_14NAD_entire_ligand": {
        "function": calculate_rmsd_and_save,  # your analysis function
        "csv_file": f'rmsd_to_reference_PYR_14NAD_entire_ligand.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_14NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    "rmsd_to_reference_14NAD_middle_carbons": {
        "function": calculate_rmsd_and_save,  # your analysis function
        "csv_file": f'rmsd_to_reference_PYR_14NAD_middle_carbons.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_14NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    "rmsd_to_reference_14NAD_reaction": {
        "function": calculate_rmsd_and_save,  # your analysis function
        "csv_file": f'rmsd_to_reference_PYR_14NAD_reaction.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_14NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    "rmsd_to_reference_16NAD_entire_ligand": {
        "function": calculate_rmsd_and_save, 
        "csv_file": f'rmsd_to_reference_PYR_16NAD_entire_ligand.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_16NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    "rmsd_to_reference_16NAD_middle_carbons": {
        "function": calculate_rmsd_and_save, 
        "csv_file": f'rmsd_to_reference_PYR_16NAD_middle_carbons.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_16NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    "rmsd_to_reference_16NAD_reaction": {
        "function": calculate_rmsd_and_save, 
        "csv_file": f'rmsd_to_reference_PYR_16NAD_reaction.csv',
        "csv_header": ['gene_name', 'model', 'rmsd'],
        "args": {
            "reference_path": [os.path.join(project_path,'prep_work/PYR_16NAD.pdb')],
            "resname_filter": "PYR"
        }
    },
    # Add other analyses here similarly
}

def build_variant_groups_from_file(folder_list_path):
    variant_groups = {}
    with open(folder_list_path, 'r') as file:
        folders = [line.strip() for line in file if line.strip()]  # Read and strip lines

    for folder in folders:
        # Verify the folder exists
        if not os.path.isdir(folder):
            print(f"Warning: folder '{folder}' does not exist, skipping.")
            continue
        # List PDB files in folder
        pdb_files = [os.path.join(folder, f) for f in os.listdir(folder) if f.lower().endswith('.pdb')]
        if pdb_files:
            variant_groups[folder] = pdb_files
        else:
            print(f"No PDB files found in folder '{folder}', skipping.")

    return variant_groups

# Main Execution
variant_groups = build_variant_groups_from_file(folder_list_file)
# Prepare all selected csv files: write their headers once
for name in selected_analyses:
    csv_file = analysis_registry[name]["csv_file"]
    header = analysis_registry[name]["csv_header"]
    with open(csv_file, 'w', newline='') as f:
        csv.writer(f).writerow(header)

# Run selected analyses on each variant group
for variant, pdb_paths in variant_groups.items():
    for analysis_name in selected_analyses:
        analysis_def = analysis_registry[analysis_name]
        func = analysis_def["function"]
        args = analysis_def["args"].copy()
        args["comp_paths"] = pdb_paths  # set target pdbs dynamically
        func(**args)
