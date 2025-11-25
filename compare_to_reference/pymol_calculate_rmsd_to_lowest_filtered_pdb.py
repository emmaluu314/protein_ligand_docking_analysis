"""
RMSD Relative to Zero-ID Model within Each gene_name in PyMOL

This script calculates the RMSD of each model relative to the model with ID 0 within each gene_name group in a PyMOL session.
Model names must follow the format <gene_name>_<ID>, where ID ending with '0' is treated as the reference for that gene_name.
For each gene group, the script finds the model with ID '0' and calculates RMSD of every other model against it.
Results are saved in a CSV with columns: gene_name, reference_model, model, rmsd.

REQUIREMENTS:
- Run inside PyMOL with models loaded and named as <gene_name>_<ID>.
- tqdm is recommended for progress display (pip install tqdm).

USAGE:
    1. Ensure the required models are loaded and enabled in PyMOL.
    2. Adjust 'selection_substring' if RMSD is to be calculated on a subset of atoms/residues.
    3. Set 'csv_path' for output.
    4. Run this script in PyMOL: 
           run /path/to/this_script.py
    5. Call the function with no arguments to execute calculations and save results.

EXAMPLE MAIN EXECUTION:
    selection_substring = " and resn PYR"
    csv_path = "./rmsd_vs_zero_id.csv"
    calculate_rmsd_vs_zero_id(selection_substring, csv_path)

OUTPUT:
- CSV file with columns: gene_name, reference_model, model, rmsd
- Each row contains RMSD of a model relative to the zero-ID model in the same gene_name group.
- Models without a zero-ID counterpart are skipped with a warning.

"""
import csv
import os
from collections import defaultdict
import sys

#selection_substring and suffix will be global variables set by piecewise_pymol_analysis.py
# --- User configuration ---
selection_substring = " and resn PYR"  # Modify as needed or empty string for whole object
suffix = "_entire_ligand"
csv_path = f"./rmsd_vs_lowest_filtered_pdb{suffix}.csv"
# --------------------------


def calculate_rmsd_to_lowest_filtered_pdb(substring, csv_filename):
    # Get all loaded objects
    all_objects = cmd.get_names("objects", enabled_only=1)
    
    cmd.remove("hydrogens")
    # Extract gene_name and ID from each object name
    parsed = []
    for obj in all_objects:
        if not isinstance(obj, str) or '_' not in obj:
            continue
        gene_name, id_part = obj.rsplit('_', 1)
        parsed.append((gene_name, id_part, obj))
    
    # Group objects by gene_name
    groups = defaultdict(list)
    for gene_name, id_part, obj in parsed:
        groups[gene_name].append((id_part, obj))
    
    results = []
    for gene_name, models in groups.items():
        # Find zero-ID model for group
        zero_id_models = [obj for id_part, obj in models if id_part == '0']
        if len(zero_id_models) == 0:
            print(f"Warning: No zero-ID model found for gene '{gene_name}', skipping this group.")
            continue
        if len(zero_id_models) > 1:
            print(f"Warning: Multiple zero-ID models found for gene '{gene_name}', choosing first one.")
        ref_model = zero_id_models[0]
        
        # Calculate RMSD of each model relative to zero-ID ref_model
        for id_part, model_obj in models:
            if model_obj == ref_model:
                continue  # Skip RMSD of model to itself
            try:
                cmd.cealign(model_obj, ref_model)  # Align structures before RMSD calculation
                rmsd_value = cmd.rms_cur(model_obj + substring, ref_model + substring)
            except Exception as e:
                print(f"Error calculating RMSD between {model_obj} and {ref_model}: {e}")
                rmsd_value = None
            results.append([gene_name, ref_model, model_obj, rmsd_value])
    
    # Determine file mode: write header only if file does not exist
    write_header = not os.path.exists(csv_filename) or os.stat(csv_filename).st_size == 0
    
    # Write or append to CSV
    with open(csv_filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(['gene_name', 'reference_model', 'model', 'rmsd'])
        writer.writerows(results)

# Example main execution - customize before executing in PyMOL
calculate_rmsd_to_lowest_filtered_pdb(selection_substring, csv_path)
