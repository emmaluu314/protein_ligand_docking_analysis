"""
Piecewise Batch Loader and Analysis Framework for PyMOL

This script automates the loading and analysis of large collections of PDB files in manageable batches within a PyMOL session.
It is designed for workflows where loading all structures simultaneously is impractical due to performance or memory constraints.

How it works:
- Reads a provided load_pdbs.py file containing lines of PyMOL cmd.load commands.
- Loads a user-specified number of structures at a time ("batch size").
- For each batch:
    - Loads the structures.
    - Sequentially executes each analysis script specified in the 'analysis_scripts' list.
    - Deletes all loaded objects before moving to the next batch.
- Continues until all structures have been processed.

To use:
1. Set 'load_pdbs_path' to your load_pdbs.py file.
2. List all desired analysis scripts in the 'analysis_scripts' variable.
3. Adjust 'batch_size' to the desired number of structures per batch.
4. Run this script inside PyMOL (e.g., 'run piecewise_pymol_analysis.py').

This approach enables reproducible, scalable, and memory-efficient PyMOL structural analyses for large datasets.
"""

import re
import csv
import os
from tqdm import tqdm

# --- User configuration ---
load_pdbs_path = "load_pdbs.py"  # Path to your load script
analysis_script_path = "/Users/emluu/Documents/Siegel lab/Scripts/protein_ligand_docking_analysis/compare_to_reference/"
analysis_scripts = [
    "calculate_rmsd_to_lowest_filtered_pdb.py",
    "calculate_ligand_pairwise_rmsd.py",
    "calculate_rmsd_to_reference.py"
    "compile_polar_contacts.py"
    # Add as many as you like and specifying the full or relative path
]
suffix = '_entire_ligand'  # Suffix to append to output CSV filenames
selection_substring = " and resn PYR"  # Modify as needed or empty string for whole object
# --------------------------

def extract_variant_prefix(obj_name):
    """
    Extract variant prefix by stripping trailing underscore and numeric suffix.
    If no numeric suffix, return full object name.
    """
    m = re.match(r"^([^_]+)", obj_name)
    if m:
        return m.group(1)
    else:
        return obj_name

def group_load_lines_by_variant(load_pdbs_path):
    """
    Reads load_pdbs.py lines of the form: cmd.load(<pdb_path>, "<object_name>")
    Groups the lines by 'variant prefix' (object_name stripped of last _rank suffix if numeric).
    Uses standard dict and manual checks instead of defaultdict.
    """
    with open(load_pdbs_path) as f:
        load_lines = [line.strip() for line in f if line.strip().startswith('cmd.load')]

    variant_groups = {}
    for line in load_lines:
        match = re.search(r'cmd\.load\([^)]+,\s*"([^"]+)"\)', line)
        if not match:
            continue
        obj_name = match.group(1)
        variant_prefix = extract_variant_prefix(obj_name)

        if variant_prefix not in variant_groups:
            variant_groups[variant_prefix] = []
        variant_groups[variant_prefix].append(line)

    return variant_groups

def piecewise_by_variant(load_pdbs_path, analysis_scripts):
    variant_groups = group_load_lines_by_variant(load_pdbs_path)

    sorted_variants = sorted(variant_groups.keys())
    total = len(sorted_variants)
    tqdm_bar = tqdm(enumerate(sorted_variants), total=total)
    cmd.do(f'selection_substring = "{selection_substring}"')
    cmd.do(f'suffix = "{suffix}"')

    for i, variant in tqdm_bar:
        # Check if the folder exists for variant
        if not os.path.isdir(variant):
            tqdm_bar.set_description(f"Skipping missing folder: {variant}")
            continue

        tqdm_bar.set_description(f"Processing variant {variant}")
        lines = variant_groups[variant]
        # Load all objects for this variant group
        for line in lines:
            exec(line, globals())  # runs cmd.load in PyMOL context

        # Run all specified analysis scripts on the currently loaded group
        for script in analysis_scripts:
            cmd.do(f'run {os.path.join(analysis_script_path, script)}')

        # Delete all loaded objects of this group
        obj_names = [
            re.search(r'cmd\.load\([^)]+,\s*"([^"]+)"\)', l).group(1)
            for l in lines
        ]
        for obj in obj_names:
            cmd.delete(obj)
        tqdm_bar.update(1)

# --- Execute once loaded in PyMOL ---
csv_files_and_headers = {
    f'rmsd_vs_lowest_filtered_pdb{suffix}.csv': ['gene_name', 'reference_model', 'model', 'rmsd'], # for calculate_rmsd_to_lowest_filtered_pdb.py
    f'rmsd_all_to_all_pairwise{suffix}.csv': ['gene_name', 'model_1', 'model_2', 'rmsd'], # for calculate_ligand_pairwise_rmsd.py
    f'rmsd_to_reference_PYR_14NAD{suffix}.csv': ['gene_name', 'model', 'rmsd'], # for calculate_rmsd_to_reference.py
    f'rmsd_to_reference_PYR_16NAD{suffix}.csv': ['gene_name', 'model', 'rmsd'], # for calculate_rmsd_to_reference.py
    f'compiled_polar_contacts{suffix}.csv': ['gene_name', 'model', 'contact_list', 'number_of_contacts'], # for compile_polar_contacts.py
}

for csv_path, header in csv_files_and_headers.items():
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)
piecewise_by_variant(load_pdbs_path, analysis_scripts)

