"""
RMSD Calculation to Reference Models in PyMOL

This script calculates the RMSD between specified reference models and all other loaded objects in the PyMOL session.
For each comparison, the gene name is extracted from the model's name by splitting at the last underscore.
Results are saved to a CSV file with columns: object, gene_name, model, rmsd.

REQUIREMENTS:
- Runs in a PyMOL session with models already loaded and named using the format <gene_name>_<ID>.
- tqdm must be installed: pip install tqdm

USAGE:
    1. Specify a list of reference model names (reference_names) and your desired atom/residue selection substring (selection_substring).
    2. Set the output CSV file path (csv_path)
    3. Load the models and only have the ones you want included in the calculations relative to the references enabled. 
    4. Execute the script:
           run /path/to/ligand_pairwise_rmsd.py
    5. The script will save the CSV file in the specified location with RMSD values for all loaded and enabled models.

EXAMPLE OF MAIN EXECUTION SECTION:
    reference_names = ["PYR_14NAD", "PYR_16NAD"]
    selection_substring = " and resn PYR"
    csv_path = "./rmsd_results.csv"
    calculate_rmsd_and_save(reference_names, selection_substring, csv_path)

Output:
- CSV file with columns: object, gene_name, model, rmsd
- Each row contains the RMSD of one reference object against another session object, the extracted gene name, and the model names compared.

"""
import csv
import os
import sys
from tqdm import tqdm
#selection_substring and suffix will be global variables set by piecewise_pymol_analysis.py
# --- User configuration ---
reference_names = ["PYR_14NAD", "PYR_16NAD"]
selection_substring = " and resn PYR"  # Modify as needed or empty string for whole object
#suffix is defined in piecewise_pymol_analysis.py
csv_path_prefix = f"./rmsd_to_reference"
# --------------------------


def calculate_rmsd_and_save(object_names, substring, csv_filename_prefix):
    all_names = cmd.get_names("objects", 1)
    cmd.remove("hydrogens")
    prefix = csv_filename_prefix

    for ref_obj in object_names:
        results = []
        for comp_obj in tqdm(all_names, desc=f"Comparing to reference {ref_obj}", leave=False):
            ref_obj_name = ref_obj.rsplit('_', 1)[0]
            cmd.cealign(f"model {ref_obj_name}", f"model {comp_obj}")  # Align structures before RMSD calculation
            try:
                rmsd_value = cmd.rms_cur(comp_obj + substring, ref_obj + substring)
            except Exception:
                rmsd_value = None
            results.append([comp_obj.rsplit('_', 1)[0], comp_obj, rmsd_value])

        output_file = f'{prefix}_{ref_obj}{suffix}.csv'
        print( f"Writing results to {output_file}" )

        # Open file in append mode and check if it is new/empty
        write_header = not os.path.exists(output_file) or os.stat(output_file).st_size == 0

        with open(output_file, mode='a', newline='') as file:
            writer = csv.writer(file)
            if write_header:
                writer.writerow(['gene_name', 'model', 'rmsd'])  # Write header only if new/empty
            writer.writerows(results)  # Append data rows

calculate_rmsd_and_save(reference_names, selection_substring, csv_path_prefix)
