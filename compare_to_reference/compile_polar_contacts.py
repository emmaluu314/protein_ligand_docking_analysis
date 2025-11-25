import os
import csv
from collections import defaultdict
from tqdm import tqdm
import sys
#selection_substring and suffix will be global variables set by piecewise_pymol_analysis.py
# --- User configuration ---
csv_path = f"./compiled_polar_contacts{suffix}.csv"
selection_radius=3.5
include_query_in_report = False # temporary while I figure out a way to not include the exact atoms in the query for more nuanced queries
# --------------------------

def analyze_contacts_by_gene(selection_substring, csv_filename, sphere_radius=3.5,include_self_flag=False):
    """
    Analyze potential contacts for all loaded PyMOL objects,
    group by gene_name and output CSV.
    """
    all_objects = cmd.get_names("objects", enabled_only=1)
    parsed = []
    # Parse gene name and id from object names
    for obj in all_objects:
        if not isinstance(obj, str) or '_' not in obj:
            continue
        gene_name, id_part = obj.rsplit('_', 1)
        parsed.append((gene_name, id_part, obj))

    # Group by gene name
    groups = defaultdict(list)
    for gene_name, id_part, obj in parsed:
        groups[gene_name].append((id_part, obj))

    # removing hydrogens
    cmd.remove("hydrogens")
    results = []
    for gene_name, models in groups.items():
        for id_part, object_name in tqdm(models,desc="Analyzing models",disable=True):
            # Compose selection: each object + user-specified substring (e.g. " and resn PYR")
            full_selection = f"model {object_name}{selection_substring} and (donors or acceptors)"

            # Select all atoms within sphere_radius of the selection in the object
            #print("Between:",full_selection, f"model {object_name} and (donors or acceptors)")
            pairs = cmd.find_pairs(full_selection, f"model {object_name} and (donors or acceptors)", cutoff=sphere_radius,mode=1)
            atom_indices = [item[1] for pair in pairs for item in pair]
            atom_indicies_str = [str(index) for index in atom_indices]
            
            # Python/PyMOL code: make space for iterate to use
            stored.sidechain_list = []
            cmd.iterate(f"model {object_name} and index "+"+".join(atom_indicies_str), "stored.sidechain_list.append([chain,resv])")
            #cmd.iterate(f"model {object_name} and index "+"+".join(atom_indicies_str), "print(model,chain,resv,name)")
            sidechain_list = stored.sidechain_list
            unique_list = list(dict.fromkeys(tuple(x) for x in sidechain_list))
            unique_list = [item for item in unique_list if item[0] != ""] # not sure why empty entries for chain are showing up

            if not include_self_flag:
                stored.exclude_chain = ''
                cmd.iterate(full_selection,"stored.exclude_chain = chain")
                unique_list = [item for item in unique_list if item[0] != stored.exclude_chain]

            contact_list = []
            for chain, residue_position in tqdm(unique_list,desc="Double checking distances",total=len(unique_list),disable=True):
                contact_list.append(f'(chain {chain} and resi {residue_position})')

            results.append([gene_name, object_name, contact_list, len(contact_list)])

    # Write/append to CSV. Use the header if the file does not exist or is empty.
    write_header = not os.path.exists(csv_filename) or os.stat(csv_filename).st_size == 0
    with open(csv_filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(['gene_name', 'model', 'contact_list', 'number_of_contacts'])
        writer.writerows(results)

# Example usage within PyMOL standalone or a script:
analyze_contacts_by_gene(selection_substring=selection_substring,
                         csv_filename=csv_path,
                         sphere_radius=selection_radius,
                         include_self_flag = include_query_in_report)


