from Bio.PDB import PDBParser, Superimposer
import numpy as np
import csv
import os
import tqdm

def get_atoms(structure, resname_filter=None, atom_names=None):
    """
    Collect atoms from the structure, filtered by residue name and optionally atom names.

    Args:
        structure: Biopython Structure object.
        resname_filter: Single residue name string or list/set of residue names to include, or None to include all.
        atom_names: Single atom name string or list/set of atom names to include, or None to include all.

    Returns:
        List of Bio.PDB.Atom objects matching filters, excluding hydrogens.
    
    Example:
        atoms = get_atoms(structure, resname_filter=['PYR'], atom_names=['C2', 'C3'])
    """

    if isinstance(resname_filter, str):
        resname_filter = {resname_filter}
    elif resname_filter is not None:
        resname_filter = set(resname_filter)

    if isinstance(atom_names, str):
        atom_names = {atom_names}
    elif atom_names is not None:
        atom_names = set(atom_names)

    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname_filter and resname not in resname_filter:
                    continue
                for atom in residue:
                    if atom.element.strip() == 'H':  # Exclude hydrogens
                        continue
                    if atom_names and atom.get_name() not in atom_names:
                        continue
                    atoms.append(atom)
    return atoms

from Bio.PDB import PDBParser, Superimposer
import numpy as np
import csv
import os
from tqdm import tqdm

def calculate_rmsd_and_save(ref_path, comp_paths, resname_filter=None, atom_names=None, csv_filename='temp.csv'):
    """
    Calculate RMSD of comparison PDBs to a reference PDB with filtering options.

    Args:
        ref_path (str): Path to reference PDB file.
        comp_paths (list[str]): List of paths to comparison PDB files.
        resname_filter (str or list or set, optional): Residue names to include (e.g. "PYR").
        atom_names (str or list or set, optional): Atom names to include (e.g. "C1" or ["C1", "N2"]).
        csv_filename_prefix (str): Prefix path for output CSV files.
        suffix (str): Suffix to append to CSV filename.
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure('reference', ref_path)
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]

    results = []
    for comp_path in tqdm(comp_paths, desc=f"Comparing to reference {ref_name}"):
        comp_struct = parser.get_structure('comp', comp_path)
        comp_name = os.path.splitext(os.path.basename(comp_path))[0]

        # Get atom lists with filtering
        ref_atoms = get_atoms(ref_struct, resname_filter, atom_names)
        comp_atoms = get_atoms(comp_struct, resname_filter, atom_names)

        n = min(len(ref_atoms), len(comp_atoms))
        if n == 0:
            rmsd_value = None  # No atoms to compare
        else:
            atoms1 = np.array([atom.get_coord() for atom in ref_atoms[:n]])
            atoms2 = np.array([atom.get_coord() for atom in comp_atoms[:n]])

            si = Superimposer()
            try:
                si.set_atoms(atoms1, atoms2)
                si.apply(comp_struct.get_atoms())  # Apply superposition to comp structure atoms (optional)
                rmsd_value = si.rms
            except Exception:
                rmsd_value = None

        gene_name = comp_name.rsplit('_', 1)[0]  # Parse gene name from filename pattern
        results.append([gene_name, comp_name, rmsd_value])

    write_header = not os.path.exists(csv_filename) or os.stat(csv_filename).st_size == 0
    with open(csv_filename, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(['gene_name', 'model', 'rmsd'])
        writer.writerows(results)
