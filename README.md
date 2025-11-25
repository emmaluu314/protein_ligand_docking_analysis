Protein–Ligand Docking Analysis
================================

Python utilities for pulling, filtering, and analyzing Rosetta ligand docking runs. The scripts cover score extraction, constraint/geometry measurements, contact frequency analyses, RMSD comparisons to references, and publication-ready plotting helpers.

## Environment
- Python `>=3.12`
- uv for dependency management (`pipx install uv` or the official installer)
- PyMOL available on the `PATH` for scripts that import `cmd` (contact/constraint scoring, RMSD utilities)
- SFTP access to your HPC cluster if you use the pull/filter helpers

Setup:
```bash
uv sync
# run any script via
uv run python <script>.py
```

## Repository map (what each script does)
- `compile_ligand_docking_results.py` – parses Rosetta `.cst` files, measures template↔motif distances in `.pdb` outputs, compiles constraint and score CSVs, and can plot histograms.
- `compile_cst_score.py` – PyMOL-driven scoring of distance/angle/dihedral constraints across many PDBs; writes `<subfolder>_cst_scores.csv`.
- `compile_contacts.py` – PyMOL selection of residues within a sphere of a target residue/resn; writes `bonds_dataframe*.csv`.
- `analyze_bond_frequency.py` – aggregates `bonds_dataframe*.csv` to contact-frequency tables and barplots by gene/group.
- `filtering/1_extract_score_files.py` – builds an SFTP script to pull `.sc` score files from a cluster (`pull_score_files.sh`).
- `filtering/2_filter_runs.py` – ranks docked PDBs by constraint score/interface energy, produces pull scripts for PDBs, optional PyMOL loader scripts, and diagnostic plots.
- `filtering/2_translate_TC.py` – helper for translating Rosetta TC outputs.
- `filtering/score_functions.py`, `pdb_functions.py`, `user_functions.py` – helpers used by the filtering pipeline.
- `compare_to_reference/` – RMSD utilities (PyMOL-based) for comparing docked structures to reference PDBs and for piecewise analyses.
- `plotting_functions/` – standalone plotting helpers for histograms and grouped/individual box-and-whisker plots.
- `custom_analyses/` – one-off analyses (e.g., best correlation, focused energy analysis).
- `download_alphafold_structures.py` – downloads AlphaFold `.cif` files for a list of UniProt IDs.
- Notebooks: `compile_ligand_docking_results.ipynb`, `generate_polar_contact_map.ipynb`, `plotting_functions/plot_histograms.ipynb`, `translate_constraints.ipynb` for interactive exploration/figure-making.
- `terminal_pymol.sh` – convenience launcher for PyMOL in a terminal session.

## Typical workflow
1) **Pull score files from the cluster**  
   - Edit user variables at the top of `filtering/1_extract_score_files.py` (server, paths, target list).  
   - `uv run python filtering/1_extract_score_files.py` → generates `pull_score_files.sh`; run it to download `.sc` files.

2) **Filter runs and prepare downloads**  
   - Edit thresholds and flags at the top of `filtering/2_filter_runs.py` (constraint threshold, interface/total score %, whether to pull all PDBs, etc.).  
   - `uv run python filtering/2_filter_runs.py` → writes `pull_pdb_files.sh`, `filtered_pdbs_*.txt`, optional `load_pdbs.py`, and diagnostic plots per folder.

3) **Compile docking results with constraint measurements**  
   - In `compile_ligand_docking_results.py`, set `originalWorkingDirectory` (root of docking subfolders) and `suffix` for outputs.  
   - `uv run python compile_ligand_docking_results.py` → emits `[folder][suffix]_scoreFile.csv` and `[folder][suffix]_measurements.csv`; can also make histograms.

4) **Contact/residue frequency analyses**  
   - Set `main_folder_path`, `code`/`queries`, radius, and selection mode in `compile_contacts.py`; run under PyMOL (`pymol -cq compile_contacts.py`).  
   - Run `uv run python analyze_bond_frequency.py` to convert contact lists into frequency tables and barplots (configure `gene_groups`, `WT_string`, flags).

5) **Reference comparisons**  
   - Use `compare_to_reference/` scripts to compute RMSDs vs. reference PDBs or across docking subsets (edit paths and selections in each script; run through PyMOL when they import `cmd`).

6) **Plotting**  
   - `plotting_functions/*.py` can be run with `uv run python plotting_functions/<script>.py` after pointing them at the compiled CSVs from steps 3–4.

7) **Extra utilities**  
   - `download_alphafold_structures.py`: set `downloadPath` and `uniProtIDs_path`, then `uv run python download_alphafold_structures.py`.

## Data layout assumptions
- Rosetta outputs live in per-run subfolders containing `.pdb`, `.sc`, and `.cst` files; many scripts assume a `results/` subfolder (see `results_folder_name` flags).
- Paths to HPC and local directories are hard-coded at the top of most scripts—update them before running.
- PyMOL-based scripts expect to be executed in a PyMOL session (`pymol -cq script.py`) so that `cmd` is available.

## Conventions when running
- Prefer `uv run python <path/to/script>.py` so dependencies resolve via uv.
- For batch jobs on HPC, generate SFTP pull scripts locally (`pull_score_files.sh`, `pull_pdb_files.sh`) and execute them where you have network access.
- Outputs are written alongside the input folders (CSV summaries, plots, and pull scripts); check the printed paths per script.
