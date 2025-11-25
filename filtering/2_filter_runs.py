import os
import warnings
import score_functions as sf

warnings.filterwarnings("ignore")

# === User-defined global settings (edit as needed) ===
current_username = 'emmaluu'
remote_server_name = 'hive.hpc.ucdavis.edu'
main_server_path = "/quobyte/jbsiegelgrp/emmaluu/standard_tools/docking/cofactors/NAD_16/LDH/runs/" 
main_local_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/2025_10_30_16NAD_mutants/'
results_folder_name = "results"
filetag = ".pdb"
check_distribution = True  # Whether to generate diagnostic plots of score distributions


# Filtering & output options:
docking_flag = True
fetch_all_flag = True
results_folder_flag = True
filtered_number_flag = True
pymol_script = True
all_cst_flag = False

amount_requested = 5
script_name = 'pull_pdb_files.sh'
pymol_script_name = 'load_pdbs.py'
cst_name = 'all_cst'

index_order = [1,0]  # [interface energy, total_score]
constraint_threshold = 1
interface_energy_threshold_percent = 0.2
total_threshold_percent = 0.1

def create_request_script(dataframe, request_flag, request_number, server_path, local_path,
                          main_path, iteration_name, file_tag):
    """
    Generates an SFTP request script and writes out filtered results into session and summary files.
    """
    if request_flag:
        file_suffix = 'all'
    elif request_number <= len(dataframe):
        file_suffix = 'truncated'
    else:
        file_suffix = 'not_filtered'
    
    # Write SFTP commands for pulling filtered PDBs
    request_path = os.path.join(main_path, script_name)
    with open(request_path, 'a') as request_file:
        for name_index in range(len(dataframe['description'])):
            output_identifier = dataframe['description'].iloc[name_index]
            request_file.write(f"\nget {server_path}{output_identifier}{file_tag} \"{local_path}{output_identifier}{file_tag}\"")
    # Export the sorted DataFrame for record-keeping    
    filtered_txt_path = os.path.join(local_path, f'filtered_pdbs_{file_suffix}.txt')
    if docking_flag:
        dataframe.sort_values("total_score").to_csv(filtered_txt_path,
                        columns=["description", "total_score", "all_cst", "total_interface_energy"],
                        index=None, sep='\t', mode='w')
    else:
        dataframe.sort_values("total_score").to_csv(filtered_txt_path,
                        columns=["description", "total_score"],
                        index=None, sep='\t', mode='w')
    # Add summary info for PyMOL automation
    with open(filtered_txt_path, 'a') as summary_file:
        summary_file.write('order ')
        summary_file.write(' '.join(dataframe.sort_values("total_score")["description"]))
        summary_file.write(' , no')
        summary_file.write(f'\nshared_prefix = \"{iteration_name}_\"')
        summary_file.write('\nnames = cmd.get_names(\"objects\",1); for index in range(len(names)): cmd.set_name('
                           'names[index],shared_prefix+str(index))')
        summary_file.write('\ndisable all')

def filter_runs():
    """
    Executes filtering over subfolders, generates request scripts, and runs diagnostics plotting.
    """
    with open(os.path.join(main_local_path, script_name), 'w') as pull_file:
        pull_file.write('#!/bin/bash\n')
        pull_file.write(f'sftp -q {current_username}@{remote_server_name} << EOF')
    # For each subfolder, filter, script, and generate diagnostics
    for folder in os.listdir(main_local_path):
        subfolder_server_path = os.path.join(main_server_path, folder, '')
        if results_folder_flag:
            subfolder_server_path = os.path.join(subfolder_server_path, results_folder_name, '')
        current_path = os.path.join(main_local_path, folder, '')
        if os.path.isdir(current_path):
            print(folder)
            score_dataframe = sf.compile_score_data_frame(current_path, docking_flag)
            if filtered_number_flag:
                print(str(len(score_dataframe))+" total pdbs")
            filtered_dataframe = sf.filter_rosetta_pdbs(
                score_dataframe, docking_flag, index_order, all_cst_flag, cst_name, 
                constraint_threshold, interface_energy_threshold_percent, total_threshold_percent,
                folder, filtered_number_flag, current_path,check_distribution)
            create_request_script(
                filtered_dataframe, fetch_all_flag, amount_requested, subfolder_server_path,
                current_path, main_local_path, folder, filetag)
            if pymol_script:
                filtered_dataframe = filtered_dataframe.sort_values("total_score").reset_index(drop=True)
                with open(os.path.join(main_local_path, pymol_script_name), 'a') as pymol_file:
                    for index in range(len(filtered_dataframe)):
                        current_file = filtered_dataframe["description"][index]
                        pymol_file.write(f'cmd.load("{current_path}{current_file}.pdb","{folder}_{index}")\n')
    # Finish SFTP script
    with open(os.path.join(main_local_path, script_name), 'a') as pull_file:
        pull_file.write('\nEOF')
    # Save filter settings for reproducibility
    with open(os.path.join(main_local_path, 'filter_settings.txt'), 'w') as setting_file:
        constraint_version = 'all_cst' if all_cst_flag else 'individual all_cst'
        filter_order_names = ['interface energy', 'total score']
        setting_file.write(f'Filter order: {constraint_version}, {filter_order_names[index_order[0]]}, {filter_order_names[index_order[1]]}\n')
        setting_file.write('Specific settings:\n')
        setting_file.write(f'\tConstraint threshold: {constraint_threshold}\n')
        setting_file.write(f'\tInterface energy % of lowest: {interface_energy_threshold_percent*100}%\n')
        setting_file.write(f'\tTotal score % of lowest: {total_threshold_percent*100}%\n')

if __name__ == '__main__':
    # Remove old PyMOL script for a clean run
    if pymol_script:
        pymol_path = os.path.join(main_local_path, pymol_script_name)
        if os.path.isfile(pymol_path):
            os.remove(pymol_path)
    filter_runs()
