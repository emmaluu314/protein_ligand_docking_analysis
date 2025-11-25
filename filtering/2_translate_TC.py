import os
import pandas as pd
import numpy as np

# User defined variables
current_username = 'emmaluu'
remote_server_name = 'cacao.genomecenter.ucdavis.edu'
main_server_path = "/share/siegellab/bryan/Rare_Sugars/4NP_Docking/" # main location of subfolders
main_local_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/Training/Bryan/2025_07_24/'
results_folder_name = "results/" #this is what your results folder is called. If you do not use a results folder, write "".
target_list_filename = 'AllRunsTogether.csv'
filetag = ".pdb"
script_name = 'pull_pdb_files.sh'
pymol_script_name = 'load_pdbs.py'

sugar_flag = True # This is true if you need to be consistent with Tim
#Change these if sugar_flag is false
cst_filter = 25
total_percent = 0.1
interface_percent = 0.2

# Change these if sugar_flag is true
SR1_cutoff = 5 # Tim did 5
SR2_cutoff = 5 # Tim did 5
SR3_cutoff = 10 # Tim did 1
SR4_cutoff = 2 # Tim did 1

def translate_to_TC(main_path, folder):
    listdata = []

    for i in os.listdir(os.path.join(main_path,folder)): # Reads in all score files (*.sc) within current working directory (CWD).
        if i.startswith('score'):
            try:
                d = pd.read_csv(os.path.join(main_path,folder,i),sep=r'\s+')
                d['indexnum'] = (np.arange(len(d)) % 100) + 1
                listdata.append(d)
            except FileNotFoundError: # This skips missing score files.
                continue

    d = pd.concat(listdata)
    print(f'{folder} has {len(d)} nstructs')

    if sugar_flag:
        #d = d.sort_values(by="total_score")[0:60]
        d = d.sort_values(by="total_score")[0:int(len(d)*0.4)]
        print(f'{len(d)} got sorted by total_score and moved to the next filter')
        
        d = d[(d['SR_1_all_cst'] < SR1_cutoff) & (d['SR_2_all_cst'] < SR2_cutoff) & (d['SR_3_all_cst'] < SR3_cutoff) & (d['SR_4_all_cst'] < SR4_cutoff) ]
        #d = d[(d['SR_1_all_cst'] < SR1_cutoff) & (d['SR_2_all_cst'] < SR2_cutoff) & (d['SR_3_all_cst'] < SR4_cutoff) ]
        print(f'{len(d)} passed the all_cst filters')
        
        interf_column_name = "".join(d.columns[d.columns.str.contains('interf_E_1_3')])
        d = d.sort_values(by = interf_column_name)[0:10] # I forgot if this was commented out or not, but I did now
        # d = d.sort_values(by = 'SR_5_interf_E_1_3')[0:10]

        SR_column_names = list(d.columns[d.columns.str.contains('all_cst')]) + list(d.columns[d.columns.str.contains('interf')])
        filtered_df = d[['total_score']+SR_column_names+['description' , 'indexnum']]
        #filtered_df = d[['total_score','SR_1_all_cst','SR_2_all_cst','SR_3_all_cst','SR_4_all_cst','SR_5_interf_E_1_3','SR_5_all_cst','SR_6_interf_E_1_2','SR_6_all_cst','description' , 'indexnum']]
        print(f'{len(filtered_df)} passed the filters')

    else:
        cst_columns = d.columns[d.columns.str.contains('_all_cst')]
        mask = d[cst_columns].le(cst_filter).all(axis=1)
        filtered_df = d[mask]
        print(f'{len(filtered_df)} passed the filter of {cst_filter} in all the *_all_cst columns')
        
        total_score_filter = round(len(filtered_df)*total_percent)
        filtered_df = filtered_df.sort_values(by="total_score")[0:total_score_filter]
        print(f'{len(filtered_df)} passed the filter of {total_score_filter} in the total_score column')

        interface_score_filter = round(len(filtered_df)*interface_percent)
        interface_column = d.columns[d.columns.str.contains('interf')][0]
        filtered_df = filtered_df.sort_values(by=interface_column)[0:interface_score_filter]
        print(f'{len(filtered_df)} passed the filter of {interface_score_filter} in the {interface_column} column')

    filtered_df.to_csv(os.path.join(main_local_path,folder, target_list_filename))

# Read the filtered CSV file
def create_scripts(folder):
    # Process each row to generate filenames and write to scripts
    counter = 0
    for _, row in data.iterrows():
        description = row['description']
        indexnum = row['indexnum']
        suffix = f"{indexnum + 1:04d}"  # Convert indexnum+1 to 4-digit string
        output_identifier = f"{description}_{suffix}"
        filepath = f"{output_identifier}{filetag}"
        
        # Add to pull_pdb_files.sh
        with open(os.path.join(main_local_path, script_name), 'a') as pull_file:
            server_path = main_server_path  # Adjust if you need per-folder paths
            if results_folder_name:
                server_path = server_path + folder + '/' +results_folder_name
            pull_file.write(f'get {server_path}{filepath} "{main_local_path}{folder}/{filepath}"\n')
        
        # Add to load_pdbs.py
        with open(os.path.join(main_local_path, pymol_script_name), 'a') as pymol_file:
            # For pymol, you may want a unique object name; here, using the file identifier
            pymol_file.write(f'cmd.load("{main_local_path}{folder}/{filepath}","{folder}_{counter}")\n')
        counter += 1

if __name__ == '__main__':
    # Clear and initialize scripts
    with open(os.path.join(main_local_path, script_name), 'w') as pull_file:
        pull_file.write('#!/bin/bash\n')
        pull_file.write(f'sftp -q {current_username}@{remote_server_name} << EOF\n')
    
    if os.path.isfile(os.path.join(main_local_path, pymol_script_name)):
        os.remove(os.path.join(main_local_path, pymol_script_name))

    # going through the folderss    
    for subfolder in os.listdir(main_local_path):
        if os.path.isdir(os.path.join(main_local_path,subfolder)):
            translate_to_TC(main_local_path,subfolder)
            data_unsorted = pd.read_csv(os.path.join(main_local_path,subfolder, target_list_filename))
            data = data_unsorted.sort_values(by='total_score')
            create_scripts(subfolder)
    
        # Close the sftp script
    with open(os.path.join(main_local_path, script_name), 'a') as pull_file:
        pull_file.write('EOF\n')