import os
import score_functions as sf

# User defined variables:
current_username = 'emmaluu'
remote_server_name = 'cacao.genomecenter.ucdavis.edu'
main_server_path = '/share/siegellab/emmaluu/docking/single_runs/Han_Li/NMN/ALDH/' # main location of subfolders
results_folder_name = "results_in/" #this is what your results folder is called. If you do not use a results folder, write "".
main_local_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/Unnatural_cofactor/ALDH/Bovine paper/NMN_in/'
target_list_filename = 'master_list.txt' #name of file with a column of the subfolders you'd like to analyze
filetag = ".pdb"
# Boolean flags for 
docking_flag = True # If analyzing relaxation, this is false
fetch_all_flag = True #if false, edit "amount_requested"
results_folder_flag = True # Have this as true if you put all of your outputs in a results subfolder on the server
filtered_number_flag = True # print how many passed each filter
pymol_script = True # create a python script to run in pymol for automating the creation of a session
all_cst_flag = False # If true, will use the compiled cst score for filtering instead of the individual ones of each constraint.


# variables related to the boolean flags
amount_requested = 2  # max number of pdbs you want to be returned if fetch_all_flag is False
script_name = 'pull_pdb_files.sh'
pymol_script_name = 'load_pdbs.py'
cst_name = 'all_cst' # this is based on your .xml script definition

index_order = [1,0] # order of filtering after constraints; 0: interface energy, 1: total_score
constraint_threshold = 1
interface_energy_threshold_percent = 0.2 #0.2
total_threshold_percent = 0.1 #0.1


def create_request_script(dataframe, request_flag, request_number, server_path, local_path,
                          main_path,iteration_name,file_tag):
    if request_flag:
        file_suffix = 'all'
    elif request_number <= len(dataframe):
        file_suffix = 'truncated'
    else:
        file_suffix = 'not_filtered'

    with open(main_path+script_name, 'a') as request_file:
        for name_index in range(len(dataframe['description'])):
            output_identifier = dataframe['description'].iloc[name_index]
            request_file.write("\nget "+server_path+output_identifier+file_tag+' "'+local_path+output_identifier+file_tag+'"')
    if docking_flag:
        dataframe.sort_values("total_score").to_csv(local_path + '/filtered_pdbs_' + file_suffix + '.txt',
                         columns=["description", "total_score", "all_cst", "total_interface_energy"],
                         index=None, sep='\t', mode='w')
    else:
        dataframe.sort_values("total_score").to_csv(local_path + '/filtered_pdbs_' + file_suffix + '.txt',
                         columns=["description", "total_score"],
                         index=None, sep='\t', mode='w')
    with open(local_path + '/filtered_pdbs_' + file_suffix + '.txt', 'a') as summary_file:
        summary_file.write('order ')
        summary_file.write(' '.join(dataframe.sort_values("total_score")["description"]))
        summary_file.write(' , no')
        summary_file.write(f'\nshared_prefix = "{iteration_name}_"')
        summary_file.write(f'\nnames = cmd.get_names("objects",1); for index in range(len(names)): cmd.set_name('
                           f'names[index],shared_prefix+str(index))')
        summary_file.write('\ndisable all')


def filter_runs():
    with open(main_local_path + script_name, 'w') as pull_file:
        pull_file.write('#!/bin/bash\n')
        pull_file.write('sftp -q '+current_username+'@'+remote_server_name+' << EOF')

    for folder in os.listdir(main_local_path):
        subfolder_server_path = main_server_path + folder + '/'
        if results_folder_flag:
            subfolder_server_path = subfolder_server_path + results_folder_name + '/'
        current_path = main_local_path + folder + '/'
        if os.path.isdir(current_path):
            print(folder)
            score_dataframe = sf.compile_score_data_frame(current_path, docking_flag)
            if filtered_number_flag:
                print(str(len(score_dataframe))+" total pdbs")
            filtered_dataframe = sf.filter_rosetta_pdbs(score_dataframe, docking_flag, index_order, all_cst_flag, cst_name, 
                                                        constraint_threshold, interface_energy_threshold_percent,
                                                        total_threshold_percent, folder, filtered_number_flag)
            create_request_script(filtered_dataframe, fetch_all_flag, amount_requested, subfolder_server_path,
                                  current_path, main_local_path,folder,filetag)
            if pymol_script:
                filtered_dataframe = filtered_dataframe.sort_values("total_score")
                filtered_dataframe = filtered_dataframe.reset_index(drop=True)
                with open(main_local_path+pymol_script_name, 'a') as pymol_file:
                    for index in range(len(filtered_dataframe)):
                        current_file = filtered_dataframe["description"][index]
                        pymol_file.write(f'cmd.load("{current_path}{current_file}.pdb","{folder}_{index}")\n')

    with open(main_local_path + script_name, 'a') as pull_file:
        pull_file.write('\nEOF')


if __name__ == '__main__':
    if pymol_script:
        if os.path.isfile(main_local_path+pymol_script_name):
            os.remove(main_local_path+pymol_script_name)
    filter_runs()
