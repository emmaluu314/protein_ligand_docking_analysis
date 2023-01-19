import os
import score_functions as sf

# User defined variables:
main_folder_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/ARPA-E/CAR/Cousins temp/'
main_server_path = '/share/siegellab/emmaluu/docking/template_for_automation/CAR/NMN/'
current_username = 'emmaluu'
remote_server_name = 'cacao.genomecenter.ucdavis.edu'
docking_flag = True
fetch_all_flag = True
results_folder_flag = True
amount_requested = 10  # max number of pdbs you want to be returned if fetch_all_flag is False

index_order = [0, 2, 1]
constraint_threshold = 1
interface_energy_threshold_percent = 0.2
total_threshold_percent = 0.1

script_name = 'pull_pdb_files.sh'


def create_request_script(dataframe, request_flag, request_number, server_path, local_path,
                          main_path):
    if request_flag:
        file_suffix = 'all'
    elif request_number <= len(dataframe):
        file_suffix = 'truncated'
    else:
        file_suffix = 'not_filtered'

    # with open(main_path + script_name, 'a') as table_file:
    #     table_file.write("\nos.system('")
    #     table_file.write('scp -r ' + username + '@' + server_name+':"')
    #     for name_index in range(len(dataframe['description'])):
    #         output_identifier = dataframe['description'].iloc[name_index]
    #         table_file.write(server_path)
    #         table_file.write(output_identifier + '.pdb')
    #         if name_index == len(dataframe['description'])-1:
    #             pass
    #         else:
    #             table_file.write(' ')
    #     table_file.write('" "' + local_path + '"')
    #     table_file.write("')\n")

    with open(main_path+script_name, 'a') as request_file:
        for name_index in range(len(dataframe['description'])):
            output_identifier = dataframe['description'].iloc[name_index]
            request_file.write("\nget "+server_path+output_identifier+'.pdb "'+local_path+output_identifier+'.pdb"')

    dataframe.to_csv(local_path + '/filtered_pdbs_' + file_suffix + '.txt',
                     columns=["description", "total_score", "all_cst", "total_interface_energy"],
                     index=None, sep='\t', mode='w')
    with open(local_path + '/filtered_pdbs_' + file_suffix + '.txt', 'a') as summary_file:
        summary_file.write('order ')
        summary_file.write(' '.join(dataframe["description"]))
        summary_file.write(' , no')


def filter_runs():
    with open(main_folder_path + script_name, 'w') as pull_file:
        pull_file.write('#!/bin/bash\n')
        pull_file.write('sftp -q '+current_username+'@'+remote_server_name+' << EOF')

    for folder in os.listdir(main_folder_path):
        subfolder_server_path = main_server_path + folder + '/'
        if results_folder_flag:
            subfolder_server_path = subfolder_server_path + 'results/'
        current_path = main_folder_path + folder + '/'
        if os.path.isdir(current_path):
            score_dataframe = sf.compile_score_data_frame(current_path, docking_flag)
            filtered_dataframe = sf.filter_rosetta_pdbs(score_dataframe, docking_flag, index_order,
                                                        constraint_threshold, interface_energy_threshold_percent,
                                                        total_threshold_percent, folder)
            create_request_script(filtered_dataframe, fetch_all_flag, amount_requested, subfolder_server_path,
                                  current_path, main_folder_path)
    with open(main_folder_path + script_name, 'a') as pull_file:
        pull_file.write('\nEOF')


if __name__ == '__main__':
    filter_runs()
