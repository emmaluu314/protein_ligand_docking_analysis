import os
import pandas as pd
import numpy as np


def pull_from_cluster(target_list, remote_server_name, username, main_server_path, main_local_path):
    python_script_name = 'pull_score_files.py'

    with open(main_local_path + python_script_name, 'w') as file:
        file.write("import os \n")
        for target in target_list:
            if '/' in target:
                index_of_target = target.find('/')
                target_folder = target[:index_of_target] + '_' + target[index_of_target + 1:]
            else:
                target_folder = target
            # edit the line below such that this script is compatible with different operating systems
            file.write("os.system('mkdir " + target_folder + "')\n")
            command_string = 'scp "' + username + '@' + remote_server_name + ':' + main_server_path + target + \
                             '/results/score*.sc" ' + target_folder
            file.write("os.system('" + command_string + "')\n")


def compile_score_data_frame(path, docking_flag):
    """
    This function compiles the score files in the specified path into a pandas dataframe.
    Modified from Augustine Arredono

    :param path: (string) the folder that contains the score files (.sc)
    :param docking_flag: (boolean) indicator for whether the score files came from a Relax protocol (False) or
    Docking protocol (True).
    :return pd.concat(score_list): a pandas dataframe of all the score files' info compiled
    """

    filenames = [filename for filename in os.listdir(path) if '.sc' in filename]  # list comprehension to compile .sc
    score_list = []  # an empty variable to compile all the files contents

    for name in filenames:  # Reading all the .sc files
        if docking_flag:
            header_value = 0
        else:
            header_value = 1
        dataframe = pd.read_csv(path + name, header=header_value, sep='\s+')  # header=1 for relax score files
        if 'SCORE:' in dataframe.columns:
            del dataframe['SCORE:']
        score_list.append(dataframe)
    final_dataframe = pd.concat(score_list)
    final_dataframe['total_interface_energy'] = np.sum(
        final_dataframe.loc[:, final_dataframe.columns.str.contains('interf_')], axis=1)
    return final_dataframe


def filter_rosetta_pdbs(score_dataframe, docking_flag, index_order, constraint_threshold,
                        interface_energy_threshold_percentage, total_energy_threshold_percentage, iteration_name):
    minimum_pdbs = 2  # bandaid fix
    column_names = ['all_cst', 'total_interface_energy', 'total_score']
    if docking_flag:
        ordered_columns = [column_names[index] for index in index_order]
        current_dataframe = score_dataframe.copy(deep=True)
        current_dataframe['total_interface_energy'] = np.sum(
            current_dataframe.loc[:, current_dataframe.columns.str.contains('interf_')], axis=1)

        for current_feature in ordered_columns:
            if current_feature == 'all_cst':
                current_threshold = constraint_threshold
            elif current_feature == 'total_score':
                current_threshold = max(current_dataframe[current_feature].nsmallest(
                    n=round(current_dataframe.shape[0] * total_energy_threshold_percentage), keep="all"))
            else:
                current_threshold = max(current_dataframe[current_feature].nsmallest(
                    n=round(current_dataframe.shape[0] * interface_energy_threshold_percentage), keep="all"))

            boolean_filter = current_dataframe[current_feature] <= current_threshold
            passed_filter_count = sum(1 for element in boolean_filter if element)
            if passed_filter_count >= minimum_pdbs:
                filtered_dataframe = current_dataframe[current_dataframe[current_feature] <= current_threshold]
                current_dataframe = filtered_dataframe.copy(deep=True)
            else:
                print(iteration_name,
                      ': Ended filtering early due to low number of pdbs after filtering ' + current_feature)
                break
        # print(iteration_name, ':', len(score_dataframe), len(current_dataframe))
    else:
        total_energy_threshold = max(
            score_dataframe.nsmallest(n=round(score_dataframe.shape[0] * total_energy_threshold_percentage),
                                      keep="all"))
        current_dataframe = score_dataframe[score_dataframe['total_score'] <= total_energy_threshold]
    return current_dataframe
