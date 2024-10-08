import os
import pandas as pd
import numpy as np


def pull_from_cluster(target_list, remote_server_name, username, main_server_path, results_folder_name, main_local_path):
    script_name = 'pull_score_files.sh'
    with open(main_local_path+script_name, 'w') as request_file:
        request_file.write("#!/bin/bash\nsftp -q ")
        request_file.write(username+"@"+remote_server_name+" << EOF")

    with open(main_local_path+script_name, 'a') as request_file:
        for target in target_list:
            if '/' in target:
                index_of_target = target.find('/')
                target_folder = target[:index_of_target] + '_' + target[index_of_target + 1:]
            else:
                target_folder = target
            request_file.write("\nget "+main_server_path+target+'/'+results_folder_name+'score*.sc "'+main_local_path+target_folder+'"')
            if not os.path.isdir(main_local_path+target_folder):
                os.mkdir(main_local_path+target_folder)


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
        else: # header=1 for relax score files
            header_value = 1
        dataframe = pd.read_csv(path + name, header=header_value, sep='\s+')  
        if 'SCORE:' in dataframe.columns:
            del dataframe['SCORE:']
        score_list.append(dataframe)
    final_dataframe = pd.concat(score_list)
    final_dataframe['total_interface_energy'] = np.sum(
        final_dataframe.loc[:, final_dataframe.columns.str.contains('interf_')], axis=1)
    return final_dataframe


def filter_rosetta_pdbs(score_dataframe, docking_flag, index_order, all_cst_flag, cst_name, constraint_threshold,
                        interface_energy_threshold_percentage, total_energy_threshold_percentage, iteration_name,
                        display_flag):
    minimum_pdbs = 0  # bandaid fix
    column_names = ['total_interface_energy', 'total_score']  
    if docking_flag:
        cst_passed_flag = True
        ordered_columns = [column_names[index] for index in index_order]
        current_dataframe = score_dataframe.copy(deep=True)
        current_dataframe['total_interface_energy'] = np.sum(
            current_dataframe.loc[:, current_dataframe.columns.str.contains('interf_')], axis=1)

        if all_cst_flag:  # using all_cst score
            boolean_filter = current_dataframe[cst_name] <= constraint_threshold
            print(np.median(current_dataframe[cst_name]))
            passed_filter_count = sum(constraint_threshold for element in boolean_filter if element)
            if passed_filter_count > minimum_pdbs:
                filtered_dataframe = current_dataframe[current_dataframe[cst_name] <= constraint_threshold]
                current_dataframe = filtered_dataframe.copy(deep=True)
            else:
                print(iteration_name, ': Ended filtering early due to low number of pdbs after filtering all_cst')
                cst_passed_flag = False
        else:  # using individual cst scores
            cst_columns = current_dataframe.columns[current_dataframe.columns.str.contains('_' + cst_name)]
            passed_filter_count = sum(
                (current_dataframe[cst_columns] <= constraint_threshold).sum(axis=1) == len(cst_columns))
            if passed_filter_count > minimum_pdbs:
                filtered_dataframe = current_dataframe[
                    (current_dataframe[cst_columns] < constraint_threshold).sum(axis=1) == len(cst_columns)]
                current_dataframe = filtered_dataframe.copy(deep=True)
            else:
                print(iteration_name, ': Ended filtering early due to low number of pdbs after filtering all_cst')
                cst_passed_flag = False

        if display_flag:
            print(f"{passed_filter_count} passed the cst filter.")

        if cst_passed_flag:
            for current_feature in ordered_columns:
                if current_feature == 'total_score':
                    current_threshold = max(current_dataframe[current_feature].nsmallest(
                        n=round(current_dataframe.shape[0] * total_energy_threshold_percentage), keep="all"))
                else:
                    current_threshold = max(current_dataframe[current_feature].nsmallest(
                        n=round(current_dataframe.shape[0] * interface_energy_threshold_percentage), keep="all"))

                boolean_filter = current_dataframe[current_feature] <= current_threshold
                passed_filter_count = sum(1 for element in boolean_filter if element)
                if display_flag:
                    print(f"{passed_filter_count} passed the {current_feature} filter of {current_threshold}")

                if passed_filter_count >= minimum_pdbs:
                    filtered_dataframe = current_dataframe[current_dataframe[current_feature] <= current_threshold]
                    current_dataframe = filtered_dataframe.copy(deep=True)
                else:
                    print(iteration_name,
                          ': Ended filtering early due to low number of pdbs after filtering ' + current_feature)
                    print(f'Median of {current_feature} is {np.median(current_dataframe[current_feature])}')
                    break
        else:
            total_energy_threshold = max(
                score_dataframe['total_score'].nsmallest(
                    n=round(score_dataframe.shape[0] * total_energy_threshold_percentage),
                    keep="all"))
            current_dataframe = score_dataframe[score_dataframe['total_score'] <= total_energy_threshold]
    else:
        current_dataframe = score_dataframe
    return current_dataframe
