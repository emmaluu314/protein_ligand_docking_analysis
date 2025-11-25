import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import score_functions as sf


def get_template_atoms(file_path, alignment_dataframe, template_name):
    template_atoms = pd.read_csv(file_path, sep=" ", names=["Residue", "Position"])
    template_atoms["Position"] = template_atoms["Position"].astype(int)

    alignment_series = alignment_dataframe.loc[alignment_dataframe[template_name] != '-', template_name]
    template_atoms = alignment_series.index[template_atoms["Position"].to_list()]
    return template_atoms


def compile_shell_residues(target_name, target_path, template_atoms, alignment_dataframe):
    temp = alignment_dataframe[target_name]
    position = alignment_dataframe[target_name].iloc[template_atoms]
    numeric_position = []
    gap_count = 0
    counter = 0

    for index in range(len(position)):
        if '-' in position.iloc[index]:
            gap_count += 1
        else:
            numeric_position.append(int(position.iloc[index][1:]))

    current_score_dataframe = sf.compile_score_data_frame(target_path, docking_flag=True)
    interface_list = []
    subscore_list = []
    total_score_list = []
    for filename in os.listdir(target_path):
        if filename.endswith('.pdb'):
            row_boolean = current_score_dataframe['description'] == filename[:-4]
            row_index = [i for i, value in enumerate(row_boolean) if value]
            interface_score = current_score_dataframe['total_interface_energy'].iloc[row_index[0]]
            total_score = current_score_dataframe['total_score'].iloc[row_index[0]]
            interface_list.append(interface_score)
            total_score_list.append(total_score_list)
            # need to make an array for calculating the average interface_score

            pdb_file = open(target_path + '/' + filename)
            pdb = pdb_file.read()
            pdb_file.close()
            if '.gz' in pdb:
                gz = '.gz'
            else:
                gz = ''

            pdb = pdb.split("\n")
            search_string_beginning = "#BEGIN_POSE_ENERGIES_TABLE " + filename + gz
            search_string_ending = "#END_POSE_ENERGIES_TABLE " + filename + gz

            start_index = pdb.index(search_string_beginning)
            end_index = pdb.index(search_string_ending)

            new_list = []
            for index in range(start_index + 2, end_index):
                new_list.append(pdb[index].split(" "))
            dataframe = pd.DataFrame(new_list, columns=pdb[start_index + 1].split(' '))

            subtotal_score = 0

            if counter == 0:  # bandaid fix
                numeric_position = numeric_position + [dataframe.shape[0] - 1]
                counter += 1
            else:
                pass

            for index in numeric_position:
                try:
                    subtotal_score += float(dataframe["total"].iloc[index])
                except:
                    print(target_name + "'s structure is probably truncated.")

            subscore_list.append(subtotal_score)
            position = pd.concat([position, dataframe["label"].tail(1)], ignore_index=True)
            # need to make an array for calculating the average of subtotal_score

            if 'contribution_dataframe' in locals():
                contribution_dataframe = pd.concat([contribution_dataframe, dataframe.iloc[numeric_position, [0, -1]]],
                                                   ignore_index=True)
            else:
                contribution_dataframe = dataframe.iloc[numeric_position, [0, -1]]

    return [subscore_list], [interface_list], [position], gap_count, [total_score_list], contribution_dataframe


def residue_contribution(gene, filename_tag, check_dataframe, current_path):
    check_dataframe.sort_values(by="label")
    new_dataframe = pd.DataFrame(columns=["residue", "contribution"])
    residue_means = []
    for residue in check_dataframe.label.unique():
        current_list = check_dataframe.loc[check_dataframe["label"] == residue, "total"].to_list()
        current_subdataframe = pd.DataFrame(data={"residue": residue, "contribution": [np.float_(current_list)]})
        new_dataframe = pd.concat([new_dataframe, current_subdataframe], ignore_index=True)
        residue_means.append(np.mean(np.float_(current_list)))
    index_sort_by_means = np.argsort(residue_means)[::-1]
    figure = plt.figure(figsize=(8, 8))
    axs = plt.gca()
    box_plot = axs.boxplot(new_dataframe.contribution[index_sort_by_means],
                           vert=False,  # vertical box alignment
                           patch_artist=True,  # fill with color
                           labels=new_dataframe.residue[index_sort_by_means])  # will be used to label x-ticks
    axs.set_title("Shell score residue breakdown for " + gene)
    axs.set_xlabel("Weighted total score for individual residue")
    plt.close(figure)
    figure.savefig(current_path + "residue_contribution_" + gene + filename_tag + ".png")
