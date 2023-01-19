import pandas as pd
import pdb_functions as pf
import user_functions as uf
import numpy as np
import matplotlib.pyplot as plt

# User defined variables:
main_folder_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/ARPA-E/CAR/Cousins temp/'
target_list_path = main_folder_path + 'target_list.txt'
alignment_path = main_folder_path + 'VerticalAlignment_aldh_cousins.csv'
sphere_path = main_folder_path + '8_angstrom_NMN.txt'
amount_requested = 10  # max number of pdbs you want to be returned if fetch_all_flag is False
template_name = 'P51977'
filename_tag = ''


def calculate_shell_score():
    vertical_alignment_dataframe = pd.read_csv(alignment_path)
    template_atoms = pf.get_template_atoms(sphere_path, vertical_alignment_dataframe, template_name)
    target_genes = uf.compile_target_list(target_list_path)
    shell_dataframe = pd.DataFrame(
        columns=["gene", "sub_score", "interf_score", "total_score", "position", "numFailed"])

    for gene in target_genes:
        current_path = main_folder_path + gene + '/'
        score, interf_score, positions, failed_count, total_score, check_dataframe = \
            pf.compile_shell_residues(gene, current_path, template_atoms, vertical_alignment_dataframe)
        if len(positions):
            temporary_dataframe = pd.DataFrame(
                {"gene": gene, "sub_score": score, "interf_score": interf_score, "total_score": total_score,
                 "position": positions, "numFailed": failed_count})
            shell_dataframe = pd.concat([shell_dataframe, temporary_dataframe], ignore_index=True)
            # pf.residue_contribution(gene, filename_tag, check_dataframe, current_path)

    for row in range(shell_dataframe.shape[0]):
        shell_dataframe.sub_score[row] = [-1 * element for element in shell_dataframe.sub_score[row]]
        shell_dataframe.interf_score[row] = [-1 * element for element in shell_dataframe.interf_score[row]]
    return shell_dataframe


def analyze_shell_score():
    shell_dataframe = calculate_shell_score()
    genes_dataframe = pd.DataFrame(columns=["gene", "median"])
    for index in range(len(shell_dataframe['gene'])):
        current_gene = shell_dataframe.gene[index]
        current_median = np.median(shell_dataframe.sub_score[index])
        current_gene = pd.DataFrame(data=[[current_gene, current_median]], columns=["gene", "median"])
        genes_dataframe = pd.concat([genes_dataframe, current_gene], ignore_index=True)
    fig = plt.figure()
    ax = plt.gca()
    ax.hist(genes_dataframe["median"].to_list(), bins=100)
    ax.set_xlim([50, 150])
    fig.savefig(main_folder_path + 'hist.png')
    genes_dataframe.to_csv(main_folder_path + 'median.csv')


if __name__ == '__main__':
    analyze_shell_score()
