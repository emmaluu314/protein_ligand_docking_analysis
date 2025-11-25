import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from tqdm import tqdm

main_folder_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/Unnatural_cofactor/PTDH_main/PTDH_new/Multiple_sanity_check/'

code = 'P3N'
main_suffix = f""#"_out_analysis"
additional_suffix = ""
queries = [code]#[186,285,363,366,369,373,416]
resi_flag = False

WT_string = '' # if this is not applicable, define the variable to be ''. Else write the corresponding gene name.
gene_groups=[['P3NA'],['P3NA_OMe']]

#gene_groups=[['WT','P05091_K369R','P05091_T365R'],
#             ['WT','P05091_T365R_Q366H','P05091_T365R_K369R'],
#             ['WT','P05091_T365R','P05091_T365R_K369R','P05091_T365R_Q366H','P05091_RHR'],
#             ['WT','P05091_Q366H','P05091_T365R_Q366H', 'P05091_Q366H_K369R','P05091_RHR'],
#             ['WT','P05091_K369R','P05091_T365R_K369R', 'P05091_Q366H_K369R','P05091_RHR']]

figure_size_individual_positions = (15,5) # this must be a tuple
figure_size_compilation = (10,8) # this must be a tuple
flag_plot_residue_breakdown = False
flag_specify_gene_groups = True

if resi_flag:
    sub_suffix = "_resi_"
else:
    sub_suffix = "_"


for index in range(len(queries)):
    file_suffix = main_suffix+sub_suffix+str(queries[index])+additional_suffix
    file_path = f'{main_folder_path}bonds_dataframe{file_suffix}.csv'
    
    bonds_dataframe = pd.read_csv(file_path,converters={'bond list': pd.eval})
    gene_list = bonds_dataframe.gene.unique().tolist()
    all_contacts = bonds_dataframe['bond list'].sum()
    all_unique_positions = list(set(all_contacts))
    all_unique_positions.sort()
    frequency_dataframe = pd.DataFrame(columns = ['gene']+[str(position) for position in all_unique_positions])
    for gene in tqdm(gene_list):
        current_subset = bonds_dataframe[bonds_dataframe.gene == gene]
        total_structures = len(current_subset)
        current_compiled_contacts = current_subset['bond list'].sum()
        compiled_frequency = []
        for position in all_unique_positions:
            compiled_frequency.append(sum([item == position for item in current_compiled_contacts])/total_structures)
        frequency_dataframe.loc[len(frequency_dataframe)] = [gene]+compiled_frequency

    if flag_plot_residue_breakdown:
        if WT_string:
            WT_index = frequency_dataframe.gene.to_list().index(WT_string)
        for position in tqdm(all_unique_positions):
            fig = plt.figure(figsize = figure_size_individual_positions)
            barplot = plt.bar(frequency_dataframe.gene.values.tolist(),frequency_dataframe[str(position)].values.tolist(),color = 'cornflowerblue')
            ax = plt.gca()
            if WT_string:
                ax.get_children()[WT_index].set_color('gold') 
            plt.xticks(rotation=45, ha='right')
            ax.grid(axis='y')
            fig.savefig(f'{main_folder_path}frequency_{position}{file_suffix}.png',bbox_inches='tight', dpi=300)
        plt.close('all')

    if flag_specify_gene_groups:
        mutation_points = gene_groups # this is a list in a list
    else:
       mutation_points = [gene_list] # this is a list in a list
    counter = 0
    for genes in tqdm(mutation_points):
        if WT_string:
            current_subset = frequency_dataframe[frequency_dataframe.gene.isin(genes+[WT_string])]
        else:
            current_subset = frequency_dataframe[frequency_dataframe.gene.isin(genes)]
        summed_subset = current_subset.sum().to_list()
        index_of_nonzero = [index for index, value in enumerate(summed_subset) if value != 0]
        ax = current_subset.iloc[:,index_of_nonzero].plot(x='gene',kind='bar', align='center', stacked=False,figsize= figure_size_compilation)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xticks(rotation=0)
        ax.grid(axis='y')
        ax.set_axisbelow(True)
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))
        if 'single' in file_suffix:
            ax.figure.savefig(f'{main_folder_path}frequency_{genes[0][:-1]}{file_suffix}.png',bbox_inches='tight', dpi=300)
        else:
            ax.figure.savefig(f'{main_folder_path}frequency_all{file_suffix}_{str(counter)}.png',bbox_inches='tight', dpi=300)
            counter += 1
    
    frequency_dataframe.to_csv(main_folder_path+'frequency_dataframe'+file_suffix+'.csv')

