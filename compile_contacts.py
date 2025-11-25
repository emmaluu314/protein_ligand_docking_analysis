import os
import pandas as pd
from tqdm import tqdm
main_folder_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/Unnatural_cofactor/PTDH_main/PTDH_new/Multiple_sanity_check/'
sphere_radius = 6

code = 'P3N'
queries = [code]#[186,285,363,366,369,373,416]
resi_flag = False # turn off if you want to use the resn selector

main_suffix = f"" #_out_analysis"

if resi_flag:
    sub_suffix = "_resi_"
    chain_selection = "chain A and (not name N+C+O) and resi "
else:
    sub_suffix = "_"
    chain_selection = "resn "

for index in range(len(queries)):
    file_suffix = main_suffix+sub_suffix+str(queries[index])
    selection = chain_selection+str(queries[index])
    print("Checking out "+selection)
    residue_dataframe = pd.DataFrame()
    bonds_dataframe = pd.DataFrame(columns=['gene','structure', 'number of residues considered', 'bond list'])
    folders = [folder for folder in os.listdir(main_folder_path) if os.path.isdir(main_folder_path + folder + '/')]
    for folder in folders:
        current_path = main_folder_path + folder + '/'
        names = [filename for filename in os.listdir(current_path) if '.pdb' in filename]
        counter = 0
        for name in tqdm(names):
            cmd.load(f"{current_path}{name}")
            cmd.select("br. all within " + str(sphere_radius) + " of " + selection + f" and {name[:-4]}")
            sidechain_list = [];
            cmd.iterate("sele and (name C)", "sidechain_list.append(resv)")
            contact_list = []
            for residue_position in sidechain_list:
                current_measurement = cmd.distance(selection, "resi " + str(residue_position), mode="2")
                if current_measurement > 0:
                    contact_list.append(residue_position)
            cmd.delete('all')
            bonds_dataframe.loc[len(bonds_dataframe)] = [folder,name[:-4], len(sidechain_list), contact_list]
    bonds_dataframe.to_csv(main_folder_path+'bonds_dataframe'+file_suffix+'.csv')

