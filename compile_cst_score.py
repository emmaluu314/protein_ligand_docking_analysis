# something to be run in pymol

import pandas as pd
import os 
from tqdm import tqdm

### USER DEFINED VARIABLES

main_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/Unnatural_cofactor/NfsA/2025_04_16/'
constraint_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/Unnatural_cofactor/Constraints/Constraints_NfsA/'
output_extension = '.pdb'
# max_folders = 10 not doing this anymore. all the outputs should be compiled into one folder
# max_files = 200 # this is to limit the calculation time if needed. I think in the future we need to paralize it
# The structure of the folders should be : 
# main_folder > subfolders with their own cst
# main_folder > subfolders with their own pdbs

### FUNCTIONS
def list_folders(path):
    return [f for f in os.listdir(path) if os.path.isdir(os.path.join(path,f))]

def get_file_of_type(path, extension):
    return [f for f in os.listdir(path) if f.endswith(extension)]

def load_prediction_files(path,folder):
    current_path = os.path.join(path,folder)
    files = get_file_of_type(current_path, output_extension)

    for file in tqdm(files, desc="Loading files"):
        print(os.path.join(current_path,file))
        cmd.load(os.path.join(current_path,file))

        cmd.set_name(file[:-4],folder+'_'+file[:-4])

def parse_constraint_files(path):
    remark_file = "".join(get_file_of_type(path,'.pdb'))
    cst_file = "".join(get_file_of_type(path,'.cst'))

    # Dictionary to store constraints for each residue
    constraints = {}

    # Parse the REMARK file
    with open(os.path.join(path,remark_file), 'r') as f: 
        block_index = -1

        for line in f:
            if line.startswith('REMARK 666 MATCH TEMPLATE'):
                block_index += 1
                parts = line.split()
                residue_indices = [parts[6],parts[11]]
                chains = [parts[4],parts[9]]
                constraints[block_index] = {
                    'first_selection': f" and chain {chains[0]} and resi {residue_indices[0]}",
                    'second_selection': f" and chain {chains[1]} and resi {residue_indices[1]}"
                }

    # Parse the CST file
    with open(os.path.join(path,cst_file), 'r') as f:
        block_index = -1

        for line in f:
            line = line.strip()
            if line.startswith('CST::BEGIN'):
                block_index += 1
                constraints[block_index]['first_atoms'] = []
                constraints[block_index]['second_atoms'] = []
                constraints[block_index]['measure_type'] = []
                constraints[block_index]['measure_value'] = []
            elif line.startswith('TEMPLATE::') and 'atom_name:' in line:
                parts = line.split('atom_name:')
                atom_map = int(parts[0].split(':')[-1])
                atoms = parts[1].strip().split()

                if atom_map == 1:
                    constraints[block_index]['first_atoms'] = atoms
                elif atom_map == 2:
                    constraints[block_index]['second_atoms'] = atoms
            elif line.startswith('CONSTRAINT::'):
                parts = line.split('::', 1)[1].split()
                constraint_type = parts[0].rstrip(':')
                constraints[block_index]['measure_type'].append(constraint_type)
                constraints[block_index]['measure_value'].append([float(number) for number in parts[1:-1]])

    return pd.DataFrame.from_dict(constraints, orient='index')

    
def calculate_distance(object_name, selection_one, selection_two, atoms_one, atoms_two):
    # return the distance
    first_atom = "model "+ object_name + selection_one + ' and name ' + atoms_one[0]
    second_atom =  "model "+ object_name + selection_two + ' and name ' + atoms_two[0]
    return cmd.get_distance(first_atom, second_atom)
    

def calculate_angle(metric_name,object_name, selection_one, selection_two, atoms_one, atoms_two):
    residue_one =  "model "+ object_name + selection_one + ' and name ' 
    residue_two =  "model "+ object_name + selection_two + ' and name ' 
    if 'A' in metric_name:
        first_atom = residue_two + atoms_two[0]
        second_atom =  residue_one + atoms_one[0]
        third_atom =  residue_one + atoms_one[1]
    elif 'B' in metric_name:
        first_atom = residue_one + atoms_one[0]
        second_atom =  residue_two + atoms_two[0]
        third_atom =  residue_two + atoms_two[1]
    
    return cmd.get_angle(first_atom,second_atom,third_atom)

def calculate_torsion(metric_name,object_name, selection_one, selection_two, atoms_one, atoms_two):
    residue_one =  "model "+ object_name + selection_one + ' and name ' 
    residue_two =  "model "+ object_name + selection_two + ' and name ' 

    if "AB" in metric_name:
        first_atom = residue_one + atoms_one[1]
        second_atom =  residue_one + atoms_one[0]
        third_atom =  residue_two + atoms_two[0]
        fourth_atom = residue_two + atoms_two[1]
    elif "A" in metric_name:
        first_atom = residue_one + atoms_one[2]
        second_atom =  residue_one + atoms_one[1]
        third_atom =  residue_one + atoms_one[0]
        fourth_atom = residue_two + atoms_two[0]
    elif "B" in metric_name:
        first_atom = residue_one + atoms_one[0]
        second_atom =  residue_two + atoms_two[0]
        third_atom =  residue_two + atoms_two[1]
        fourth_atom = residue_two + atoms_two[2]

    return cmd.get_dihedral(first_atom,second_atom,third_atom,fourth_atom)

def calculate_penalty(measurement,tolerance):
    difference = abs(measurement - tolerance[0])
    if difference <= tolerance[1]:
        penalty = 0
    else:
        penalty = tolerance[2]*(difference-tolerance[1])
    return penalty

def calculate_corresponding_constraint_score(path,current_constraints):
    names = cmd.get_names("objects",1)
    print(names)
    score_dataframe = pd.DataFrame(columns=['filename','measurement','score','constraint','type','block'])
    current_measurement = -1


    for block_index in tqdm(range(len(current_constraints))):
        current_first_atoms = current_constraints.loc[block_index,'first_atoms']
        current_second_atoms =  current_constraints.loc[block_index,'second_atoms']
        current_first_selection = current_constraints.loc[block_index,'first_selection']
        current_second_selection = current_constraints.loc[block_index,'second_selection']

        for metric_index in range(len(current_constraints.loc[block_index,'measure_type'])):
                current_metric = current_constraints.loc[block_index,'measure_type'][metric_index]
                current_tolerance = current_constraints.loc[block_index,'measure_value'][metric_index]
                for name in tqdm(names, desc="going through block "+str(block_index)+" "+current_metric):
                    if 'distance' in current_metric:
                        current_measurement = calculate_distance(name, 
                                                                 current_first_selection, current_second_selection,
                                                                 current_first_atoms, current_second_atoms)
                    elif 'angle' in current_metric:
                        current_measurement = calculate_angle(current_metric, name, 
                                                              current_first_selection, current_second_selection,
                                                              current_first_atoms, current_second_atoms)
                    elif 'torsion' in current_metric:
                        current_measurement = calculate_torsion(current_metric, name, 
                                                              current_first_selection, current_second_selection,
                                                              current_first_atoms, current_second_atoms)
                        calculate_penalty(current_measurement,current_tolerance)
                    current_score = calculate_penalty(current_measurement,current_tolerance)
                    score_dataframe = pd.concat([score_dataframe, pd.DataFrame({
                            'filename': name,
                            'measurement': current_measurement,
                            'score': current_score,
                            'constraint': [current_tolerance],
                            'type': current_metric,
                            'block': block_index})], ignore_index = True)
                    #print(name,current_measurement,current_tolerance,current_score)

    score_dataframe.to_csv(path+'_cst_scores.csv')


### START OF MAIN CODE
subfolders = list_folders(main_path)

for subfolder in tqdm(subfolders,desc="Going through enzymes"): 
    # prepping for analysis 

    current_path = os.path.join(main_path,subfolder)
    load_prediction_files(main_path,subfolder) # all of the files you want to analyze should be loaded into pymol now
    current_constraints = parse_constraint_files(os.path.join(constraint_path,subfolder))

    calculate_corresponding_constraint_score(current_path,current_constraints)
    #delete all of the objects loaded before moving on to the next folder
    cmd.delete('all')

print('Finished')




