"""
compile_ligand_docking_results.py

Date of last update: 8/22/25

Description:
------------
This script processes Rosetta ligand docking results. It extracts constraint information 
from constraint (`.cst`) files, computes geometric measurements (currently distances) 
between template and motif molecules in docked complexes, and compiles the results 
into structured `.csv` files for downstream analysis. Additionally, score files 
(`.sc`) are compiled, and measurement summaries can be used for histogram/figure generation.

Key Features:
-------------
1. Reads multiple docking output folders containing `.pdb`, `.sc`, and `.cst` files.
2. Extracts and parses constraint definitions (distances, angles, dihedrals).
3. Computes geometric distances between ligand (motif) atoms and protein (template) atoms.
4. Aggregates results across runs into summary `.csv` files.
5. Provides functions for visualizing results (histograms).
6. User only needs to update:
   - `originalWorkingDirectory` → root directory where docking subfolders reside.
   - `suffix` → naming suffix for output result files.

Main Outputs:
-------------
- `[folder][suffix]_scoreFile.csv` : compiled constraint scores from `.sc` files
- `[folder][suffix]_measurements.csv` : table of measured constraints per `.pdb`
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from biopandas.pdb import PandasPdb
from tqdm import tqdm
import warnings

# Only suppress specific warnings if needed
warnings.filterwarnings("ignore", category=UserWarning, module="biopandas")

# ========== PLOTTING STYLE ==========
def set_pub():
    """
    Update matplotlib plotting defaults for publication-quality output.

    Modifies plt.rcParams to:
    - Place axes grid lines below plotted data
    - Save figures with 300 dpi resolution

    Inputs:
        None
    Outputs:
        None (modifies global matplotlib state)
    """
    plt.rcParams.update({
        "axes.axisbelow": True,
        "savefig.dpi": 300,
    })

# ========== COMPILATION FUNCTIONS ==========

def compile_scoreFiles(filePathway, desiredFileName):
    """
    Compile Rosetta score files (.sc) into a single CSV.

    This function looks for `.sc` files in the given path, merges them into a
    single DataFrame, and saves the compiled results.

    Inputs:
        filePathway : str
            Directory path containing .sc files.
        desiredFileName : str
            Base prefix for output compiled CSV file.

    Outputs:
        score_dataFrame : pd.DataFrame
            Combined DataFrame of all score files with constraints summed into 
            a single 'Constraint' column (if present).
        Also saves: [desiredFileName]_scoreFile.csv
    """
    listdata = []
    for fname in os.listdir(filePathway):
        if fname.endswith('.sc'):
            d = pd.read_csv(os.path.join(filePathway, fname), sep=r'\s+')
            listdata.append(d)
    if not listdata:
        print("No .sc files found in", filePathway)
        return pd.DataFrame()
    score_dataFrame = pd.concat(listdata, ignore_index=True)
    cst_cols = [col for col in score_dataFrame if col.endswith('_all_cst')]
    if cst_cols:
        score_dataFrame['Constraint'] = score_dataFrame[cst_cols].sum(axis=1)
    score_dataFrame = score_dataFrame.reset_index(drop=True)
    score_dataFrame.to_csv(os.path.join(filePathway, f"{desiredFileName}_scoreFile.csv"), index=False)
    return score_dataFrame

def compile_constraintParameters(pdb_DataFrame, cstFile_DataFrame, blockIndices_DataFrame):
    """
    Extract atom-specific constraint information from a CST file and integrate with PDB remarks.

    Inputs:
        pdb_DataFrame : PandasPdb
            Parsed PDB file in BioPandas format.
        cstFile_DataFrame : pd.DataFrame
            Raw CST file lines stored under column 'Original line'.
        blockIndices_DataFrame : pd.DataFrame
            Indices marking CST::BEGIN and CST::END blocks.

    Outputs:
        constraints_DataFrame : pd.DataFrame
            Each row corresponds to one constraint block, with fields for:
            - Template/motif atom names
            - Target distances/angles/dihedrals
            - Molecule/residue identifiers parsed from REMARK lines
    """
    columns = [
        'TEMPLATE Atoms', 'MOTIF Atoms', 'Distance AB',
        'Angle A', 'Angle B', 'Dihedral A', 'Dihedral B', 'Dihedral AB'
    ]
    constraintsAtom_DataFrame = pd.DataFrame(columns=columns)
    atomIdentifier = 'atom_name:'
    constraintIdentifier = 'CONSTRAINT::'

    for k in range(blockIndices_DataFrame.shape[0]):
        lineIndices = blockIndices_DataFrame.iloc[k, 0:blockIndices_DataFrame.shape[1]]
        block_DataFrame = pd.DataFrame(
            data=[[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']],
            columns=columns,
            index=[k+1]
        )
        for l in np.arange(lineIndices['start']+1, lineIndices['end']):
            currentLine = cstFile_DataFrame['Original line'][l].split()
            if atomIdentifier in currentLine:
                if len(currentLine) > 5:
                    currentSubset = [currentLine[-3:]]
                else:
                    currentSubset = currentLine[-1]
                if '1' in currentLine:
                    block_DataFrame['TEMPLATE Atoms'] = [currentSubset]
                elif '2' in currentLine:
                    block_DataFrame['MOTIF Atoms'] = [currentSubset]
            elif constraintIdentifier in currentLine:
                if 'distanceAB:' in currentLine:
                    block_DataFrame['Distance AB'] = [currentLine[-4:]]
                elif 'angle_A:' in currentLine:
                    block_DataFrame['Angle A'] = [currentLine[-4:]]
                elif 'angle_B:' in currentLine:
                    block_DataFrame['Angle B'] = [currentLine[-4:]]
                elif 'torsion_AB:' in currentLine:
                    block_DataFrame['Dihedral AB'] = [currentLine[-4:]]
                elif 'torsion_A:' in currentLine:
                    block_DataFrame['Dihedral A'] = [currentLine[-4:]]
                elif 'torsion_B:' in currentLine:
                    block_DataFrame['Dihedral B'] = [currentLine[-4:]]
        constraintsAtom_DataFrame = pd.concat([constraintsAtom_DataFrame, block_DataFrame])

    constraints_DataFrame = pd.DataFrame(
        pdb_DataFrame.df['OTHERS'][pdb_DataFrame.df['OTHERS']['record_name'] == 'REMARK']
    )
    placeholder_DataFrame = pd.DataFrame()
    for i in range(constraints_DataFrame.shape[0]):
        tempEntry = constraints_DataFrame.iloc[i, 1].split()
        temp_DataFrame = pd.DataFrame({
            'TEMPLATE Molecule': tempEntry[3],
            'TEMPLATE Res Name': tempEntry[4],
            'TEMPLATE Res Number': tempEntry[5],
            'MOTIF Molecule': tempEntry[8],
            'MOTIF Res Name': tempEntry[9],
            'MOTIF Res Number': tempEntry[10]
        }, index=[i+1])
        placeholder_DataFrame = pd.concat([placeholder_DataFrame, temp_DataFrame])
    constraints_DataFrame = pd.concat([constraints_DataFrame, placeholder_DataFrame], axis=1)
    constraints_DataFrame = pd.concat([constraints_DataFrame, constraintsAtom_DataFrame], axis=1)
    return constraints_DataFrame

def compile_moleculeInfo(moleculeType, pdb_DataFrame, constraints_DataFrame):
    """
    Retrieve molecule atom coordinates for either template or motif residues 
    referenced in constraint definitions.

    Inputs:
        moleculeType : str
            Either 'template' or 'motif'.
        pdb_DataFrame : PandasPdb
            Parsed PDB structure.
        constraints_DataFrame : pd.DataFrame
            Constraint metadata table with residue numbers and chain IDs.

    Outputs:
        molecule_DataFrame : pd.DataFrame
            Subset of the PDB DataFrame containing only atoms from relevant residues.
    """
    if moleculeType == 'template':
        moleculeString = 'TEMPLATE Molecule'
        residueString = 'TEMPLATE Res Number'
    elif moleculeType == 'motif':
        moleculeString = 'MOTIF Molecule'
        residueString = 'MOTIF Res Number'
    else:
        print('Please specify molecule type (template or motif).')
        return pd.DataFrame()

    molecule_DataFrame = pd.DataFrame()
    previousResNumber = []
    previousChainID = []
    for k in range(constraints_DataFrame.shape[0]):
        index = k + 1
        chainID = constraints_DataFrame[moleculeString][index]
        residueNumber = int(constraints_DataFrame[residueString][index])
        for j in pdb_DataFrame.df.keys():
            if j != 'OTHERS':
                subset = pdb_DataFrame.df[j][pdb_DataFrame.df[j]['chain_id'] == chainID]
                if subset.shape[0] > 0:
                    if (residueNumber in previousResNumber) and (chainID in previousChainID):
                        pass
                    else:
                        pdb_DataFrameSubset = pd.DataFrame(subset)
                        finalSubset = pdb_DataFrameSubset[pdb_DataFrameSubset['residue_number'] == residueNumber]
                        molecule_DataFrame = pd.concat([molecule_DataFrame, finalSubset], ignore_index=True)
                        previousResNumber.append(residueNumber)
                        previousChainID.append(chainID)
    return molecule_DataFrame

# ========== PLACEHOLDER: You must implement these functions ==========

def identify_constraintBlocks(cstFile_DataFrame):
    """
    Identify the line indices of constraint blocks in a CST file.

    CST files use 'CST::BEGIN' and 'CST::END' tags to separate blocks.

    Inputs:
        cstFile_DataFrame : pd.DataFrame
            Lines of the CST file as a column 'Original line'.

    Outputs:
        blockIndices_DataFrame : pd.DataFrame
            DataFrame with 'start' and 'end' indices for each block.
    """
    blockStart = 'CST::BEGIN'
    blockEnd = 'CST::END'
    blockIndices_DataFrame = pd.DataFrame(columns = ['start', 'end'])
    blockCounter = -1
    # Getting the indices of the blocks
    for i in range(cstFile_DataFrame.shape[0]):
        currentLine = cstFile_DataFrame['Original line'][i]
        #print(currentLine)
        if blockStart in currentLine:
            blockCounter = blockCounter + 1
            blockIndices_DataFrame = pd.concat([blockIndices_DataFrame,pd.DataFrame([[i,0]],columns = ['start','end'])], ignore_index= True)
            #print('start')
        elif blockEnd in currentLine:
            #blockIndices_DataFrame['end'][blockCounter] = i
            blockIndices_DataFrame.loc[blockCounter, "end"] = i
            #print('end')
        else:
            pass
    return blockIndices_DataFrame

def calculate_distanceAB(templateMolecule_DataFrame, motifMolecule_DataFrame, constraints_DataFrame):
    """
    Calculate Euclidean distances between template and motif atom pairs defined in constraints.

    Inputs:
        templateMolecule_DataFrame : pd.DataFrame
            Template molecule coordinates (subset of PDB atoms).
        motifMolecule_DataFrame : pd.DataFrame
            Motif molecule coordinates (subset of PDB atoms).
        constraints_DataFrame : pd.DataFrame
            Defines expected atom pairs and target distances.

    Outputs:
        dist_temp_array : list[float]
            Measured distances (Å) for each constraint.
        ssr_temp_array : list[float]
            Normalized squared residuals [(measured-expected)/expected]^2 for each constraint.
    """
    dist_temp_array = []
    for m in range(constraints_DataFrame.shape[0]):
        #print(constraints_DataFrame['TEMPLATE Atoms'].iloc[m], constraints_DataFrame['MOTIF Atoms'].iloc[m])
        #print(type(constraints_DataFrame['TEMPLATE Atoms'].iloc[m]), type(constraints_DataFrame['MOTIF Atoms'].iloc[m]))
        if  isinstance(constraints_DataFrame['TEMPLATE Atoms'].iloc[m],str): # checking if one or more atoms were specified in constraint block
            #print('just one atom')
            template_atom_condition = templateMolecule_DataFrame['atom_name']==constraints_DataFrame['TEMPLATE Atoms'].iloc[m]
            template_residue_condition = templateMolecule_DataFrame['residue_number']==int(constraints_DataFrame['TEMPLATE Res Number'].iloc[m])
        else:
            #print('more than one atom')
            template_atom_condition = templateMolecule_DataFrame['atom_name']==constraints_DataFrame['TEMPLATE Atoms'].iloc[m][0][0]
            template_residue_condition = templateMolecule_DataFrame['residue_number']==int(constraints_DataFrame['TEMPLATE Res Number'].iloc[m])
        #print(templateMolecule_DataFrame[['atom_name','residue_number']][template_atom_condition & template_residue_condition])

        if isinstance(constraints_DataFrame['MOTIF Atoms'].iloc[m],str): # checking if one or more atoms were specified in constraint block
            #print('just one atom')
            motif_atom_condition = motifMolecule_DataFrame['atom_name']==constraints_DataFrame['MOTIF Atoms'].iloc[m]
            motif_residue_condition = motifMolecule_DataFrame['residue_number']==int(constraints_DataFrame['MOTIF Res Number'].iloc[m])
        else:
            #print('more than one atom')
            motif_atom_condition = motifMolecule_DataFrame['atom_name']==constraints_DataFrame['MOTIF Atoms'].iloc[m][0][0] #bandaid fix
            motif_residue_condition = motifMolecule_DataFrame['residue_number']==int(constraints_DataFrame['MOTIF Res Number'].iloc[m])
        #print(motifMolecule_DataFrame[['atom_name','residue_number']][motif_atom_condition & motif_residue_condition])

        template_coordinates = templateMolecule_DataFrame[['x_coord','y_coord','z_coord']][template_atom_condition & template_residue_condition].to_numpy()
        motif_coordinates = motifMolecule_DataFrame[['x_coord','y_coord','z_coord']][motif_atom_condition & motif_residue_condition].to_numpy()
        #print(constraints_DataFrame['TEMPLATE Atoms'].iloc[m], template_coordinates,motif_coordinates)
        sq_dist = np.sum(np.square(template_coordinates - motif_coordinates))
        #print(np.sqrt(sq_dist))
        dist_temp_array.append(np.sqrt(sq_dist))
    # calculating the ssr
    user_distAB_array = [float(constraints_DataFrame['Distance AB'].iloc[n][0]) for n in range(len(constraints_DataFrame))]
    ssr_temp_array = [np.square((dist_temp_array[p] - user_distAB_array[p])/user_distAB_array[p]) for p in range(len(user_distAB_array))]
    return dist_temp_array, ssr_temp_array

# ========== MAIN MEASUREMENT COMPILATION FUNCTION ==========

def compile_measurements(filePathway, cstFileName, userInputConstraintNames, desiredFileName):
    """
    Main function: compile geometric constraint measurements across all PDBs.

    Steps:
    1. Read the CST file and identify constraint blocks.
    2. For each PDB file in the folder:
        - Parse PDB structure
        - Extract constraint atom/residue mappings
        - Compute distances between key atoms
    3. Save results as a CSV with one row per PDB.

    Inputs:
        filePathway : str
            Directory containing PDB and CST files.
        cstFileName : str
            Name of CST file (must match PDB dataset).
        userInputConstraintNames : list[str]
            User-defined labels for constraints, e.g., ["Cu-His420", "AFB1-His498"].
        desiredFileName : str
            Prefix for output CSV file.

    Outputs:
        measurements_DF : pd.DataFrame
            Table of constraint measurements per PDB file.
        Also saves: [desiredFileName]_measurements.csv
    """
    # Predefine all column names
    distance_cols = [f"Distance AB {name}" for name in userInputConstraintNames]
    angleA_cols = [f"Angle A {name}" for name in userInputConstraintNames]
    angleB_cols = [f"Angle B {name}" for name in userInputConstraintNames]
    dihedralA_cols = [f"Dihedral A {name}" for name in userInputConstraintNames]
    dihedralB_cols = [f"Dihedral B {name}" for name in userInputConstraintNames]
    dihedralAB_cols = [f"Dihedral AB {name}" for name in userInputConstraintNames]
    
    # Single list to collect all measurement data
    all_data = []
    file_names = []

    # Load constraint file
    cst_path = os.path.join(filePathway, cstFileName)
    cstFile_DataFrame = pd.read_csv(cst_path, header=None)
    cstFile_DataFrame.columns = ['Original line']
    
    # Identify constraint blocks
    blockIndices_DataFrame = identify_constraintBlocks(cstFile_DataFrame)

    for file in os.listdir(filePathway):
        if file.endswith('.pdb'):
            # Process PDB file
            pdb_path = os.path.join(filePathway, file)
            pdb_DataFrame = PandasPdb().read_pdb(pdb_path)
            file_names.append(file.replace('.pdb', ''))
            
            # Get constraint parameters
            constraints_DataFrame = compile_constraintParameters(
                pdb_DataFrame, cstFile_DataFrame, blockIndices_DataFrame
            )
            template_DF = compile_moleculeInfo('template', pdb_DataFrame, constraints_DataFrame)
            motif_DF = compile_moleculeInfo('motif', pdb_DataFrame, constraints_DataFrame)

            # Calculate measurements
            dist_values, _ = calculate_distanceAB(template_DF, motif_DF, constraints_DataFrame)
            # Create measurement dictionary
            measurement_dict = {}

            for i, name in enumerate(userInputConstraintNames):
                measurement_dict[f"Distance AB {name}"] = dist_values[i]
                # Add placeholders for other measurements
                measurement_dict[f"Angle A {name}"] = np.nan
                measurement_dict[f"Angle B {name}"] = np.nan
                measurement_dict[f"Dihedral A {name}"] = np.nan
                measurement_dict[f"Dihedral B {name}"] = np.nan
                measurement_dict[f"Dihedral AB {name}"] = np.nan
            
            all_data.append(measurement_dict)

    # Create final DataFrame
    measurements_DF = pd.DataFrame(all_data)
    measurements_DF['fileName'] = file_names
    
    # Save results
    output_path = os.path.join(filePathway, f"{desiredFileName}_measurements.csv")
    measurements_DF.to_csv(output_path, index=False)
    return measurements_DF


def extract_constraint_names(file_path):
    """
    Extract human-readable names for constraints from REMARK lines in a PDB file.

    Logic:
        - For Cu ligands, names like 'Cu-His420'
        - For LG1 ligands (e.g., AFB1), names like 'AFB1-His498'
        - Default: [Ligand]-[Residue][Number]

    Inputs:
        file_path : str
            Path to any docked PDB file containing REMARKs.

    Outputs:
        names : list[str]
            Annotated constraint names for user labeling.
    """
    names = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('REMARK'):
                parts = line.split()
                # Example line:
                # REMARK 666 MATCH TEMPLATE X CU 512 MATCH MOTIF A HIS 420 1 1
                # parts indices:  0     1   2     3        4 5  6   7   8  9  10 11 12
                # For the last line: REMARK 666 MATCH TEMPLATE Y LG1 1 MATCH MOTIF A HIS 498 5 1
                # parts:           [0]    [1]   [2]   [3]      [4] [5] [6]  [7] [8] [9] [10] [11] [12]

                # Determine ligand name
                ligand = parts[5]
                motif_resname = parts[10]
                motif_resnum = parts[11]

                if 'CU' in line:
                    name = f'Cu-{motif_resname.capitalize()}{motif_resnum}'
                elif 'LG1' in line:
                    name = f'AFB1-{motif_resname.capitalize()}{motif_resnum}'
                else:
                    # Default fallback, can be customized
                    name = f'{ligand}-{motif_resname.capitalize()}{motif_resnum}'

                names.append(name)
    return names


# # User defined elements and execution

# ========== USER DEFINED VARIABLES ==========
# Change these before running the script
originalWorkingDirectory = "/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/Laccases/C12 docking/"
suffix = "_C1" # Output file suffix for labeling results
structure_extension = '.cif' #'.pdb'

# ===== Start of the execution of the script === 
items = os.listdir(originalWorkingDirectory)
folders = [item for item in items if os.path.isdir(os.path.join(originalWorkingDirectory,item))]

counter = 0
for folder in tqdm(folders,desc=f"Measuring constraints"): 
    # ========== FOLDER LOCATIONS ==========
    #print("Original working directory:", originalWorkingDirectory)
    currentWorkingDirectory = os.path.join(originalWorkingDirectory, folder)
    #print("Current working directory:", currentWorkingDirectory)

    iterationFileName = f'{folder}{suffix}'  # Prefix for output files
    cstFileName = f'{folder}.enzdes.cst'
    items = os.listdir(currentWorkingDirectory)
    first_filename = [item for item in items if item.endswith(structure_extension)][0]
    userInputConstraintNames = extract_constraint_names(os.path.join(currentWorkingDirectory,first_filename))
    
    # ========== USAGE EXAMPLE ==========
    # Uncomment and implement the required functions before running!
    scoreFile = compile_scoreFiles(currentWorkingDirectory, iterationFileName)
    measurements_DataFrame = compile_measurements(currentWorkingDirectory, cstFileName, userInputConstraintNames, iterationFileName)
    counter +=1

