import os
import csv
import pandas as pd
from Bio import SeqIO

def extract_pose_energies_table(lines):
    in_table = False
    table_lines = []
    for line in lines:
        if '#BEGIN_POSE_ENERGIES_TABLE' in line:
            in_table = True
            continue
        if '#END_POSE_ENERGIES_TABLE' in line:
            break
        if in_table:
            table_lines.append(line.strip())
    return table_lines

def parse_table(table_lines):
    headers = []
    data = []
    for line in table_lines:
        if line.startswith('label'):
            headers = line.split()
        elif line and not line.startswith('#'):
            data.append(line.split())
    # Convert to DataFrame
    df = pd.DataFrame(data, columns=headers)
    # Optionally, convert all numeric columns to floats
    for col in headers:
        try:
            df[col] = pd.to_numeric(df[col])#df[col].astype(float)
        except ValueError:
            pass  # Ignore conversion error for non-numeric columns
    return df

def sum_total_for_indices(df, indices):
    df_valid = df.iloc[2:].copy()  # skip first two lines
    df_valid['total'] = pd.to_numeric(df_valid['total'], errors='coerce')
    # Extract the part after '_' safely, converting errors to NaN
    label_str_indices = df_valid['label'].str.split('_').str[1]
    label_indices = pd.to_numeric(label_str_indices, errors='coerce')
    # Filter only rows with valid numeric indices
    valid_rows = df_valid[label_indices.notna()]
    
    # Convert to integers for matching
    valid_label_indices = label_indices[label_indices.notna()].astype(int)

    # Select rows whose index in label is in provided indices
    selected_rows = valid_rows[valid_label_indices.isin(indices)]
    
    sum_total = selected_rows['total'].sum()
    return sum_total


def infer_protocol(filename):
    if 'fastrelaxlabmate' in filename:
        return 'FastRelaxCustom'
    elif 'fastrelaxstd' in filename:
        return 'FastRelaxStandard'
    else:
        return 'RelaxClassic'

def map_ref_to_query_position(ungapped_ref_pos, aligned_ref_seq, aligned_query_seq):
    ungapped_count = 0
    alignment_index = None
    for i, ch in enumerate(aligned_ref_seq):
        if ch != '-':
            if ungapped_count == ungapped_ref_pos:
                alignment_index = i
                break
            ungapped_count += 1
    if alignment_index is None:
        raise ValueError(f"Reference position {ungapped_ref_pos} not found in aligned sequence.")

    ungapped_count_query = 0
    for i in range(alignment_index + 1):
        if aligned_query_seq[i] != '-':
            ungapped_count_query += 1
    if aligned_query_seq[alignment_index] == '-':
        # Gap in query at this alignment position, no corresponding ungapped index
        return None
    return ungapped_count_query - 1

def identify_indices_by_sequence(filename, ref_index_list, reference_sequence_name, fasta_path):
    target_id_substring = filename.split('_')[0]  # adjust as needed for filename format
    
    with open(fasta_path, 'r') as fasta_file:
        seqs = list(SeqIO.parse(fasta_file, 'fasta'))

    aligned_ref_seq = None
    aligned_query_seq = None

    # Find reference sequence and query sequence by header
    for record in seqs:
        if reference_sequence_name in record.id:
            aligned_ref_seq = str(record.seq)
        if target_id_substring in record.id:
            aligned_query_seq = str(record.seq)

    if aligned_ref_seq is None:
        raise ValueError(f"Reference sequence '{reference_sequence_name}' not found in MSA.")
    if aligned_query_seq is None:
        raise ValueError(f"Query sequence matching '{target_id_substring}' not found in MSA.")

    index_list = []
    for ref_pos in ref_index_list:
        query_pos = map_ref_to_query_position(ref_pos, aligned_ref_seq, aligned_query_seq)
        if query_pos is not None:
            index_list.append(query_pos)

    return index_list

def fetch_experimental_result(gene_name,data_csv_path):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(data_csv_path)
    
    # Find the row where ALDH_# equals gene_name
    row = df[df['ALDH_#'] == gene_name]
    
    if not row.empty:
        # Extract and return the log(Ratio) value
        return row.iloc[0]['log(Ratio)']
    else:
        # Handle the case where the gene_name is not found
        return None

def process_folder(output_csv,folder_path, index_list,ref_name,aligned_file_path, expdata_csv_path):
    rows = []
    for fname in os.listdir(folder_path):
        if fname.endswith('.pdb'):
            mapped_indices = identify_indices_by_sequence(fname, index_list, 
                                                          ref_name,aligned_file_path)
            with open(os.path.join(folder_path, fname), 'r') as f:
                pdb_lines = f.readlines()
            table_lines = extract_pose_energies_table(pdb_lines)
            dataframe = parse_table(table_lines)
            score = sum_total_for_indices(dataframe, mapped_indices)
            protocol = infer_protocol(fname)
            full_score = dataframe[dataframe.label == 'pose']['total'].values[0]
            gene_name = fname.split('_')[0]
            experimental_result = fetch_experimental_result(gene_name,expdata_csv_path)  
            rows.append({'protocol': protocol,'gene': gene_name, 'log(Ratio)': experimental_result,
                         'filename': fname, 'whole_total': full_score,
                         'sub_total': score, 'indices':mapped_indices})

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['protocol', 'gene','log(Ratio)','filename',
                      'whole_total', 'sub_total','indices']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

# Usage
main_path = '/Users/emluu/Documents/Siegel lab/Manuscript work/ALDH/ALDH_phase_2/'
folder = os.path.join(main_path,'relaxed_pdbs')
indices = [165,166,167,168,169,170,171,172,173,174,175]  # specify the residue indices 
reference_sequence_name = 'ALDH-002'
output_csv = os.path.join(main_path,'compiled_pose_energies.csv')
experimental_csv = os.path.join(main_path,'files_from_Li Lab','aldh_exp_results.csv')
aligned_fasta_path = os.path.join(main_path,'muscle_aldh_flexibility.aln-fasta')
process_folder(output_csv, folder, indices, reference_sequence_name, 
               aligned_fasta_path, experimental_csv)
