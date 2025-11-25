import os
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from itertools import combinations
from sklearn.metrics import r2_score
from scipy.stats import spearmanr


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
    df = pd.DataFrame(data, columns=headers)
    for col in headers:
        try:
            df[col] = pd.to_numeric(df[col])
        except ValueError:
            pass
    return df


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
        return None
    return ungapped_count_query - 1


def identify_indices_by_sequence(filename, ref_index_list, reference_sequence_name, fasta_path):
    target_id_substring = filename.split('_')[0]
    with open(fasta_path, 'r') as fasta_file:
        seqs = list(SeqIO.parse(fasta_file, 'fasta'))

    aligned_ref_seq = None
    aligned_query_seq = None

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


def fetch_experimental_result(gene_name, data_csv_path):
    df = pd.read_csv(data_csv_path)
    row = df[df['ALDH_#'] == gene_name]
    if not row.empty:
        return row.iloc[0]['log(Ratio)']
    else:
        return None


def process_folder_to_prediction_df(folder_path, index_list, ref_name, aligned_file_path, expdata_csv_path):
    pred_rows = []
    gene_names = []
    protocols = []
    filenames = []

    # Load experimental data once outside the loop, indexed by gene name
    exp_df = pd.read_csv(expdata_csv_path).set_index('ALDH_#')

    for fname in tqdm(os.listdir(folder_path), desc="Processing PDB files", leave=False):
        if fname.endswith('.pdb'):
            mapped_indices = identify_indices_by_sequence(fname, index_list, ref_name, aligned_file_path)
            with open(os.path.join(folder_path, fname), 'r') as f:
                pdb_lines = f.readlines()
            table_lines = extract_pose_energies_table(pdb_lines)
            df = parse_table(table_lines)

            df_valid = df.iloc[2:].copy()
            label_str_indices = df_valid['label'].str.split('_').str[1]
            label_indices = pd.to_numeric(label_str_indices, errors='coerce')

            pos_scores = {}
            for i, idx in enumerate(mapped_indices):
                rows_at_idx = df_valid[label_indices == idx]
                col_name = f"resi_{i+1}"
                if not rows_at_idx.empty:
                    pos_scores[col_name] = rows_at_idx['total'].values[0]
                else:
                    pos_scores[col_name] = np.nan

            gene_name = fname.split('_')[0]
            gene_names.append(gene_name)
            protocols.append(infer_protocol(fname))
            filenames.append(fname)

            pred_rows.append(pos_scores)

    pred_df = pd.DataFrame(pred_rows)
    pred_df['gene_name'] = gene_names
    pred_df['filename'] = filenames
    pred_df['protocol'] = protocols

    # Add experimental data column aligned by gene_name
    pred_df['experimental_log_ratio'] = pred_df['gene_name'].map(lambda g: exp_df.loc[g]['log(Ratio)'] if g in exp_df.index else np.nan)

    return pred_df


def correlation_combo_analysis(pred_df, max_combo_size=3, output_csv='combo_correlations.csv', combination_list=None):
    results = []
    unique_protocols = pred_df['protocol'].unique()

    for protocol in unique_protocols:
        pred_subset = pred_df[pred_df['protocol'] == protocol]
        pos_columns = pred_subset.columns.difference(['protocol', 'experimental_log_ratio', 'gene_name', 'filename'])

        if combination_list is None:
            combo_iter = []
            for k in range(1, max_combo_size + 1):
                combo_iter.extend(combinations(pos_columns, k))
        else:
            # Filter provided combinations for protocol columns - must be tuples of column names
            combo_iter = [combo for combo in combination_list if all(col in pos_columns for col in combo)]

        for cols in tqdm(combo_iter, desc=f"Analyzing combos {protocol}", leave=False):
            combo_df = pred_subset[list(cols)].apply(pd.to_numeric, errors='coerce')
            combo_pred = combo_df.sum(axis=1, min_count=1)

            valid_samples = combo_pred.dropna().index
            experimental_values = pred_subset.loc[valid_samples, 'experimental_log_ratio']
            experimental_values = pd.to_numeric(experimental_values, errors='coerce').dropna()

            idx_intersection = combo_pred.index.intersection(experimental_values.index)
            y_true = experimental_values.loc[idx_intersection]
            y_pred = combo_pred.loc[idx_intersection]

            if len(y_true) > 1:
                r2 = r2_score(y_true, y_pred)
                spearman_corr, _ = spearmanr(y_true, y_pred)
            else:
                r2, spearman_corr = np.nan, np.nan

            results.append({
                "protocol": protocol,
                "positions": ";".join(cols),
                "combo_size": len(cols),
                "R2": r2,
                "Spearman": spearman_corr
            })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by=['protocol', 'R2'], ascending=[True, False])
    results_df.to_csv(output_csv, index=False)
    return results_df


def main():
    main_path = '/Users/emluu/Documents/Siegel lab/Manuscript work/ALDH/ALDH_phase_2/'
    folder = os.path.join(main_path, 'relaxed_pdbs')
    indices = list(range(160, 198))  # residues 160 to 197
    reference_sequence_name = 'ALDH-002'
    experimental_csv = os.path.join(main_path, 'files_from_Li Lab', 'aldh_exp_results.csv')
    aligned_fasta_path = os.path.join(main_path, 'muscle_aldh_flexibility.aln-fasta')

    # Process folder and build prediction DataFrame
    pred_df = process_folder_to_prediction_df(folder, indices, reference_sequence_name, aligned_fasta_path, experimental_csv)
    # Load experimental data and align with predictions
    # experimental_df = pd.read_csv(experimental_csv).set_index('ALDH_#')
    # experimental_values = experimental_df.loc[pred_df.index, 'log(Ratio)']
    # Run correlation analysis on combinations of positions

    my_combos = [('resi_1',), ('resi_3','resi_5'), ('resi_2','resi_4','resi_6')]
    results_df = correlation_combo_analysis(pred_df, combination_list=my_combos,
                                            output_csv=os.path.join(main_path, 'combo_correlations.csv'))
    # results_df = correlation_combo_analysis(pred_df, max_combo_size=len(indices),
    #                                        output_csv=os.path.join(main_path, 'combo_correlations.csv'))

    print("Correlation combination analysis complete.")
    print(results_df.head(5))


if __name__ == "__main__":
    main()
