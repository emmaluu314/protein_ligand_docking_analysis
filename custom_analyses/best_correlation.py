import os
import csv
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from itertools import combinations
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso, LassoCV, LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
import matplotlib.pyplot as plt


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
    pred_df['experimental_log_ratio'] = pred_df['gene_name'].map(
        lambda g: exp_df.loc[g]['log(Ratio)'] if g in exp_df.index else np.nan)

    return pred_df


def correlation_combo_analysis(pred_df, max_combo_size=3, output_csv='combo_correlations.csv',
                               allowed_positions_by_protocol=None, plot_dir=None):
    results = []
    unique_protocols = pred_df['protocol'].unique()
    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)

    for protocol in unique_protocols:
        # Filter by protocol
        pred_subset = pred_df[pred_df['protocol'] == protocol]
        # All columns except metadata columns
        pos_columns = pred_subset.columns.difference(['protocol', 'experimental_log_ratio', 'gene_name', 'filename'])
        if allowed_positions_by_protocol and protocol in allowed_positions_by_protocol:
            pos_columns = pos_columns.intersection(allowed_positions_by_protocol[protocol])
        pos_columns = list(pos_columns)
        for k in tqdm(range(1, max_combo_size + 1), desc=f"Analyzing combos {protocol}", leave=False):
            for cols in tqdm(combinations(pos_columns, k), desc=f"Combo size {k}", leave=False):
                # Select columns and convert all to numeric, coerce errors to NaN
                combo_df = pred_subset[list(cols)].apply(pd.to_numeric, errors='coerce')
                # Sum across columns with min_count=1: if all NaN, sum is NaN, not 0
                combo_pred = combo_df.sum(axis=1, min_count=1)
                # Build aligned frame so gene_name stays attached to true/pred values
                aligned = pd.DataFrame({
                    "gene_name": pred_subset["gene_name"],
                    "experimental_log_ratio": pd.to_numeric(
                        pred_subset["experimental_log_ratio"], errors="coerce"
                    ),
                    "combo_pred": combo_pred,
                }).dropna(subset=["experimental_log_ratio", "combo_pred"])
                y_true = aligned["experimental_log_ratio"]
                y_pred = aligned["combo_pred"]

                if len(y_true) > 1:
                    r2 = r2_score(y_true, y_pred)
                    spearman_corr, _ = spearmanr(y_true, y_pred)
                    # Fit a simple linear regression (y_pred ~ experimental) and capture its R^2
                    lr_model = LinearRegression()
                    X = y_true.values.reshape(-1, 1)
                    lr_model.fit(X, y_pred.values)
                    r2_fit = lr_model.score(X, y_pred.values)
                    if plot_dir:
                        positions_label = "_".join(cols)
                        fig, ax = plt.subplots(figsize=(4.5, 4.5))
                        ax.scatter(y_true, y_pred, alpha=0.7, s=25)
                        diag_min = min(y_true.min(), y_pred.min())
                        diag_max = max(y_true.max(), y_pred.max())
                        ax.plot([diag_min, diag_max], [diag_min, diag_max], color="gray", linestyle="--", linewidth=1)
                        ax.set_xlabel("Experimental log ratio")
                        ax.set_ylabel("Predicted sum")
                        ax.set_title(f"{protocol} | {positions_label}")
                        ax.grid(alpha=0.2)
                        fname = f"{protocol}_k{k}_{positions_label}.png".replace(os.sep, "_")
                        fig.tight_layout()
                        fig.savefig(os.path.join(plot_dir, fname), dpi=150)
                        plt.close(fig)
                else:
                    r2, spearman_corr, r2_fit = np.nan, np.nan, np.nan

                results.append({
                    "protocol": protocol,
                    "positions": ";".join(cols),
                    "combo_size": k,
                    "R2": r2,
                    "Spearman": spearman_corr,
                    "R2_fit": r2_fit
                })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by=['protocol', 'R2'], ascending=[True, False])
    results_df.to_csv(output_csv, index=False)
    return results_df


def lasso_preselection(pred_df, output_csv='lasso_coefficients.csv', alphas=None, cv=5, max_iter=20000,
                       bootstrap_reps=0, random_state=0):
    """
    Fit LASSO per protocol and write residue coefficients to CSV.
    Rows with any NaNs in the position columns or target are dropped.
    Optionally estimate stability via bootstrap frequency of non-zero coefficients.
    """
    rows = []
    for protocol in pred_df['protocol'].unique():
        subset = pred_df[pred_df['protocol'] == protocol]
        pos_cols = subset.columns.difference(['protocol', 'experimental_log_ratio', 'gene_name', 'filename'])
        X = subset[pos_cols].apply(pd.to_numeric, errors='coerce')
        # Drop columns that are entirely NaN to avoid singular features
        X = X.loc[:, ~X.isna().all()]
        if X.empty:
            continue

        y = pd.to_numeric(subset['experimental_log_ratio'], errors='coerce')
        valid_mask = (~X.isna().any(axis=1)) & (~y.isna())
        X_valid = X.loc[valid_mask]
        y_valid = y.loc[valid_mask]

        if len(X_valid) < 2:
            continue

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_valid)

        actual_cv = min(cv, len(y_valid))
        lasso = LassoCV(alphas=alphas, cv=actual_cv, max_iter=max_iter).fit(X_scaled, y_valid)
        coefs = pd.Series(lasso.coef_, index=X.columns)

        # Report validation MSE (CV) at chosen alpha and train-set R2 for visibility
        alpha_idx = list(lasso.alphas_).index(lasso.alpha_)
        cv_mse_mean = lasso.mse_path_[alpha_idx].mean()
        cv_mse_std = lasso.mse_path_[alpha_idx].std()
        train_r2 = r2_score(y_valid, lasso.predict(X_scaled))

        # Optional bootstrap stability: frequency of non-zero selection at fixed alpha
        selection_freq = pd.Series(0, index=X.columns, dtype=float)
        if bootstrap_reps > 0:
            rng = np.random.default_rng(random_state)
            for _ in range(bootstrap_reps):
                idx = rng.integers(0, len(y_valid), size=len(y_valid))
                X_boot = X_scaled[idx]
                y_boot = y_valid.iloc[idx]
                lasso_boot = Lasso(alpha=lasso.alpha_, max_iter=max_iter).fit(X_boot, y_boot)
                selection_freq += (lasso_boot.coef_ != 0)
            selection_freq /= bootstrap_reps

        for position, coef in coefs.items():
            rows.append({
                "protocol": protocol,
                "position": position,
                "coefficient": coef,
                "selection_freq": selection_freq[position] if bootstrap_reps > 0 else np.nan
            })

    coef_df = pd.DataFrame(rows)
    coef_df.to_csv(output_csv, index=False)
    return coef_df


def build_allowed_positions_from_lasso(coef_df, coef_threshold=0.0, selection_freq_min=0.0, top_n=None):
    """
    Build a dict of protocol -> list of positions to keep based on LASSO output.
    - coef_threshold: keep residues with |coef| > threshold
    - selection_freq_min: if present, keep residues with selection_freq >= this
    - top_n: if set, cap to top N by |coef| per protocol after thresholds
    """
    allowed = {}
    for protocol, group in coef_df.groupby('protocol'):
        filtered = group[np.abs(group['coefficient']) > coef_threshold]
        if 'selection_freq' in group.columns:
            filtered = filtered[(filtered['selection_freq'].isna()) | (filtered['selection_freq'] >= selection_freq_min)]
        filtered = filtered.reindex(filtered['coefficient'].abs().sort_values(ascending=False).index)
        if top_n is not None:
            filtered = filtered.head(top_n)
        allowed[protocol] = filtered['position'].tolist()
    return allowed


def randomized_combo_search(pred_df, allowed_positions_by_protocol, samples_per_k=1000, max_combo_size=4,
                            random_state=0, output_csv='random_combo_correlations.csv'):
    """
    Sample random combinations of allowed positions per protocol and score correlations.
    samples_per_k controls how many combos per size; capped by total possible combos.
    """
    rng = np.random.default_rng(random_state)
    results = []

    for protocol in pred_df['protocol'].unique():
        subset = pred_df[pred_df['protocol'] == protocol]
        pos_columns = subset.columns.difference(['protocol', 'experimental_log_ratio', 'gene_name', 'filename'])
        if allowed_positions_by_protocol and protocol in allowed_positions_by_protocol:
            pos_columns = pos_columns.intersection(allowed_positions_by_protocol[protocol])
        pos_columns = list(pos_columns)
        if not pos_columns:
            continue

        X = subset[pos_columns].apply(pd.to_numeric, errors='coerce').to_numpy()
        y = pd.to_numeric(subset['experimental_log_ratio'], errors='coerce').to_numpy()
        n_cols = len(pos_columns)

        for k in range(1, min(max_combo_size, n_cols) + 1):
            total_combos = math.comb(n_cols, k)
            sample_count = min(samples_per_k, total_combos)
            if sample_count == 0:
                continue

            if sample_count == total_combos:
                combo_indices = list(combinations(range(n_cols), k))
            else:
                seen = set()
                combo_indices = []
                while len(combo_indices) < sample_count:
                    idxs = tuple(sorted(rng.choice(n_cols, size=k, replace=False)))
                    if idxs not in seen:
                        seen.add(idxs)
                        combo_indices.append(idxs)

            for idxs in tqdm(combo_indices, desc=f"Random combos k={k} {protocol}", leave=False):
                arr = X[:, idxs]
                combo_sum = np.nansum(arr, axis=1)
                all_nan_mask = np.isnan(arr).all(axis=1)
                combo_sum[all_nan_mask] = np.nan

                valid_mask = (~np.isnan(combo_sum)) & (~np.isnan(y))
                y_true = y[valid_mask]
                y_pred = combo_sum[valid_mask]

                if len(y_true) > 1:
                    r2 = r2_score(y_true, y_pred)
                    spearman_corr, _ = spearmanr(y_true, y_pred)
                else:
                    r2, spearman_corr = np.nan, np.nan

                positions = [pos_columns[i] for i in idxs]
                results.append({
                    "protocol": protocol,
                    "positions": ";".join(positions),
                    "combo_size": k,
                    "R2": r2,
                    "Spearman": spearman_corr
                })

    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df = results_df.sort_values(by=['protocol', 'R2'], ascending=[True, False])
        results_df.to_csv(output_csv, index=False)
    return results_df


def main():
    main_path = '/Users/emluu/Documents/Siegel lab/Manuscript work/ALDH/ALDH_phase_2/'
    folder = os.path.join(main_path, 'relaxed_pdbs_curated')
    indices = list(range(160, 198))  # residues 160 to 197
    reference_sequence_name = 'ALDH-002'
    experimental_csv = os.path.join(main_path, 'files_from_Li Lab', 'aldh_exp_results.csv')
    aligned_fasta_path = os.path.join(main_path, 'muscle_aldh_flexibility.aln-fasta')

    # Process folder and build prediction DataFrame
    pred_df = process_folder_to_prediction_df(folder, indices, reference_sequence_name, aligned_fasta_path, experimental_csv)

    # Run LASSO preselection and save coefficients
    k = 5
    lasso_output_csv = os.path.join(main_path, f'lasso_coefficients_k{k}.csv')
    coef_df = lasso_preselection(pred_df, output_csv=lasso_output_csv, cv=k, bootstrap_reps=50, random_state=0)

    # Build allowed positions per protocol from LASSO (tunable thresholds)
    allowed_positions = build_allowed_positions_from_lasso(
        coef_df,
        coef_threshold=0.0,      # keep non-zero coefficients
        selection_freq_min=0.5,  # require selection in at least 50% of bootstraps
        top_n=None               # limit to top N positions by |coef|; set None to keep all passing thresholds
    )  

    # Set max combo size to the largest allowed list (fallback to full index length if empty)
    max_combo_size = max((len(v) for v in allowed_positions.values()), default=len(indices))

    # Randomized search over allowed positions (controls runtime)
    random_results = randomized_combo_search(
        pred_df,
        allowed_positions_by_protocol=allowed_positions,
        samples_per_k=2000,
        max_combo_size=max_combo_size,
        random_state=0,
        output_csv=os.path.join(main_path, 'random_combo_correlations.csv')
    )
    if not random_results.empty:
        print("Top randomized combos (by R2):")
        print(random_results.head(10))

    # Optional exhaustive search (may be expensive)
    run_exhaustive = False
    if run_exhaustive:
        results_df = correlation_combo_analysis(
            pred_df,
            max_combo_size=max_combo_size,
            output_csv=os.path.join(main_path, 'combo_correlations.csv'),
            allowed_positions_by_protocol=allowed_positions
        )
        print("Correlation combination analysis complete.")
        print(results_df.head(5))

if __name__ == "__main__":
    main()
