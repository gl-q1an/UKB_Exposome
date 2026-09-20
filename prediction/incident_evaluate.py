import glob
import re
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
from sklearn.metrics import roc_auc_score, roc_curve
from joblib import Parallel, delayed

# Limit thread usage per process
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

warnings.filterwarnings('ignore')

BASE_ANALYSIS_DIR = os.path.join(DPATH, 'Results/IncidentPrediction/')

NB_CPUS = 10
NB_ITERS = 1000

def sort_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c.replace("_", "")) for c in re.split('([0-9]+)', key)]
    l.sort(key=alphanum_key)
    return l

def Find_Optimal_Cutoff(target, predicted):
    fpr, tpr, thresholds = roc_curve(target, predicted)
    specificity = 1 - fpr
    ix = np.argmin(np.abs(tpr - specificity))
    return thresholds[ix]

def get_metrics_numpy(y_test, pred_prob, cutoff):
    pred_binary = (pred_prob >= cutoff).astype(int)

    tp = np.sum((y_test == 1) & (pred_binary == 1))
    tn = np.sum((y_test == 0) & (pred_binary == 0))
    fp = np.sum((y_test == 0) & (pred_binary == 1))
    fn = np.sum((y_test == 1) & (pred_binary == 0))

    epsilon = 1e-10

    acc = (tp + tn) / (tp + tn + fp + fn + epsilon)
    sens = tp / (tp + fn + epsilon)
    spec = tn / (tn + fp + epsilon)
    prec = tp / (tp + fp + epsilon)
    Youden = sens + spec - 1
    f1 = 2 * prec * sens / (prec + sens + epsilon)

    try:
        auc = roc_auc_score(y_test, pred_prob)
    except:
        auc = 0.5

    return [auc, acc, sens, spec, prec, Youden, f1]

def get_avg_output_optimized(mydf, gt_col, pred_col, cutoff, nb_iters):
    y_true = mydf[gt_col].values
    y_pred = mydf[pred_col].values
    n_samples = len(y_true)

    results = []

    # Bootstrap Loop
    for i in range(nb_iters):
        idx = np.random.randint(0, n_samples, n_samples)
        y_t_bt = y_true[idx]
        y_p_bt = y_pred[idx]

        if len(np.unique(y_t_bt)) < 2:
            res = [np.nan] * 7
        else:
            res = get_metrics_numpy(y_t_bt, y_p_bt, cutoff)
        results.append(res)

    result_matrix = np.array(results)
    result_matrix = result_matrix[~np.isnan(result_matrix).any(axis=1)]

    if len(result_matrix) == 0:
        return pd.Series([np.nan] * 7)

    medians = np.median(result_matrix, axis=0)
    lbds = np.percentile(result_matrix, 2.5, axis=0)
    ubds = np.percentile(result_matrix, 97.5, axis=0)

    output_lst = []
    for i in range(7):
        output_lst.append('{:.3f}'.format(medians[i]) + ' [' +
                          '{:.3f}'.format(lbds[i]) + ' - ' +
                          '{:.3f}'.format(ubds[i]) + ']')

    cols = ['AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Youden-index', 'F1-score']
    return pd.Series(output_lst, index=cols)

def process_single_file(file_path, omics_name, output_dir):
    tgt = os.path.basename(file_path)[:-4]

    try:
        tgt_pred_df = pd.read_csv(file_path)

        if len(tgt_pred_df) < 10:
            return None
        if tgt_pred_df['target_y'].nunique() < 2:
            return None

        # Detect which prediction columns are available
        col_omics = f'y_pred_{omics_name}'
        col_cov = 'y_pred_cov'
        col_combined = f'y_pred_{omics_name}_cov'

        results_list = []
        index_list = []

        # Omics only
        if col_omics in tgt_pred_df.columns:
            ct = Find_Optimal_Cutoff(tgt_pred_df.target_y, tgt_pred_df[col_omics])
            res = get_avg_output_optimized(tgt_pred_df, 'target_y', col_omics, ct, nb_iters=NB_ITERS)
            results_list.append(res)
            index_list.append(omics_name)

        # Covariates only
        if col_cov in tgt_pred_df.columns:
            ct = Find_Optimal_Cutoff(tgt_pred_df.target_y, tgt_pred_df[col_cov])
            res = get_avg_output_optimized(tgt_pred_df, 'target_y', col_cov, ct, nb_iters=NB_ITERS)
            results_list.append(res)
            index_list.append('Covariates')

        # Combined
        if col_combined in tgt_pred_df.columns:
            ct = Find_Optimal_Cutoff(tgt_pred_df.target_y, tgt_pred_df[col_combined])
            res = get_avg_output_optimized(tgt_pred_df, 'target_y', col_combined, ct, nb_iters=NB_ITERS)
            results_list.append(res)
            index_list.append(f'{omics_name}+Covariates')

        if not results_list:
            print(f"[Skip] {tgt}: no recognized prediction columns found")
            return None

        res_df = pd.concat(results_list, axis=1).T
        res_df.index = index_list

        out_path = os.path.join(output_dir, tgt + '.csv')
        res_df.to_csv(out_path, index=True)
        return tgt

    except Exception as e:
        print(f"[Error] {tgt}: {e}")
        return None

if __name__ == "__main__":
    if not os.path.exists(BASE_ANALYSIS_DIR):
        sys.exit(0)

    all_pred_dirs = sorted([d for d in os.listdir(BASE_ANALYSIS_DIR)
                            if d.endswith('_Predictions')
                            and os.path.isdir(os.path.join(BASE_ANALYSIS_DIR, d))])

    for pred_folder_name in all_pred_dirs:
        omics_name = pred_folder_name.replace('_Predictions', '')
        input_dir = os.path.join(BASE_ANALYSIS_DIR, pred_folder_name)
        output_dir = os.path.join(BASE_ANALYSIS_DIR, f'{omics_name}_Metrics/')

        os.makedirs(output_dir, exist_ok=True)

        search_pattern = os.path.join(input_dir, '*.csv')
        tgt_dir_lst = sort_nicely(glob.glob(search_pattern))

        finished_files = set([f[:-4] for f in os.listdir(output_dir) if f.endswith('.csv')])
        todo_files = [f for f in tgt_dir_lst if os.path.basename(f)[:-4] not in finished_files]

        if not todo_files:
            continue

        if len(sys.argv) >= 3:
            chunk_id = int(sys.argv[1])
            total_chunks = int(sys.argv[2])

            chunk_size = (len(todo_files) + total_chunks - 1) // total_chunks
            start_idx = chunk_id * chunk_size
            end_idx = min((chunk_id + 1) * chunk_size, len(todo_files))

            current_batch = todo_files[start_idx:end_idx]
        else:
            current_batch = todo_files

        if not current_batch:
            continue

        print(f"Processing {omics_name}: {len(current_batch)} files")

        Parallel(n_jobs=NB_CPUS, backend='loky')(
            delayed(process_single_file)(f, omics_name, output_dir) for f in tqdm(current_batch)
        )