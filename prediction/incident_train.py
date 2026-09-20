import sys
import os
import gc
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from lightgbm import LGBMClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import random
from itertools import product

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

NB_CPUS       = 10
NB_PARAMS     = 100
MY_SEED       = 2026
TOP_K_FEATURES = 50

params_dict = {
    'n_estimators':    [100, 200, 300, 400, 500],
    'max_depth':       np.linspace(5, 30, 6).astype('int32').tolist(),
    'num_leaves':      np.linspace(5, 30, 6).astype('int32').tolist(),
    'subsample':       np.linspace(0.6, 1, 9).tolist(),
    'learning_rate':   [0.1, 0.05, 0.01, 0.001],
    'colsample_bytree': np.linspace(0.6, 1, 9).tolist()
}
PARAM_KEYS = list(params_dict.keys())


def select_params_combo(my_dict, nb_items, my_seed):
    combo_list = [dict(zip(my_dict.keys(), v)) for v in product(*my_dict.values())]
    random.seed(my_seed)
    return random.sample(combo_list, nb_items)


def normal_imp(mydict):
    total = sum(mydict.values())
    return {k: v / total for k, v in mydict.items()} if total > 0 else mydict


def get_best_params(mydf, my_f_lst, my_params_lst, my_seed):
    X, y = mydf[my_f_lst], mydf.target_y
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=my_seed)
    my_params_res_lst = []
    for my_params in my_params_lst:
        auc_cv_lst = []
        my_params0 = my_params.copy()
        for train_idx, val_idx in skf.split(X, y):
            in_cv_X_train, in_cv_y_train = X.iloc[train_idx], y.iloc[train_idx]
            in_cv_X_test,  in_cv_y_test  = X.iloc[val_idx],   y.iloc[val_idx]
            my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False,
                                    verbosity=-1, seed=my_seed, n_jobs=1)
            my_lgb.set_params(**my_params)
            my_lgb.fit(in_cv_X_train, in_cv_y_train)
            y_pred_prob = my_lgb.predict_proba(in_cv_X_test)[:, 1]
            try:
                auc_cv_lst.append(roc_auc_score(in_cv_y_test, y_pred_prob))
            except ValueError:
                auc_cv_lst.append(0.5)
        my_params0['AUC_cv_MEAN'] = np.round(np.mean(auc_cv_lst), 5)
        my_params_res_lst.append(my_params0)
    my_params_res_df = pd.DataFrame(my_params_res_lst)
    my_params_res_df.sort_values(by='AUC_cv_MEAN', ascending=False, inplace=True)
    best_param = dict(my_params_res_df.iloc[0][PARAM_KEYS])
    best_param['n_estimators'] = int(best_param['n_estimators'])
    best_param['max_depth']    = int(best_param['max_depth'])
    best_param['num_leaves']   = int(best_param['num_leaves'])
    return best_param


def get_expo_f_lst(mydf, f_lst, my_seed):
    X_train, y_train = mydf[f_lst], mydf.target_y
    lgb_full = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False,
                               verbosity=-1, seed=my_seed, n_jobs=1)
    lgb_full.fit(X_train, y_train)
    gain_imp  = dict(zip(lgb_full.feature_name_, lgb_full.booster_.feature_importance(importance_type='gain')))
    cover_imp = dict(zip(lgb_full.feature_name_, lgb_full.booster_.feature_importance(importance_type='split')))
    top_features = sorted(gain_imp, key=gain_imp.get, reverse=True)[:TOP_K_FEATURES]
    return top_features, gain_imp, cover_imp


def model_train_pred(traindf, testdf, f_lst, my_params, my_seed):
    X_train, y_train = traindf[f_lst], traindf.target_y
    X_test = testdf[f_lst]
    my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False,
                             verbosity=-1, seed=my_seed, n_jobs=1)
    my_lgb.set_params(**my_params)
    calibrate = CalibratedClassifierCV(my_lgb, method='isotonic', cv=5)
    calibrate.fit(X_train, y_train)
    y_pred = calibrate.predict_proba(X_test)[:, 1].tolist()
    return y_pred


def get_iter_predictions(mydf, full_expo_f_lst, cov_f_lst, fold_id, candidate_params_lst, my_seed):
    traindf = mydf.loc[mydf.Custom_Region != fold_id].reset_index(drop=True)
    testdf  = mydf.loc[mydf.Custom_Region == fold_id].reset_index(drop=True)

    top_expo_f_lst, gain_imp, cover_imp = get_expo_f_lst(traindf, full_expo_f_lst, my_seed)

    params_expo     = get_best_params(traindf, top_expo_f_lst,              candidate_params_lst, my_seed)
    params_cov      = get_best_params(traindf, cov_f_lst,                   candidate_params_lst, my_seed)
    params_expo_cov = get_best_params(traindf, top_expo_f_lst + cov_f_lst,  candidate_params_lst, my_seed)

    y_pred_expo     = model_train_pred(traindf, testdf, top_expo_f_lst,             params_expo,     my_seed)
    y_pred_cov      = model_train_pred(traindf, testdf, cov_f_lst,                  params_cov,      my_seed)
    y_pred_expo_cov = model_train_pred(traindf, testdf, top_expo_f_lst + cov_f_lst, params_expo_cov, my_seed)

    return (fold_id, gain_imp, cover_imp,
            testdf.eid.tolist(), testdf.target_y.tolist(),
            y_pred_expo, y_pred_cov, y_pred_expo_cov)


if __name__ == "__main__":
    print("--- Starting Incident Analysis Pipeline ---")
    candidate_params_lst = select_params_combo(params_dict, NB_PARAMS, MY_SEED)

    cov_cols = ['eid', 'Age_p53_i0', 'Sex', 'AccCeni0']
    cov_df = pd.read_csv(COV_FILE, sep='\t', usecols=cov_cols)
    cov_df['eid']       = cov_df['eid'].astype(int)
    cov_df['Sex']       = cov_df['Sex'].map({'Male': 1, 'Female': 0})
    cov_df['Ethnicity'] = (cov_df['Ethnicity'] == 'White').astype(int)
    cov_df['AccCeni0']  = LabelEncoder().fit_transform(cov_df['AccCeni0'].astype(str))
    cov_f_lst = [c for c in cov_cols if c != 'eid']

    time_df = pd.read_csv(TIME_FILE, sep='\t', usecols=['eid', 'p53_i0'])
    fold_df = pd.read_csv(FOLD_FILE, usecols=['eid', 'Custom_Region'])
    base_static_df = pd.merge(cov_df, time_df, on='eid', how='left').merge(fold_df, on='eid', how='left')
    del cov_df, time_df, fold_df
    gc.collect()

    summary_df      = pd.read_csv(SUMMARY_FILE)
    new_disease_df  = pd.read_excel(NEW_FILTER_FILE, usecols=['NAME', 'SEX'])
    valid_tgts_new  = set(new_disease_df['NAME'].astype(str))
    sex_map         = dict(zip(new_disease_df['NAME'].astype(str), new_disease_df['SEX']))

    valid_tgts_summary = summary_df[summary_df['case_new_p53_i0'] >= 200]['file_name'].astype(str).tolist()
    valid_tgts  = [t for t in valid_tgts_summary if t in valid_tgts_new]
    all_files   = [os.path.join(DIAG_DIR, f + '.csv') for f in valid_tgts
                   if os.path.exists(os.path.join(DIAG_DIR, f + '.csv'))]

    chunk_id      = int(sys.argv[1]) if len(sys.argv) >= 3 else 0
    total_chunks  = int(sys.argv[2]) if len(sys.argv) >= 3 else 1
    chunk_size    = (len(all_files) + total_chunks - 1) // total_chunks
    diag_files_chunk = all_files[chunk_id * chunk_size : min((chunk_id + 1) * chunk_size, len(all_files))]

    expo_file_list = sorted([f for f in os.listdir(EXPO_DIR) if f.endswith('.csv')])

    for expo_file_name in expo_file_list:
        expo_name = os.path.splitext(expo_file_name)[0].replace('Category_', '')
        print(f"\nProcessing Expo: {expo_name}")

        expo_df = pd.read_csv(os.path.join(EXPO_DIR, expo_file_name))
        expo_df['eid'] = expo_df['eid'].astype(int)
        full_expo_f_lst = [c for c in expo_df.columns if c != 'eid']
        expo_df[full_expo_f_lst] = expo_df[full_expo_f_lst].astype(np.float32)

        base_df = pd.merge(base_static_df, expo_df, how='inner', on='eid')
        del expo_df
        gc.collect()

        fold_id_lst = sorted(base_df['Custom_Region'].dropna().unique().astype(int))

        out_imp_dir  = os.path.join(DPATH, f'Results/IncidentPrediction/{expo_name}_Importance/')
        out_pred_dir = os.path.join(DPATH, f'Results/IncidentPrediction/{expo_name}_Predictions/')
        for d in [out_imp_dir, out_pred_dir]:
            os.makedirs(d, exist_ok=True)

        finished_files    = set(os.listdir(out_pred_dir))
        current_diag_files = [f for f in diag_files_chunk if os.path.basename(f) not in finished_files]

        for diag_file in tqdm(current_diag_files, desc=f"Targets ({expo_name})"):
            tgt    = os.path.basename(diag_file)[:-4]
            tgt_df = pd.read_csv(diag_file)
            tgt_df['eid'] = tgt_df['eid'].astype(int)

            tmp_df = pd.merge(base_df, tgt_df, how='inner', on='eid')
            if len(tmp_df) < 200:
                continue

            tmp_df['date']   = pd.to_datetime(tmp_df['date'],   errors='coerce')
            tmp_df['p53_i0'] = pd.to_datetime(tmp_df['p53_i0'], errors='coerce')
            tmp_df['target_y'] = np.where(
                (tmp_df['status'] == 1) & (tmp_df['date'] > tmp_df['p53_i0']), 1,
                np.where(tmp_df['status'] == 0, 0, np.nan)
            )
            tmp_df.dropna(subset=['target_y', 'Custom_Region'], inplace=True)
            if tmp_df['target_y'].sum() < 200:
                continue

            curr_sex_id  = sex_map.get(tgt, 0)
            run_cov_f_lst = [c for c in cov_f_lst if c != 'Sex'] if curr_sex_id in [1, 2] else list(cov_f_lst)

            fold_results = Parallel(n_jobs=NB_CPUS)(
                delayed(get_iter_predictions)(
                    tmp_df.reset_index(drop=True), full_expo_f_lst, run_cov_f_lst,
                    fold_id, candidate_params_lst, MY_SEED
                ) for fold_id in fold_id_lst
            )

            imp_df_gain  = pd.DataFrame({'Biomarker_code': full_expo_f_lst})
            imp_df_cover = pd.DataFrame({'Biomarker_code': full_expo_f_lst})
            all_res = {
                'eid': [], 'target_y': [],
                f'y_pred_{expo_name}': [], 'y_pred_cov': [], f'y_pred_{expo_name}_cov': []
            }

            for res in fold_results:
                fold_id, gain_imp, cover_imp = res[0], res[1], res[2]

                gain_norm  = normal_imp(gain_imp)
                cover_norm = normal_imp(cover_imp)

                tg_fold_df = pd.DataFrame({'Biomarker_code': list(gain_norm.keys()),
                                           f'Gain_fold{fold_id}': list(gain_norm.values())})
                tc_fold_df = pd.DataFrame({'Biomarker_code': list(cover_norm.keys()),
                                           f'Cover_fold{fold_id}': list(cover_norm.values())})

                imp_df_gain  = pd.merge(imp_df_gain,  tg_fold_df, on='Biomarker_code', how='left')
                imp_df_cover = pd.merge(imp_df_cover, tc_fold_df, on='Biomarker_code', how='left')

                all_res['eid'].extend(res[3])
                all_res['target_y'].extend(res[4])
                all_res[f'y_pred_{expo_name}'].extend(res[5])
                all_res['y_pred_cov'].extend(res[6])
                all_res[f'y_pred_{expo_name}_cov'].extend(res[7])

            gain_cols  = [c for c in imp_df_gain.columns  if c.startswith('Gain_fold')]
            cover_cols = [c for c in imp_df_cover.columns if c.startswith('Cover_fold')]
            imp_df_gain['TotalGain_cv']   = imp_df_gain[gain_cols].mean(axis=1)
            imp_df_cover['TotalCover_cv'] = imp_df_cover[cover_cols].mean(axis=1)

            imp_save_df = pd.merge(
                imp_df_gain[['Biomarker_code'] + gain_cols + ['TotalGain_cv']],
                imp_df_cover[['Biomarker_code'] + cover_cols + ['TotalCover_cv']],
                on='Biomarker_code', how='outer'
            ).fillna(0).sort_values('TotalGain_cv', ascending=False)

            imp_save_df.to_csv(os.path.join(out_imp_dir, tgt + '.csv'), index=False)
            pd.DataFrame(all_res).to_csv(os.path.join(out_pred_dir, tgt + '.csv'), index=False)

            del tmp_df, fold_results, all_res, imp_df_gain, imp_df_cover
            gc.collect()

        del base_df
        gc.collect()

    print("--- All Processing Completed ---")