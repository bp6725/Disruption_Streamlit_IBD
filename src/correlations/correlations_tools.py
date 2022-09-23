from data_api import DataAPI
import pandas as pd
import scipy
import numpy as np
from statsmodels.stats import multitest,weightstats
from numba import njit



class CorrelationsTools():
    def __init__(self):
        self.data_api = DataAPI()

    def build_correlation_dfs(self ,med_sig_df ,scd_sig_df ,union_ref_df ,med_with_nutriens ,scd_with_nutriens, \
                              before_med_sig_df_with_nutriens ,before_scd_sig_df_with_nutriens ,all_refrence_with_nuts, alpha, threshold):
        med_corr, med_asymptotic_p = self.calc_spearman(med_sig_df)
        scd_corr, scd_asymptotic_p = self.calc_spearman(scd_sig_df)
        backbone_corr, backbone_asymptotic_p = self.calc_spearman(union_ref_df)

        with_nut_med_corr, with_nut_med_asymptotic_p = self.calc_spearman(med_with_nutriens)
        with_nut_scd_corr, with_nut_scd_asymptotic_p = self.calc_spearman(scd_with_nutriens)

        with_nut_before_med_corr, with_nut_before_med_asymptotic_p = self.calc_spearman(before_med_sig_df_with_nutriens)
        with_nut_before_scd_corr, with_nut_before_scd_asymptotic_p = self.calc_spearman(before_scd_sig_df_with_nutriens)
        with_nut_reference_corr, with_nut_reference_asymptotic_p = self.calc_spearman(all_refrence_with_nuts)

        rejected_med, _, _, _ = multitest.multipletests(med_asymptotic_p[0], alpha=alpha, method='fdr_tsbh')
        rejected_scd, _, _, _ = multitest.multipletests(scd_asymptotic_p[0], alpha=alpha, method='fdr_tsbh')
        rejected_backbone, _, _, _ = multitest.multipletests(backbone_asymptotic_p[0], alpha=alpha, method='fdr_tsbh')
        rejected_with_nut_med, _, _, _ = multitest.multipletests(with_nut_med_asymptotic_p[0], alpha=alpha,
                                                                 method='fdr_tsbh')
        rejected_with_nut_scd, _, _, _ = multitest.multipletests(with_nut_scd_asymptotic_p[0], alpha=alpha,
                                                                 method='fdr_tsbh')

        rejected_with_nut_before_med, _, _, _ = multitest.multipletests(with_nut_before_med_asymptotic_p[0],
                                                                        alpha=alpha, method='fdr_tsbh')
        rejected_with_nut_before_scd, _, _, _ = multitest.multipletests(with_nut_before_scd_asymptotic_p[0],
                                                                        alpha=alpha, method='fdr_tsbh')
        rejected_with_nut_refrence, _, _, _ = multitest.multipletests(with_nut_reference_asymptotic_p[0], alpha=alpha,
                                                                      method='fdr_tsbh')

        med_corr = med_corr[rejected_med]
        scd_corr = scd_corr[rejected_scd]
        ref_corr = backbone_corr[np.logical_or(rejected_med, rejected_scd)]
        backbone_corr = backbone_corr[rejected_backbone]
        with_nut_med_corr = with_nut_med_corr[rejected_with_nut_med]
        with_nut_scd_corr = with_nut_scd_corr[rejected_with_nut_scd]

        with_nut_before_med_corr = with_nut_before_med_corr[rejected_with_nut_before_med]
        with_nut_before_scd_corr = with_nut_before_scd_corr[rejected_with_nut_before_scd]
        with_nut_reference_corr = with_nut_reference_corr[rejected_with_nut_refrence]

        med_corr = self.filter_correlation_by_threshold(med_corr, threshold)
        scd_corr = self.filter_correlation_by_threshold(scd_corr, threshold)
        backbone_corr = self.filter_correlation_by_threshold(backbone_corr, threshold)
        ref_corr = self.filter_correlation_by_threshold(ref_corr, threshold)
        with_nut_med_corr = self.filter_correlation_by_threshold(with_nut_med_corr, threshold)
        with_nut_scd_corr = self.filter_correlation_by_threshold(with_nut_scd_corr, threshold)

        with_nut_before_med_corr = self.filter_correlation_by_threshold(with_nut_before_med_corr, threshold)
        with_nut_before_scd_corr = self.filter_correlation_by_threshold(with_nut_before_scd_corr, threshold)
        with_nut_reference_corr = self.filter_correlation_by_threshold(with_nut_reference_corr, threshold)

        return ref_corr, med_corr, scd_corr, with_nut_med_corr, with_nut_scd_corr, with_nut_before_med_corr, with_nut_before_scd_corr, with_nut_reference_corr

    def filter_correlation_by_threshold(self ,correlation_df, threshold):
        inter_type_thr, cross_type_thr, key_type_thr = threshold

        correlation_df['is_same_type'] = correlation_df.index.map(
            lambda x: self.data_api.feature_to_type_dict[x[0]] == self.data_api.feature_to_type_dict[x[1]])
        correlation_df['is_key'] = correlation_df.index.map(
            lambda x: (self.data_api.feature_to_type_dict[x[0]] == 'key') |
                        (self.data_api.feature_to_type_dict[x[1]] == 'key'))

        inter_filtered = correlation_df['is_same_type'] & (correlation_df['is_key'] != True) & (
                correlation_df[correlation_df.columns[0]] > inter_type_thr)
        cross_filtered = (correlation_df['is_same_type'] != True) & (correlation_df['is_key'] != True) & (
                correlation_df[correlation_df.columns[0]] > cross_type_thr)
        cross_filtered = correlation_df['is_key'] & (correlation_df[correlation_df.columns[0]] > key_type_thr)

        new_correlation_df = correlation_df[inter_filtered | cross_filtered | cross_filtered]

        return new_correlation_df.drop(columns=['is_same_type'])

    def copy_df(self, df):
        if (type(df) is pd.DataFrame):
            return pd.DataFrame(columns=df.columns.copy(deep=True), index=df.index.copy(deep=True),
                                data=df.values.copy())
        else:
            return pd.Series(data=df.values.copy(), index=df.index.copy(deep=True))

    def remove_duplicate_correlation(self, corr_df):
        temp = corr_df
        temp['temp_index'] = temp.apply(
            lambda row: (min(str(row['level_0']), str(row['level_1'])), max(str(row['level_0']), str(row['level_1']))),
            axis=1)
        temp = temp.drop_duplicates(subset='temp_index').drop(columns=['temp_index'])
        return temp

    def calc_spearman(self, df, remove_duplicates=True):
        X, p_values = scipy.stats.spearmanr(df)

        if (type(df.columns) is pd.core.indexes.multi.MultiIndex):
            p_df = pd.DataFrame(index=df.columns.droplevel(0), columns=df.columns.droplevel(0), data=p_values)
            x_df = pd.DataFrame(index=df.columns.droplevel(0), columns=df.columns.droplevel(0), data=X)
        else:
            p_df = pd.DataFrame(index=df.columns, columns=df.columns, data=p_values)
            x_df = pd.DataFrame(index=df.columns, columns=df.columns, data=X)

        flatten_x = x_df.unstack()
        flatten_p = p_df.unstack()

        if (remove_duplicates):
            self_corr_mask = np.array(list(map(lambda tup: tup[0] == tup[1], zip(flatten_x.index.get_level_values(0),
                                                                                 flatten_x.index.get_level_values(1)))))
            flatten_x = flatten_x.drop(labels=flatten_x.index[self_corr_mask])
            flatten_p = flatten_p.drop(labels=flatten_p.index[self_corr_mask])

        flatten_x = self.copy_df(flatten_x).reset_index()
        flatten_p = self.copy_df(flatten_p).reset_index()

        if (remove_duplicates):
            flatten_x = self.remove_duplicate_correlation(flatten_x)
            flatten_p = self.remove_duplicate_correlation(flatten_p)

        #     flatten_x

        return flatten_x.set_index(['level_0', 'level_1']), flatten_p.set_index(['level_0', 'level_1'])

    def fisher_transformation_over_flatt_pd(self,flatt_df, n):
        return flatt_df.applymap(lambda r: self.fisher_transformation(r, n))

    # @njit
    def fisher_transformation(self,rho, n):
        if (abs(rho - 1) < 1e-5):
            rho = np.sign(rho) * (1 - 1e-5)
        return np.arctanh(rho)
