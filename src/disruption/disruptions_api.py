from correlations_tools import CorrelationsTools
import pandas as pd
import numpy as np
from global_utils import GlobalUtils
from functools import reduce
from sklearn.metrics import matthews_corrcoef


class DisruptionAPI():
    BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET = ['reference', 'med', 'scd', 'med with nut', 'scd with nut',
                                               'before med with nuts', 'before scd with nuts', 'refrence with nuts']
    DIET_WIDG_FOR_DISRUPTION_MODULE_NET = ['reference', 'med', 'scd', 'med with nut', 'scd with nut',
                                           'before med with nuts', 'before scd with nuts', 'refrence with nuts']

    def __init__(self):
        self.corr_tools = CorrelationsTools()

    # region private

    @GlobalUtils.cache_me
    def leave_one_out_reference_z_disruption_network(self, df, correlation_to_take, all_in_z_corr=None):
        unique_index_df = self.corr_tools.copy_df(df)
        unique_index_df.index = pd.Index([str(pn) + '_' + str(it) for it, pn in enumerate(unique_index_df.index)])

        if (all_in_z_corr is None):
            all_in_z_corr = self.corr_tools.fisher_transformation_over_flatt_pd(
                self.corr_tools.calc_spearman(unique_index_df, False)[0].loc[correlation_to_take],
                unique_index_df.shape[0])

        res_dfs = []
        for pn in unique_index_df.index:
            res_df = pd.DataFrame(index=all_in_z_corr.index, columns=[pn.split('_')[0]])
            one_out_df = unique_index_df.loc[unique_index_df.index.difference([pn])]
            one_out_z_coor = self.corr_tools.fisher_transformation_over_flatt_pd(
                self.corr_tools.calc_spearman(one_out_df, False)[0].loc[correlation_to_take], one_out_df.shape[0])
            # -1* is for positive disruption
            z_disruption = -1 * (abs(all_in_z_corr) - abs(one_out_z_coor))
            res_df[pn.split('_')[0]] = z_disruption
            res_dfs.append(res_df)
        res = pd.concat(res_dfs, axis=1)
        res[res == -1 * np.inf] = -1
        return res

    def normalize_pn(self,series, reference):
        diff_from_ref = reference.apply(lambda col: col - series).dropna(axis=0)
        return ((diff_from_ref > 0).sum(axis=1) / ((reference > 0).sum(axis=1))).fillna(0)

    def calc_explanation_score(self,disruption, feature):
        m = disruption.to_frame().merge(feature.to_frame(), left_index=True, right_index=True)
        # score is pvalue.
        _f_vec = m[m.columns[0]]
        _d_vec = m[m.columns[1]]

        original_distance = matthews_corrcoef(_f_vec, _d_vec)

        permut_results = []
        for i in range(99):
            _f_permu_vec = _f_vec.sample(frac=1).reset_index(drop=True)
            _d_permu_vec = _d_vec.sample(frac=1).reset_index(drop=True)

            _permut = matthews_corrcoef(_f_permu_vec, _d_permu_vec)
            _is_smaller = original_distance > _permut
            permut_results.append(_is_smaller)

        p_val_score = 1 - sum(permut_results) / len(permut_results)

        return p_val_score

    # endregion

    #region public

    def calculate_ref_z_disruption(self,backbone_df,relevent_correlations,reference_diet_idx,treatment_diet_idx,alpah_fdr,corr_th):
        ref_z_disruption = self.leave_one_out_reference_z_disruption_network(backbone_df, relevent_correlations.index,
                                                                        for_cache={'fdr_alpha': alpah_fdr,'corr_th': corr_th,'refernce_diet':
                                                                                       DisruptionAPI.BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET[reference_diet_idx],
                                                                                   'for_diet':DisruptionAPI.DIET_WIDG_FOR_DISRUPTION_MODULE_NET[treatment_diet_idx]})
        return ref_z_disruption

    @GlobalUtils.cache_me
    def calculate_z_disruption_networks(self,ref_df , df,correlation_to_take):
        pn_to_disruption_dic = {}
        base_corr = self.corr_tools.fisher_transformation_over_flatt_pd(
            self.corr_tools.calc_spearman(ref_df, remove_duplicates=False)[0].loc[correlation_to_take], ref_df.shape[0])

        for pn in df.index:
            if type(df.loc[pn]) is pd.DataFrame:
                ref_with_pn = ref_df.append(df.loc[pn].mean().to_frame().rename(columns={0: pn}).T)
            else:
                ref_with_pn = ref_df.append(df.loc[pn])
            corr_with_pn = self.corr_tools.fisher_transformation_over_flatt_pd(self.corr_tools.calc_spearman(ref_with_pn)[0].loc[correlation_to_take],
                                                               ref_with_pn.shape[0])
            # -1* is for positive disruption
            disruption = -1 * (abs(corr_with_pn) - abs(base_corr))
            pn_to_disruption_dic[pn] = disruption.rename(columns={0: pn})

        return reduce(lambda left_df, right_df: left_df.merge(right_df, left_index=True, right_index=True),
                      pn_to_disruption_dic.values())

    def get_patient_in_correlation(self,ref_disruption):
        # negtive disruption : someone that help to the correlation
        # positive : someone thet disrupt
        is_sample_in_corr = ref_disruption < 0
        pn_in_corr = pd.DataFrame(index=is_sample_in_corr.index, columns=is_sample_in_corr.columns.unique())
        for pn in pn_in_corr.columns:
            if is_sample_in_corr.shape[1] == 2:
                pn_in_corr[pn] = is_sample_in_corr[pn].all(axis=1)
            else:
                pn_in_corr[pn] = is_sample_in_corr[pn]
        return pn_in_corr


        #endregion

    @GlobalUtils.cache_me
    def calc_disruption_score_and_remove_redundant(self,df, reference, _pn_in_corr):
        df = df[_pn_in_corr].fillna(0)
        res = df.apply(lambda pn: self.normalize_pn(pn, reference))
        res[res > 0] = 1 - res[res > 0]
        return res

    def get_enrichment_scores(self,disruption_df,diet_df,explanatory_features,feature_trh):
        _disruption_df = disruption_df
        _disruption_df[_disruption_df>0]=1
        explanation_per_feature = {}

        for exp_feature in explanatory_features:
            exp_feature_values = diet_df[exp_feature]
            exp_feature_bin = pd.Series(index=exp_feature_values.index,data=0)

            #label all pn with the feature trh condition
            if (feature_trh >= 0):
                exp_feature_bin[exp_feature_values>=feature_trh] = 1
            else :
                exp_feature_bin[exp_feature_values<=feature_trh] = 1

            explanation_per_feature[exp_feature] =  _disruption_df.apply(lambda row:self.calc_explanation_score(row,exp_feature_bin),axis=1).fillna(0)
        return explanation_per_feature

    def get_enriched_features(self,enrichment_score_df, explanatory_features, trh):
        enriched_features_df = {}
        for exp_feature in explanatory_features:
            score = enrichment_score_df[exp_feature]
            enriched_features_df[exp_feature] = score[score < trh]
        return enriched_features_df