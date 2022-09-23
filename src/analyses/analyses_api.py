from statsmodels.stats import multitest
import pandas as pd
from scipy import stats
import numpy as np
from global_utils import GlobalUtils

class AnalysesApi() :
    def __init__(self):
        pass

    @staticmethod
    def extract_disruptions_per__lihi_module(disruption_df, module):
        disruption_df = GlobalUtils.copy_df(disruption_df)
        relvent_idx = disruption_df.index.map(lambda x: ((str(x[0]) in module) or (str(x[1]) in module)))
        return disruption_df.loc[relvent_idx]

    @staticmethod
    def extract_indexes_from_list(df, idx_list):
        idx_as_index = pd.Index(idx_list)
        level0_map = df.index.isin(idx_as_index, level=0)
        level1_map = df.index.isin(idx_as_index, level=1)
        _edges = np.logical_or(level0_map, level1_map)

        idx_list_int = [int(idx) for idx in idx_list if idx.isnumeric()]
        idx_as_index_int = pd.Index(idx_list_int)
        level0_map_int = df.index.isin(idx_as_index_int, level=0)
        level1_map_int = df.index.isin(idx_as_index_int, level=1)
        _edges_int = np.logical_or(level0_map_int, level1_map_int)

        return df.loc[np.logical_or(_edges, _edges_int)]

    @staticmethod
    def keep_only_lihi_picked_features(enriched_features_df,picked_features_by_lihi):
        _new = {}
        for fe, enr_df in enriched_features_df.items():
            feature_enrichment_in_comm = AnalysesApi().extract_indexes_from_list(enr_df, picked_features_by_lihi)
            _new[fe] = feature_enrichment_in_comm
        return _new

    @staticmethod
    def calculate_correlation_to_disruption(disruption_df,data_df,pn_in_corr,pn_in_corr_limit=15) :
        _pn_in_corr = pn_in_corr.copy(deep=True)
        _pn_in_corr.index = _pn_in_corr.index.map(lambda x: f"{x[0]}-{x[1]}")

        correlations_to_disruption_matrix = []
        p_value_of_correlations_to_disruption_matrix = []

        for dis_corr in disruption_df.index:
            dis_corr_data = disruption_df.loc[dis_corr][_pn_in_corr.loc[dis_corr]]

            if _pn_in_corr.loc[dis_corr].sum() < pn_in_corr_limit:
                continue

            for imp in data_df.columns:
                imp_data = data_df[imp][_pn_in_corr.loc[dis_corr]]
                rho, p_val = stats.spearmanr(pd.merge(dis_corr_data, imp_data, left_index=True, right_index=True))
                if np.isnan(p_val):
                    correlations_to_disruption_matrix.append([dis_corr, imp, 0])
                    p_value_of_correlations_to_disruption_matrix.append([dis_corr, imp, 0.99])
                    continue

                correlations_to_disruption_matrix.append([dis_corr, imp, rho])
                p_value_of_correlations_to_disruption_matrix.append([dis_corr, imp, p_val])

        if len(correlations_to_disruption_matrix) == 0:
            return None, None

        correlations_to_disruption = pd.DataFrame(data=correlations_to_disruption_matrix,
                                                  columns=["disruption", "ser", "corr"])
        correlations_to_disruption = correlations_to_disruption.pivot_table(index="disruption", columns="ser",
                                                                            values="corr")

        p_value_of_correlations_to_disruption = pd.DataFrame(data=p_value_of_correlations_to_disruption_matrix,
                                                             columns=["disruption", "ser", "p_val"])
        p_value_of_correlations_to_disruption = p_value_of_correlations_to_disruption.pivot_table(index="disruption",
                                                                                                  columns="ser",
                                                                                                  values="p_val")

        return correlations_to_disruption, p_value_of_correlations_to_disruption

    @staticmethod
    def calculate_corrected_p_value_of_correlations_to_disruption(disruption_of_sig_modules, seg_imp_data, pn_in_corr):
        # consider only pn thet take place in the correaltion - disruption is high for pn not from the correaltion with "no good reason"
        correlations_to_disruption, p_value_of_correlations_to_disruption = AnalysesApi().calculate_correlation_to_disruption(disruption_of_sig_modules, seg_imp_data, pn_in_corr)

        if correlations_to_disruption is None:
            return None

        # FDR for the correlation to baseline abundance
        all_corrected_pval = []
        for disruption in p_value_of_correlations_to_disruption.index:
            pvals_to_dis = p_value_of_correlations_to_disruption.loc[disruption]
            rejected, pvals_corrected, _, _ = multitest.multipletests(pvals_to_dis, 1.0, method="fdr_bh")
            series = pd.DataFrame(index=pvals_to_dis.index, columns=[disruption], data=pvals_corrected, copy=True)
            all_corrected_pval.append(series)
        corrcted_p_value_of_correlations_to_disruption = pd.concat(all_corrected_pval, axis=1).T

        return corrcted_p_value_of_correlations_to_disruption

    @staticmethod
    def calculate_joined_correlation_to_disruption(disruption_of_sig_modules, seg_imp_data, pn_in_corr,pn_in_corr_limit = 15):
        correlations_to_disruption,p_value_of_correlations_to_disruption = AnalysesApi().calculate_correlation_to_disruption(disruption_of_sig_modules,seg_imp_data,pn_in_corr,pn_in_corr_limit=pn_in_corr_limit)
        corrcted_p_value_of_correlations_to_disruption = AnalysesApi().calculate_corrected_p_value_of_correlations_to_disruption(disruption_of_sig_modules, seg_imp_data, pn_in_corr)

        if corrcted_p_value_of_correlations_to_disruption is None :
            return None

        # reshape tables - FDR for the correlation to baseline abundance
        melted_corr = correlations_to_disruption.reset_index().melt(id_vars=("disruption")).set_index(
            ["disruption", "ser"])
        melted_corr = melted_corr.fillna(0)
        melted_pval = p_value_of_correlations_to_disruption.reset_index().melt(id_vars=("disruption")).set_index(
            ["disruption", "ser"])
        melted_corrcted_pval = corrcted_p_value_of_correlations_to_disruption.reset_index().rename(
            columns={"index": "disruption"}).melt(id_vars=("disruption")).set_index(["disruption", "ser"])

        melted_corr_and_disruption = pd.merge(melted_corr, melted_pval, left_index=True, right_index=True,
                                              suffixes=("_corr", "_pval"))
        melted_corr_and_disruption = pd.merge(melted_corr_and_disruption, melted_corrcted_pval, left_index=True,
                                              right_index=True).rename(columns={"value": "corrected_pval"})

        return melted_corr_and_disruption

