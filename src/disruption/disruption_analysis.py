from correlations_api import CorrelationsAPI
from data_api import DataAPI
from correlations_tools import CorrelationsTools
import pandas as pd
from global_utils import GlobalUtils
from disruptions_api import DisruptionAPI
from analyses_api import AnalysesApi
import config
from itertools import permutations

class DisruptionAnalysis() :
    def __init__(self, reference_diet_idx = 0, treatment_diet_idx = 1, alpah_fdr = 0.15, corr_th = (0.5, 0.4, 0.3), edge_cut = 3, enrichment_trh = 0.15,
                  disruption_cutoff = 0.05, explanatory_features = ['CRP']):
        self.reference_diet_idx = reference_diet_idx
        self.treatment_diet_idx = treatment_diet_idx
        self.alpah_fdr  = alpah_fdr
        self.corr_th = corr_th
        self.edge_cut = edge_cut
        self.enrichment_trh = enrichment_trh
        self.disruption_cutoff = disruption_cutoff
        self.explanatory_features = explanatory_features

        self.corr_tools = CorrelationsTools()

        self.disruption_api = DisruptionAPI()
        self.corr_api = CorrelationsAPI(self.alpah_fdr,self.corr_th)
        self.data_api = DataAPI()
        self.analyses_api = AnalysesApi()

    #region private

    def cast_to_int(self,st):
        if (st.isdigit()):
            return int(st)
        return st

    def update_class_properties(self, kwargs):
        self.reference_diet_idx = self.update_property(self.reference_diet_idx, kwargs, "reference_diet_idx")
        self.treatment_diet_idx = self.update_property(self.treatment_diet_idx, kwargs, "treatment_diet_idx")
        self.alpah_fdr = self.update_property(self.alpah_fdr, kwargs, "alpah_fdr")
        self.corr_th = self.update_property(self.corr_th, kwargs, "corr_th")
        self.edge_cut = self.update_property(self.edge_cut, kwargs, "edge_cut")
        self.enrichment_trh = self.update_property(self.enrichment_trh, kwargs, "enrichment_trh")
        self.disruption_cutoff = self.update_property(self.disruption_cutoff, kwargs, "disruption_cutoff")
        self.explanatory_features = self.update_property(self.explanatory_features, kwargs, "explanatory_features")

    @staticmethod
    def update_property(org_val,kwargs,arg_name):
        if arg_name in kwargs.keys() :
            return kwargs[arg_name]
        return org_val

    #endregion

    #region public

    def calculate_ref_z_disruption(self):

        # relevent_correlations here is the correlations exists in the backbone .. there are the "exist" correlations we want the the diet to disurpt
        relevent_correlations = self.corr_api.return_correlation_by_index(self.reference_diet_idx)
        backbone_df = self.data_api.return_unagg_by_index(self.reference_diet_idx)

        return self.disruption_api.calculate_ref_z_disruption(backbone_df,relevent_correlations,
                                                              self.reference_diet_idx,self.treatment_diet_idx,
                                                              self.alpah_fdr,self.corr_th)

    def calculate_z_disruption_networks(self):
        relevent_correlations = self.corr_api.return_correlation_by_index(self.reference_diet_idx)
        backbone_df = self.data_api.return_unagg_by_index(self.reference_diet_idx)
        diet_df = self.data_api.return_unagg_by_index(self.treatment_diet_idx)

        for_cache_args = {'fdr_alpha': self.alpah_fdr, 'corr_th': self.corr_th,
                          'refernce_diet': GlobalUtils.BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET[self.reference_diet_idx],
                          'for_diet': GlobalUtils.DIET_WIDG_FOR_DISRUPTION_MODULE_NET[self.treatment_diet_idx]}
        return self.disruption_api.calculate_z_disruption_networks(backbone_df,diet_df,relevent_correlations.index,for_cache=for_cache_args)

    def get_patient_in_correlation(self):
        return self.disruption_api.get_patient_in_correlation(self.calculate_ref_z_disruption())

    def get_disruption_matrix(self):
        '''
        you probably need get_disruption_score ..
        :return:
        '''
        for_cache_args = {'fdr_alpha': self.alpah_fdr, 'corr_th': self.corr_th,
                          'refernce_diet': GlobalUtils.BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET[self.reference_diet_idx],
                          'for_diet': GlobalUtils.DIET_WIDG_FOR_DISRUPTION_MODULE_NET[self.treatment_diet_idx]}

        return self.disruption_api.calc_disruption_score_and_remove_redundant(
            self.calculate_z_disruption_networks(), self.calculate_ref_z_disruption(),
            self.get_patient_in_correlation(), for_cache=for_cache_args)

    def get_disruption_score(self):
        for_cache_args = {'fdr_alpha': self.alpah_fdr, 'corr_th': self.corr_th,
                          'refernce_diet': GlobalUtils.BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET[self.reference_diet_idx],
                          'for_diet': GlobalUtils.DIET_WIDG_FOR_DISRUPTION_MODULE_NET[self.treatment_diet_idx]}

        disruption_matrix = self.disruption_api.calc_disruption_score_and_remove_redundant(self.calculate_z_disruption_networks(), self.calculate_ref_z_disruption(),
                                                                                           self.get_patient_in_correlation(),for_cache=for_cache_args)

        clean_disruption_df = disruption_matrix.copy(deep=True)
        clean_disruption_df[clean_disruption_df < self.disruption_cutoff] = 0
        clean_disruption_df = clean_disruption_df[(clean_disruption_df > 0).sum(axis=1) > self.edge_cut]

        return clean_disruption_df

    def get_diet_df_with_delta(self):
        diet_df = self.data_api.return_unagg_by_index(self.treatment_diet_idx)

        if self.treatment_diet_idx in [1, 2, 3, 4]:
            explanatory_features_origion = {1: 'MED', 2: 'SCD', 3: 'MED', 4: 'SCD'}[self.treatment_diet_idx]
            diet_df_with_delta = self.data_api.get_diet_df_with_delta(diet_df, self.explanatory_features,
                                                                      explanatory_features_origion)
        else:
            explanatory_features_origion = {1: 'MED', 2: 'SCD', 3: 'MED', 4: 'SCD'}[self.reference_diet_idx]
            diet_df_with_delta = self.data_api.get_diet_df_with_delta(diet_df, self.explanatory_features,
                                                                      explanatory_features_origion)

        return diet_df_with_delta

    def get_enrichment_scores(self,feature_trh):
        clean_disruption_df = self.get_disruption_score()
        diet_df_with_delta = self.get_diet_df_with_delta()

        return self.disruption_api.get_enrichment_scores(clean_disruption_df,diet_df_with_delta,[f + '_delta' for f in self.explanatory_features],feature_trh)

    def get_enriched_features(self,feature_trh,only_pre_picked_features = True):
        enriched_features_df = self.disruption_api.get_enriched_features(self.get_enrichment_scores(feature_trh),
                                                                         [f + '_delta' for f in
                                                                          self.explanatory_features],
                                                                         self.enrichment_trh)
        if only_pre_picked_features :
            picked_features_by_lihi = [str(mat) for mat in
                                       pd.read_excel(config.disruption_analysis["features_by_lihi"],header=None)[1].tolist()]
            enriched_features_df = self.analyses_api.keep_only_lihi_picked_features(enriched_features_df, picked_features_by_lihi)

        return enriched_features_df

    def get_disruption_of_significant_modules(self,feature_trh,only_pre_picked_features = True):
        # only disruption of sig features - same indexes as enriched_features_df. its only for reshapeing the df
        disruption_of_sig_modules = self.calculate_z_disruption_networks().loc[[v for v in self.get_enriched_features(feature_trh,only_pre_picked_features).values()][0].index]
        disruption_of_sig_modules.index = disruption_of_sig_modules.index.map(lambda x: f"{x[0]}-{x[1]}")

        return disruption_of_sig_modules

    #endregion


    #region Tools

    def extract_disruptions_per_module(self,disruption_df, modules):
        disruption_df = GlobalUtils.copy_df(disruption_df)
        modules_dfs = []
        for module_idx, module in enumerate(modules):
            l = list(map(lambda tup: (self.cast_to_int(tup[0]), self.cast_to_int(tup[1])), permutations(module, 2)))
            i = pd.MultiIndex.from_tuples(l)
            module_df = disruption_df.loc[i]
            module_df = module_df.dropna()
            module_df['module_idx'] = module_idx
            modules_dfs += [module_df]

        return pd.concat(modules_dfs)

    #endregion