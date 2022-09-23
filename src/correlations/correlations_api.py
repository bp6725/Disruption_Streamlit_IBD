from data_api import DataAPI
from correlations_tools import CorrelationsTools
from networkx.algorithms import community as nx_community
from global_utils import GlobalUtils
import networkx as nx


class CorrelationsAPI():

    CORRELATION_NAME_TO_IDX = {"ref_corr":0, "med_corr":1, "scd_corr":2, "with_nut_med_corr":3, "with_nut_scd_corr":4, "with_nut_before_med_corr":5, "with_nut_before_scd_corr":6, "with_nut_reference_corr":7}


    def __init__(self,alpah_fdr,corr_th):
        self.data_api = DataAPI()
        self.correlation_class = CorrelationsTools()

        ed_sig_df, scd_sig_df, union_ref_df, med_with_nutriens, scd_with_nutriens, \
        before_med_sig_df_with_nutriens, before_scd_sig_df_with_nutriens, all_refrence_with_nuts = self.data_api.return_unagg_biom_data()
        self.corr_dfs = self.correlation_class.build_correlation_dfs(ed_sig_df, scd_sig_df, union_ref_df, med_with_nutriens, scd_with_nutriens, \
                                                   before_med_sig_df_with_nutriens, before_scd_sig_df_with_nutriens, all_refrence_with_nuts,
                                                   alpah_fdr, corr_th)

    def return_correlation_by_index(self,idx):
        return self.corr_dfs[idx]

    def return_correlation_by_name(self,name):
        return self.corr_dfs[CorrelationsAPI.CORRELATION_NAME_TO_IDX[name]]

    #region network analysis

    def extract_communites_from_g(self,G, modularity_base=True):
        if (modularity_base):
            communities_generator = nx_community.modularity_max.greedy_modularity_communities(G)
            return [[mem for mem in comm] for comm in communities_generator]
        else:
            communities_generator = nx_community.girvan_newman(G)
        top_level_communities = next(communities_generator)
        next_level_communities = next(communities_generator)
        return sorted(map(sorted, next_level_communities))

    def generete_nx_network(self,series):
        series = GlobalUtils.copy_df(series)
        series = series.reset_index()

        series[['level_0', 'level_1']] = series[['level_0', 'level_1']].astype(str)
        return nx.from_pandas_edgelist(series, source='level_0', target='level_1')

    def extract_spaciel_features_base_communites(self,G, features=['PDAI', 'Calprotectin', 'CRP', 'Shannon', 'MDI']):
        _communites = []
        for node in features:
            if not G.has_node(node):
                continue
            _community = [node]
            neighbors = [n for n in G.neighbors(node)]
            neighbors_of_neighbors = [[sn for sn in G.neighbors(n)] for n in neighbors]
            neighbor_list = neighbors + neighbors_of_neighbors[0]
            for n in neighbor_list:
                _community.append(n)
            _communites.append(_community)
        return _communites

    #endregion