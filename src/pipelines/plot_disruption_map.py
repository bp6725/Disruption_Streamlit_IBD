from correlations_api import CorrelationsAPI
from data_api import DataAPI
import pandas as pd
from analyses_api import AnalysesApi
from disruption_analysis import DisruptionAnalysis
from functools import reduce
import holoviews as hv
from global_utils import GlobalUtils
import networkx as nx
from bokeh.plotting import show

hv.extension('bokeh')

class PlotDisruptionMap():
    def __init__(self,stage_for_disruption_prediction = "Before_MED",
                 explanatory_features = ['CRP','Cholesterol','HDL','LDL','Calprotectin','Shannon','Zonulin','PDAI','IL6'],
                 disruption_analysis = None,corr_api = None,data_api = None,analyses_api = None):

        self.stage_for_disruption_prediction = stage_for_disruption_prediction
        self.explanatory_features = explanatory_features

        if disruption_analysis is None:
            self.disruption_analysis = DisruptionAnalysis()
        else :
            self.disruption_analysis = disruption_analysis

        if corr_api is None:
            self.corr_api = CorrelationsAPI(self.disruption_analysis.alpah_fdr, self.disruption_analysis.corr_th)
        else :
            self.corr_api = corr_api

        if data_api is None:
            self.data_api = DataAPI()
        else :
            self.data_api = data_api

        if analyses_api is None:
            self.analyses_api = AnalysesApi()
        else :
            self.analyses_api = analyses_api


    def plot_correlation_network(self,network, communites=None,layout_graph = "circular"):
        if ((type(network) is pd.Series) or (type(network) is pd.DataFrame)):
            G = self.corr_api.generete_nx_network(network)
        else:
            G = network

        if (communites is not None):
            node_to_color = self.node_to_color_from_communites(communites)
            nx.set_node_attributes(G, node_to_color, 'color_att')
        nx.set_node_attributes(G, self.data_api.feature_name_dict, 'name')

        padding = dict(x=(-1.1, 1.1), y=(-1.1, 1.1))

        if layout_graph == "circular" :
            return hv.Graph.from_networkx(G, nx.layout.circular_layout) \
            .redim.range(**padding) \
            .options(color_index='color_att', cmap='Category20', fontsize={'label': 1})
        else :
            return hv.Graph.from_networkx(G, nx.layout.spring_layout) \
                .redim.range(**padding) \
                .options(color_index='color_att', cmap='Category20', fontsize={'label': 1})

    def plot_enriched_disruption_network(self,reference_diet_idx,treatment_diet_idx,feature_trh=0.1,only_spaciel_features_based_comunites = False):
        relevent_correlations = self.corr_api.return_correlation_by_index(reference_diet_idx)

        kwargs = {"reference_diet_idx":reference_diet_idx,"treatment_diet_idx":treatment_diet_idx,"explanatory_features":self.explanatory_features}
        self.disruption_analysis.update_class_properties(kwargs)

        enriched_features_df = self.disruption_analysis.get_enriched_features(feature_trh,False)

        if only_spaciel_features_based_comunites:
            _communites = self.corr_api.extract_spaciel_features_base_communites(self.corr_api.generete_nx_network(relevent_correlations),self.explanatory_features)
        else:
            _communites = relevent_correlations.index.get_level_values(0).unique().union(
                relevent_correlations.index.get_level_values(1).unique()).unique()
            _communites = [[str(i) for i in _communites.tolist()]]

        _communites = [list(set(reduce(lambda x, y: x + y, _communites)))]

        # disruption_matrix = self.disruption_analysis.get_disruption_matrix()
        # disruption_matrix_with_comm = self.disruption_analysis.extract_disruptions_per_module(disruption_matrix, _communites)

        modules = _communites
        feature_enrichments = enriched_features_df

        results = {}
        for enriched_feature, _feature_enrichment in feature_enrichments.items():
            feature_enrichment = GlobalUtils.copy_df(_feature_enrichment)
            comm_plots = []
            for comm in modules:
                feature_enrichment_in_comm = AnalysesApi.extract_indexes_from_list(feature_enrichment, comm)
                if (feature_enrichment_in_comm.shape[0] == 0):
                    continue
                G = self.plot_correlation_network(feature_enrichment_in_comm)
                labels = hv.Labels(G.nodes, ['x', 'y'], 'name')

                padding = dict(x=(-1.1, 1.1), y=(-1.1, 1.1))
                comm_plots.append(G * labels)
            if len(comm_plots) > 0:
                comm_final_plot = reduce(lambda x, y: x + y, comm_plots)

                results[enriched_feature] =  comm_final_plot

        return results
