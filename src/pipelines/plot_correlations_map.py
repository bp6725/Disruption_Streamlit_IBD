from correlations_api import CorrelationsAPI
from data_api import DataAPI
import pandas as pd
from analyses_api import AnalysesApi
from global_utils import GlobalUtils

import networkx as nx
import numpy as np
import holoviews as hv
import matplotlib
from functools import reduce
from bokeh.plotting import show
from bokeh.palettes import Category20


# hv.extension('bokeh')

class PlotCorrelationsMap():
    def __init__(self, alpah_fdr = 0.15, corr_th = (0.5, 0.4, 0.3),corr_api = None,data_api = None,analyses_api = None):
        self.alpah_fdr = alpah_fdr
        self.corr_th = corr_th

        if corr_api is None :
            self.corr_api = CorrelationsAPI(alpah_fdr, corr_th)

        if data_api is None :
            self.data_api = DataAPI()

        if analyses_api is None :
            self.analyses_api = AnalysesApi()

    #region plot tools

    def plot_ordered_colored_network(self,df, title,communites_order,communites_color_hex):
        correlation_df_for_plot = df.copy(deep=True)

        correlation_df_for_plot = correlation_df_for_plot.rename(index=self.rename_index).rename(index=self.reshape_str_for_circos)
        correlation_df_for_plot = correlation_df_for_plot.to_frame()
        correlation_df_for_plot["order"] = correlation_df_for_plot.index.get_level_values(0).map(communites_order)
        correlation_df_for_plot = correlation_df_for_plot.sort_values("order")
        correlation_df_for_plot = correlation_df_for_plot[correlation_df_for_plot.columns[0]]

        series = GlobalUtils.copy_df(correlation_df_for_plot)
        series = series.reset_index()

        series[['level_0', 'level_1']] = series[['level_0', 'level_1']].astype(str)
        G = nx.from_pandas_edgelist(series, source='level_0', target='level_1')

        nx.set_node_attributes(G, communites_color_hex, 'color_att')
        nx.set_edge_attributes(G, correlation_df_for_plot.to_dict(), "Weight")

        padding = dict(x=(-1.1, 1.1), y=(-1.1, 1.1))

        return hv.Graph.from_networkx(G, nx.layout.circular_layout) \
            .redim.range(**padding) \
            .options(color_index='color_att', cmap='Category20', fontsize={'label': 1}, edge_alpha='Weight',
                     title=title)

    def node_to_color_from_communites(self,communites):
        #     colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','#ff81c0','#ff81c0','#ff81c0']

        node_to_color = {}
        for idx, com in enumerate(communites):
            if idx > 19:
                idx = np.random.choice(range(19))
            for node in com:
                node_to_color[node] = Category20[20][idx]
        return node_to_color

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

    def plot_correlation_graph_with_communites(self,corr_df):
        G = self.corr_api.generete_nx_network(corr_df)
        _communites = self.corr_api.extract_communites_from_g(G)
        return self.plot_correlation_network(corr_df, _communites)

    #endregion

    #region unique tools

    def rename_index(self,ind):
        if type(ind) is int:
            return f"{self.data_api.imp_id_to_path[ind]}_{ind}"

        if str.isdigit(ind):
            if (int(ind) not in self.data_api.imp_id_to_path.keys()) and (str(ind) not in self.data_api.imp_id_to_path.keys()):
                return f"metabolite_{ind}"
            else:
                return f"{self.data_api.imp_id_to_path[int(ind)]}_{ind}"

        if type(ind) is str:
            if "-" in ind:
                _f_dis_ind, _s_dis_ind = ind.split("-")
                if _f_dis_ind.isdigit():
                    f_dis_ind = f"{_f_dis_ind} {self.data_api.imp_id_to_path[int(_f_dis_ind)]}"
                else:
                    f_dis_ind = _f_dis_ind

                if _s_dis_ind.isdigit():
                    s_dis_ind = f"{_s_dis_ind} {self.data_api.imp_id_to_path[int(_s_dis_ind)]}"
                else:
                    s_dis_ind = _s_dis_ind

                return f"({f_dis_ind} - {s_dis_ind})"
            else:
                return ind

        dis_ind = ind[0]
        ser_ind = ind[1]

        _f_dis_ind, _s_dis_ind = dis_ind.split("-")

        if _f_dis_ind.isdigit():
            f_dis_ind = f"{_f_dis_ind} {self.data_api.imp_id_to_path[int(_f_dis_ind)]}"
        else:
            f_dis_ind = _f_dis_ind

        if _s_dis_ind.isdigit():
            s_dis_ind = f"{_s_dis_ind} {self.data_api.imp_id_to_path[int(_s_dis_ind)]}"
        else:
            s_dis_ind = _s_dis_ind

        new_dis_ind = f"({f_dis_ind} - {s_dis_ind})"
        new_ser_ind = f"{ser_ind} {self.data_api.imp_id_to_path[ser_ind]}"

        new_ind = (new_dis_ind, new_ser_ind)

        return new_ind

    def reshape_str_for_circos(self,ind):
        return ind.replace(' ', '_').replace(',', '_').replace('-', '_').replace('/', '_').replace('(', '_').replace(
            ')', '_').replace('__', '_')

    def rename_bacteria(self,ind):
        x_0 = ind[0]
        x_1 = ind[1]

        if str.isdigit(str(x_0)) and str.isdigit(str(x_1)):
            return ind

        if "(" in str(x_0):
            x_0 = x_0.split(", '")[-1].split("')")[0]
            x_0 = x_0.replace('_', '')

        if "(" in str(x_1):
            x_1 = x_1.split(", '")[-1].split("')")[0]
            x_1 = x_1.replace('_', '')

        return (x_0, x_1)

    #endregion

    def plot_all_correlation_circus_plot(self):
        '''
        plot the 5 circus plot : per stage,disruption and creation. each one separately
        :return:
        '''
        ref_corr, med_corr, scd_corr = self.corr_api.return_correlation_by_index(0).copy(deep=True),\
                                       self.corr_api.return_correlation_by_index(1).copy(deep=True),\
                                       self.corr_api.return_correlation_by_index(2).copy(deep=True)

        ref_corr.index = ref_corr.index.map(lambda x: self.rename_bacteria(x))
        med_corr.index = med_corr.index.map(lambda x: self.rename_bacteria(x))
        scd_corr.index = scd_corr.index.map(lambda x: self.rename_bacteria(x))

        ref_corr["unique_index"] = ref_corr.index.map(
            lambda x: f"{sorted([str(x[0]), str(x[1])])[0]}_{sorted([str(x[0]), str(x[1])])[1]}")
        med_corr["unique_index"] = med_corr.index.map(
            lambda x: f"{sorted([str(x[0]), str(x[1])])[0]}_{sorted([str(x[0]), str(x[1])])[1]}")
        scd_corr["unique_index"] = scd_corr.index.map(
            lambda x: f"{sorted([str(x[0]), str(x[1])])[0]}_{sorted([str(x[0]), str(x[1])])[1]}")

        ref_corr = ref_corr.rename(columns={0: "ref"})
        med_corr = med_corr.rename(columns={0: "med"})
        scd_corr = scd_corr.rename(columns={0: "scd"})

        merged_data = ref_corr.merge(med_corr, left_on=["unique_index"], right_on=["unique_index"], how='outer')
        merged_data = merged_data.merge(scd_corr, left_on=["unique_index"], right_on=["unique_index"], how='outer')
        merged_data = merged_data.drop(columns=[col for col in merged_data.columns if "is_key" in str(col)])
        merged_data["level_0"] = merged_data["unique_index"].map(lambda x: x.split("_")[0])
        merged_data["level_1"] = merged_data["unique_index"].map(lambda x: x.split("_")[1])

        merged_data = merged_data.drop(columns="unique_index").set_index(["level_0", "level_1"])
        merged_data = merged_data.fillna(0)

        ref_final = merged_data["ref"]
        med_final = merged_data["med"]
        scd_final = merged_data["scd"]

        med_disrupted_final = pd.Series(index=merged_data.index, data=0.0)
        med_disrupted_final[((merged_data["ref"] > 0) & (merged_data["med"] == 0))] = \
        merged_data[((merged_data["ref"] > 0) & (merged_data["med"] == 0))]["ref"]

        scd_disrupted_final = pd.Series(index=merged_data.index, data=0.0)
        scd_disrupted_final[((merged_data["ref"] > 0) & (merged_data["scd"] == 0))] = \
        merged_data[((merged_data["ref"] > 0) & (merged_data["scd"] == 0))]["ref"]

        med_created_final = pd.Series(index=merged_data.index, data=0.0)
        med_created_final[((merged_data["ref"] == 0) & (merged_data["med"] > 0))] = \
        merged_data[((merged_data["ref"] == 0) & (merged_data["med"] > 0))]["med"]

        scd_created_final = pd.Series(index=merged_data.index, data=0.0)
        scd_created_final[((merged_data["ref"] == 0) & (merged_data["scd"] > 0))] = \
        merged_data[((merged_data["ref"] == 0) & (merged_data["scd"] > 0))]["scd"]

        ref_matrix_final_for_circos = ref_final.reset_index().pivot_table(index="level_0", columns="level_1").fillna(0)
        ref_matrix_final_for_circos.columns = ref_matrix_final_for_circos.columns.droplevel(0)

        ref_matrix_final_for_circos.columns = ref_matrix_final_for_circos.columns.map(lambda x: self.rename_index(str(x)))
        ref_matrix_final_for_circos.index = ref_matrix_final_for_circos.index.map(lambda x: self.rename_index(str(x)))

        ref_matrix_final_for_circos.columns = ref_matrix_final_for_circos.columns.map(
            lambda x: self.reshape_str_for_circos(str(x)))
        ref_matrix_final_for_circos.index = ref_matrix_final_for_circos.index.map(
            lambda x: self.reshape_str_for_circos(str(x)))

        melted_ref_matrix_final = pd.melt(ref_matrix_final_for_circos.reset_index(), id_vars='level_0',
                                          value_vars=ref_matrix_final_for_circos.columns)
        G = nx.from_pandas_edgelist(melted_ref_matrix_final, source='level_0', target='level_1')
        _communites = self.corr_api.extract_communites_from_g(G)

        _communites_ordered = reduce(lambda x, y: x + y, _communites)
        communites_order = {fea: i for i, fea in enumerate(_communites_ordered)}

        communites_color = {}
        communites_color_hex = {}
        color_to_color_name = {}
        for commun, color_name in zip(_communites, matplotlib.colors.cnames):
            _hax = matplotlib.colors.cnames[color_name]
            for fea in commun:
                _rgb = matplotlib.colors.to_rgb(_hax)
                r = int(_rgb[0] * 255)
                g = int(_rgb[1] * 255)
                b = int(_rgb[2] * 255)
                communites_color[fea] = f"{r},{g},{b}"
                communites_color_hex[fea] = matplotlib.colors.to_hex(_hax)
                color_to_color_name[f"{r},{g},{b}"] = color_name.capitalize()

        hvg1 = self.plot_ordered_colored_network(ref_final, "ref",communites_order,communites_color_hex)
        hvg2 = self.plot_ordered_colored_network(med_final, "med",communites_order,communites_color_hex)
        hvg3 = self.plot_ordered_colored_network(scd_final, "scd",communites_order,communites_color_hex)
        hvg4 = self.plot_ordered_colored_network(med_disrupted_final, "med_disrupted",communites_order,communites_color_hex)
        hvg5 = self.plot_ordered_colored_network(scd_disrupted_final, "scd_disrupted",communites_order,communites_color_hex)
        hvg6 = self.plot_ordered_colored_network(med_created_final, "med_created",communites_order,communites_color_hex)
        hvg7 = self.plot_ordered_colored_network(scd_created_final, "scd_created",communites_order,communites_color_hex)

        org_plot = hvg1 + hvg2 + hvg3
        show(hv.render(org_plot))

        dis_plot = hvg4 + hvg5
        show(hv.render(dis_plot))

        cre_plot = hvg6 + hvg7
        show(hv.render(cre_plot))

    def plot_single_correlation_graph(self,stage_idx):
        diet_df = self.corr_api.return_correlation_by_index(stage_idx).copy(deep=True)
        hv_G = self.plot_correlation_graph_with_communites(diet_df)

        labels = hv.Labels(hv_G.nodes, ['x', 'y'], 'name')

        return (hv_G * labels , hv_G)

