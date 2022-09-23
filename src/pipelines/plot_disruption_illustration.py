from correlations_api import CorrelationsAPI
from data_api import DataAPI
import pandas as pd
from analyses.analyses_api import AnalysesApi
from disruption.disruption_analysis import DisruptionAnalysis
import seaborn as sns
from IPython.display import display
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st
import streamlit as st
from global_utils import GlobalUtils


class PlotDisruptionIllustration():
    def __init__(self,feature_trhs=[0.1,0.3,-0.1,-0.3],diet_combinations = [(0,1),(1,0),(0,2),(2,0)],stage_for_disruption_prediction = "Before_MED",
                 corrected_pval_for_prediction = 0.25,disruption_analysis = None,corr_api = None,data_api = None,analyses_api = None):
        self.feature_trhs=feature_trhs
        self.diet_combinations = diet_combinations
        self.stage_for_disruption_prediction = stage_for_disruption_prediction
        self.corrected_pval_for_prediction = corrected_pval_for_prediction

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

    #region ploting tools

    def _genSankey(self, df, cat_cols=[], value_cols='', title='Sankey Diagram'):
        # maximum of 6 value cols -> 6 colors
        colorPalette = ['#4B8BBE', '#306998', '#FFE873', '#FFD43B', '#646464']
        labelList = []
        colorNumList = []
        for catCol in cat_cols:
            labelListTemp = list(set(df[catCol].values))
            colorNumList.append(len(labelListTemp))
            labelList = labelList + labelListTemp

        # remove duplicates from labelList
        labelList = list(dict.fromkeys(labelList))

        # define colors based on number of levels
        colorList = []
        for idx, colorNum in enumerate(colorNumList):
            colorList = colorList + [colorPalette[idx]] * colorNum

        # transform df into a source-target pair
        for i in range(len(cat_cols) - 1):
            if i == 0:
                sourceTargetDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
                sourceTargetDf.columns = ['source', 'target', 'count']
            else:
                tempDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
                tempDf.columns = ['source', 'target', 'count']
                sourceTargetDf = pd.concat([sourceTargetDf, tempDf])
            sourceTargetDf = sourceTargetDf.groupby(['source', 'target']).agg({'count': 'sum'}).reset_index()

        # add index for source-target pair
        sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
        sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))

        # creating the sankey diagram
        data = dict(
            type='sankey',
            node=dict(
                pad=15,
                thickness=20,
                line=dict(
                    color="black",
                    width=0.5
                ),
                label=labelList,
                color=colorList
            ),
            link=dict(
                source=sourceTargetDf['sourceID'],
                target=sourceTargetDf['targetID'],
                value=sourceTargetDf['count']
            )
        )

        layout = dict(
            title=title,
            font=dict(
                size=10
            )
        )

        fig = dict(data=[data], layout=layout)
        return fig

    def _plot_correlation(self, chosen_idx, _df, _df_for_regline, to_color, title, color_tag):
        _x = str(chosen_idx[0])
        _y = str(chosen_idx[1])

        x = _df[chosen_idx[0]]
        y = _df[chosen_idx[1]]

        h = pd.Series(data=f"non-{color_tag}", index=_df.index)

        h.loc[to_color] = color_tag

        _data = pd.DataFrame({_x: x, _y: y, "dis": h})

        #     sns.regplot(x=_x, y=_y, data=_data,ax=ax)

        if (len(_df.index) > len(_df.index.unique())):
            new_data = pd.DataFrame(index=(_data.index).unique(), columns=_data.columns)

            new_data[_x] = _data.groupby(level=0).mean()[_x]
            new_data[_y] = _data.groupby(level=0).mean()[_y]
            new_data["dis"] = _data.groupby(level=0)["dis"].first()

            _data = new_data.copy(deep=True)

        if color_tag == "disruptor":
            #         graph = sns.lmplot(x=_x, y=_y, hue='dis', data=_data, fit_reg=False)
            graph = sns.lmplot(x=_x, y=_y, hue='dis', palette=['purple', 'yellow'], data=_data, fit_reg=False)
        else:
            graph = sns.lmplot(x=_x, y=_y, hue='dis', palette=[ 'gray','black'], data=_data, fit_reg=False)

        x_for_rl = _df_for_regline[chosen_idx[0]]
        y_for_rl = _df_for_regline[chosen_idx[1]]

        h_for_rl = pd.Series(data=f"non-{color_tag}", index=_df_for_regline.index)

        _data_for_rl = pd.DataFrame({_x: x_for_rl, _y: y_for_rl, "dis": h_for_rl})

        sns.regplot(x=_x, y=_y, data=_data_for_rl, scatter=False, ax=graph.axes[0, 0])
        graph.fig.suptitle(title)
        return graph

    def _illustrate_disruption(self,exp_feature, relevent_correlations, seg_disruption_df, pn_in_corr, ref_z_disruption, z_dis_net,
                               diet_df, backbone_df, melted_corr_and_disruption, seg_imp_data, feature_trh,
                               is_creation_correlation=False):

        predictions = melted_corr_and_disruption[melted_corr_and_disruption["corrected_pval"] < self.corrected_pval_for_prediction]
        for idx in range(seg_disruption_df.shape[0]):

            chosen_idx = seg_disruption_df.index[idx]

            if not self._is_corr_in_pre_picked_disruption(chosen_idx):
                continue

            st.text(f"------{self.data_api.feature_name_dict[str(chosen_idx[0])]},{self.data_api.feature_name_dict[str(chosen_idx[1])]}------")

            disruption_values = z_dis_net.loc[chosen_idx]

            _pn_in_corr = pn_in_corr.loc[chosen_idx][pn_in_corr.loc[chosen_idx]].index.to_list()

            disruption_values[disruption_values.index.difference(pd.Index(_pn_in_corr))] = 0
            positive_disrupted = disruption_values[disruption_values > 0].index.to_list()

            # disrupted pn relevent to disruption score - clean_disruption_df
            _clean_disruption = seg_disruption_df.loc[chosen_idx].copy(deep=True)
            disrupted_pn = _clean_disruption[_clean_disruption >= 0.05].index.to_list()
            non_disrupted_pn = self.intersection(_clean_disruption[_clean_disruption < 0.05].index.to_list(),
                                            _pn_in_corr)

            try:
                is_created_str = "created" if is_creation_correlation else ""
                graph = self._plot_correlation(chosen_idx, backbone_df.rank(), backbone_df.rank(), _pn_in_corr,
                                 f"reference rank {is_created_str}", "contributors")
                st.pyplot(graph.fig)

                graph = self._plot_correlation(chosen_idx, self.diet_rank_in_refrence(diet_df.copy(deep=True), backbone_df.copy(deep=True),
                                                                              _pn_in_corr), backbone_df.rank(), positive_disrupted,
                                 f"diet rank {is_created_str}", "disruptor")
                st.pyplot(graph.fig)

                #dist plot:
                kwargs = dict(hist_kws={'alpha': .6}, kde_kws={'linewidth': 2, 'bw': 0.15})
                fig = plt.figure(figsize=(10, 7), dpi=80)
                sns.distplot(diet_df[exp_feature].loc[disrupted_pn].dropna(), color="dodgerblue", label="disrupted pn",
                             **kwargs)
                sns.distplot(diet_df[exp_feature].loc[~diet_df.index.isin(disrupted_pn)].dropna(), color="orange",
                             label="non disrupted pn", **kwargs)
                plt.legend()
                st.pyplot(fig)


                fig = plt.figure(figsize=(10, 7), dpi=80)
                sns.distplot(diet_df[exp_feature].loc[set(disrupted_pn).intersection(_pn_in_corr)].dropna(), hist=False,
                             color="purple", label="disrupted pn", **kwargs)
                sns.distplot(diet_df[exp_feature].loc[set(non_disrupted_pn).intersection(_pn_in_corr)].dropna(),
                             hist=False, color="yellow", label="non disrupted pn", **kwargs)
                plt.legend()
                st.pyplot(fig)

                if feature_trh >= 0:
                    pn_over_crp = diet_df[exp_feature][diet_df[exp_feature] >= feature_trh].index
                else:
                    pn_over_crp = diet_df[exp_feature][diet_df[exp_feature] <= feature_trh].index

                # sankey plot
                all_pns = set(diet_df.index.to_list())

                n_pn_in_corr = pd.Index(_pn_in_corr)

                disrupted = len(set(disrupted_pn))
                not_disrupted = len(set(non_disrupted_pn))

                reponders_not_in_corr = len(self.intersection(self.Diff(n_pn_in_corr.to_list(), all_pns), pn_over_crp.to_list()))
                not_in_corr_not_responders = len(
                    self.intersection(self.Diff(n_pn_in_corr.to_list(), all_pns), self.Diff(pn_over_crp.to_list(), all_pns)))

                disrupted_responders = len(self.intersection(disrupted_pn, pn_over_crp.to_list()))
                disrupted_non_responders = len(self.intersection(disrupted_pn, self.Diff(pn_over_crp.to_list(), all_pns)))

                not_disrupted_responders = len(self.intersection(non_disrupted_pn, pn_over_crp.to_list()))
                not_disrupted_not_responders = len(self.intersection(non_disrupted_pn, self.Diff(pn_over_crp.to_list(), all_pns)))

                # print(f"disrupted_pn : {disrupted_pn}")
                # print(f"non_disrupted_pn : {non_disrupted_pn}")
                # print(f"n_pn_in_corr : {n_pn_in_corr}")
                # print(f"not n_pn_in_corr : {self.Diff(n_pn_in_corr.to_list(), all_pns)}")
                # print(f"incresed : {pn_over_crp}")
                # print(f"decresed : {self.Diff(pn_over_crp.to_list(), all_pns)}")

                labels_dielation = ["contributors to ref corr", "non contributors to ref corr", "disrupted patients",
                                    "non disrupted patients", "increased CRP", "decreased CRP"]
                labels_creation = ["contributors to diet corr", "non contributors to diet corr", "wasnt in corr in ref",
                                   "was in corr in ref", "increased CRP", "decreased CRP"]

                fig = go.Figure(data=[go.Sankey(
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="black", width=0.5),
                        label=["correlation contributors", "not in correlation", "disrupted patients",
                               "non disrupted patients", "increased CRP", "decreased CRP"],
                        color=['orange', "blue", 'purple', 'yellow', 'red', 'green']
                    ),
                    link=dict(
                        source=[0, 0, 1, 1, 2, 2, 3, 3],
                        target=[2, 3, 4, 5, 4, 5, 4, 5],
                        value=[disrupted, not_disrupted, reponders_not_in_corr, not_in_corr_not_responders,
                               disrupted_responders, disrupted_non_responders, not_disrupted_responders,
                               not_disrupted_not_responders]
                    ))])

                fig.update_layout(title_text="Basic Sankey Diagram", font_size=15)
                st.plotly_chart(fig)


                corr_str = f"{chosen_idx[0]}-{chosen_idx[1]}"
                try:
                    if corr_str in predictions.index.get_level_values(0):
                        #             predictions
                        rel_dis = disruption_values[_pn_in_corr].sort_values(ascending=True)
                        predictor = predictions.loc[corr_str].idxmin()["corrected_pval"]
                        mat = seg_imp_data[predictor].loc[rel_dis.index]

                        fig, axes = plt.subplots(ncols=2, sharey=True)
                        axes[0].barh(rel_dis.index, rel_dis, align='center')
                        axes[0].set_xlabel('disruption score')
                        axes[1].barh(mat.index, mat, align='center')
                        axes[1].set_xlabel(predictor)
                        axes[0].invert_xaxis()
                        fig.suptitle("disruption score vs metabolism at reference")
                        st.pyplot(fig)

                except:
                    print("hello")
            except Exception as  e:
                print(f"error : {e}")

    #endregion

    #region tools

    def is_in_list_of_disruptions(self,idx, list_of_disruptions):
        return any([((str(lod[0]) == str(idx[0])) and (str(lod[1]) == str(idx[1]))) or (
                    (str(lod[0]) == str(idx[1])) and (str(lod[1]) == str(idx[0]))) for lod in list_of_disruptions])

    def _is_corr_in_pre_picked_disruption(self,chosen_idx):
        return True
        return     ((61587 in chosen_idx) and (35879 in chosen_idx) or (61587 in chosen_idx) and (
                     35879 in chosen_idx) or (61857 in chosen_idx) and (54890 in chosen_idx) or (
                     61857 in chosen_idx) and (32398 in chosen_idx) or (61857 in chosen_idx) and (
                     62374 in chosen_idx) or (62798 in chosen_idx) and (36600 in chosen_idx) or (
                     62798 in chosen_idx) and (48733 in chosen_idx) or (40007 in chosen_idx) and (
                     38102 in chosen_idx) or (40007 in chosen_idx) and (52608 in chosen_idx) or (
                     40007 in chosen_idx) and (62746 in chosen_idx) or (40007 in chosen_idx) and (
                     62746 in chosen_idx) or (62746 in chosen_idx) and (38102 in chosen_idx) or (
                     62746 in chosen_idx) and (52608 in chosen_idx) or (1561 in chosen_idx) and (
                     33418 in chosen_idx) or (1561 in chosen_idx) and (33418 in chosen_idx) or (
                     33418 in chosen_idx) and (52473 in chosen_idx) or (33418 in chosen_idx) and (
                     52473 in chosen_idx) or (47666 in chosen_idx) and (44876 in chosen_idx) or (
                     47666 in chosen_idx) and (62280 in chosen_idx) or (47666 in chosen_idx) and (
                     62280 in chosen_idx) or (15747 in chosen_idx) and (43256 in chosen_idx) or (
                     38652 in chosen_idx) and (43256 in chosen_idx) or (38652 in chosen_idx) and (
                     57651 in chosen_idx) or (57513 in chosen_idx) and (57651 in chosen_idx) or (
                     57635 in chosen_idx) and (57804 in chosen_idx) or (57637 in chosen_idx) and (
                     57651 in chosen_idx) or (1105 in chosen_idx) and (33447 in chosen_idx) or (
                     1105 in chosen_idx) and (34035 in chosen_idx) or (1105 in chosen_idx) and (
                     21232 in chosen_idx) or (1105 in chosen_idx) and (32506 in chosen_idx) or (
                     1105 in chosen_idx) and (32506 in chosen_idx) or (1105 in chosen_idx) and (
                     33587 in chosen_idx) or (1105 in chosen_idx) and (33971 in chosen_idx) or (
                     1105 in chosen_idx) and (54920 in chosen_idx) or (32980 in chosen_idx) and (
                     32415 in chosen_idx) or (32980 in chosen_idx) and (57443 in chosen_idx) or (
                     32980 in chosen_idx) and ("Calprotectin" in chosen_idx) or (32980 in chosen_idx) and (
                     47123 in chosen_idx) or (32415 in chosen_idx) and (15747 in chosen_idx) or (
                     32415 in chosen_idx) and (33587 in chosen_idx) or (32415 in chosen_idx) and (
                     33968 in chosen_idx) or (32415 in chosen_idx) and (44877 in chosen_idx) or (
                     44877 in chosen_idx) and (57432 in chosen_idx) or (44877 in chosen_idx) and (
                     57517 in chosen_idx) or (44877 in chosen_idx) and (57517 in chosen_idx) or (
                     57432 in chosen_idx) and (57517 in chosen_idx) or (57432 in chosen_idx) and (
                     57525 in chosen_idx))

    def Diff(self,li1, li2):
        return (list(list(set(li1) - set(li2)) + list(set(li2) - set(li1))))

    def intersection(self,lst1, lst2):
        return list(set(lst1) & set(lst2))

    def cast_to_int(self,st):
        if (st.isdigit()):
            return int(st)
        return st

    def diet_rank_in_refrence(self,diet_df, backbone_df, _pn_in_corr):
        _added_df = diet_df.copy(deep=True)
        value_of_diet_in_backbone = pd.DataFrame(columns=backbone_df.columns)

        for pn in _pn_in_corr:
            _pn_dis = f"{pn}_dis"
            _ref_df = backbone_df.copy(deep=True)
            if type(_added_df.loc[pn]) is pd.Series:
                new_line = pd.Series(index=_added_df.loc[pn].index, data=_added_df.loc[pn], name=_pn_dis)
                new_ref_df = _ref_df.append(new_line)
            else:
                new_df = _added_df.loc[pn].copy(deep=True)
                new_df.index = new_df.index.map(lambda x: _pn_dis)
                new_ref_df = _ref_df.append(new_df)

            value_of_diet_in_backbone = value_of_diet_in_backbone.append(new_ref_df.rank().loc[_pn_dis])

        value_of_diet_in_backbone.index = value_of_diet_in_backbone.index.map(lambda x: x.split("_")[0])
        return value_of_diet_in_backbone

    @GlobalUtils.cache_me
    def _calculate_all_for_illustration(self,backbone_index,diet_index,feature_trh):
        self.disruption_analysis.update_class_properties(
            {"reference_diet_idx": backbone_index, "treatment_diet_idx": diet_index})

        relevent_correlations = self.corr_api.return_correlation_by_index(self.disruption_analysis.reference_diet_idx)
        ref_z_disruption = self.disruption_analysis.calculate_ref_z_disruption()
        z_dis_net = self.disruption_analysis.calculate_z_disruption_networks()
        clean_disruption_df = self.disruption_analysis.get_disruption_score()
        pn_in_corr = self.disruption_analysis.get_patient_in_correlation()

        diet_df = self.disruption_analysis.get_diet_df_with_delta()
        backbone_df = self.data_api.return_unagg_by_index(backbone_index)

        # only disruption of sig features - same indexes as enriched_features_df. its only for reshapeing the df
        enriched_features_df = self.disruption_analysis.get_enriched_features(feature_trh,
                                                                              only_pre_picked_features=True)

        results_per_exp_feature = {}
        for exp_feature,features_df in enriched_features_df.items() :

            disruption_of_sig_modules = z_dis_net.loc[features_df.index]
            disruption_of_sig_modules.index = disruption_of_sig_modules.index.map(lambda x: f"{x[0]}-{x[1]}")

            if disruption_of_sig_modules.empty:
                return None

            sc_imp_data = self.data_api.return_relvent_sc_imp_data(stage=self.stage_for_disruption_prediction)
            seg_imp_data = sc_imp_data[self.data_api.union_ref_df.columns.intersection(sc_imp_data.columns)]

            melted_corr_and_disruption = self.analyses_api.calculate_joined_correlation_to_disruption(
                disruption_of_sig_modules, seg_imp_data, pn_in_corr)

            if melted_corr_and_disruption is None:
                return None

            all_melted_corr_and_disruption = melted_corr_and_disruption[
                melted_corr_and_disruption["corrected_pval"] < self.corrected_pval_for_prediction]
            all_melted_corr_and_disruption = all_melted_corr_and_disruption.copy(deep=True)

            list_of_disruptions = all_melted_corr_and_disruption.index.get_level_values(0).drop_duplicates().map(
                lambda x: (str(self.cast_to_int(x.split("-")[0])), str(self.cast_to_int(x.split("-")[1])))).tolist()

            all_melted_corr_and_disruption.index = all_melted_corr_and_disruption.index.map(self.data_api.rename_index)

            seg_disruption_df = clean_disruption_df.loc[
                clean_disruption_df.index.map(lambda x: self.is_in_list_of_disruptions(x, list_of_disruptions))]

            results_per_exp_feature[exp_feature] = {"relevent_correlations":relevent_correlations.copy(deep=True),
                                                    "seg_disruption_df":seg_disruption_df.copy(deep=True),
                                                    "pn_in_corr":pn_in_corr.copy(deep=True),
                                                    "ref_z_disruption":ref_z_disruption.copy(deep=True),
                                                    "z_dis_net":z_dis_net.copy(deep=True),
                                                    "diet_df":diet_df.copy(deep=True),
                                                    "backbone_df":backbone_df.copy(deep=True),
                                                    "melted_corr_and_disruption":melted_corr_and_disruption.copy(deep=True),
                                                    "seg_imp_data":seg_imp_data.copy(deep=True),
                                                    "diet_index":diet_index}

        return results_per_exp_feature

    #endregion

    def illustrate_disruption(self,backbone_index,diet_index,feature_trh,for_cache = None):
        if for_cache is None :
            results_per_exp_feature = self._calculate_all_for_illustration(backbone_index,diet_index, feature_trh)
        else :
            results_per_exp_feature = self._calculate_all_for_illustration(backbone_index,
                                                                                                        diet_index,
                                                                                                        feature_trh,for_cache = for_cache)

        # print("start illustrate")
        for exp_feature, disruption_data in results_per_exp_feature.items() :
            relevent_correlations = disruption_data["relevent_correlations"]
            seg_disruption_df = disruption_data["seg_disruption_df"]
            pn_in_corr = disruption_data["pn_in_corr"]
            ref_z_disruption = disruption_data["ref_z_disruption"]
            z_dis_net = disruption_data["z_dis_net"]
            diet_df = disruption_data["diet_df"]
            backbone_df = disruption_data["backbone_df"]
            melted_corr_and_disruption = disruption_data["melted_corr_and_disruption"]
            seg_imp_data = disruption_data["seg_imp_data"]
            diet_index = disruption_data["diet_index"]

            st.text(f"------------------------------------{exp_feature}-------------------------------------------")
            self._illustrate_disruption(exp_feature,relevent_correlations, seg_disruption_df, pn_in_corr,
                                    ref_z_disruption, z_dis_net, diet_df, backbone_df,
                                    melted_corr_and_disruption, seg_imp_data, feature_trh, diet_index == 0)

    def get_all_possible_disruptions(self,diet_index,backbone_index,feature_trh):
        elevent_correlations, seg_disruption_df, pn_in_corr, \
        ref_z_disruption, z_dis_net, diet_df, backbone_df, \
        melted_corr_and_disruption, seg_imp_data, diet_index = self._calculate_all_for_illustration(diet_index,
                                                                                                    backbone_index,
                                                                                                    feature_trh)

        return seg_disruption_df.index.to_list()