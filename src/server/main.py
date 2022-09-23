import sys
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis\analyses")
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis\correlations")
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis\data_process")
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis\disruption")
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis\pipelines")
sys.path.append(r"C:\Repos\MBA_beilinson\disruption_analysis")

import asyncio
loop = asyncio.new_event_loop()
asyncio.set_event_loop(loop)

from plot_correlations_map import PlotCorrelationsMap
from plot_disruption_map import PlotDisruptionMap
from plot_disruption_illustration import PlotDisruptionIllustration

from disruption_analysis import DisruptionAnalysis
from correlations_api import CorrelationsAPI

import streamlit as st
import holoviews as hv
from bokeh.plotting import figure
import bokeh
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import matplotlib.pyplot as plt


mode_selectbox = st.sidebar.selectbox("analysis", ["correlation plot", "disruption map", "disruption illustration"],
                     index=2, key = "mode_selectbox")

diet_opt = ['reference', 'med', 'scd', 'med with nut', 'scd with nut', 'before med with nuts'
                        ,'before scd with nuts','refrence with nuts']

baseline_diet_selectbox = st.sidebar.selectbox("baseline", diet_opt, index=1, key="baseline_diet_selectbox")

if mode_selectbox in ["disruption map", "disruption illustration"] :
    treatment_diet_selectbox = st.sidebar.selectbox("Diet", diet_opt, index=0, key="treatment_diet_selectbox")

alpha_fdr_selection = st.sidebar.slider("correlation alpha FDR", min_value=0.05, max_value=0.2,
                                        value=0.15, step=0.05,key="alpha_fdr_selection")
corr_thr_selection = st.sidebar.slider("correlation threshold", min_value=0.3, max_value=0.6,
                                       value=0.4, step=0.1,key="corr_thr_selection")
cross_omics_corr_thr_selection = st.sidebar.slider("correlation threshold - cross omics", min_value=0.3, max_value=0.6,
                                                   value=0.4, step=0.1,key="cross_omics_corr_thr_selection")
key_features_corr_thr_selection = st.sidebar.slider("correlation threshold - to key features", min_value=0.3, max_value=0.6,
                                                    value=0.3, step=0.1,key="key_features_corr_thr_selection")

if mode_selectbox in ["disruption map", "disruption illustration"] :
    explanatory_features_opts = ['CRP', 'Cholesterol', 'HDL', 'LDL', 'Calprotectin', 'Shannon', 'Zonulin', 'PDAI', 'IL6']
    explanatory_features_multiselect = st.sidebar.multiselect("response measure",explanatory_features_opts,["CRP","Calprotectin"] )
    feature_trh_selection = st.sidebar.slider("fold change for responder", min_value=-0.3, max_value=0.3,
                                                    value=0.1, step=0.1,key="feature_trh_selection")
    enrichment_trh_selector = st.sidebar.slider("p-value of disruption enrichment to respond", min_value= 0.05 , max_value=0.25,
                                                    value=0.15, step=0.05,key="enrichment_trh_selector")

    edge_cut_selector = st.sidebar.slider("minimum number of PN in disruption ", min_value= 1 , max_value=10,
                                                    value=3, step=1,key="edge_cut_selector")

if mode_selectbox == "disruption illustration":
    corrected_pval_for_prediction = st.sidebar.slider("corrected alpha for disruption prediction", min_value= 0.05 , max_value=0.3,
                                                    value=0.25, step=0.05,key="corrected_pval_for_prediction")

# plor_correlation_plot(alpha_fdr_selection,(corr_thr_selection, cross_omics_corr_thr_selection, key_features_corr_thr_selection),baseline_diet_selectbox)

if mode_selectbox == "correlation plot" :
    pcm = PlotCorrelationsMap(alpah_fdr = alpha_fdr_selection, corr_th = (corr_thr_selection, cross_omics_corr_thr_selection, key_features_corr_thr_selection))
    with_labels_figure,no_labels_figure =  pcm.plot_single_correlation_graph(diet_opt.index(baseline_diet_selectbox))

    plot = hv.render(with_labels_figure.opts(fontscale=2, width=550, height=700, title='correlation plot'), backend='bokeh')
    st.bokeh_chart(plot,use_container_width=True)

    plot2 = hv.render(with_labels_figure.opts(fontscale=2, width=550, height=700, title='correlation plot'),
                     backend='bokeh')
    st.bokeh_chart(plot2, use_container_width=True)

if mode_selectbox in ["disruption map", "disruption illustration"] :
    reference_diet_idx = diet_opt.index(baseline_diet_selectbox)
    treatment_diet_idx = diet_opt.index(treatment_diet_selectbox)
    alpah_fdr = alpha_fdr_selection
    feature_trh = feature_trh_selection
    corr_th = (corr_thr_selection, cross_omics_corr_thr_selection, key_features_corr_thr_selection)
    edge_cut = edge_cut_selector
    enrichment_trh = enrichment_trh_selector
    explanatory_features = [ef for ef in explanatory_features_multiselect]

    da = DisruptionAnalysis(reference_diet_idx = reference_diet_idx, treatment_diet_idx = treatment_diet_idx,
                            alpah_fdr = alpah_fdr, corr_th = corr_th, edge_cut = edge_cut,
                            enrichment_trh = enrichment_trh, explanatory_features = explanatory_features)
    ca = CorrelationsAPI(alpah_fdr = alpah_fdr,corr_th=corr_th)


if mode_selectbox == "disruption map" :
    pdm = PlotDisruptionMap(disruption_analysis = da,corr_api = ca,explanatory_features = explanatory_features)

    figs_to_plot = pdm.plot_enriched_disruption_network(reference_diet_idx = reference_diet_idx,
                                                        treatment_diet_idx = treatment_diet_idx,
                                                        feature_trh=feature_trh)

    for enriched_feature,hvplot in figs_to_plot.items() :
        plot = hv.render(hvplot.opts(fontscale=2, width=550, height=700, title=enriched_feature), backend='bokeh')
        st.bokeh_chart(plot, use_container_width=True)


if mode_selectbox == "disruption illustration":
    feature_trhs = [feature_trh]
    diet_combinations = [(reference_diet_idx, treatment_diet_idx)]
    corrected_pval_for_prediction =corrected_pval_for_prediction

    for_cache = {"ref_idx": reference_diet_idx,
                 "tre_idx": treatment_diet_idx,
                 "alpah": alpah_fdr,
                 "corr": corr_thr_selection,
                 "cross_corr": cross_omics_corr_thr_selection,
                 "key_corr": key_features_corr_thr_selection,
                 "feature": feature_trh,
                 "enrichment": enrichment_trh_selector,
                 "edge_cut": edge_cut_selector,
                 "corrected_p": corrected_pval_for_prediction,
                 "expl_features": explanatory_features}


    pdi = PlotDisruptionIllustration(feature_trhs=feature_trhs,diet_combinations = diet_combinations,
                                     corrected_pval_for_prediction = corrected_pval_for_prediction,
                                     disruption_analysis=da, corr_api=ca)
    # all_possible_disruptions = pdi.get_all_possible_disruptions(reference_diet_idx,treatment_diet_idx,feature_trh)
    # st.text([f"{pdi.data_api.feature_name_dict[str(dist[0])]},{pdi.data_api.feature_name_dict[str(dist[1])]}" for dist in all_possible_disruptions])
    pdi.illustrate_disruption(reference_diet_idx,treatment_diet_idx,feature_trh,for_cache = for_cache)





