from unittest import TestCase
from pipelines.plot_correlations_map import PlotCorrelationsMap


class TestPlotCorrelationsMap(TestCase):
    def test_plot_circus_plot_per_correlation(self):
        pcm = PlotCorrelationsMap()

        pcm.plot_all_correlation_circus_plot()

        self.assert_(True)


class TestPlotCorrelationsMap(TestCase):
    def test_plot_single_correlation_graph(self):
        pcm = PlotCorrelationsMap()

        pcm.plot_single_correlation_graph(0)

        self.assert_(True)

