from unittest import TestCase
from pipelines.plot_disruption_map import PlotDisruptionMap

class TestPlotDisruptionMap(TestCase):
    def test_plot_enriched_disruption_network(self):
        pdm = PlotDisruptionMap()
        pdm.plot_enriched_disruption_network(0,1)
        self.assert_(True)
