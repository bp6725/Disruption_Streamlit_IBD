from unittest import TestCase
from pipelines.plot_disruption_illustration import PlotDisruptionIllustration


class TestDisruptionIllustration(TestCase):
    def test_illustrate_all_permutation(self):
        da = PlotDisruptionIllustration()
        da.illustrate_disruption()
        self.assert_(True)
