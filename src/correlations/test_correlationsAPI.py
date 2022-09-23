from unittest import TestCase
from correlations_api import CorrelationsAPI

class TestCorrelationsAPI(TestCase):
    def test_build_correlation_dfs(self):
        alpah_fdr = 0.15
        corr_th = (0.5, 0.4, 0.3)
        res = CorrelationsAPI(alpah_fdr,corr_th)
        print('OK')
