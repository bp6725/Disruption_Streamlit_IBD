from unittest import TestCase
from data_process.data_api import DataAPI


class TestDataAPI(TestCase):
    def test_return_unagg_global_diet_to_idx(self):
        data_api = DataAPI()
        res = data_api.return_unagg_global_diet_to_idx()
        # self.assert_(True)
