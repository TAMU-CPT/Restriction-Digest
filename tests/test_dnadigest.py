import unittest
import dnadigest

class TestDnaDigest(unittest.TestCase):

    def test_process_data_no_enzymes(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
        assert True
    def test_no_restriction_site(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
    def test_one_restriction_site(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
        break
    def test_multi_restriction_sites(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
        break
    def test_one_site_requiring_loop(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
        break
    def test_multi_sites_requiring_loop(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')

