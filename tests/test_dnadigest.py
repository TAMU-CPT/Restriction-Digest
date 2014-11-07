import unittest
import dnadigest

class TestDnaDigest(unittest.TestCase):

    def test_process_data_no_enzymes(self):
        dd = dnadigest.Dnadigest()
        enzyme_dict = dd.get_dict('enzyme_data.yaml')
        assert True

