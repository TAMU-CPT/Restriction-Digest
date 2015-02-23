import unittest
import dnadigest

class TestDnaDigest(unittest.TestCase):


    def test_expansion(self):
        dd = dnadigest.Dnadigest()
        self.assertEqual(dd.expand_multiple('A3N2B9'), 'AAANNBBBBBBBBB')

    def test_dict_filter(self):
        orig = {'A': [1,2,3], 'B': ['c', 'd', 'e'], 'B2': ['0']}
        filtered_1 = {'B': ['c', 'd', 'e'], 'B2': ['0']}
        filtered_2 = {'A': [1,2,3]}

        dd = dnadigest.Dnadigest()
        self.assertEqual(dd.enzyme_dict_filter(orig, ['B', 'B2']), filtered_1)
        self.assertEqual(dd.enzyme_dict_filter(orig, ['A']), filtered_2)
