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

    def test_string_cutter(self):
        dd = dnadigest.Dnadigest()
        self.assertEqual(dd.string_cutter('qqqqqqACTGnnnnnn', 'ACTG', 2, 'linear'), ['qqqqqqAC', 'TGnnnnnn'])
        self.assertEqual(dd.string_cutter('qqqqqqACTGnnnnnn', 'ACTG', 1, 'linear'), ['qqqqqqA', 'CTGnnnnnn'])

        # Heavy duty testing :)
        cut_concatamer = ('n' * 10 + 'ACC', 'CCG' + 'q' * 10)
        # Join the two halves
        concatamer = ''.join(cut_concatamer)
        # A bunch of copies
        full_concatamer = concatamer * 5
        # Generate a set, chopping at each location in the original concatamer
        concatamers = [full_concatamer[x:] + full_concatamer[0:x] for x in
                       range(len(concatamer))]
        for c in concatamers:
            self.assertEqual(dd.string_cutter(c, 'ANNNNG', 3, 'circular'),
                             [cut_concatamer[1] + cut_concatamer[0] for x in range(5)])
            self.assertEqual(dd.string_cutter(c, 'CNNNNT', 3, 'circular'),
                             [cut_concatamer[1] + cut_concatamer[0] for x in range(5)])
