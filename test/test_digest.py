import unittest
import os
import dnadigest as dd


class TestDnaDigest(unittest.TestCase):

    def test_enzyme(self):
        e = dd.Enzyme(name='EcoRI', forward=('G', 'A2T2C'), reverse=('CTTAA', 'G'))
        self.assertEqual(e.expand_multiple('A3N2B9'), 'AAANNBBBBBBBBB')
        fr, rv = e.get_regex()
        self.assertTrue(fr.match('GAATTC'))
        self.assertFalse(rv.match('GAATTC'))

        self.assertFalse(fr.match('CTTAAG'))
        self.assertTrue(rv.match('CTTAAG'))

    def test_enzymelibrary(self):
        bsp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'bsp.yaml')
        el = dd.EnzymeLibrary(bsp)
        print el
        # sequence = 'qqqCGnCGwwwwwwTCCGGAeeeeeAGGCCTrrrrrGACNNNNNNGTCoooGCnGCooo'
