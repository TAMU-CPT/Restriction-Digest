#!/usr/bin/env python
"""Test cases for DNA Digest library"""
import unittest
import os
import dnadigest as dd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestDnaDigest(unittest.TestCase):
    """Test DnaDigest classes"""

    def test_enzyme(self):
        """test enzyme class"""
        enzyme = dd.Enzyme(name='EcoRI', forward=('G', 'A2T2C'), reverse=('CTTAA', 'G'))
        self.assertEqual(enzyme.expand_multiple('A3N2B9'), 'AAANNBBBBBBBBB')
        forward, reverse = enzyme.get_regex()
        self.assertTrue(forward.match('GAATTC'))
        self.assertFalse(reverse.match('GAATTC'))

        self.assertFalse(forward.match('CTTAAG'))
        self.assertTrue(reverse.match('CTTAAG'))

        self.assertEqual(enzyme.iupac_to_regex('KAYAK'),
                         '[GT]A[CT]A[GT]')

    def test_dnadigest(self):
        """test dnadigest"""
        bsp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'bsp.yaml')
        digester = dd.DnaDigest(bsp)
        sequence = SeqRecord(
            Seq('nnnCGnCGnnnnnnTCCGGAnnnnnAGGCCTnnnnnGACNNNNNNGTCnnnGCnGCnnn'), id='test')
        fragments, cut_sites, did_cut = digester.digest_sequence(sequence, ['Bsp6I'])
        self.assertTrue(did_cut)
        self.assertDictEqual(
            {65: ['Bsp6I']},
            cut_sites
        )

        self.assertEqual(
            'nGCnnnnnnCGn',
            str(fragments[0].seq)
        )

        self.assertEqual(
            'CGnnnnnnTCCGGAnnnnnAGGCCTnnnnnGACNNNNNNGTCnnnGC',
            str(fragments[1].seq)
        )
