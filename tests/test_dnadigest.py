import unittest
import dnadigest

class TestDnaDigest(unittest.TestCase):


    def test_expansion(self):
        dd = dnadigest.Dnadigest()
        self.assertEqual(dd.expand_multiple('A3N2B9'), 'AAANNBBBBBBBBB')
