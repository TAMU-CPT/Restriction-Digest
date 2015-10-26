#!/usr/bin/env python
import sys
import argparse
import dnadigest
from Bio import SeqIO


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool',
                                     epilog="")
    # Input filename, required
    parser.add_argument('file', help='Input fasta genome(s)')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--data', help='Enzyme cut site dataset')
    # A list of one or more strings, at the end
    parser.add_argument('enzyme', help='Comma separated list of enzymes')
    args = parser.parse_args()

    dd = dnadigest.DnaDigest(enzyme_data_file=args.data)
    enzymes = args.enzyme.split(',')

    for record in SeqIO.parse(args.file, 'fasta'):
        fragments, cut_sites, did_cut = dd.digest_sequence(record, enzymes)

        for i, fragment in enumerate(fragments):
            fragment.description = '[orig={};status={};cut_with={}]'.format(
                fragment.id, did_cut, enzymes)
            fragment.id = '%s_%s' % (fragment.id, i)

            SeqIO.write([fragment], sys.stdout, 'fasta')
