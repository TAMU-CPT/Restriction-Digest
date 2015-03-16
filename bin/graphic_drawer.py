#!/usr/bin/env python
from Bio import SeqIO
import argparse
import dnadigest


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool',
                                     epilog="")
    # Input filename, required
    parser.add_argument('file', help='Input fasta genome(s)')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--data', help='Enzyme cut site dataset',
                        default='enzyme_data.yaml')
    # A list of one or more strings, at the end
    parser.add_argument('enzyme', help='Comma separated list of enzymes')
    args = parser.parse_args()

    dd = dnadigest.Dnadigest()

    print '<html><head></head><body>'
    for record in SeqIO.parse(args.file, 'fasta'):
        cut_sites = dd.find_cut_sites(str(record.seq), args.enzyme.split(','))
        print '<h1>%s</h1>' % record.id
        print dd.drawer(len(record.seq), cut_sites, sequence_id=record.id)
    print '</body></html>'
