#!/usr/bin/env python
from Bio import SeqIO, Seq
import argparse
import dnadigest


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool')
    # Input filename, required
    parser.add_argument('file', help='Input fasta genome(s)')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--data', help='Enzyme cut site dataset', default='enzyme_data.yaml')
    # A list of one or more strings, at the end
    parser.add_argument('enzyme', metavar='E', type=str, nargs='+', help='Space separated list of enzymes')
    args = parser.parse_args()

    dd = dnadigest.Dnadigest()
    enzyme_dict = dd.get_dict(args.data)

    for record in SeqIO.parse(args.file,'fasta'):
        fragments, status = dd.process_data(str(record.seq), enzyme_dict,
                                    cut_with=args.enzyme)

        for i, fragment in enumerate(fragments):
            fragseq = Seq.Seq(fragment)
            print '>%s_%s [%s;status=%s]\n%s\n' % (record.id, i,
                                                   record.description, status,
                                                   fragseq)
