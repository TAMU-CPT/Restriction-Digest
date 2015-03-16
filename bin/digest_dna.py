#!/usr/bin/env python
from Bio import SeqIO, Seq
import argparse
import dnadigest


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

    dd = dnadigest.Dnadigest(enzyme_data_file=args.data)
    template = '>%s_%s [orig=%s;status=%s;cut_with=%s]\n%s\n'

    for record in SeqIO.parse(args.file, 'fasta'):
        processed_results = dd.process_data(str(record.seq), cut_with=args.enzyme.split(','))

        for i, fragment in enumerate(processed_results['fragment_list']):
            fragseq = Seq.Seq(fragment)
            print template % (record.id, i,
                              record.description,
                              processed_results['status'],
                              ','.join(processed_results['cut_with']),
                              fragseq)
