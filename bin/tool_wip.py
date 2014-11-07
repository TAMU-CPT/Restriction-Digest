from Bio import SeqIO
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

    seqs = [str(record.seq) for record in SeqIO.parse(args.file,'fasta')]
    dd = dnadigest.Dnadigest()
    enzyme_dict = dd.get_dict('enzyme_data.yaml')

    dd.process_data(seqs, enzyme_dict, cut_with=args.enzyme)
