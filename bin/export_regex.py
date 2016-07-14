#!/usr/bin/env python
import argparse
import dnadigest
import json


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool',
                                     epilog="")
    # Input filename, required
    parser.add_argument('--data', help='Enzyme cut site dataset')
    args = parser.parse_args()

    dd = dnadigest.DnaDigest(data_file=args.data)
    output_data = {}
    for enzyme in dd.library:
        e = dd.library[enzyme]

        f, r = e._gen_regex_str()
        reg_f = e.iupac_to_regex(f)
        reg_r = e.iupac_to_regex(r)
        output_data[enzyme] = [
            reg_f, reg_r
        ]

    print json.dumps(output_data, indent=2)
