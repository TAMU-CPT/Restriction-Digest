#!/usr/bin/env python
import re
import yaml
from filecache import filecache
import logging
logging.basicConfig()
log = logging.getLogger('tmp')
import urllib2
mapping = {
    '1': 'enzyme',
    '2': 'isoscizomers',
    '3': 'cut',
    '4': 'methylation_site',
    '5': 'microorganism',
    '6': 'source',
    '7': 'commercial_availability',
    '8': 'references',
}

@filecache(24 * 60 * 60)
def fetchfile():
    URL = "ftp://ftp.neb.com/pub/rebase/withrefm.txt"
    response = urllib2.urlopen(URL)
    data = response.read()
    return data

rep = {
    'A': 'T',
    'B': 'V',
    'C': 'G',
    'D': 'H',
    'G': 'C',
    'H': 'D',
    'K': 'M',
    'M': 'K',
    'N': 'N',
    'N': 'N',
    'R': 'Y',
    'S': 'W',
    'T': 'A',
    'V': 'B',
    'W': 'S',
    'Y': 'R',
}
# stackoverflow.com/questions/6116978
rep = dict((re.escape(k), v) for k, v in rep.iteritems())
pattern = re.compile('|'.join(rep.keys()))


def com(sequence):
    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], sequence)
    return text

def parse_data(data):
    header = data[0:data.index('<1>')]
    data = data[data.index('<1>'):]
    enzymes = data.split('\n\n')
    fixed_data = {}

    gb = {'+': 0, '-': 0}
    for enzyme in enzymes:
        es = enzyme[1:].split('\n<')
        ed = {}
        # Parse
        for x in es:
            try:
                idx, data = x.split('>')
                ed[mapping[idx]] = data
            except:
                continue

        if len(ed.keys()) > 0:
            # Reformat
            # ed['isoscizomers'] = ed['isoscizomers'].split(',')
            if 'isoscizomers' in ed:
                del ed['isoscizomers']
            # Build new data
            edcut = ed['cut']

            if edcut.count('^') == 1:
                print '+', ed['enzyme']
                cut_1, cut_2 = edcut.split('^')
                cut_3 = com(edcut).replace('^', '')
                cut_4 = cut_3[0: len(cut_3) - len(cut_1)]
                cut_5 = cut_3[len(cut_3) - len(cut_1):]
                ed['cut'] = [
                    "5' ---%s  %s--- 3'" % (cut_1, cut_2),
                    "3' ---%s  %s--- 5'" % (cut_4, cut_5)
                ]

                ed['recognition_sequence'] = [
                    "5' %s" % edcut.replace('^', ''),
                    "3' %s" % com(edcut).replace('^', ''),
                ]
                fixed_data[ed['enzyme']] = ed
                gb['+'] += 1
            else:
                gb['-'] += 1
    print gb

    with open('rebase.yaml', 'w') as handle:
        handle.write(yaml.dump(fixed_data, default_flow_style=False))

if __name__ == '__main__':
    parse_data(fetchfile())
