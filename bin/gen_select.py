#!/usr/bin/env python
import yaml
import sys

with open(sys.argv[1], 'r') as handle:
    data = yaml.load(handle)

for x in data:
    for q in ([data[x]['enzyme']] + data[x]['isoscizomers']):
    #if len(data[x]['isoscizomers']) > 0:
        #q += " [%s]" % ' '.join(data[x]['isoscizomers'])
        print '<option value="%s">%s</option>' % (q, q)
