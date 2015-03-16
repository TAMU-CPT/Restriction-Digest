#!/usr/bin/env python
import yaml
import sys

with open(sys.argv[1], 'r') as handle:
    data = yaml.load(handle)

full_list = []
for x in data:
    for q in ([data[x]['enzyme']] + data[x]['isoscizomers']):
        if q not in full_list:
            full_list.append(q)

for q in sorted(full_list):
    print '<option value="%s">%s</option>' % (q, q)
