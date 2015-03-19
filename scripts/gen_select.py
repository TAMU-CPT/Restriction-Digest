#!/usr/bin/env python
import yaml
import sys

with open(sys.argv[1], 'r') as handle:
    data = yaml.load(handle)

for x in data:
    if len(data[x]['isoscizomers']) > 0:
        print '<option value="%s">%s [%s]</option>' % (data[x]['enzyme'], data[x]['enzyme'], ', '.join( data[x]['isoscizomers']))
    else:
        print '<option value="%s">%s</option>' % (data[x]['enzyme'], data[x]['enzyme'])
