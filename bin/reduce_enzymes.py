#!/usr/bin/env python
import yaml
import sys

if len(sys.argv) != 3:
    print "Run with %s input.yaml output.yaml" % sys.argv[0]
    sys.exit(1)

with open(sys.argv[1], 'r') as handle:
    data = yaml.load(handle)

recog_site_map = {}

for e in data:
    en = e['enzyme']
    key = ':'.join(e['cut'])
    if key in recog_site_map:
        # Clean out empty item in list
        iz = recog_site_map[key]['isoscizomers']
        if len(iz) == 1 and iz[0] == '':
            recog_site_map[key]['isoscizomers'] = []

        if en not in recog_site_map[key]['isoscizomers']:
            recog_site_map[key]['isoscizomers'].append(en)
    else:
        recog_site_map[key] = e

        # Clean out empty item in list
        iz = recog_site_map[key]['isoscizomers']
        if len(iz) == 1 and iz[0] == '':
            recog_site_map[key]['isoscizomers'] = []

print "Reduced from %s to %s" % (len(data), len(recog_site_map.keys()))

restructured = {}
for e in recog_site_map:
    restructured[recog_site_map[e]['enzyme']] = recog_site_map[e]

with open(sys.argv[2], 'w') as handle:
    yaml.dump(restructured, handle)
