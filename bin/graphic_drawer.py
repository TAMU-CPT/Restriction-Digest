#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse
import dnadigest
import math
import logging
from pkg_resources import resource_stream
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

html_tpl = u"""
<!DOCTYPE html>
<html>
    <head>
        <script type="text/javascript">
        {}
        </script>
    </head>
    <body>
        {}
        {}
    </body>
</html>
"""

style_static = u"""
<style>
    body {font-family: 'Lato';font-weight:400;}
    .boundary {stroke-dasharray:2,2;stroke-width:2px}
    .mdlabel {font-size:12px}
    .smlabel {font-size:8px}
</style>
"""

region_tpl = u"""
<trackmarker start="{start}" end="{end}" markerstyle="fill:{color}" arrowendlength="{strand}4" arrowstartlength="{notstrand}4">
    <markerlabel type="path" class="mdlabel white" text="{label}"></markerlabel>
</trackmarker>

<trackmarker start="{start}" markerstyle="stroke:{color}" class="boundary" wadjust="20">
    <markerlabel class="smlabel" text="{start}" vadjust="30"></markerlabel>
</trackmarker>
<trackmarker start="{end}" markerstyle="stroke:{color}" class="boundary" wadjust="20">
    <markerlabel class="smlabel" text="{end}" vadjust="30"></markerlabel>
</trackmarker>
"""

cut_site_tpl = u"""
<trackmarker start="{position}" markerstyle="stroke:rgba(0, 0, 0, 0.7)" wadjust="20" vadjust="-6" class="marker">
    <markerlabel text="{label}" labelclass="smlabel" vadjust="30"></markerlabel>
</trackmarker>
"""

plasmid_tpl = u"""
<plasmid sequencelength="{sequence_length}" plasmidheight="{size_full}" plasmidwidth="{size_full}">
    <plasmidtrack trackstyle="fill:#ccc" width="5" radius="{size_half}"></plasmidtrack>

    <plasmidtrack trackstyle="fill:rgba(225,225,225,0.5)" radius="{size_half_inner}">
        <tracklabel text="{sequence_id}" labelstyle='font-size:20px;font-weight:400'></tracklabel>
        <tracklabel text="{sequence_length} bp" labelstyle='font-size:10px' vadjust="20"></tracklabel>

        <!-- draw the main markers and labels -->
        {regions}

        <!-- draw the scales -->
        <trackscale interval="{minor_tick_interval}" style='stroke:#999' direction="in" ticksize="3"></trackscale>
        <trackscale interval="{minor_tick_interval}" style='stroke:#999' ticksize="3"></trackscale>
        <trackscale interval="{major_tick_interval}" style="stroke:#f00" direction="in" showlabels="1" labelstyle="fill:#999;stroke:none;text-anchor:middle;alignment-baseline:middle;font-size:10px"></trackscale>
    </plasmidtrack>
</plasmid>
"""


def BestTick(largest, mostticks):
    # http://stackoverflow.com/questions/361681/algorithm-for-nice-grid-line-intervals-on-a-graph/361687#361687
    minimum = largest / mostticks
    magnitude = 10 ** math.floor(math.log(minimum) / math.log(10))
    residual = minimum / magnitude
    if residual > 5:
        tick = 10 * magnitude
    elif residual > 2:
        tick = 5 * magnitude
    elif residual > 1:
        tick = 2 * magnitude
    else:
        tick = magnitude
    return tick

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restriction Digest Tool',
                                     epilog="")
    # Input filename, required
    parser.add_argument('file', help='Input fasta genome(s)')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--data', help='Enzyme cut site dataset')
    # A single string with default value of 'enzyme_data.yaml'
    parser.add_argument('--regions', type=file, help='Regions to label (must be bed6+)')
    # A list of one or more strings, at the end
    parser.add_argument('enzyme', help='Comma separated list of enzymes')
    args = parser.parse_args()

    annotated_regions = {}
    if args.regions:
        for line in args.regions:
            data = line.strip().split('\t')
            (chrom, start, end, name, score, strand) = data[0:6]
            if chrom not in annotated_regions:
                annotated_regions[chrom] = []
            annotated_regions[chrom].append(data[1:6])

    dd = dnadigest.DnaDigest(data_file=args.data)
    enzymes = args.enzyme.split(',')

    plasmids = []
    for record in SeqIO.parse(args.file, 'fasta'):
        fragments, cut_sites, did_cut = dd.digest_sequence(record, enzymes)
        regions = []
        for annotated_region in annotated_regions.get(record.id, []):
            regions.append(region_tpl.format(
                start=annotated_region[0],
                end=annotated_region[1],
                label=annotated_region[2],
                color='rgba(255, 0, 0, %s)' % (float(annotated_region[3]) / 1000),
                strand='' if annotated_region[4] == '+' else '-',
                notstrand='' if annotated_region[4] != '+' else '-',
            ))

        for site in cut_sites:
            regions.append(cut_site_tpl.format(position=site, label=', '.join(cut_sites[site])))

        size = 500
        major_tick_interval = BestTick(len(record.seq), 20)
        minor_tick_interval = major_tick_interval / 5
        plasmids.append(
            plasmid_tpl.format(
                sequence_length=len(record.seq),
                sequence_id=record.id,
                regions='\n'.join(regions),
                size_full=size,
                size_half=.40 * size,
                size_half_inner=.40 * size - 10,
                major_tick_interval=major_tick_interval,
                minor_tick_interval=minor_tick_interval,
            )
        )

    data = html_tpl.format(
        resource_stream('dnadigest', 'angularplasmid.complete.min.js').read().decode("utf-8"),
        style_static,
        u'\n'.join(plasmids)
    )
    print data.encode('utf-8')
