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
    <markerlabel type="path" class="smlabel" text="{label}"></markerlabel>
</trackmarker>

<!--
<trackmarker start="{start}" markerstyle="stroke:rgba(0,0,0,0.3)" class="boundary" wadjust="10">
    <markerlabel class="smlabel" text="{start}" vadjust="10"></markerlabel>
</trackmarker>
<trackmarker start="{end}" markerstyle="stroke:rgba(0,0,0,0.3)" class="boundary" wadjust="10">
    <markerlabel class="smlabel" text="{end}" vadjust="10"></markerlabel>
</trackmarker>
-->
"""

cut_site_tpl = u"""
<trackmarker start="{position}" markerstyle="stroke:rgba(0, 0, 0, 0.3);" wadjust="{rad}" vadjust="-6" class="marker">
    <markerlabel text="{label}" labelclass="smlabel"
    vadjust="{rad_label}"></markerlabel>
</trackmarker>
"""
# labelstyle="
# -webkit-transform: rotate({rotation}deg);
# -moz-transform: rotate({rotation}deg);
# -ms-transform: rotate({rotation}deg);
# -o-transform: rotate({rotation}deg);
# -webkit-transform-origin: {pos_x}px {pos_y}px;
# -moz-transform-origin: {pos_x}px {pos_y}px;
# -ms-transform-origin: {pos_x}px {pos_y}px;
# -o-transform-origin: {pos_x}px {pos_y}px;
# "

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

    enzyme_distance_adjustment_map = {enzyme: 20 + 15 * idx for idx, enzyme in enumerate(enzymes)}
    length_adjustments = {}

    plasmids = []
    for record in SeqIO.parse(args.file, 'fasta'):
        cut_sites = dd.multidigest_sites(record, enzymes)
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

        import pprint
        log.debug(pprint.pformat(cut_sites))

        cut_site_keys = sorted(cut_sites.keys())
        for idx, site in enumerate(cut_site_keys):
            prev_site = cut_site_keys[idx-1] if idx > 1 else None
            # next_site = cut_sites[idx+1] if idx < len(cut_sites) else None
            local_length_adjustment = 20
            if prev_site is not None and site - prev_site < 20:
                if prev_site in length_adjustments:
                    local_length_adjustment = length_adjustments[prev_site] + 15
            length_adjustments[site] = local_length_adjustment

            for cut_enz in cut_sites[site]:
                regions.append(cut_site_tpl.format(
                    position=site,
                    label=cut_enz,
                    rad=local_length_adjustment + 10,
                    rad_label=local_length_adjustment + 15,
                ))

        size = 500
        major_tick_interval = BestTick(len(record.seq), 20)
        minor_tick_interval = major_tick_interval / 5
        radial_scaling = .35
        plasmids.append(
            plasmid_tpl.format(
                sequence_length=len(record.seq),
                sequence_id=record.id,
                regions='\n'.join(regions),
                size_full=size,
                size_half=radial_scaling * size,
                size_half_inner=radial_scaling * size - 10,
                major_tick_interval=major_tick_interval,
                minor_tick_interval=minor_tick_interval,
            )
        )

    print html_tpl.format(
        resource_stream('dnadigest', 'angularplasmid.complete.min.js').read().decode("utf-8"),
        style_static,
        u'\n'.join(plasmids)
    ).encode('utf-8')
