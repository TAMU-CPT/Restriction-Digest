Restriction Digest Tool
=======================

Galaxy tool with the following aims:

-  digest DNA into fragments, produce FASTA output of that
-  plot restriction digests in easy to use SVG format
-  run "hypothetical gels" to help biologists understand the output of
   these fancy computer things
-  optimise restrictions with multiple genomes in order to easily
   determine enzymes for unique resetriction patterns

Testing
-------

::

    python setup.py test

Running
-------

::

    python bin/tool_wip.py  test.fa Bsp6I Bsp13I DseDI

