Restriction Digest Tool
=======================

Python/Galaxy tool with the following aims:

-  digest DNA into fragments, produce FASTA output of that
-  plot restriction digests in easy to use SVG format
-  run "hypothetical gels" to help biologists understand the output of
   these fancy computer things
-  optimise restrictions with multiple genomes in order to easily
   determine enzymes for unique resetriction patterns

Running
-------

::

    python bin/digest_dna.py test.fa Bsp6I,Bsp13I,DseDI

Authors
-------

* Stephen Crosby
* Eric Rasche

Acknowledgements
----------------

-  Thanks to NEB's REBASE for providing the data backing
   `dnadigest/rebase.yaml <dnadigest/rebase.yaml>`__.
-  Thanks to Wikipedia for their enzyme data backing
   `dnadigest/wikipedia.yaml <dnadigest/wikiepdia.yaml>`__
