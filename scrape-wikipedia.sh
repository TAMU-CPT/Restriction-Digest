#!/bin/bash
for i in {'A','C','D','E','F','G','H','I','K','L','M','N','O','P','R','S','T','U','V','X','Y','Z','Ba-Bc','Bd-Bp','Bsa-Bso','Bsp','Bsr-Bss','Bst','Bsu-Bv'};
do
    wget http://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:$i -O $i.html;
    sleep 1;
done;
