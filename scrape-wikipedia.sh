#!/bin/bash
for i in A C-D E-F G-K L-N O-R S T-Z Ba-Bc Bd-Bp Bsa-Bso Bsp-Bss Bst-Bv
do
   wget  http://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_$i  -O $i.html;
 sleep 1;
done;
