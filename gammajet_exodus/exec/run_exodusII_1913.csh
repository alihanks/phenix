#!/bin/tcsh
echo ' '
echo 'START: '`date`
cd /phenix/hhj/ahanks/gammajet_exodus
hostname 
/phenix/hhj/ahanks/gammajet_exodus/exodus<<EOF 
6 
25000 
/phenix/hhj/ahanks/gammajet_exodus/fe_sharks/cent0_II_1913.root
/phenix/hhj/ahanks/gammajet_exodus/fe_sharks/cent0_all_1913.txt
88 
6 
EOF 
echo ' ' 
echo 'END: '`date` 
echo ' ' 
