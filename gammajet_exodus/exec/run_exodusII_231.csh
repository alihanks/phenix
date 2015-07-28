#!/bin/tcsh
echo ' '
echo 'START: '`date`
cd /phenix/hhj/ahanks/gammajet_exodus
hostname 
/phenix/hhj/ahanks/gammajet_exodus/exodus<<EOF 
6 
25000 
/phenix/hhj/ahanks/gammajet_exodus/fe_sharks/cent0_II_231.root
/phenix/hhj/ahanks/gammajet_exodus/fe_sharks/cent0_all_231.txt
88 
6 
EOF 
echo ' ' 
echo 'END: '`date` 
echo ' ' 
