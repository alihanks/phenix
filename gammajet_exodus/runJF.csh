#!/bin/tcsh
echo ' '
echo 'START: '`date`
cd /direct/phenix+hhj/jfrantz/gammajet_exodus/
hostname 
/phenix/hhj/jfrantz/gammajet_exodus/exodus<<EOF 
6 
20000
/phenix/hhj/jfrantz/gammajet_exodus/output/seccent0_II_20k.root
/phenix/hhj/jfrantz/gammajet_exodus/output/cent0_all_0.txt
88 
6 
EOF 
echo ' ' 
echo 'END: '`date` 
echo ' ' 
