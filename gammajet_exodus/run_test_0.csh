#!/bin/tcsh
echo ' '
echo 'START: '`date`
cd /direct/phenix+u/workarea/mjuszkie/gammajet_exodus/
hostname 
/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF 
6 
20000 
/direct/phenix+hp/data10/mjuszkie/exodus/cent0/new/cent0_test_0.root
/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_all_0.txt
88 
6 
EOF 
echo ' ' 
echo 'END: '`date` 
echo ' ' 
