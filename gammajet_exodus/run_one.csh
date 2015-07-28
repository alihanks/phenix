#!/bin/tcsh
echo ' '
echo 'START: '`date`
echo ' '
hostname
/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF
6
100
/phenix/scratch/mjuszkie/gammajet.exodus/testI.root
/phenix/scratch/mjuszkie/gammajet.exodus/testI.txt
500
1
EOF
echo ' '
echo 'END: '`date`
echo ' '
