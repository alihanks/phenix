#!/bin/tcsh
echo ' '
echo 'START: '`date`
echo ' '
hostname
/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF
6
100000000
/phenix/scratch/mjuszkie/gammajet.exodus/lowpt100M.root
/phenix/scratch/mjuszkie/gammajet.exodus/lowpt100M.txt
500
1
EOF
echo ' '
echo 'END: '`date`
echo ' '
