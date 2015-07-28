
#!/bin/tcsh
echo ' '
echo 'START: '`date`
echo ' '
hostname
/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF
3
1000000
/phenix/scratch/mjuszkie/gammajet.exodus/test_sharkfin1.root
46.671
890.7
1.0
5
EOF
echo ' '
echo 'END: '`date`
echo ' '
