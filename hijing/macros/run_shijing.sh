#! /bin/sh
date

IT=$1
NEVT=$2
CONE=$3

echo "Running job ${IT} with ${NEVT} and cone cut = ${CONE}"
cd /phenix/scratch/ahanks/run8/shijing/

rm -rf ${IT}
mkdir ${IT}
cd ${IT}

SEED=`expr $IT + 10040`
echo "seed set to ${SEED}"
# now I have to create the xml file in a new directory
echo "<?xml version=\"1.0\"?>" > sHijing.xml
echo "<HIJING>" >> sHijing.xml
echo " <EFRM>200.0</EFRM>" >> sHijing.xml
echo " <FRAME>CMS</FRAME>" >> sHijing.xml
echo " <PROJ>A</PROJ>" >> sHijing.xml
echo " <TARG>A</TARG>" >> sHijing.xml
echo " <IAP>2</IAP>" >> sHijing.xml
echo " <IZP>1</IZP>" >> sHijing.xml
echo " <IAT>197</IAT>" >> sHijing.xml
echo " <IZT>79</IZT>" >> sHijing.xml
echo " <N>${NEVT}</N>" >> sHijing.xml
echo " <BMIN>0.0</BMIN>" >> sHijing.xml
echo " <BMAX>8.0</BMAX>" >> sHijing.xml
echo " <KEEP_SPECTATORS>1</KEEP_SPECTATORS>" >> sHijing.xml
echo " <OUTPUT>sHijing${IT}_${CONE}.dat</OUTPUT>" >> sHijing.xml
echo " <RANDOM>" >> sHijing.xml
echo "  <SEED>${SEED}</SEED>" >> sHijing.xml
echo " </RANDOM>" >> sHijing.xml
echo " <HIPR1>" >> sHijing.xml
echo "  <5>2</5>" >> sHijing.xml
echo "  <10>-10.0</10>" >> sHijing.xml
echo " </HIPR1>" >> sHijing.xml
echo " <IHPR2>" >> sHijing.xml
echo "  <2>3</2>" >> sHijing.xml
# 1 for standard hard process, 2 to force a direct photon in the final state
echo "  <3>1</3>" >> sHijing.xml
#echo "  <3>2</3>" >> sHijing.xml
# 1 for only direct photons, 0 for no decay photons but pi0s
echo "  <12>0</12>" >> sHijing.xml
#echo "  <12>1</12>" >> sHijing.xml
echo "  <21>1</21>" >> sHijing.xml
echo " </IHPR2>" >> sHijing.xml
echo "</HIJING>" >> sHijing.xml

sHijing


output=/phenix/hhj/ahanks/shijing/run8/output

date

cd /phenix/hhj/ahanks/shijing/run8
root.exe -b <<EOF
.x /phenix/u/workarea/ahanks/run8/hijing/macros/shijing_ana_macro.C("/phenix/scratch/ahanks/run8/shijing/${IT}/sHijing${IT}_${CONE}.dat","${output}/hijing_analysis_pi0_${IT}_${CONE}.root",0,${CONE});
.q
EOF

date

rm -rf /phenix/scratch/ahanks/run8/shijing/${IT}
