#! /bin/sh                                                                                               
# directory

if [ $# -eq 4 ]
then

NEVT=$1
NIT=$2
OFFSET=$3
CONE=$4

echo "submitting $NIT jobs ..."

WORKDIR=/phenix/u/workarea/ahanks/run8/hijing/macros
OUTDIR=/phenix/hhj/ahanks/shijing/run8

IT=$OFFSET
MAX=$(( $NIT+$OFFSET ))
echo "running over iteration $IT to $MAX"

cd ${OUTDIR}
echo "Universe       = vanilla" > run_hijing.job
echo "Notification   = Error" >> run_hijing.job
echo "Notify_user     = (USER)@rcf2.rhic.bnl.gov" >> run_hijing.job
echo "GetEnv         = True" >> run_hijing.job
echo "Initialdir     = ${WORKDIR}" >> run_hijing.job
echo "Executable     = ${WORKDIR}/run_shijing.sh" >> run_hijing.job

while [ $IT -lt $MAX ]; do
 
echo "Arguments      = ${IT} ${NEVT} ${CONE}" >> run_hijing.job
#echo "Arguments      = /direct/phenix+scratch/sahlmul/dAuMC/\$(Process).csh" >> run_hijing.job
echo "Log            = ${OUTDIR}/logfiles/log_hijing_${IT}_${CONE}.log" >> run_hijing.job
echo "Output         = ${OUTDIR}/logfiles/run_hijing_${IT}_${CONE}.out" >> run_hijing.job
echo "Error          = ${OUTDIR}/logfiles/run_hijing_${IT}_${CONE}.out" >> run_hijing.job
echo "+Experiment    = \"phenix\"" >> run_hijing.job
echo "+Job_Type      = \"cas\"" >> run_hijing.job
echo "Queue" >> run_hijing.job
IT=`expr $IT + 1`

done

condor_submit run_hijing.job

else

echo " Usage: NEVT(per-job) NJOBS FIRST_RUN# CONE_SIZE"

fi
