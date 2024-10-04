magfile=$1
time_step=$2
nx=$3
ny=$4
nt=$5
amlife=$6
act=$7
det_radius=$8
det_height=$9
xfov=${10}
yfov=${11}
act_file=${12}
ndetectors=${13}
a2=${14}
if_reco=1
gyro=${15}
B0=${16}
gammaMRI_exe=gammaMRI_det_v3_opt.x
#gammaMRI_exe=gammaMRI_det_v2_opt.x
rm list_gammaMRI.raw
echo ./$gammaMRI_exe $magfile $time_step $nx $ny $nt $amlife $act $det_radius $det_height $xfov $yfov $act_file $ndetectors $a2 $if_reco $gyro $B0
./$gammaMRI_exe $magfile $time_step $nx $ny $nt $amlife $act $det_radius $det_height $xfov $yfov $act_file $ndetectors $a2 $if_reco $gyro $B0
ls -lrht list_gammaMRI.raw
