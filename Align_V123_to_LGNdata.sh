#v3
surfIDs=(161108PSY 161108K 170314JBH 170225HJH 170703KB 161205SHY 180409LSY 180724HJH 180907CES 181018YYH 181210JHY 190305KYJ )
SNs=(01 03 04 05 06 07 08 09 10 11 13 14)
retIDs=(05 04 47 44 69 08 98 123 128 132 144 153)
cpROIs=( lh.V1 lh.V2d lh.V3d lh.V2v lh.V3v lh.V4v rh.V1 rh.V2d rh.V3d rh.V2v rh.V3v rh.V4v )

#RSVP

surfIDs=(161108PSY 161205SHY 161108K 161027C_Tra_Local 170222JYJ 170110KDH 170225HJH 170314JBH 170106LSC)
SNs=(01 02 03 05 06 07 08 09 10)
retIDs=(05 08 04 07 41 25 44 47 21)
cpROIs=( lh.V1 lh.V2d lh.V3d lh.V2v lh.V3v lh.V4v rh.V1 rh.V2d rh.V3d rh.V2v rh.V3v rh.V4v )


#Ecc
surfIDs=(170628JYS 170515PJH 161108PSY 170222JYJ 170110KDH 170912LYR 161108K 171009SSK 171120HSH 170804HJY 180326BKD)
SNs=(05 06 07 08 09 10 12 13 14 15 16)
retIDs=(65 58 05 41 25 75 04 78 82 73 95)
cpROIs=( lh.V1 lh.V2d lh.V3d lh.V2v lh.V3v lh.V4v rh.V1 rh.V2d rh.V3d rh.V2v rh.V3v rh.V4v )

rm -rf *h.V4v.1D.roi
mv lh.V4v_new.1D.roi lh.V4v.1D.roi
mv rh.V4v_new.1D.roi rh.V4v.1D.roi
exp=v3
exp=Ecc
exp=RSVP

SURF_DIR=/sas2/PECON/freesurfer/subjects
LOC_DIR=/sas2/PECON/HJY/CrM/${SN}/LocalizerAnalysis/
RET_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/Img_data
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh


for ((i=0;i<=10;i++));
do
surfID=${surfIDs[$i]}
SN=${SNs[$i]}

SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data
LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results
master=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results/final_epi_external_volreg_base+orig
cd ${SN_DIR}

rm -f SVol@*
#align surface volume to T1@main (SVol@main)
@SUMA_AlignToExperiment -exp_anat ${SN}_loc.results/T1@main+orig. \
	-surf_anat $SURF_DIR/${surfID}/SUMA/${surfID}_SurfVol_SS+orig. \
	-wd \
	-prefix SVol@main -align_centers #sometimes it's better to leave -wd option out 

done
afni SVol@main+orig. T1@main+orig.

#---------------Copy ROI from Retinotopy data ------------------# 

for ((k=0;k<=13;k++));
do
surfID=${surfIDs[$k]}
retID=${retIDs[$k]}
SN=${SNs[$k]}

thresmask=thresmask_${exp}_q.001+orig
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data
LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results
master=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results/final_epi_external_volreg_base+orig
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh
cd ${LOC_DIR}

mkdir roi
cd roi

#cp roi from Retinotopy data or draw on your own
for i in "${cpROIs[@]}"
do
cp ${ROI_DIR}/${i}.1D.roi .
done

#run ROI2Vol.sh lh and rh -> this process aligns ROIs to our basegrid! 
cp /home/jiyeongha/Script/ROI2Vol_LGN.sh .; 
./ROI2Vol_LGN.sh lh ${surfID} ${retID} ${SN}
./ROI2Vol_LGN.sh rh ${surfID} ${retID} ${SN}

ROIs=( V1 V4v center1dg center3dg )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r lh.${r}+orig. rh.${r}+orig.
done

ROIs=( V2 V3 )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix ${r} lh.${r}d+orig. rh.${r}d+orig. lh.${r}v+orig. rh.${r}v+orig.
done

# V1 & V2 -> exclude from V1 mask
mv V1+orig.HEAD tmp.V1+orig.HEAD; mv V1+orig.BRIK tmp.V1+orig.BRIK
3dcalc -a tmp.V1+orig -b V2+orig -expr 'a+b' -prefix tmp.V1+V2+orig
3dcalc -a tmp.V1+V2+orig -expr 'ispositive(a-1.5)' -prefix mask.V1+V2+orig
3dcalc -a tmp.V1+orig -b mask.V1+V2+orig -expr 'ispositive(a-b)' \
	-prefix V1+orig

# V2 & V3 -> exclude from V2 mask
mv V2+orig.HEAD tmp.V2+orig.HEAD
mv V2+orig.BRIK tmp.V2+orig.BRIK
3dcalc -a tmp.V2+orig -b V3+orig -expr 'a+b' -prefix tmp.V2+V3+orig
3dcalc -a tmp.V2+V3+orig -expr 'ispositive(a-1.5)' -prefix mask.V2+V3+orig
3dcalc -a tmp.V2+orig -b mask.V2+V3+orig -expr 'ispositive(a-b)' \
	-prefix V2+orig

# V3 & V4(V4v+VO) -> exclude from V3 mask
mv V3+orig.HEAD tmp.V3+orig.HEAD
mv V3+orig.BRIK tmp.V3+orig.BRIK
3dcalc -a tmp.V3+orig -b V4v+orig -expr 'a+b' -prefix tmp.V3+V4v+orig
3dcalc -a tmp.V3+V4v+orig -expr 'ispositive(a-1.5)' -prefix mask.V3+V4v+orig
3dcalc -a tmp.V3+orig -b mask.V3+V4v+orig -expr 'ispositive(a-b)' \
	-prefix V3+orig


rm -f tmp.* mask.V*


#---------------ThresMask X ROI (Final Masks) ------------------# 

ROIs=( V1 V2 V3 V4v )
# threshold with functional contrast map
for i in "${ROIs[@]}"
do
#3dresample -dxyz 2.0 2.0 2.0 -prefix rsmp_${r} -input ${r}+orig. #->match ROI voxel size - thesmask voxel size 
3dcalc -a ${i}+orig.HEAD -b ../$thresmask -expr 'step(a*b)' \
-prefix ${i}_masked_q.001 -overwrite #
3dcalc -a ${i}_masked_q.001+orig.HEAD -b ../full_mask.${SN}_loc+orig.HEAD \
-expr 'step(a*b)' -prefix  ${SN_DIR}/forwardmodel/roi/${i}_fmasked_q.001.nii -overwrite
done

rm -f *_masked+orig.*
done

#---------------fmasked X center3dg, center1dg  ------------------# 

ROIs=( V1 V2 V3 V4v )
# threshold with functional contrast map
for i in "${ROIs[@]}"
do
3dcalc -a ${SN_DIR}/forwardmodel/roi/${i}_fmasked_q.001.nii -b center3dg+orig. -expr 'step(a*b)' \
-prefix ${SN_DIR}/forwardmodel/roi/${i}_3dg_fmasked_q.001 -overwrite 
3dcalc -a ${SN_DIR}/forwardmodel/roi/${i}_fmasked_q.001.nii -b center1dg+orig. -expr 'step(a*b)' \
-prefix ${SN_DIR}/forwardmodel/roi/${i}_1dg_fmasked_q.001 -overwrite 
3dcalc -a ${SN_DIR}/forwardmodel/roi/${i}_fmasked_q.001.nii -b center3dg+orig. -c center1dg+orig. -expr 'step(a*(step(b)-ispositive(b*c)))' \
-prefix ${SN_DIR}/forwardmodel/roi/${i}_3-1dg_fmasked_q.001 -overwrite 
done


#---------------fmasked X center3dg, center1dg  ------------------# 

LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data/${SN}_loc.results
cd ${LOC_DIR}
3dcalc -a thresmask_v3_q.001+orig -b ../forwardmodel/roi/center3dg.nii -expr 'step(a*b)' \
-prefix ../forwardmodel/roi/center3dg_fmasked_q.001.nii -overwrite 
3dcalc -a thresmask_v3_q.001+orig  -b ../forwardmodel/roi/center1dg.nii -expr 'step(a*b)' \
-prefix ../forwardmodel/roi/center1dg_fmasked_q.001 -overwrite 



#---------------ThresMask X ROI (Final Masks) ------------------# 


for ((k=0;k<=12;k++));
do
surfID=${surfIDs[$k]}
retID=${retIDs[$k]}
SN=${SNs[$k]}
exp=v3
thresmask=thresmask_${exp}_c-a_p.05+orig
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data
LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results
master=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results/pb01.${SN}_loc.r01.volreg+orig
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh
cd ${LOC_DIR}/roi

ROIs=( V1 V2 V3 V4v )
for i in "${ROIs[@]}"
do
#3dresample -dxyz 2.0 2.0 2.0 -prefix rsmp_${r} -input ${r}+orig. #->match ROI voxel size - thesmask voxel size 

3dcalc -a ${i}+orig.HEAD -b ${SN_DIR}/forwardmodel/roi/${i}_fmasked.nii \
-expr 'ispositive(a-b)' -prefix  tmp.${i}_notoverlapped.nii
3dcalc -a tmp.${i}_notoverlapped.nii -b ../full_mask.${SN}_loc+orig.HEAD \
-expr 'ispositive(a*b)' -prefix  ${SN_DIR}/forwardmodel/roi/${i}_notoverlapped.nii
done
done


ROIs=( V2 V3 )
for i in "${ROIs[@]}"
do
#3dresample -dxyz 2.0 2.0 2.0 -prefix rsmp_${r} -input ${r}+orig. #->match ROI voxel size - thesmask voxel size 
3dmerge -gmax -prefix ${i}d.nii lh.${i}d+orig. rh.${i}d+orig.
3dcalc -a ${i}d.nii -b ${SN_DIR}/forwardmodel/roi/${i}_fmasked.nii \
-expr 'ispositive(a-b)' -prefix  tmp.${i}d.nii
3dcalc -a ${i}d.nii -b ${SN_DIR}/forwardmodel/roi/${i}_fmasked.nii \
-expr 'ispositive(a*b)' -prefix  ${SN_DIR}/forwardmodel/roi/${i}d_fmasked.nii
3dcalc -a tmp.${i}d.nii -b ../full_mask.${SN}_loc+orig.HEAD \
-expr 'ispositive(a*b)' -prefix  ${SN_DIR}/forwardmodel/roi/${i}d_notoverlapped.nii
done


done

#---------------ThresMask X ROI (Final Masks) ------------------# 


for ((k=0;k<=11;k++));
do
surfID=${surfIDs[$k]}
retID=${retIDs[$k]}
SN=${SNs[$k]}
exp=RSVP
thresmask=thresmask_${exp}_q.001
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data
LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results
master=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results/pb01.${SN}_loc.r01.volreg+orig
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh
cd ${LOC_DIR}/roi

ROIs=( V1 V2 V3 V4v )
for r in "${ROIs[@]}"
do
3dcalc -a ${r}+orig. -b ../${thresmask}+orig. \
-expr 'ispositive(a*b)' -prefix ${r}_masked -overwrite

3dcalc -a ${r}_masked+orig. -b ../full_mask.${SN}_loc+orig. \
-expr 'ispositive(a*b)' -prefix ${SN_DIR}/forwardmodel/roi/${r}_fmasked_q.001.nii -overwrite
done
done


#---------------ThresMask X ROI (Final Masks) ------------------# 


for ((k=0;k<=14;k++));
do
surfID=${surfIDs[$k]}
retID=${retIDs[$k]}
SN=${SNs[$k]}
exp=Ecc
thresmask=thresmask_${exp}_q.001
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data
LOC_DIR=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results
master=/group_hpc/WMShimLab2/PSY_Color/Color${exp}/${SN}/Img_data/${SN}_loc.results/pb01.${SN}_loc.r01.volreg+orig
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh
cd ${LOC_DIR}/roi

ROIs=( V1 V2 V3 V4v )
for r in "${ROIs[@]}"
do
3dcalc -a ${r}+orig. -b ../${thresmask}+orig. \
-expr 'ispositive(a*b)' -prefix ${r}_masked -overwrite

3dcalc -a ${r}_masked+orig. -b ../full_mask.${SN}_loc+orig. \
-expr 'ispositive(a*b)' -prefix ${SN_DIR}/forwardmodel/roi/${r}_fmasked_q.001.nii -overwrite
done
done








