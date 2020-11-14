###Align everything to T1@main
###T1@main is anatomy aligned to main epi data. (from ?_main.results)
###

#0429 2020
#Align LGN_?_hk2 roi to main
#LGN_?_hk2 is a bigger roi than ?_hk

SN=14
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LGN_PD_DIR=${SN_DIR}/LGN_PD

PD_id=190314KYJ
PD_DIR=/group_hpc/WMShimLab/PD/${PD_id}

#make a directory for HK2 roi 
mkdir -p ${LGN_PD_DIR}/HK2
chmod -R 777 ${LGN_PD_DIR}/HK2

cd ${LGN_PD_DIR}/HK2

#copy PD_mean_al (HK2 roi background), T1_PD (which PD_mean_al is aligned to) & HK2
#make sure that LGN_?_hk2+orig. is drawn on PD_mean_al!
3dcopy ${PD_DIR}/T1_SS+orig. ${LGN_PD_DIR}/HK2/T1_PD+orig.
3dcopy ${PD_DIR}/PD_mean_al+orig. ${LGN_PD_DIR}/HK2/meanPD@T1_PD
cp ${PD_DIR}/LGN_*_hk2+orig.* ${LGN_PD_DIR}/HK2/

#copy T1@main here
3dcopy ../../T1@main+orig. ${LGN_PD_DIR}/HK2/ -overwrite

#align T1_PD @ T1@main 
align_epi_anat.py -dset1 meanPD+orig -dset2 T1@main+orig -dset1to2 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no -overwrite

#merge LGN 
3dmerge -gmax -prefix LGN_hk2 LGN_l_hk2+orig. LGN_r_hk2+orig.

#align LGN @ T1@main using T1_PD_al_mat 
3dAllineate -cubic -1Dmatrix_apply meanPD_al_mat.aff12.1D \
-master ../../${SN}_main.results/pb01.${SN}_main.r01.volreg+orig \
-prefix LGN_hk@main+orig.HEAD \
LGN_hk2+orig.HEAD -overwrite

cd ${SN_DIR}


#LGN_hk2@main X LGN_thresmask
run=( 1 2 3 4 )
for r in "${run[@]}"
do
@Align_Centers -no_cp -base T1@main+orig. -dset loc_${r}.nii
done
#align epi to T1 (LGN)
subj=${SN}_loc
afni_proc.py -subj_id ${subj}     \
-dsets loc_?.nii \
-blocks align volreg mask regress 	  \
-volreg_base_dset loc_1.nii['0']      \
-volreg_align_e2a                         \
-align_opts_aea -giant_move               \
-copy_anat T1@main+orig.     	  \
-regress_censor_motion 0.5                \
-regress_censor_outliers 0.1 		  

tcsh -xef proc.${subj} |& tee output.proc.${subj}
cd ${subj}.results/
run=( 01 02 03 04 )
for r in "${run[@]}"
do
3dClipLevel pb01.${subj}.r${r}.volreg+orig.HEAD >> clip.txt
3dTstat -mean -prefix r.${r}.base pb01.${subj}.r${r}.volreg+orig.HEAD'[0..$]'

done

clip=$(sort -n clip.txt | sed -n '1p')
run=( 01 02 03 04 )
for r in "${run[@]}"
do
3dcalc -a pb01.${subj}.r${r}.volreg+orig. -b r.${r}.base+orig. \
-expr "(100 * a/b) * step(b-$clip)" -prefix pb01.${subj}.r${r}.scaled
3dTstat -mean -prefix r.${r}.sc_base pb01.${subj}.r${r}.scaled+orig.HEAD'[0..$]'

# detrend linear trend
3dDetrend -polort 1 -prefix pb01.${subj}.r${r}.sc_dt pb01.${subj}.r${r}.scaled+orig
# add mean after detrend
3dcalc -a pb01.${subj}.r${r}.sc_dt+orig.HEAD -b r.${r}.sc_base+orig.HEAD \
-expr 'a+b' -prefix pb01.${subj}.r${r}.sc_dt_am

# calculate mean image
3dTstat -mean -prefix r.${r}.sc_dt_base pb01.${subj}.r${r}.sc_dt_am+orig.HEAD'[0..$]'

# high-pass filter
3dBandpass -prefix pb01.${subj}.r${r}.sc_dt_hp 0.01 99999 pb01.${subj}.r${r}.sc_dt_am+orig
# add mean after hp filter
3dcalc -a pb01.${subj}.r${r}.sc_dt_hp+orig.HEAD -b r.${r}.sc_dt_base+orig.HEAD \
-expr 'a+b' -prefix pb01.${subj}.r${r}.sc_dt_hp_am

#blur
3dmerge -1blur_fwhm 3 -doall -prefix rm.pb01.$subj.r${r}.sc_dt_hp_am_blur    \
pb01.${subj}.r${r}.sc_dt_hp_am+orig

# and set boundaries using anat mask 
3dcalc -a rm.pb01.${subj}.r${r}.sc_dt_hp_am_blur+orig -b full_mask.${subj}+orig. \
-expr 'a*b' -prefix pb01.${subj}.r${r}.sc_dt_hp_am_blur+orig

rm -f rm.pb01*

done

onsetdir=/group_hpc/WMShimLab/PSY_ColorStudy/
3dDeconvolve -input pb01.${subj}.r*.sc_dt_hp_am_blur+orig.HEAD                     \
    -censor censor_${subj}_combined_2.1D                                    \
    -polort 3                                                               \
    -num_stimts 9                                                          \
    -stim_times 1 ${onsetdir}loc_color.txt 'BLOCK(12,1)'           \
    -stim_label 1 Color                                                    \
    -stim_times 2 ${onsetdir}loc_gray.txt 'BLOCK(12,1)'           \
    -stim_label 2 Grayscale                                                 \
    -stim_times 3 ${onsetdir}loc_fix.txt 'BLOCK(12,1)'         \
    -stim_label 3 Fixation                                              \
    -stim_file 4 motion_demean.1D'[0]' -stim_base 4 -stim_label 4 roll  \
    -stim_file 5 motion_demean.1D'[1]' -stim_base 5 -stim_label 5 pitch \
    -stim_file 6 motion_demean.1D'[2]' -stim_base 6 -stim_label 6 yaw   \
    -stim_file 7 motion_demean.1D'[3]' -stim_base 7 -stim_label 7 dS    \
    -stim_file 8 motion_demean.1D'[4]' -stim_base 8 -stim_label 8 dL    \
    -stim_file 9 motion_demean.1D'[5]' -stim_base 9 -stim_label 9 dP    \
    -local_times                                                        \
    -gltsym 'SYM: +1*Color -1*Grayscale'  				\
    -glt_label 1 Color-Gray                                             \
    -gltsym 'SYM: +1*Grayscale -1*Color'  				\
    -glt_label 2 Gray-Color                                             \
    -gltsym 'SYM: +2*Grayscale -1*Color -1*Fixation'  			\
    -glt_label 3 Gray-all                                               \
    -gltsym 'SYM: +1*Color -1*Grayscale -1*Fixation'  			\
    -glt_label 4 Color-all                                               \
    -gltsym 'SYM: +1*Color -1*Fixation'  			\
    -glt_label 5 Color-Fix                                               \
    -gltsym 'SYM: +1*Grayscale -1*Fixation'  			\
    -glt_label 6 Grayscale-Fix                                               \
    -gltsym 'SYM: +1*Color +1*Grayscale -2*Fixation'  			\
    -glt_label 7 Color+Gray-Fix                                               \
    -float                                                              \
    -jobs 8                                                             \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                             \
    -x1D_uncensored X.nocensor.xmat.1D                                  \
    -bucket stats.all

3dbucket stats.all+orig.HEAD[23] -prefix stat.Color-Fix
3dcalc -a stat.Color-Fix+orig. -expr 'ispositive(a-3.936)' -prefix thresmask_v3_q.001
3dcalc -a stat.Color-Fix+orig. -expr 'ispositive(a-1.648)' -prefix thresmask_v3_p.1
###

SNs=(01 03 04 05 06 07 08 09 10 11 12 13 14)
for SN in "${SNs[@]}"
do
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LGN_PD_DIR=${SN_DIR}/LGN_PD
ROI_DIR=${SN_DIR}/forwardmodel/roi

3dcalc -a ${LGN_PD_DIR}/HK2/LGN_hk@main+orig. -b ${SN_DIR}/${SN}_loc.results/thresmask_v3_p.1+orig. \
-expr 'ispositive(a*b)' -prefix ${ROI_DIR}/LGN_hk2_p.1.nii -overwrite
done


for SN in "${SNs[@]}"
do

SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data/${SN}_main.results
LGN_PD_DIR=${SN_DIR}/LGN_PD
NII_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data/forwardmodel/nii


mv ${NII_DIR}/${SN}_main_combined_sc_dt_hp_am.nii  ${NII_DIR}/${SN}_main_combined_sc_dt_hp_am_old.nii
cd ${SN_DIR}
3dTcat -prefix ${NII_DIR}/${SN}_main_combined_sc_dt_hp_am.nii \
	pb01.${SN}_main.r0?.sc_dt_hp_am+orig.HEAD 

done







###

SN=13
rm -rf ${SN}_main.results
rm -rf ${SN}_loc.results
rm -rf proc*
rm -rf T1@main+orig.*


epi_name=main
run=( 01 02 03 04 05 06 07 08 )
for r in "${run[@]}"
do
@Align_Centers -no_cp -base T1_${epi_name}_SS+orig. -dset ${epi_name}_${r}.nii
done

afni_proc.py -subj_id ${SN}_main -dsets main_01.nii main_02.nii main_03.nii  \
     main_04.nii main_05.nii main_06.nii main_07.nii main_08.nii           \
     -volreg_base_dset 'main_05.nii[0]' -copy_anat T1_${epi_name}_SS+orig.             \
     -regress_censor_motion 0.5 -regress_censor_outliers 0.1 -blocks align \
     volreg mask regress

tcsh -xef proc.${SN}_main |& tee output.proc.${SN}_main



cd ${SN}_main.results/

align_epi_anat.py -anat2epi -anat T1_main_SS+orig -suffix _al_do -epi  external_volreg_base+orig -epi_base 0 -epi_strip 3dAutomask -volreg off  -tshift off -deoblique off -anat_has_skull no -giant_move -overwrite
3drename T1_main_SS_al_do+orig. T1@main
cp T1@main+orig.* ../




#align T1_main_SS to main_05.nii[0]
align_epi_anat.py -dset1 T1_main_SS+orig -dset2 main_05.nii'[0]' -dset1to2 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no
3drename inplane_main_SS_al+orig. inplane@main

#align T1_main_SS to main_05.nii[0]
align_epi_anat.py -dset1 ../T1_main_SS+orig -dset2 inplane@main+orig -dset1to2 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no
3drename T1_main_SS_al+orig. T1@main




#-------------------make thresmask_v3_p.05
SNs=(01 03 04 05 06 07 08 09 10 11 12 13 14 )

for SN in "${SNs[@]}"
do

SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LOC_DIR=${SN_DIR}/${SN}_loc.results
LGN_PD_DIR=${SN_DIR}/LGN_PD
ROI_DIR=${SN_DIR}/forwardmodel/roi
MAIN_DIR=${SN_DIR}/${SN}_main.results

cd ${LOC_DIR}
3dcalc -a stat.Color-Fix+orig. -expr 'ispositive(a-1.964)' -prefix thresmask_v3_p.05


3dcalc -a ${LGN_PD_DIR}/HK2/LGN_hk@main+orig. -b ${SN_DIR}/${SN}_loc.results/thresmask_v3_p.05+orig. \
-expr 'ispositive(a*b)' -prefix ${ROI_DIR}/LGN_hk2_p.05.nii -overwrite

done


#align hk3
#-------------------make thresmask_hk3_p.05

SN=13
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LGN_PD_DIR=${SN_DIR}/LGN_PD

PD_id=190107JHY
PD_DIR=/group_hpc/WMShimLab/PD/${PD_id}

mkdir -p ${LGN_PD_DIR}/HK3
cd ${LGN_PD_DIR}/HK3
chmod -R 777  ${LGN_PD_DIR}/HK3


#3dcopy ${PD_DIR}/T1_SS+orig. ${LGN_PD_DIR}/HK3/T1_PD+orig.
3dcopy ${PD_DIR}/PD_mean_al+orig. ${LGN_PD_DIR}/HK3/meanPD@T1_PD
cp ${PD_DIR}/LGN_l_hk3.nii.gz ${LGN_PD_DIR}/HK3/
cp ${PD_DIR}/LGN_r_hk3.nii.gz ${LGN_PD_DIR}/HK3/
3dcopy ../HK2/T1@main+orig. ${LGN_PD_DIR}/HK3/ -overwrite


cp ../HK2/*_al_mat.aff12.1D .

#align PD @ T1@main
#align_epi_anat.py -dset1 meanPD@T1_PD+orig -dset2 T1@main+orig -dset1to2 \
# -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no -overwrite



#align T1_PD @ T1@main 
align_epi_anat.py -dset1 PD_mean_SS_al+orig. -dset2 T1@main+orig -dset1to2 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no -overwrite


lr=(l r)
for r in "${lr[@]}"
do
#align LGN @ T1@main using T1_PD_al_mat 
3dAllineate  -final NN -1Dmatrix_apply PD_mean_SS_al_al_mat.aff12.1D \
-master ../../${SN}_main.results/pb01.${SN}_main.r01.volreg+orig \
-prefix LGN_${r}_hk3@main+orig.HEAD \
LGN_${r}_hk3.nii -overwrite
done

#merge LGN 
3dmerge -gmax -prefix LGN_hk3@main LGN_l_hk3@main+orig. LGN_r_hk3@main+orig. -overwrite

SNs=(01 02 03 05 06 07 08 09 10)

for SN in "${SNs[@]}"
do

SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LOC_DIR=${SN_DIR}/${SN}_loc.results
LGN_PD_DIR=${SN_DIR}/LGN_PD
ROI_DIR=${SN_DIR}/forwardmodel/roi

3dcalc -a ${LGN_PD_DIR}/HK3/LGN_hk3@main+orig. -b ${SN_DIR}/${SN}_loc.results/thresmask_v3_p.05+orig. \
-expr 'ispositive(a*b)' -prefix ${ROI_DIR}/LGN_hk3_p.05.nii 
done

##############################################

#align hk3
#-------------------make thresmask_hk3_p.05

SN=10CES
SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LGN_PD_DIR=${SN_DIR}/LGN_PD

PD_id=181022CES
PD_DIR=/group_hpc/WMShimLab/PD/${PD_id}

mkdir -p ${LGN_PD_DIR}/HK4
cd ${LGN_PD_DIR}/HK4
chmod -R 777  ${LGN_PD_DIR}/HK4


#3dcopy ${PD_DIR}/T1_SS+orig. ${LGN_PD_DIR}/HK3/T1_PD+orig.
3dcopy ${PD_DIR}/PD_mean_al+orig. ${LGN_PD_DIR}/HK4/meanPD@T1_PD
cp ${PD_DIR}/LGN_l_hk4_rsm.nii ${LGN_PD_DIR}/HK4/
cp ${PD_DIR}/LGN_r_hk4_rsm.nii ${LGN_PD_DIR}/HK4/
3dcopy ../HK2/T1@main+orig. ${LGN_PD_DIR}/HK4/ -overwrite

cp ../HK2/*_al_mat.aff12.1D .

#align PD @ T1@main
#align_epi_anat.py -dset1 meanPD@T1_PD+orig -dset2 T1@main+orig -dset1to2 \
# -cost mi -deoblique off -feature_size 0.5 -ginormous_move -anat_has_skull no -overwrite



lr=(l r)
for r in "${lr[@]}"
do
#align LGN @ T1@main using T1_PD_al_mat 
3dAllineate  -final NN -1Dmatrix_apply PD_mean_al_al_mat.aff12.1D \
-master ../../${SN}_main.results/pb01.${SN}_main.r01.volreg+orig \
-prefix LGN_${r}_hk4@main+orig.HEAD \
LGN_${r}_hk4_rsm.nii -overwrite
done

#merge LGN 
3dmerge -gmax -prefix LGN_hk4@main LGN_l_hk4@main+orig. LGN_r_hk4@main+orig. -overwrite

SNs=(01 02 03 05 06 07 08 09 10)

for SN in "${SNs[@]}"
do

SN_DIR=/group_hpc/WMShimLab2/PSY_Color/Colorv3/${SN}/Img_data
LOC_DIR=${SN_DIR}/${SN}_loc.results
LGN_PD_DIR=${SN_DIR}/LGN_PD
ROI_DIR=${SN_DIR}/forwardmodel/roi

3dcalc -a ${LGN_PD_DIR}/HK4/LGN_hk4@main+orig. -b ${SN_DIR}/${SN}_loc.results/thresmask_v3_p.05+orig. \
-expr 'ispositive(a*b)' -prefix ${ROI_DIR}/LGN_hk4_p.05.nii 
done
###############33
