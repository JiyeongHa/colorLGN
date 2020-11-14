#!/bin/sh


SN=( 11 )
run=( 01 02 03 04 05 06 07 08 )
PDid=181224YYH
exp=main
 main_02.nii main_03.nii  \
     main_04.nii main_05.nii main_06.nii main_07.nii main_08.nii 

afni_proc.py -subj_id ${SN}_${exp} -dsets ${exp}_0?.nii          \
     -volreg_base_dset 'main_01.nii[0]' -copy_anat inplane_${exp}.nii             \
     -regress_censor_motion 0.5 -regress_censor_outliers 0.1 -blocks align \
     volreg mask regress

tcsh -ef proc.${SN}_${exp}

#root_dir=/sas2/PECON/PSY/Colorv3/
SNdir=/group_hpc/WMShimLab/PSY_ColorStudy/v3

	##------------------------02 Scale & remove skull in EPI -----------
for xsn in "${SNs[@]}"
do

subj=${SNs}_${exp}
cd ${SNdir}/${xsn}/Img_data/${subj}.results
for r in "${run[@]}"
do
3dClipLevel pb01.${subj}.r${r}.volreg+orig.HEAD >> clip.txt
3dTstat -mean -prefix r.${r}.base pb01.${subj}.r${r}.volreg+orig.HEAD'[0..$]'
done

more clip.txt # check the smallest clip value across all runs
clip=$(cut -f1 -d"," clip.txt | sort -n | head -1) # assign the clip value

for r in "${run[@]}"
do
3dcalc -a pb01.${subj}.r${r}.volreg+orig. -b r.${r}.base+orig. -expr "(100 * a/b) * step(b-$clip)" -prefix pb01.${subj}.r${r}.scaled
3dTstat -mean -prefix r.${r}.sc_base pb01.${subj}.r${r}.scaled+orig.HEAD'[0..$]'
# detrend linear trend
3dDetrend -polort 1 -prefix pb01.${subj}.r${r}.sc_dt pb01.${subj}.r${r}.scaled+orig
# add mean after detrend
3dcalc -a pb01.${subj}.r${r}.sc_dt+orig.HEAD -b r.${r}.sc_base+orig.HEAD -expr 'a+b' -prefix pb01.${subj}.r${r}.sc_dt_am

# calculate mean image
3dTstat -mean -prefix r.${r}.sc_dt_base pb01.${subj}.r${r}.sc_dt_am+orig.HEAD'[0..$]'

# high-pass filter
3dBandpass -prefix pb01.${subj}.r${r}.sc_dt_hp 0.01 99999 pb01.${subj}.r${r}.sc_dt_am+orig
# add mean after hp filter
3dcalc -a pb01.${subj}.r${r}.sc_dt_hp+orig.HEAD -b r.${r}.sc_dt_base+orig.HEAD \
-expr 'a+b' -prefix pb01.${subj}.r${r}.sc_dt_hp_am	

done
done

##Finished. Let's combine all data-------------------------------------
	3dTcat -prefix FNL.nii \
	pb01.${ww}.r0?.sc_dt_hp_am+orig.HEAD 

	# (tip) when calling variables(like $clip) within 3dcalc -expr, be sure to use "", and not ''

============aligning LGN image to your experiment===============#
# align T1 from PD scan to your T1, and then align PD with your T1
# double-check by aligning mean PD image to your T1 image!

#align T1_main to inplane_main
3dSkullStrip -input T1_main.nii -prefix T1_main_SS.nii
3dSkullStrip -input inplane_main.nii -prefix inplane_main_SS.nii
@Align_Centers -base inplane_main_SS.nii -dset T1_main_SS.nii -no_cp
align_epi_anat.py -dset1 inplane_main_SS.nii -dset2 T1_main_SS.nii -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move
3drename T1_SS_main_al+orig. T1@main


PDdir=/group_hpc/WMShimLab/PD/${PDid}
mkdir -p ${expdir}
cp ${PDdir}/LGN_*masked+orig.* .
cp ${PDdir}/T1_SS_al+orig.HEAD T1_PD_SS+orig.HEAD
cp ${PDdir}/T1_SS_al+orig.BRIK T1_PD_SS+orig.BRIK

align_epi_anat.py -dset1 T1_PD_SS+orig -dset2 T1@main+orig -dset1to2 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move

3dAllineate -cubic -1Dmatrix_apply T1_PD_SS_al_mat.aff12.1D \
-master T1@main+orig.HEAD \
-prefix LGN_l_masked+orig.HEAD \
LGN_l+orig.HEAD

3dAllineate -cubic -1Dmatrix_apply T1_PD_SS_al_mat.aff12.1D \
-master T1@main+orig.HEAD \
-prefix LGN_r_masked+orig.HEAD \
LGN_r+orig.HEAD

# to double-check with PD image...
cp ${PDdir}/PD_mean_al+orig.* .
3dAllineate -cubic -1Dmatrix_apply T1_PD_SS_al_mat.aff12.1D \
-master T1@main+orig.HEAD \
-prefix PD@main+orig.HEAD \
PD_mean_al+orig.HEAD

LGN=(l r)
for r in "${LGN[@]}"
do
3dresample -dxyz 2 2 2 -inset LGN_${r}_al+orig. -prefix LGN_${r}_rsm+orig. -overwrite 
3dcalc -a LGN_${r}_rsm+orig. -b loc.results/thresmask_p.1+orig. -expr 'ispositive(a*b)' -prefix LGN_${r}_masked
3dBrickStat -sum LGN_${r}_masked+orig.
done



























