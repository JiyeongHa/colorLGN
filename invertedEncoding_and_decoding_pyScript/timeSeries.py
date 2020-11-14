def load_attributes(attr_file):
    x = os.path.join(attr_file)
    attr = ColumnData(x, header=True)
    # attr = SampleAttributes(x)
    return attr
def load_nii(nii_file, mask_file, attr):
    """load experiment dataset"""

    fds = fmri_dataset(samples=os.path.join(nii_file),
                       targets=attr.xMD, chunks=attr.run,
                       mask=os.path.join(mask_file))
    return fds
def lag_correction(fds, runTRs, lagTRs):
    """correct dataset for hemodynamic lag"""

    # split dataset into runs
    nRuns = len(fds) / float(runTRs)
    if int(nRuns) != nRuns:
        print 'Error! number of TRs per run must be a factor of total TRs'
        raise SystemExit

    nRuns = int(nRuns)

    split_fds = []

    for i in range(nRuns):  # split dataset into separate runs
        split_fds.append(fds[i * runTRs:(i + 1) * runTRs])

    # do the shift for each run
    for i in range(len(split_fds)):
        split_fds[i].sa.color[lagTRs:] = \
            split_fds[i].sa.color[:-lagTRs]
        split_fds[i].sa.censor[lagTRs:] = \
            (split_fds[i].sa.censor[:-lagTRs])
        split_fds[i].sa.trials[lagTRs:] = \
            split_fds[i].sa.trials[:-lagTRs]
        split_fds[i].sa.ecc[lagTRs:] = \
            split_fds[i].sa.ecc[:-lagTRs]
        split_fds[i].sa.chunks[lagTRs:] = \
            split_fds[i].sa.chunks[:-lagTRs]
        split_fds[i].sa.carInt[lagTRs:] = \
            split_fds[i].sa.carInt[:-lagTRs]
        split_fds[i].sa.LMS[lagTRs:] = \
            split_fds[i].sa.LMS[:-lagTRs]

        split_fds[i] = (split_fds[i])[lagTRs:]

    ##  merge back datasets
    fds = split_fds[0]
    for i in range(1, len(split_fds)):
        fds.append(split_fds[i])

    return fds

# libraries needed by pymvpa
import os
import sys
from mvpa2.suite import *
import numpy as np
import matplotlib.pyplot as plt

# a few more settings for searchlight
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
mvpa2.debug.active = ['APERM']


# subject info.
sName = ['01', '03', '04', '05', '06', '07', '08', '09', '10', '11', '13', '14','12']
# sName = ['01', '03', '04', '05', '06', '07', '08', '09', '10', '11', '13', '14']

# experiment info.
exp = 'v3'
cond = ["In", "Out"]
ROI = ["LGN_hk", "V1_3dg_fmasked_q.001"]
# ROI = ["V1_fmasked_q.001", "V2_fmasked_q.001", "V3_fmasked_q.001", "V4v_fmasked_q.001","center3dg_fmasked_q.001", "center1dg_fmasked_q.001"]
# ROI = ["center3dg_fmasked_q.001", "center1dg_fmasked_q.001"]

nRun = 4
nTR = 197
lagTR = 2
nTargets = 4
roi_ts_results = np.zeros((len(ROI), 2))
roi_ts_sd_results = np.zeros((len(ROI), 2))

roi_ts_results_new = np.zeros((len(ROI), 2))

for xROI in range(0, len(ROI)):
group_ts_results = np.zeros((len(sName), 2))
print "ROI: %s" % (ROI[xROI])
for xSN in range(0, len(sName)):
print "SN no.%s" % (sName[xSN])
# path
basedir = '/group_hpc/WMShimLab2/PSY_Color/Color%(Exp)s' % {"Exp": exp}
os.chdir(basedir)
subj_path = os.path.abspath('%(subj)s/Img_data/forwardmodel' % {"subj": sName[xSN]})
fds_path = '%(sjpath)s/nii' % {"sjpath": subj_path}
#fds_path = '/group_hpc/WMShimLab2/PSY_Color/Color%(Exp)s/%(subj)s/Img_data/%(subj)s_loc.results/'
attr_path = '%(sjpath)s/attr' % {"sjpath": subj_path}
roi_path = '%(sjpath)s/roi' % {"sjpath": subj_path}
result_path = '/home/jiyeongha/colorLGN/%(Exp)s' % {"Exp": exp}
if not os.path.exists(result_path):
    os.makedirs(result_path)

nii_file = '%(fds_path)s/%(subj)s_main_combined_sc_dt_hp_am.nii' % {"fds_path": fds_path ,"subj": sName[xSN]}
#nii_file = '%(fds_path)s/%(subj)s_loc_combined_sc_dt_hp_am.nii' % {"fds_path": fds_path ,"subj": sName[xSN]}
roi_file = '%(roi_path)s/%(roi)s' % {"roi_path": roi_path, "roi": ROI[xROI] + '.nii'}
attr_file = '%(attr_path)s/%(subj)s_ColorTuningATTEcc_8con_Att.txt' % {"attr_path": attr_path, "subj": sName[xSN]}
mask_file = '%(base_path)s/%(subj)s/Img_data/%(subj)s_main.results/full_mask.%(subj)s_main.nii' % {"base_path": basedir, "subj": sName[xSN]}
result_acc_file = '%(result_path)s/%(Exp)s_%(roi)s_LMS4_decoding_2TR.csv' % {"Exp": exp, "result_path": result_path, "roi": ROI[xROI]}
result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_LMS4_cfm_2TR.csv' % {"Exp": exp, "result_path": result_path, "roi": ROI[xROI]}
result_ts_file = '%(result_path)s/%(Exp)s_%(roi)s_ts_LMS.csv' % {"Exp": exp, "result_path": result_path, "roi": ROI[xROI]}

attr = load_attributes(attr_file=attr_file)
#fds = fmri_dataset(samples=nii_file,
#               chunks=np.asarray([np.full(nTR,1),np.full(nTR,2), np.full(nTR,3), np.full(nTR,4)]).reshape(-1,), #attr.chunks,
#               mask=roi_file)
fds = fmri_dataset(samples=nii_file, chunks=attr.chunks, mask=roi_file)
#raw_attr = np.full(len(attr.censor),1)*7

fds.sa['color'] = attr.targets
fds.sa['censor'] = attr.censor
fds.sa['trials'] = attr.trials
fds.sa['ecc'] = attr.ecc
fds.sa['TRs'] = attr.TRs
fds.sa['carInt'] = np.full(len(attr.censor), 1)*0
fds.sa['LMS'] = np.full(len(attr.censor), 1)*0
#fds.sa['TR'] = np.tile(np.arange(1, nTR + 1), nRun)
fds.sa.carInt[fds.sa.color != 0] = np.mod(fds.sa.color[fds.sa.color != 0], 2)+1 #cardinal= 1, intercardinal = 2
fds.samples = asarray(fds.samples)
fds = lag_correction(fds=fds, runTRs=nTR, lagTRs=lagTR)  # another custom subfunction
#fds.sa['TR'] = np.tile(np.arange(1, 198), nRun)
#fds = fds[~np.logical_or(fds.sa.TR == 1, fds.sa.TR == 2)]
#fds = fds[~np.logical_or(fds.sa.TR == 3, fds.sa.TR == 196)]
#fds = fds[fds.sa.TR != 197]

averager = mean_group_sample(['chunks'])
fds_cond = fds.get_mapped(averager)

fds.sa.LMS[np.logical_or(fds.sa.color == 1, fds.sa.color == 5)] = 1 # L-M color:1
fds.sa.LMS[np.logical_or(fds.sa.color == 3, fds.sa.color == 7)] = 2 # S color:2
fds.sa.LMS[np.logical_or(fds.sa.color == 2, fds.sa.color == 6)] = 3 # int 1
fds.sa.LMS[np.logical_or(fds.sa.color == 4, fds.sa.color == 8)] = 4  # int 2
fds_color = fds[fds.sa.TRs != 5]

# ===prep===
baseline_mean = np.mean(fds_cond.samples, axis=1).reshape(1,-1)
# fds_car = fds_color[fds_color.sa.carInt == 1]
# fds_int = fds_color[fds_color.sa.carInt == 2]
# fds_car = fds_car.get_mapped(averager)
# fds_int = fds_int.get_mapped(averager)
# cardinal_mean = np.mean(fds_car.samples, axis=1).reshape(1,-1)
# intercardinal_mean = np.mean(fds_int.samples, axis=1).reshape(1, -1)
#
# percent_car = np.true_divide(cardinal_mean, baseline_mean)
# percent_car = np.mean(percent_car)
# percent_int = np.true_divide(intercardinal_mean, baseline_mean)
# percent_int = np.mean(percent_int)
fds_LM = fds_color[fds_color.sa.LMS == 1]
fds_S = fds_color[fds_color.sa.LMS == 2]
fds_LM = fds_LM.get_mapped(averager)
fds_S = fds_S.get_mapped(averager)
LM_mean = np.mean(fds_LM.samples, axis=1).reshape(1, -1)
S_mean = np.mean(fds_S.samples, axis=1).reshape(1, -1)

percent_LM = np.true_divide(LM_mean, baseline_mean)
percent_LM = np.mean(percent_LM)
percent_S = np.true_divide(S_mean, baseline_mean)
percent_S = np.mean(percent_S)

group_ts_results[xSN, 0] = percent_LM
group_ts_results[xSN, 1] = percent_S
#averager = mean_group_sample(['TRinTrial'])
#averager = mean_group_sample(['TR'])
#fds = fds.get_mapped(averager)
# group_ts_results[xSN, 0] = percent_car
# group_ts_results[xSN, 1] = percent_int
np.savetxt(result_ts_file, group_ts_results, delimiter=',')
roi_ts_sd_results[xROI,:] = np.true_divide(np.std(group_ts_results, axis=0), np.sqrt(len(group_ts_results)))
roi_ts_results[xROI, :] = np.mean(group_ts_results, axis=0)
roi_ts_results[xROI, :] = np.mean(group_ts_results[0:11, :], axis=0)
roi_ts_sd_results[xROI,:] = np.true_divide(np.std(group_ts_results[0:11, :], axis=0), np.sqrt(len(group_ts_results)))
roi_ts_results_new[:,0:5] = roi_ts_results[:,1:6]
roi_ts_results_new[:,5] = roi_ts_results[:,0]


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#x = np.arange(1,7)
x = np.arange(1,nTR+1)
y1 = roi_ts_results[0,:]
y2 = roi_ts_results[1,:]

function1 = ax1.plot(x, y1, 'r', label='LGN', lw=3)
ax1.set_xlabel('TR', color='k')
ax1.set_ylabel('LGN', color='r')
ax1.set_ylim((95,103))
function2 = ax2.plot(x, y2, 'b', label='V1', lw=3)
ax2.set_ylabel('V1', color='b')
ax2.set_ylim((95,103))
ax2.set_xticks(np.arange(1,nTR+1))

#ax1.legend(handles=[function1, function2])
stimOffset=np.arange(0,nTR,12)
stimOnset = np.arange(6,nTR,6)
ax1.vlines(stimOnset, 0, 1, transform=ax1.get_xaxis_transform(), color='k')
ax1.vlines(stimOffset, 0, 1, transform=ax1.get_xaxis_transform(), color='y')

plt.tight_layout()
plt.show()



ax1.axvline(x=1, color='k', linestyle='--', label='stim onset')
plt.tight_layout()
plt.show()





fig = plt.figure()
plt.plot(x, y1, 'r', linestyle='-', label='LGN', linewidth=5)
plt.plot(x, y2, 'b', linestyle='-', label='V1', linewidth=5)
plt.xlabel('TR', fontsize=25)
plt.legend(loc='best', fontsize=25)
plt.tight_layout()
plt.show()
