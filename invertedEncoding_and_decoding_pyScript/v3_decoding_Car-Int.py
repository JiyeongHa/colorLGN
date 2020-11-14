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
        split_fds[i].sa.TRs[lagTRs:] = \
            split_fds[i].sa.TRs[:-lagTRs]

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
import matplotlib as plt

# a few more settings for searchlight
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
mvpa2.debug.active = ['APERM']


# subject info.
sName = ['01', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14']
# sName = ['01', '03', '04', '05', '06', '07', '08', '09', '10', '11', '13', '14']
# experiment info.
exp = 'v3'
cond = ["In", "Out"]
ROI = ["LGN_hk"]
# ROI = ["V1_3dg_fmasked_q.001", "V2_3dg_fmasked_q.001", "V3_3dg_fmasked_q.001", "V4v_3dg_fmasked_q.001", "center3dg_fmasked_q.001", "center1dg_fmasked_q.001"]
#ROI = ["center3dg_fmasked_q.001", "center1dg_fmasked_q.001"]

nRun = 8
nTR = 197
lagTR =2
nTargets = 2
for xROI in range(0, len(ROI)):
    group_cond_results = np.zeros((len(sName), len(cond)))
    confusion_matrix = np.zeros((len(cond),nTargets,nTargets))
    print "ROI: %s" % (ROI[xROI])
    for xSN in range(0, len(sName)):
        print "SN no.%s" % (sName[xSN])
        # path
        basedir = '/group_hpc/WMShimLab2/PSY_Color/Color%(Exp)s' % {"Exp": exp}
        os.chdir(basedir)
        subj_path = os.path.abspath('%(subj)s/Img_data/forwardmodel' % {"subj": sName[xSN]})
        fds_path = '%(sjpath)s/nii' % {"sjpath": subj_path}
        attr_path = '%(sjpath)s/attr' % {"sjpath": subj_path}
        roi_path = '%(sjpath)s/roi' % {"sjpath": subj_path}
        result_path = '/home/jiyeongha/colorLGN/%(Exp)s' % {"Exp": exp}
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        nii_file = '%(fds_path)s/%(subj)s_main_combined_sc_dt_hp_am.nii' % {"fds_path": fds_path ,"subj": sName[xSN]}
        roi_file = '%(roi_path)s/%(roi)s' % {"roi_path": roi_path, "roi": ROI[xROI] + '.nii'}
        attr_file = '%(attr_path)s/%(subj)s_ColorTuningATTEcc_8con_Att.txt' % {"attr_path": attr_path, "subj": sName[xSN]}
        result_acc_file = '%(result_path)s/%(Exp)s_%(roi)s_CarInt_decoding_used4TRs.csv' % {"Exp": exp, "result_path": result_path, "roi": ROI[xROI]}
        result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_CarInt_cfm_used4TRs.csv' % {"Exp": exp, "result_path": result_path, "roi": ROI[xROI]}

        attr = load_attributes(attr_file=attr_file)
        fds = fmri_dataset(samples=nii_file,
                           chunks=attr.chunks,
                           mask=roi_file)
        raw_attr = np.full(len(attr.censor),1)*7

        fds.sa['color'] = attr.targets
        fds.sa['censor'] = attr.censor
        fds.sa['trials'] = attr.trials
        fds.sa['ecc'] = attr.ecc
        fds.sa['TRs'] = attr.TRs
        fds.sa['carInt'] = np.full(len(attr.censor), 1)*0
        fds.sa['LMS'] = np.full(len(attr.censor), 1)*0

        fds.sa.carInt[fds.sa.color != 0] = np.mod(fds.sa.color[fds.sa.color != 0], 2)+1 #cardinal= 1, intercardinal = 2
        fds.sa.LMS[np.logical_or(fds.sa.color == 1, fds.sa.color == 5)] = 1 # L-M color:1
        fds.sa.LMS[np.logical_or(fds.sa.color == 3, fds.sa.color == 7)] = 2 # S color:2
        fds.sa.LMS[np.mod(fds.sa.color,2)==0] = 3 # S color:2


        # ===prep===
        fds.samples = asarray(fds.samples)
        fds = lag_correction(fds=fds, runTRs=nTR, lagTRs=lagTR)  # another custom subfunction
        fds = fds[fds.sa.censor ==1]
        fds = fds[fds.sa.color != 0]
        fds.sa['targets'] = fds.sa.carInt

        for xCond in range(0, len(cond)): # 1- in, 2 - out
            print "%s condition" % (cond[xCond])
            fds_cond = fds[fds.sa.ecc == xCond+1]
            fds_cond = fds_cond[fds_cond.sa.TRs != 5]
            zscore(fds_cond, chunks_attr = 'chunks')
            averager = mean_group_sample(['targets', 'chunks'])
            fds_cond = fds_cond.get_mapped(averager)


            ######## SVM  ########
            nfeatures = ceil(fds_cond.nfeatures * 1)
            clf = LinearCSVMC()
            fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
                                                    FixedNElementTailSelector(nfeatures, tail='upper',
                                                                              mode='select', sort=False))
            fclf = FeatureSelectionClassifier(clf, fsel)
            partitioner = ChainNode([NFoldPartitioner(cvtype=1),
                                     Balancer(attr='targets',
                                              count=1,  # for real data > 1
                                              limit='partitions',
                                              apply_selection=True
                                              )],
                                    space='partitions')
            cvte = CrossValidation(fclf, partitioner,
                                   errorfx=mean_mismatch_error,
                                   postproc=mean_sample(),
                                   enable_ca=['stats'])

            err_result = cvte(fds_cond)
            group_cond_results[xSN,xCond] = np.round(cvte.ca.stats.stats['ACC%'],4)
            confusion_matrix[xCond,:,:] = confusion_matrix[xCond,:,:] + np.true_divide(cvte.ca.stats.matrix, nRun)

    confusion_matrix = confusion_matrix = confusion_matrix.reshape(nTargets*len(cond), nTargets)
    confusion_matrix_mean = np.true_divide(confusion_matrix, len(sName))

    np.savetxt(result_acc_file, group_cond_results, delimiter=',')
    np.savetxt(result_cfm_file, confusion_matrix_mean, delimiter=',')

