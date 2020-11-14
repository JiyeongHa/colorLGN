def makeSin_power(p, powpow):
    x = linspace(0, pi, 180)
    pred = sin(x) ** powpow
    pred = pred / max(pred)
    peak = where(pred == 1)[0]

    tmp = abs(x - p['u']).min()
    shiftInd = abs(x - p['u']).argmin()

    pred = roll(pred, shiftInd - peak[0])  # direction of shift is opposite from 'wshift' in MATLAB

    # here - the x axis should be derived based on teh p.cCenters, so that each
    # function peaks properly...then go back to just shifting N times over the
    # space, not 180 times...
    resampInd = ceil(float(len(x)) / float(nChans))
    pred = pred[0:-1:int(resampInd)]

    return pred
def imagesc(X, scale=None, cmap=None, **kwds):
    """ implements imagesc as in matlab"""
    import pylab

    # pylab.figure()

    if scale is None:
        vmin = None
        vmax = None
    else:
        vmin = scale[0]
        vmax = scale[1]

    if not kwds.has_key('extent'):
        kwds['extent'] = [0, 1, 0, 1]

    if not kwds.has_key('interpolation'):
        kwds['interpolation'] = 'nearest'

    a = pylab.imshow(X, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap, **kwds)
    show()
    return a
def load_attributes(attr_file, attr_dir):
    x = os.path.join(attr_dir, attr_file)
    attr = ColumnData(x, header=True)
    # attr = SampleAttributes(x)
    return attr
def load_nii(nii_file, mask_file, attr, nii_dir, mask_dir):
    """load experiment dataset"""

    fds = fmri_dataset(samples=os.path.join(nii_dir, nii_file),
                       targets=attr.targets, chunks=attr.chunks,
                       mask=os.path.join(mask_dir, mask_file))

    return fds
def lag_correction(fds, runTRs, lagTRs):
    """correct dataset for hemodynamic lag"""

    # split dataset into runs
    nRuns = len(fds) / int(runTRs)
    if int(nRuns) != nRuns:
        print "Error! number of TRs per run must be a factor of total TRs"
        raise SystemExit

    split_fds = []

    for i in range(nRuns):  # split dataset into separate runs
        split_fds.append(fds[i * runTRs:(i + 1) * runTRs])

    # do the shift for each run
    for i in range(len(split_fds)):
        split_fds[i].sa.targets[lagTRs:] = (split_fds[i]
                                                .sa.targets[:-lagTRs])  # need to shift target labels too

        split_fds[i].sa.censor[lagTRs:] = (split_fds[i]
                                               .sa.censor[:-lagTRs])  # and censor labels

        split_fds[i].sa.trials[lagTRs:] = (split_fds[i]
                                               .sa.trials[:-lagTRs])

        split_fds[i].sa.ecc[lagTRs:] = (split_fds[i]
                                            .sa.ecc[:-lagTRs])  # shift ecc condition labels

        split_fds[i].sa.TRs[lagTRs:] = (split_fds[i]
                                            .sa.TRs[:-lagTRs])  # shift ecc condition labels

        split_fds[i] = split_fds[i][lagTRs:]

    # merge back datasets
    fds = split_fds[0]
    for i in range(1, len(split_fds)):
        fds.append(split_fds[i])

    return fds

import os
from mvpa2.suite import *
from numpy import *
from pylab import *

##Subject info
subs = ['01', '03', '04', '05', '06', '07', '08', '09', '10', '11', '13', '14']

##Experiment info
nColors = 8  # number of actual orientations in your study
nChans = 8  # number of channels that you want to use as your 'basis set' to model the response in each voxel.
nTrials = 32  # number trials per orientation for the simulation (must be evenly divisible by nColors)
nCond = 2  # In or Out
nRun = 8
lagTR = 2
totalTR = 197

##ROI info
# ROIs = ['LGN_hk']
ROIs = ['V1_3dg_fmasked_q.001', 'V2_3dg_fmasked_q.001', 'V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001']
newdir = '/Volumes/Duri/data/Colorv3/'
ROIdir = 'roi'
result_dir = 'sc_dt_hp_am'

##IEM parameters
powpow = 6  # sin power for basis function
p = {'chan_num_colors': linspace(0, pi - pi / (nColors), nColors),
     'chan_num_chans': linspace(0, pi - pi / (nChans), nChans)}

for subjnum in subs:
    os.chdir(newdir + subjnum + '/Img_data/forwardmodel/')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    attr_file = subjnum + '_ColorTuningATTEcc_8con_Att.txt'
    for i in ROIs:
        print i
        csvname = subjnum + '_tuning_' + i + '2TR_used4TRs_neutralModel'
        print csvname
        channel_response_unshifted = np.zeros((nRun, nCond, nColors, nChans))
        channel_response_shifted = np.zeros((nRun, nCond, nColors, nChans))
        attr = load_attributes(attr_file=attr_file, attr_dir='attr')
        fds = load_nii(nii_file=subjnum + '_main_combined_sc_dt_hp_am.nii', mask_file=i + '.nii',
                       attr=attr, nii_dir='nii', mask_dir=ROIdir)
        fds.sa['censor'] = attr.censor
        fds.sa['trials'] = attr.trials
        fds.sa['ecc'] = attr.ecc
        fds.sa['TRs'] = attr.TRs
        print 'Voxel Number Orignial: %s' % (fds.nfeatures)
        fds.samples = asarray(fds.samples)
        # ===lag correction===
        fds = lag_correction(fds=fds, runTRs=totalTR, lagTRs=lagTR)  # another custom subfunction
        # ===censoring===
        fds = fds[fds.sa.censor == 1]  # remove censored points
        # ===detrend===
        #   poly_detrend(fds, polyord=1, chunks_attr='chunks')
        # ===zscore per run===
        zscore(fds, chunks_attr='chunks')
        # ===remove 'rest' TRs===
        fds = fds[fds.sa.targets != 0]
        fds = fds[fds.sa.TRs != 5]
        fds = fds[fds.sa.TRs != 4]  # use only 4 TRs for analysis
        averager = mean_group_sample(['ecc', 'targets', 'chunks'])
        fds = fds.get_mapped(averager)

        for rr in np.arange(1, nRun + 1):
            print "Computing iteration %d out of %d\n" % (rr, nRun)
            # ### feature selection ###
            # clf = LinearCSVMC()
            # fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
            #                                         FixedNElementTailSelector(ceil(fds.nfeatures), tail='upper',
            #
            trn = fds[fds.sa.chunks != rr]  # data from training scans (all but one scan)
            tst = fds[fds.sa.chunks == rr]
            trn.samples = trn.samples + 1  # add an arbitrary number
            tst.samples = tst.samples + 1

            # Extract fMRI data for training
            uRuns = unique(trn.sa.chunks)
            x_trn = np.zeros((nColors * len(uRuns), trn.samples.shape[1]))
            for xRun in np.arange(0, len(uRuns)):
                print "Training run no. %d\n" % (uRuns[xRun])
                for xColor in np.arange(0, nColors):
                    x_trn[nColors * xRun + xColor, :] = nanmean(
                        trn.samples[logical_and(trn.sa.chunks == uRuns[xRun], trn.sa.targets == xColor+1), :], 0)
            # Build arbitrary channel for training
            c_trn = np.zeros((shape(x_trn)[0], nChans))
            for xChan in arange(nChans):
                p['u'] = p['chan_num_chans'][xChan]
                c_trn[:, xChan] = tile(makeSin_power(p, powpow), [1, shape(x_trn)[0] / nChans])

            x_trn[isnan(x_trn)] = 0
            c_trn[isnan(c_trn)] = 0
            # Training: get weight for each color of each voxels
            w = dot(dot(linalg.inv(dot(c_trn.T, c_trn)), c_trn.T), x_trn)
            # Testing: get channel responses
            c_tst = dot(dot(linalg.inv(dot(w, w.T)), w), tst.samples.T)
            c_tst = c_tst.T
            # All done. stack both unshifted & shifted responses...
            for xEcc in np.arange(0, nCond):
                channel_response_unshifted[rr-1, xEcc, :, :] = c_tst[tst.sa.ecc == xEcc + 1]
                for xColor in np.arange(0, nColors):
                    channel_response_shifted[rr-1, xEcc, xColor, :] = roll(
                        channel_response_unshifted[rr-1, xEcc, xColor, :], int(ceil(nChans / 2) - (xColor + 1)))

        mean_channel_response_unshifted = np.mean(channel_response_unshifted, 0)
        mean_channel_response_shifted = np.mean(channel_response_shifted, 0)
        stacked_channel_response_unshifted = np.vstack((mean_channel_response_unshifted[0, :, :],
                                                        mean_channel_response_unshifted[1, :, :]))
        stacked_channel_response_shifted = np.vstack((mean_channel_response_shifted[0, :, :],
                                                      mean_channel_response_shifted[1, :, :]))
        savetxt(os.path.join(result_dir, csvname + '_unshift.txt'), stacked_channel_response_unshifted, '%f', '\t')
        savetxt(os.path.join(result_dir, csvname + '_shift.txt'), stacked_channel_response_shifted, '%f', '\t')


