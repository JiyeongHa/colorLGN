import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# experiment info.
exp = 'v3'
cond = ["In", "Out"]
ROI = ["LGN_hk", "V1_fmasked_q.001", "V2_fmasked_q.001", "V3_fmasked_q.001", "V4v_fmasked_q.001"]
ROI = ["center3dg_fmasked_q.001", "center1dg_fmasked_q.001"]

#ROI = ["V1_fmasked_q.001", "V2_fmasked_q.001", "V3_fmasked_q.001", "V4v_fmasked_q.001"]

nRun = 8
nTR = 197
lagTR = 1
nTargets = 8

result_path = '/home/jiyeongha/colorLGN/%(Exp)s' % {"Exp": exp}

xROI=0
text = '8colors'

fig, axes = plt.subplots(nrows=len(ROI), ncols=len(cond), sharex=True, sharey=True, constrained_layout=True)

for i in range(0, len(ROI)):
    # result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_8colors_cfm.csv' % {"Exp": exp, "result_path": result_path,
    #                                                                       "roi": ROI[i]}
    # result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_CarInt_cfm.csv' % {"Exp": exp, "result_path": result_path,
    #                                                                       "roi": ROI[i]}
    result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_LMS4_cfm.csv' % {"Exp": exp, "result_path": result_path,
                                                                        "roi": ROI[i]}
    # result_cfm_file = '%(result_path)s/%(Exp)s_%(roi)s_LMS2_cfm.csv' % {"Exp": exp, "result_path": result_path,
    #                                                                     "roi": ROI[i]}
    cfm = np.genfromtxt(result_cfm_file, delimiter=',')
    cfm = np.round(cfm, 4)
    nTargets = cfm.shape[1]
    in_cfm = cfm[0:nTargets, :]
    out_cfm = cfm[0 + nTargets:nTargets * 2, :]
    min_val = np.true_divide(1, nTargets) - 0.05
    sns.heatmap(in_cfm, cmap="YlGnBu", cbar = False, ax=axes[i,0], linewidths=.4, annot=False, vmin=min_val, vmax=1)
    sns.heatmap(out_cfm, cmap="YlGnBu", ax=axes[i,1], linewidths=.4, annot=False, vmin=min_val, vmax=1)
    axes[i,0].set_ylabel('%(roi)s' % {"roi": ROI[i]})

axes[0,0].set_title('In')
axes[0,1].set_title('Out')
#plt.setp(axes, xticks=np.arange(1,9), yticks=np.arange(1,9))
plt.show()


#axes[1,0].set_title('%(roi)s %(txt)s in cfm' % {"roi": ROI[i], "txt": text})
#axes[0,1].set_title('%(roi)s %(txt)s out cfm' % {"roi": ROI[i], "txt": text})
