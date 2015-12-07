# This script uses block bootstrap to randomize coral data and use different sampling time length to generate distribution plot of seasonal cycle amplitude

import cdms2
import cPickle as pickle
import math
import matplotlib.colorbar as cbar
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from os import listdir, chdir
from os.path import isfile, join
import pcoral_bootstrap
import bootstrap
from scipy import stats

dirs = ['/home/scec-00/julieneg/pmip3/historical', '../past1000', '../midHolocene', '../piControl']
#dirs = '/home/scec-00/julieneg/pmip3/piControl'  # DEBUG ONLY
pcoral_boot = {}  # define output dictionary

Nb = 1000    # number of bootstrap samples     
Lb = 2*12  #  needed to sample multiples of 12 years
windows = [25,50,75,99] # observation windows
nw = windows.__len__()
 
for dir in dirs:                # I'm a loop over experiments
    experiment_name = dir.split('/')[-1]
    # # ---test!!!
    # if experiment_name != 'historical':
    #     continue
    # # !!!test ends---
    print "experiment: %s" %experiment_name
    chdir(dir)
    # get a list of all the coral data filenames in the folder
    names = [ a for a in listdir(dir) if isfile]
    corals = [name for name in names if name.split('_')[-1]=='total.nc']
    corals.sort()

    # change in a list of the index of change in model names in the corals file name list
    change = [0]
    for i in range(len(corals)-1):
        if corals[i].split('_')[0] != corals[i+1].split('_')[0]:
            change.append(i+1)
    change.append(len(corals))

    # pcoral_boot_exp stores all the randomized series of all the models of that experiment
    pcoral_boot_exp = {}; variance = {}; seasonal_amp = {}
    for i in range(len(change)-1): # I'm a loop over models
        model_name = corals[change[i]].split('_')[0]
        print "    model: %s" %model_name

        #for name in corals[change[i]:change[i+1]]: # I'm a loop over runs of a single model
       # JEG: this loop is unenecessary since we lonly keep the last element
        name = corals[change[i]]
        f = cdms2.open(name,'r')
        st = f.getAxis('time').asComponentTime()[0].year
        et = f.getAxis('time').asComponentTime()[-1].year
        if (et-st) < max(windows):
            print " experiment too short"
        print "        working on %s" %name
        # compute bootsrapped climate statistics on the three regions of interest
        # WESTERN PACIFIC
        variance_w, seasonal_amp_w = pcoral_bootstrap.computer(name, 120, 180, -20, 0, Nb, Lb, windows)
        # CENTRAL PACIFIC
        variance_c, seasonal_amp_c = pcoral_bootstrap.computer(name, 190, 240, -5, 5, Nb, Lb, windows)
        # EASTERN PACIFIC
        variance_e, seasonal_amp_e = pcoral_bootstrap.computer(name, 270, 280, -10, 0, Nb, Lb, windows)
            
        # store the results in a temperary dictionary
        model_r = model_name.split('-')[0]
        # store variance results
        variance[model_r] = np.empty((3*nw,Nb))
        variance[model_r][0:nw,:]   = variance_w
        variance[model_r][nw:2*nw,:]  = variance_c
        variance[model_r][2*nw:3*nw,:] = variance_e
        # store seasonal amplitude results
        seasonal_amp[model_r] = np.empty((3*nw,Nb))
        seasonal_amp[model_r][0:nw,:]        = seasonal_amp_w
        seasonal_amp[model_r][nw:2*nw,:]   = seasonal_amp_c
        seasonal_amp[model_r][2*nw:3*nw,:] = seasonal_amp_e
        
    pcoral_boot_exp['var'] = variance
    pcoral_boot_exp['seas'] = seasonal_amp   
        
    pcoral_boot[experiment_name] = pcoral_boot_exp

print "Done!"

chdir('../outputData/combined')

# save the dictionary to a pickle file
# based on http://stackoverflow.com/questions/4893689/save-a-dictionary-to-a-file-alternative-to-pickle-in-python
with open('pmip3_pcoral_bootstrap.p', 'wb') as f:
    pickle.dump(pcoral_boot, f)

## save .mat
import scipy.io as io
io.savemat('pmip3_pcoral_bootstrap.mat', pcoral_boot)
#
print "saved!"

# # three lists to store ratios
# ratio_MH, ratio_LM, ratio_PI = np.empty((12,1000)), np.empty((12,1000)), np.empty((12,1000))

# # compute the ratios
# for i in range(12):
#     numerator_MH = np.array(pcoral_boot['midHolocene']['CCSM4'][i])
#     numerator_LM = np.array(pcoral_boot['past1000']['CCSM4'][i])
#     numerator_PI = np.array(pcoral_boot['piControl']['CCSM4'][i])
#     denominator = np.array(pcoral_boot['historical']['CCSM4'][i])
#     ratio_MH[i] = numerator_MH / denominator
#     ratio_LM[i] = numerator_LM / denominator
#     ratio_PI[i] = numerator_PI / denominator

# # exclude data outside of 2.5 percentile
# top = 10
# bottom = 100 - top
# ratio_MH_ma = np.empty((12,Nb*(bottom-top)/100))
# ratio_LM_ma = np.empty((12,Nb*(bottom-top)/100))
# ratio_PI_ma = np.empty((12,Nb*(bottom-top)/100))
# for i in range(len(ratio_MH)):
#     ratio_MH_bottom = np.percentile(ratio_MH[i], top)
#     ratio_MH_top = np.percentile(ratio_MH[i], bottom)
#     ratio_MH_ma[i,:] = [j for j in ratio_MH[i] if j>ratio_MH_bottom and j<ratio_MH_top]
# for i in range(len(ratio_LM)):
#     ratio_LM_bottom = np.percentile(ratio_LM[i], top)
#     ratio_LM_top = np.percentile(ratio_LM[i], bottom)
#     ratio_LM_ma[i,:] = [j for j in ratio_LM[i] if j>ratio_LM_bottom and j<ratio_LM_top]
# for i in range(len(ratio_PI)):
#     ratio_PI_bottom = np.percentile(ratio_PI[i], top)
#     ratio_PI_top = np.percentile(ratio_PI[i], bottom)
#     ratio_PI_ma[i,:] = [j for j in ratio_PI[i] if j>ratio_PI_bottom and j<ratio_PI_top]

# # transpose
# ratio_MH_T = ratio_MH.T
# ratio_LM_T = ratio_LM.T
# ratio_PI_T = ratio_PI.T

# # assign color according to the p value returned by the KS test 
# MH_colors = []
# LM_colors = []
# norm = colors.Normalize(vmin=-3,vmax=3)
# cmap = cm.RdBu_r
# m = cm.ScalarMappable(norm=norm, cmap=cmap)
# for i in range(12):
#     MH_test = stats.ks_2samp(ratio_MH_ma[i], ratio_PI_ma[i])[1]
#     LM_test = stats.ks_2samp(ratio_LM_ma[i], ratio_PI_ma[i])[1]
#     MH_colors.append(m.to_rgba(np.log10(0.05/MH_test)))
#     print np.log10(0.05/MH_test)
#     LM_colors.append(m.to_rgba(np.log10(0.05/LM_test)))
#     print np.log10(0.05/LM_test)

# # plot
# plt.clf()
# fig = plt.figure()
# plt.subplot(3,1,1)
# plt.axhline(1, color='k', linewidth=0.5,zorder=0)
# box1 = plt.boxplot(ratio_MH_T, notch=True, patch_artist=True, sym='')
# for patch, color in zip(box1['boxes'], MH_colors):
#     patch.set_facecolor(color)
# plt.xticks(range(1,13), ['25','50','75','100','25','50','75','100','25','50','75','100','25','50','75','100'])
# plt.ylim(0,4)
# plt.ylabel('MH/HT')
# plt.title("pacific pseudocoral $\delta {}^{18}O$ 2-7 yr bandpass filtered variance ratio\n(model=CCSM4)")
# #
# plt.subplot(3,1,2)
# plt.axhline(1, color='k', linewidth=0.5,zorder=0)
# box2 = plt.boxplot(ratio_LM_T, notch=True, patch_artist=True, sym='')
# for patch, color in zip(box2['boxes'], LM_colors):
#     patch.set_facecolor(color)
# plt.xticks(range(1,13), ['25','50','75','100','25','50','75','100','25','50','75','100','25','50','75','100'])
# plt.ylim(0,4)
# plt.ylabel('LM/HT')
# #
# plt.subplot(3,1,3)
# plt.axhline(1, color='k', linewidth=0.5,zorder=0)
# box3 = plt.boxplot(ratio_PI_T, notch=True, patch_artist=True, sym='')
# for patch in box3['boxes']:
#     patch.set_facecolor(m.to_rgba(-1000))
# plt.xticks(range(1,13), ['25','50','75','100','25','50','75','100','25','50','75','100','25','50','75','100'])
# plt.ylim(0,4)
# plt.ylabel('PI/HT')
# plt.xlabel('Sampling length')

# # colorbar
# fig.subplots_adjust(right=0.8, top=0.9)
# cb_ax = fig.add_axes([0.83, 0.05,0.03, 0.85])
# cb = cbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, ticks=np.arange(-3,3.5,0.5), orientation='vertical')
# cb.set_label(r'$log_{10}(0.05/p_{val})$')

# # save plot
# # plt.tight_layout()
# chdir('../historical')
# plt.savefig('var_ratio.pdf')
