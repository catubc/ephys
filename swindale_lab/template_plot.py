from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *
from random import randrange

from scipy import stats
import numpy as np
import time, math
import sys
import os.path
from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab

from util import xcorrm
import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import pandas as pd  
import matplotlib.pyplot as plt  
from scipy.stats import sem  

colors = np.array(['blue', 'red', 'green', 'orange', 'magenta', 'cyan'])


#***************************** START FUNCTION **********************************

sim_dir = '/media/cat/4TB/in_vivo/nick/ptc21/67-tr5c-blankscreen/'
tsf_name = '67-tr5c-blankscreen_scaled_notch_multiunit'
sorted_file = '67-tr5c-blankscreen_scaled_notch_multiunit_5clusters'


#********************* READ TSF FILES *******************
f = open(sim_dir + tsf_name + '.tsf', "rb")
print "Loading ", sim_dir + tsf_name + '.tsf'
tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
tsf.sim_dir = sim_dir
tsf.tsf_name = tsf_name
f.close()

ptcs_flag = 1

Sort1 = Loadptcs(sorted_file, sim_dir, ptcs_flag, save_timestamps=False)
Sort1.name=sorted_file
Sort1.filename=sorted_file
Sort1.rawfile= tsf_name
Sort1.directory=sim_dir



#********************** PLOT TEMPLATES WITH STANDARD DEVIATIONS************************

if False:
    ax = plt.subplot(1,1,1)

    width = 120
    x = np.arange(0,width,1)
    unit_color=0
    for unit in Sort1.units:
        spike_traces = np.zeros((len(unit), tsf.n_electrodes, width), dtype=np.float64)
        counter=0
        for k in unit:
            for i in range(tsf.n_electrodes):  
                #plt.plot(x, tsf.ec_traces[i][k-width/2:k+width/2]-i*40000, color=colors[counter])
                spike_traces[counter][i]=tsf.ec_traces[i][k-width/2:k+width/2]
            counter+=1
        
        spike_sd = np.zeros((tsf.n_electrodes, width), dtype=np.float64)
        spike_mean = np.zeros((tsf.n_electrodes, width), dtype=np.float64)
        for i in range(tsf.n_electrodes):
            for t in range(width):
                spike_sd[i][t] = np.std(spike_traces[:,i,t])
                spike_mean[i][t] = np.mean(spike_traces[:,i,t])
                
        time = np.arange(0, width,1) #X-axis has lable in ms.

        # Use matplotlib's fill_between() call to create error bars.  
        # Use the dark blue "#3F5D7D" as a nice fill color.  
        for i in range(tsf.n_electrodes):
            plt.fill_between(time+width*1.5*unit_color, (spike_mean[i] - spike_sd[i])-i*20000,  
            (spike_mean[i]+spike_sd[i])-i*20000, color=colors[unit_color])  
            plt.plot(time+width*1.5*unit_color,spike_mean[i]-i*20000, color="black", lw=2)  
        unit_color+=1

    ax.set_xlim(left=-50)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()


#********************** PLOT SINGLE MULTI CHANNEL SPIKES ************************


width = 120
x = np.arange(0,width,1)
unit_color=0
for unit in Sort1.units:
    counter=0
    for k in unit:
        ax = plt.subplot(1,1,1)
        for i in range(tsf.n_electrodes):  
            plt.plot(x, tsf.ec_traces[i][k-width/2:k+width/2]-i*40000, color=colors[counter], linewidth=2, alpha=.1)

        ax.set_xlim(-width*3, width*4)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        #mng = plt.get_current_fig_manager()
        #mng.resize(*mng.window.maxsize())
    plt.show()
    counter+=1

