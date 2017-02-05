from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

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
import re
import filter
import glob

global colors

colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']
state_colors = ['violet','lightseagreen','blue','green','red','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']
states = ['synchronized', 'mixed', 'desynchronized']

track_name = 'ptc21_tr5c'

if track_name == 'ptc17_tr2b':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'

if track_name == 'ptc18_tr1':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc18/'

if track_name == 'ptc20_tr1':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc20/'

if track_name == 'ptc20_tr3':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc20/'

if track_name == 'ptc21_tr5c':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

if track_name == 'ptc22_tr1':
    sim_dir = '/media/cat/4TB/in_vivo/nick/ptc22/'


print "LOADING SUA SPIKES"

#******** Load particular recordings
if track_name == 'ptc17_tr2b':
    recordings = ['44','45','46','47','48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73']
    state =      [ 1,   1,   1,   1,   0.5,  0.5,  1,    0.5,  1,    1,    1,    1,    1,    1,    0.5,  0,    .5,   1,    1,    1,    1,    1,    1,    0.,   0.5,  1,    1,    1,    1,    1]

if track_name == 'ptc18_tr1':
    recordings = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44']
    state =      [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0.5, 0,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    0.5,  0,    0,    1,    0,    0,    0  ]

if track_name == 'ptc20_tr1':
    recordings = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33']
    state =      [ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,    1,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1  ]

if track_name == 'ptc20_tr3':
    recordings = ['58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72']
    state =      [ 0,    1,    1,    1,    1,    1,    0.5,  1,    1,    1,    0.5,  0.5,  1,    1,    1  ]

if track_name == 'ptc21_tr5c':
    recordings = ['59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75']
    state =      [0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  ]

if track_name == 'ptc22_tr1':
    recordings = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21', '22']
    state =      [ 1,   1,   0,   0.5, 1,   1,   1,   0.5, 0,   0.5, 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1  ]



#******** Concatenate recordings 
#SUA_sorts, rasters = Concatenate_ptcs(sim_dir, track_name, recordings)

#******** Vectorize data into Schneidman vectors
bin_width = 0.020 #Bin width in seconds
rec = 1
actual_state = states[int(state[rec]*2)]

print "Recording: ", recordings[rec]
temp_name = glob.glob(sim_dir+recordings[rec]+"*")[0]
file_name = temp_name.replace(sim_dir,'')

vectors, Sort, last_spike = Schneidman_vectors(sim_dir, file_name, bin_width)
print vectors.shape
rec_length = last_spike/Sort.samplerate
print "Length rec: ", rec_length , " sec  = ", rec_length/60., " mins."

for spike_match in [2,3,4]:
    for k in range(10):
        print "# spikes match: ", spike_match, "  vector: ", k
        #Find a random time_bin w. required # of spikes 
        while True: 
            ind = np.random.randint(0, len(vectors[0]))
            v_sum = sum(vectors[:, ind])
            if v_sum == spike_match: 
                #print "looking for index: ", ind
                #print vectors[:, ind]
                n_matches = 0
                for p in range(len(vectors[0])):
                    if np.array_equal(vectors[:,p], vectors[:, ind]):
                        n_matches+=1
                
                if n_matches>1: break

        #Identify units that fired during time_bin
        indexes = np.where(vectors[:,ind]>0)[0]
        indexes_values = vectors[:,ind][indexes]

        #Compute firing rates for all units:
        rasters = []
        fire_rate = []
        for k in range(len(Sort.units)):
            temp_index = np.where(np.array(Sort.units[k])<1E10)[0] #Exclude very large spikes due to SS .ptcs spike time export error
            rasters.append(np.array(Sort.units[k])[temp_index])
            fire_rate.append(len(rasters[k])/(last_spike/Sort.samplerate))

        #Find indep_prob of firing
        indep_prob = 1
        for p in range(len(indexes)):
            for i in range(indexes_values[p]):
                f_rate = fire_rate[indexes[p]]
                p_rate = f_rate *bin_width  #Prob of a spike in a bin
                indep_prob = indep_prob * p_rate
                
        #print "Indep prob rate: ", indep_prob

        #Look for matches in all time bins
        #n_matches = 0
        #for p in range(len(vectors[0])):
            #if np.array_equal(vectors[:,p], vectors[:, ind]):
                #n_matches+=1

        actual_prob = float(n_matches)/len(vectors[0])
        #print "Actual occurence rate: ", actual_prob
        
        plt.scatter(actual_prob, indep_prob, color=state_colors[spike_match])

plt.plot([0,1],[0,1],  color='black', linewidth=4)

plt.xscale('log')
plt.yscale('log')
plt.ylim(1E-13,1E0)
plt.xlim(1E-5, 1E0)

plt.ylabel("Independent Assumption Probability", fontsize = 20)
plt.xlabel("Actual Occurence Rate", fontsize = 20)
plt.suptitle(track_name + " rec: " + recordings[rec] + "     cortical state: " + actual_state, fontsize = 30)
plt.show()

quit()

































#******** Compute #Spikes in each bin pattern 
ax1 = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)


bin_distribution = []
for k in range(len(vectors[0])):
    bin_distribution.append(np.sum(vectors[:,k]))

yy = np.histogram(bin_distribution, bins = np.arange(0,20,1))

ax1.plot(yy[1][:-1], np.float32(yy[0])/max(yy[0]), color=state_colors[int(state[r]*2)], linewidth = 3)
plt.xlim(0,19)

ax2.plot(yy[1][1:-1], np.float32(yy[0][1:])/max(yy[0][1:]), color=state_colors[int(state[r]*2)], linewidth = 3)
plt.xlim(0,19)

ax3.plot(yy[1][2:-1], np.float32(yy[0][2:])/max(yy[0][2:]), color=state_colors[int(state[r]*2)], linewidth = 3)
plt.xlim(0,19)
        
plt.suptitle(track_name + ",   "  + " # recs: "+ str(len(recordings)) + "\nDistribution of number of simultaneous spikes (i.e. length of spike chain in bin)", fontsize=25)
plt.xlabel("Number of spikes in time bin (20ms)", fontsize=20)
plt.show()










    #Plot distribution of correlation coefficients - TRY TWO DIFF CORTICAL STATES - ALSO CHECK MOUSE DATA










