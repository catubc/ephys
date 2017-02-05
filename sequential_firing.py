import struct, array, csv
import numpy as np
import os
import math, operator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab as pl
import time
from tsf_ptcs_classes import *
import random
import copy
import glob
from libtfr import *
import itertools


from scipy import stats
from scipy.optimize import curve_fit
import scipy.stats as sstats
from scipy.stats import chisquare

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from numpy import exp

from math import sin, pi
from matplotlib.patches import Circle # for simplified usage, import this patch

#from sklearn.decomposition import PCA
from sklearn import decomposition
from sklearn import datasets


def Plot_rasters(Sort1, Sort2):

    #****************** PLOT RASTERS BY SEQUENCE*********************
    if True:
        ax = plt.subplot(1, 1, 1)

        bin_width = 0.250 #Bin width in seconds
        y = []
        for i in range(Sort1.n_units):
            x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5
    
            ymin=np.zeros(len(Sort1.units[i]))
            ymax=np.zeros(len(Sort1.units[i]))
            ymin+=i+0.4
            ymax+=i-0.4
    
            plt.vlines(x, ymin, ymax, linewidth=2, color='black') #colors[mod(counter,7)])

            y.append(x)

        #Plot MUA
        y = np.hstack(y)
        y1 = np.histogram(y, bins = np.arange(0,max(y),bin_width))
        plt.plot(y1[1][:-1],y1[0]/10.-4, color='red', linewidth=2)

        #Plot LFP spike
        x = np.array(Sort2.units[0],dtype=float32)/float(Sort2.samplerate)*100#float(Sort1.samplerate)*2.5
        ymin=np.zeros(len(Sort2.units[0]))
        ymax=np.zeros(len(Sort2.units[0]))
        ymin+=-5.
        ymax+=-8
        plt.vlines(x, ymin, ymax, linewidth=2, color='blue') #colors[mod(counter,7)])
        

        
        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Single Unit ID',multialignment='center', fontsize=15)

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
        plt.xlim(left=-0.03)

        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Single Unit ID',multialignment='center', fontsize=15)
    
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()


    #****************** PLOT RASTERS BY DEPTH*********************
    if True:
        ax = plt.subplot(1, 1, 1)
        title("Rasters for: "+ Sort1.name, multialignment='center',fontsize=15)
    
        colors = ['blue','red', 'green', 'cyan', 'black', 'magenta', 'yellow']
    
        counter = 0 #These 2 variables used to track and color units that land on the same max channels
        depth = 0   #by cycling through the colors array above
        for i in range(Sort1.n_units):
            x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5
    
            ymin=np.zeros(len(Sort1.units[i]))
            ymax=np.zeros(len(Sort1.units[i]))
            ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*5+10.0
            ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*5-10.0
    
            #Check to see if previous cell at same depth
            if depth == (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*5):
                counter+=1
            else:
                counter=0
            plt.vlines(x, ymin, ymax, linewidth=2, color='black') #colors[mod(counter,7)])
    
            depth = (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*5)
    
        #ax.xaxis.grid()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
        plt.xlim(left=-0.03)
    
        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)
    
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()
    
    
    #**************** PLOT MUA FIRING RATES ********************    
    if True:
        ax = plt.subplot(1, 1, 1)
        bin_width = 0.250 #Bin width in seconds
    
        title("MUA Firing rates over entire recording (unsorted; bin width = " + str(bin_width) + " sec)", multialignment='center',fontsize=15)
    
        for i in range(Sort1.n_units):
            y = np.array(Sort1.units[i])/float(Sort1.samplerate)
            y = np.histogram(y, bins = np.arange(0,250.0,bin_width))[0]
            y_offset = -Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*40.0
            plt.plot(y+y_offset, color='blue')
        
        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)
    
        ax.xaxis.grid()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
        plt.xlim(left=-0.03)
        
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()
    
    #**************** PLOT TRIAL AVERAGED MUA FIRING RATES ********************    
    if True:
        ax = plt.subplot(1, 1, 1)

        #Load stimulus times - SEV DATA
        file_name = Sort1.directory + 'DeepCSD4.SweepTime' #ofname+exp_name+'/'+str(i)+'/'+str(i)+'.SweepTime.mat'
        mat = scipy.io.loadmat(file_name)
        stim_times = np.array(mat[mat.keys()[1]][0][0][4][0][1], dtype=np.float32)/Sort1.samplerate  
        print stim_times
        
        triggered_spikes=[]
        firing_rate = []
        #Find all spikes in 1.0 second window after stimulus
        for i in range(Sort1.n_units):
            triggered_spikes.append([])
            temp_array = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #WORK IN SECONDS TIMESCALE
            #Loop over all stimulus_periods
            for j in range(len(stim_times)):
                #Find all unit spikes that fall w/in +/- window of pop spike
                temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][0]+1.0))[0]
                
                triggered_spikes[i].extend(temp_array[temp2]-stim_times[j][0]) #Offset time of spike to t=t_0 
                #save_response_array_individual[i].append(temp_array[temp2]-stim_times[j][0])

        #Compute ave. firing rate for each channel
        bin_width = 0.010
        title("MUA Firing rates trial averaged (unsorted; bin width = " + str(bin_width) + " sec)", multialignment='center',fontsize=15)

        for i in range(len(Sort1.units)):
            x = np.arange(0,1.0,bin_width)
            y = np.histogram(triggered_spikes[i], bins = x)[0]/2.
            y_offset = -Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/22*5.0
            plt.plot(x[:-1],y+y_offset, color='blue', linewidth=1.5)
        
        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)

        ax.xaxis.grid()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

        plt.xlim(left=0.0) #, right = 0.250)
        
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()
    

def Plot_rasters_SUA_vs_LFP(Sorts, Sort2):
    
    if True:
        ax = plt.subplot(1, 1, 1)

        bin_width = 0.250 #Bin width in seconds
        y = []
        Sort1=Sorts[1]
        for i in range(len(Sorts)):
            x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5
    
            ymin=np.zeros(len(Sort1.units[i]))
            ymax=np.zeros(len(Sort1.units[i]))
            ymin+=i+0.4
            ymax+=i-0.4
    
            plt.vlines(x, ymin, ymax, linewidth=4, color='black') #colors[mod(counter,7)])

            y.append(x)

        #Plot LFP spike
        x = np.array(Sort2.units[0],dtype=float32)/float(Sort2.samplerate)*100#float(Sort1.samplerate)*2.5
        ymin=np.zeros(len(Sort2.units[0]))
        ymax=np.zeros(len(Sort2.units[0]))
        ymin+=-5.
        ymax+=-8
        plt.vlines(x, ymin, ymax, linewidth=4, color='blue') #colors[mod(counter,7)])
        
        plt.xlabel('Time (seconds)',fontsize=35, weight='bold')
        plt.ylabel('Single Unit ID',multialignment='center', fontsize=35)

        #ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
        plt.xlim(0,300)

        plt.ylabel('LFP Cluster Raster         Single Unit IDs',multialignment='center', fontsize=35, weight='bold')
        
        ax.tick_params(axis='both', which='major', labelsize=30)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

def Compute_binned_lfp(lfp_data, stim_times,spont_activity_length, bin_width, lowcut, highcut):
    '''Assume sample rate is 1000Hz, i.e. 1 sample / ms.; SHOULD BE GENERALIZED'''
    #Compute_binned_lfp(lfp_data, stim_times[i], spont_activity_length, bin_width, lowcut, highcut))
    
    print "Computing binned lfp vectors"
    
    #Filtering info:
    fs = 1000 #Assume this for filter
    lowcut_array=[lowcut]
    highcut_array=[highcut]
    
    #Look at cumulative firing rate for all trials
    #vector_out = []
    print "Band: ", lowcut, highcut

    for lowcut, highcut in zip(lowcut_array, highcut_array):

        #Filter each frame/segment of data
        lfp_data_filtered=[]
        for e in range(len(lfp_data)):
            lfp_data_filtered.append(butter_bandpass_filter(lfp_data[e], lowcut, highcut, fs, order = 2))

        #Loop over data in bins
        temp = []
        for i in range(int(stim_times[0]*1E3), int((stim_times[1]+spont_activity_length)*1E3), bin_width):  
            lfp_val = []

            for e in range(len(lfp_data)):  #Loop over channels
                chunk = lfp_data_filtered[e][i:(i+1)*bin_width]
                lfp_val.append(np.average(chunk))

            temp.append(lfp_val)

        vector_out= np.array(temp)
    #vector_out = np.array(vector_out)

    print vector_out.shape

    return vector_out


def cluster_matplotlib(data, fname):
    from mpl_toolkits.mplot3d import Axes3D

    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n = 100
    c='r'
    m='^'
    xs=[]
    ys=[]
    zs=[]
    print data.shape
    for i in range(len(data)):    
        xs.append(data[i][0])
        ys.append(data[i][1])
        zs.append(data[i][2])

    #KMEANS
    if False: 
        labels = KMEANS(data, 10)

    #MEAN SHIFT
    if True:
        from sklearn.cluster import MeanShift, estimate_bandwidth
        from sklearn.datasets.samples_generator import make_blobs
        
        # The following bandwidth can be automatically detected using
        bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=500)

        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(data)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        print("number of estimated clusters : %d" % n_clusters_)

    print labels.shape

    print np.unique(labels, return_counts=True)

    clrs = []
    xs=[]
    ys=[]
    zs=[]

    import matplotlib as mpl
    cmap = mpl.cm.get_cmap('jet')

    n_clusters = 20 #n_clusters_

    clrs_rand = []
    for i in range(n_clusters):
        clrs_rand.append(np.random.rand(3,))
        
    for i in range(len(labels)):
        if labels[i]<n_clusters:
            #print cmap(float(labels[i])/len(labels))
            
            clrs.append(clrs_rand[labels[i]])

            xs.append(data[i][0])
            ys.append(data[i][1])
            zs.append(data[i][2])
        
    ax.scatter(xs, ys, zs, c=clrs, marker=m,s=300)
    plt.show()

    clusters = []
    file_name = fname+"_clusters.txt"
    for i in range(n_clusters):
        clusters.append(np.where(labels==i)[0])

    import cPickle
    cPickle.dump(clusters, open(file_name, 'wb')) 

    return clusters



def Plot_rasters_and_stimulus(Sort1, stim_times):

    #*********************** VLINE SETS OF PLOTS *****************************
    ax = plt.subplot(1, 1, 1)
    title("Original file: "+ Sort1.name, multialignment='center',fontsize=15)

    colors = ['blue','red', 'green', 'cyan', 'black']

    counter = 0 #These 2 variables used to track and color units that land on the same max channels
    depth = 0   #by cycling through the colors array above
    
    for i in range(Sort1.n_units):

        x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5

        ymin=np.zeros(len(Sort1.units[i]))
        ymax=np.zeros(len(Sort1.units[i]))
        ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5
        ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5-5.0

        if depth == (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5): 
            counter+=1
        else:
            counter=0
        plt.vlines(x, ymin, ymax, linewidth=2, color=colors[mod(counter,5)])

        depth = (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5)

    for i in range(len(stim_times)):
        p = ax.axvspan(stim_times[i][0],stim_times[i][1], facecolor='red', alpha=0.15)

    ax.xaxis.grid()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.xlim(left=-0.03)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def Plot_stimulus_triggered_rasters (Sort1, stim_times):

    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    #********************************************************
    #Make histograms -Xms .. +Xms from population spikes
    min_spikes = -1          #min no. of spikes in window
    
    save_response_array=[]
    save_response_array_individual=[]
    for i in range(Sort1.n_units):
        save_response_array.append([])
        save_response_array_individual.append([])

    
    n_repeats = len(stim_times)
    freq=[]
    all_spikes_rasters=[]
    non_locked_spikes_rasters=[]
    locked_spikes_rasters=[]

    #*********** LFP EVENT RESPONSES ***********
    if True:
        lfp_events = np.load(Sort1.directory+Sort1.filename+'/'+Sort1.filename+'.lfp.events.npy')
        locked_lfp_events=[]
        for i in range(len(lfp_events)):
            #Search non_locked_spikes
            locked_lfp_events.append([])
            temp_array = np.array(lfp_events[i],dtype=float32)*1E-3
            
            #Loop over all stimulus_periods
            for j in range(n_repeats):     #Find all unit spikes that fall w/in +/- window of pop spike
                temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][1]))[0]
                locked_lfp_events[i].append(temp_array[temp2]-stim_times[j][0])

        len_recording = 1000

        #Plot LFP raster responses:
        sig = 0.020 #sd of convolved gaussian

        for j in range(len(lfp_events)):    #Loop over LFP clusters
            ax=plt.subplot(2,3,j+1)
            for i in range(n_repeats):
                plt.scatter(locked_lfp_events[j][i],[i]*len(locked_lfp_events[j][i]), s = 10, color=colors[j])
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(n_repeats,-0.5)
            plt.title("lfp_spikes #: "+str(len(lfp_events[j])))

            xx = np.linspace(0, 5.0, 500)
            fit = np.zeros(500, dtype=np.float32)
            for i in range(len(locked_lfp_events[j])):
                for k in range(len(locked_lfp_events[j][i])):
                    mu = np.array(locked_lfp_events[j][i][k]) 
                    fit += gaussian(xx, mu, sig)
            plt.plot(xx, -fit/(np.max(fit))*n_repeats+n_repeats, linewidth=4, color=colors[j], alpha=1)
            #mua_bin_width = 0.025
            #mua = list(itertools.chain.from_iterable(locked_lfp_events[j]))
            #mua_plot=np.histogram(mua, bins = np.arange(0,len_recording,mua_bin_width))
            #plt.plot(mua_plot[1][0:-1],-mua_plot[0]*10+n_repeats, linewidth=4, color=colors[j], alpha=1)
        #plt.suptitle(Sort1.filename+"    LFP Event Rasters   -  "+str(int(mua_bin_width*1E3))+"ms bins", fontsize = 20)
        plt.suptitle(Sort1.filename+"    LFP Event Rasters   -  (20ms gaussian)", fontsize = 20)
        plt.show()
        quit()
    
    #************** SINGLE UNIT RESPONSES *************
    if False:
        #Loop over all units in sorted data - look for spikes in each 4.5 sec movie repeat
        non_locked_spikes = np.load(Sort1.directory+Sort1.filename+'/'+Sort1.filename+'_non_locked_spikes.npy')
        locked_spikes = np.load(Sort1.directory+Sort1.filename+'/'+Sort1.filename+'_locked_spikes.npy')

        for i in range(Sort1.n_units):
            
            #Search all_spikes
            all_spikes_rasters.append([])
            temp_array = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)
            #Loop over all stimulus_periods
            for j in range(n_repeats):     #Find all unit spikes that fall w/in +/- window of pop spike
                temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][1]))[0]
                all_spikes_rasters[i].append(temp_array[temp2]-stim_times[j][0])

            #Search non_locked_spikes
            non_locked_spikes_rasters.append([])
            temp_array = np.array(non_locked_spikes[i],dtype=float32)
            #Loop over all stimulus_periods
            for j in range(n_repeats):     #Find all unit spikes that fall w/in +/- window of pop spike
                temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][1]))[0]
                non_locked_spikes_rasters[i].append(temp_array[temp2]-stim_times[j][0])

            #Search non_locked_spikes
            locked_spikes_rasters.append([])
            temp_array = np.array(locked_spikes[i],dtype=float32)
            #Loop over all stimulus_periods
            for j in range(n_repeats):     #Find all unit spikes that fall w/in +/- window of pop spike
                temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][1]))[0]
                locked_spikes_rasters[i].append(temp_array[temp2]-stim_times[j][0])


        #Plot movie raster responses:
        for j in range(Sort1.n_units):
            n_repeats = len(all_spikes_rasters[j])
            
            ax=plt.subplot(2,3,1)
            for i in range(n_repeats):
                plt.scatter(all_spikes_rasters[j][i],[i]*len(all_spikes_rasters[j][i]), s = 10, color='blue')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(n_repeats,-0.5)
            plt.title("all_spikes #: "+str(len(Sort1.units[j])))
            
            mua_bin_width = 0.020
            mua = list(itertools.chain.from_iterable(all_spikes_rasters[j]))
            mua_plot_all_spikes=np.histogram(mua, bins = np.arange(0,4.6,mua_bin_width))
            
            ax=plt.subplot(2,3,2)
            for i in range(n_repeats):
                plt.scatter(non_locked_spikes_rasters[j][i],[i]*len(non_locked_spikes_rasters[j][i]), s = 10, color='red')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(n_repeats,-0.5)
            plt.title("non_locked_spikes #: "+str(len(non_locked_spikes[j])))
            
            mua = list(itertools.chain.from_iterable(non_locked_spikes_rasters[j]))
            mua_plot_non_locked=np.histogram(mua, bins = np.arange(0,4.6,mua_bin_width))
                    
            ax=plt.subplot(2,3,3)
            for i in range(n_repeats):
                plt.scatter(locked_spikes_rasters[j][i],[i]*len(locked_spikes_rasters[j][i]), s = 10, color='green')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(n_repeats,-0.5)
            plt.title("locked_spikes #: "+str(len(locked_spikes[j])))
            
            mua = list(itertools.chain.from_iterable(locked_spikes_rasters[j]))
            mua_plot_locked=np.histogram(mua, bins = np.arange(0,4.6,mua_bin_width))

            #Find max y_lim
            y_max = max(np.max(mua_plot_all_spikes[0]), np.max(mua_plot_non_locked[0]), np.max(mua_plot_locked[0]))

            ax=plt.subplot(2,3,4)
            plt.plot(mua_plot_all_spikes[1][0:-1],mua_plot_all_spikes[0], linewidth=3, color='blue')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(0, y_max)

            ax=plt.subplot(2,3,5)
            plt.plot(mua_plot_non_locked[1][0:-1],mua_plot_non_locked[0], linewidth=3, color='red')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(0, y_max)

            ax=plt.subplot(2,3,6)
            plt.plot(mua_plot_locked[1][0:-1],mua_plot_locked[0], linewidth=3, color='green')
            plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            plt.ylim(0, y_max)
            
            plt.show()
    
    return

def Compute_responses(Sort, stim_times , spont_activity_length):

    #Bin across all movie repeats
    min_spikes = -1          #min no. of spikes in window

    save_response_array_individual=[]
    for i in range(Sort.n_units):
        save_response_array_individual.append([])

    n_repeats = len(stim_times)
    freq=[]
    #Loop over all units in sorted data - look for spikes in each 4.5 sec movie repeat
    for i in range(Sort.n_units):
        temp_array = np.array(Sort.units[i],dtype=float32)/float(Sort.samplerate)

        #Loop over all stimulus_periods
        for j in range(n_repeats):
            
            #Find all unit spikes that fall w/in +/- window of pop spike; add spont activity at length of 5s nat scenes
            temp2 = np.where(np.logical_and(temp_array>=stim_times[j][0], temp_array<=stim_times[j][1]+spont_activity_length))[0]
            #print temp2
            
            if len(temp2)>min_spikes: #Min spikes in window 
                #if len(temp2)==0: print temp_array[temp2]-stim_times[j][0]
                save_response_array_individual[i].append(temp_array[temp2]-stim_times[j][0])

    return save_response_array_individual

def PSTH_plot(millisecs_load, bin_width, save_response_array_individual, work_dir, file_name, movie_frames):

    total_rate = np.zeros(millisecs_load/bin_width+1, dtype=np.float32) #Make 1ms wide raster;
    for i in range(len(save_response_array_individual)):
        print "cell: ", i
        temp = np.zeros(millisecs_load/bin_width+1, dtype=np.float32) #Initialize firing rate vector
        all_spikes = np.sort(list(itertools.chain(*save_response_array_individual[i]))) #Load spikes from repeats - flatten 2D array
        
        for s in all_spikes: #Sort.units[i]:
            s = s *1E3
            bin_loc = int(s / bin_width) #Show both division and multiplication
            temp[bin_loc]+=1

        total_rate+=temp

        xx = np.arange(0,len(temp),1)*bin_width
        plt.plot(xx,temp/-50.+i,  linewidth=2, color='black') #Plot firing rates upside down
        #plt.plot(temp/-50.+i,  linewidth=2, color='black') #Plot firing rates upside down

    for j in range(0, 5500, 500):
        plt.plot([j,j],[i+10,0],'r--', color='black')
        
    xx = np.arange(0,len(total_rate),1)*bin_width
    plt.plot(xx, total_rate/-100.+ i+10, linewidth=2, color='red')
    plt.ylim(i+10,0)
    plt.xlim(0,5500)
    plt.fill_between([4480,5500], i+10,0, facecolor='blue', alpha=0.3)

    plt.ylabel("Unit # (Depth Ordered)", fontsize=30)
    plt.xlabel("Milliseconds (bin width:"+str(bin_width)+" ms)", fontsize=30)
    plt.title(work_dir + file_name +'\n'+"All repeats"+str(movie_frames[0])+'..'+str(movie_frames[-1]), fontsize=30)
    plt.show()



def LFP_PSTH_plot(lfp_matrix, bin_width, work_dir, file_name, movie_frames, lowcut, highcut, spont_activity_length):
    print "Plotting LFP repeats"

    channel_matrix = np.swapaxes(lfp_matrix,0,2)
    channel_matrix = np.swapaxes(channel_matrix,1,2)
    
    xx = np.linspace(0,len(lfp_matrix[0]),1)*bin_width #Convert to ms

    ax = plt.subplot(111)

    total_rate = np.zeros(len(xx), dtype=np.float32) 
    #Loop over each channel & compute average values
    for i in range(10):
        print "channel: ", i
        #temp = np.zeros(millisecs_load/bin_width, dtype=np.float32) #Initialize firing rate vector
        #all_spikes = np.sort(list(itertools.chain(*save_response_array_individual[i]))) #Load spikes from repeats - flatten 2D array
        
        channel_arrays = channel_matrix[i,:,:]
        temp_ave = np.average(channel_arrays,axis=0)
        temp_std = np.std(channel_arrays, axis=0)
        #print temp_array
        #for frame in range(len(channel_matrix[i])): #Sort.units[i]:
        #    
        #    s = s *1E3
        #    bin_loc = int(s / bin_width) #Show both division and multiplication
        #    temp[bin_loc]+=1

        #total_rate+=temp
        xx = np.arange(0,len(temp_ave),1)*bin_width
        plt.plot(xx, temp_ave/50.+i,  linewidth=2, color='black') #Plot firing rates upside down
        ax.fill_between(xx, temp_ave/50.+i+temp_std/50., temp_ave/50.+i-temp_std/50., facecolor='red', alpha=0.15)
        
            
    #plt.plot(xx, total_rate/-10.+ i+10, linewidth=2, color='red')
    plt.ylim(-1,11)
    plt.xlim(0,5500)
    #plt.plot([4480,4480],[-1,11], color='blue', linewidth=3)
    ax.fill_between([4480,5500], -1,11, facecolor='blue', alpha=0.3)

    plt.ylabel("LFP Channel # (Exported Order)", fontsize=30)
    plt.xlabel("Time - Milliseconds (bin width:"+str(bin_width)+" ms)", fontsize=30)
    plt.title(work_dir + file_name +'\n'+"movie repeats: "+str(movie_frames[0])+'..'+str(movie_frames[-1])+
    ",   band: "+str(lowcut)+'..'+str(highcut)+"hz", fontsize=30)
    plt.show()
    quit()

def multi_dim_scaling(matrix_in, method):

    import sklearn 
    from sklearn import manifold
    
    methods = ['MDS - SMACOF', 't-SNE', 'PCA', 'Sammon']
    print "Computing dim reduction, size of array: ", np.array(matrix_in).shape
    
    print "... running dim reduction: ", methods[method], " ..."
    
    if method==0:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        print "... pairwise dist ..."
        dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
        adist = np.array(dists)
        amax = np.amax(adist)
        adist /= amax
        
        print "... computing MDS ..."
        mds_clf = manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
        results = mds_clf.fit(adist)
        Y = results.embedding_*1E3
        
    elif method==1:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        print "... pairwise dist ..."
        dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
        adist = np.array(dists)
        amax = np.amax(adist)
        adist /= amax
        
        print "... computing tSNE ..."
        tsne = manifold.TSNE(n_components=3, init='pca', random_state=0)
        Y = tsne.fit_transform(adist)*1E2
    
    elif method==2:
        #PCA
        Y, X = PCA(matrix_in, 3)
        Y*=1E2
    
    elif method==3:
        #Sammon port from Matlab
        print "NOT IMPLEMENTED"
        
        sammon(matrix_in)    

    return Y
    
    ##Locally Linear Embedding Methods - NOT GREAT ON THIS PROBLEM
    #methods = ['standard', 'ltsa', 'hessian', 'modified']
    #Y = manifold.LocallyLinearEmbedding(10, 3,
                                        #eigen_solver='auto',
                                        #method=methods[3]).fit_transform(adist)

def sammon(x, n = 3, display = 2, inputdist = 'raw', maxhalves = 20, maxiter = 500, tolfun = 1e-9, init = 'pca'):

    import numpy as np 

    """Perform Sammon mapping on dataset x
    y = sammon(x) applies the Sammon nonlinear mapping procedure on
    multivariate data x, where each row represents a pattern and each column
    represents a feature.  On completion, y contains the corresponding
    co-ordinates of each point on the map.  By default, a two-dimensional
    map is created.  Note if x contains any duplicated rows, SAMMON will
    fail (ungracefully). 
    [y,E] = sammon(x) also returns the value of the cost function in E (i.e.
    the stress of the mapping).
    An N-dimensional output map is generated by y = sammon(x,n) .
    A set of optimisation options can be specified using optional
    arguments, y = sammon(x,n,[OPTS]):
       maxiter        - maximum number of iterations
       tolfun         - relative tolerance on objective function
       maxhalves      - maximum number of step halvings
       input          - {'raw','distance'} if set to 'distance', X is 
                        interpreted as a matrix of pairwise distances.
       display        - 0 to 2. 0 least verbose, 2 max verbose.
       init           - {'pca', 'random'}
    The default options are retrieved by calling sammon(x) with no
    parameters.
    File        : sammon.py
    Date        : 18 April 2014
    Authors     : Tom J. Pollard (tom.pollard.11@ucl.ac.uk)
                : Ported from MATLAB implementation by 
                  Gavin C. Cawley and Nicola L. C. Talbot
    Description : Simple python implementation of Sammon's non-linear
                  mapping algorithm [1].
    References  : [1] Sammon, John W. Jr., "A Nonlinear Mapping for Data
                  Structure Analysis", IEEE Transactions on Computers,
                  vol. C-18, no. 5, pp 401-409, May 1969.
    Copyright   : (c) Dr Gavin C. Cawley, November 2007.
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    """

    def euclid(a,b):
        d = np.sqrt( ((a**2).sum(axis=1)*np.ones([1,b.shape[0]]).T).T + \
            np.ones([a.shape[0],1])*(b**2).sum(axis=1)-2*(np.dot(a,b.T)))
        return d

    # Create distance matrix unless given by parameters
    if inputdist == 'distance':
        D = x
    else:
        D = euclid(x,x)

    # Remaining initialisation
    N = x.shape[0] # hmmm, shape[1]?
    scale = 0.5 / D.sum()
    D = D + np.eye(N)
    Dinv = 1 / D # Returns inf where D = 0.
    Dinv[np.isinf(Dinv)] = 0 # Fix by replacing inf with 0 (default Matlab behaviour).
    if init == 'pca':
        [UU,DD,_] = np.linalg.svd(x)
        y = UU[:,:n]*DD[:n] 
    else:
        y = np.random.normal(0.0,1.0,[N,n])
    one = np.ones([N,n])
    d = euclid(y,y) + np.eye(N)
    dinv = 1. / d # Returns inf where d = 0. 
    dinv[np.isinf(dinv)] = 0 # Fix by replacing inf with 0 (default Matlab behaviour).
    delta = D-d 
    E = ((delta**2)*Dinv).sum() 

    # Get on with it
    for i in range(maxiter):
        print "Sammon method iter: ", i
        # Compute gradient, Hessian and search direction (note it is actually
        # 1/4 of the gradient and Hessian, but the step size is just the ratio
        # of the gradient and the diagonal of the Hessian so it doesn't
        # matter).
        delta = dinv - Dinv
        deltaone = np.dot(delta,one)
        g = np.dot(delta,y) - (y * deltaone)
        dinv3 = dinv ** 3
        y2 = y ** 2
        H = np.dot(dinv3,y2) - deltaone - np.dot(2,y) * np.dot(dinv3,y) + y2 * np.dot(dinv3,one)
        s = -g.flatten(order='F') / np.abs(H.flatten(order='F'))
        y_old    = y

        # Use step-halving procedure to ensure progress is made
        for j in range(maxhalves):
            s_reshape = s.reshape(2,len(s)/2).T
            y = y_old + s_reshape
            d = euclid(y, y) + np.eye(N)
            dinv = 1 / d # Returns inf where D = 0. 
            dinv[np.isinf(dinv)] = 0 # Fix by replacing inf with 0 (default Matlab behaviour).
            delta = D - d
            E_new = ((delta**2)*Dinv).sum()
            if E_new < E:
                break
            else:
                s = np.dot(0.5,s)

        # Bomb out if too many halving steps are required
        if j == maxhalves:
            print('Warning: maxhalves exceeded. Sammon mapping may not converge...')

        # Evaluate termination criterion
        if np.abs((E - E_new) / E) < tolfun:
            if display:
                print('TolFun exceeded: Optimisation terminated')
            break

        # Report progress
        E = E_new
        if display > 1:
            print('epoch = ' + str(i) + ': E = ' + str(E * scale))

    # Fiddle stress to match the original Sammon paper
    E = E * scale
    
    return [y,E]



def Find_duplicates(data):
    
    index_duplicates = np.zeros(len(data),dtype=np.int32)
    duplicates = [[]]*len(data)

    for i in range(len(data)):
        for j in range(len(data)):
            if (np.array_equal(data[i],data[j])) and (i!=j):
                #print i, "has duplicates", j
                index_duplicates[i]+=1
                duplicates[i]=data[i]

    duplicate_indexes = np.where(index_duplicates>0)[0]
    
    index_duplicates = np.array(index_duplicates)[duplicate_indexes]
    duplicates = np.array(duplicates)[duplicate_indexes]
    
    return index_duplicates, duplicates

def Compute_firing_rates(save_response_array_individual,millisecs_load, bin_width):

    #Look at cumulative firing rate for all trials
    fire_rate = []
    cell_rate = []
    #save respons array contains raster for each cell x movie repeats
    for i in range(len(save_response_array_individual)):  #Loops over cells
        #print "cell: ", i
        fire_rate = []
        for j in range(len(save_response_array_individual[i])): #Loops over movie repeats
            temp = np.zeros(millisecs_load/bin_width+1, dtype=np.float32) #Initialize firing rate vector to all zeros (also captures no spike repeats)
            all_spikes = np.sort(save_response_array_individual[i][j]) #Load spikes from repeats
            
            for s in all_spikes: #Sort.units[i]:
                s = s *1E3
                bin_loc = int(s / bin_width) #find bin # for spike
                temp[bin_loc]+=1

            fire_rate.append(temp)
        cell_rate.append(fire_rate)

    cell_rate=np.array(cell_rate)
    print cell_rate.shape
    cell_rate_trans = np.transpose(cell_rate, (1,2,0))
    print cell_rate_trans.shape

    return cell_rate_trans


def Load_din_5sec(work_dir,file_name):

    din_data = np.fromfile(work_dir + file_name+'/'+file_name+'.din',dtype=np.int64)
    print din_data

    times = din_data[::2]       #Extract frame times
    frames = din_data[1::2]     #Extract frame #s ranging from 0 to 299 for MVI_1400_5s stimulus; each frame presented 3 times

    #Search for max frame in movie sequence; 1st set blank values 65000+ to zero temporarily
    temp_frames = frames.copy()
    temp_frames[temp_frames>60000]=0
    final_frame= np.max(temp_frames)

    start_frames = np.where(frames==0)
    end_frames = np.where(frames==final_frame)

    start_times = times[start_frames]
    end_times = times[end_frames]

    start_times = start_times[::3]  #Take every 3rd time which marks beginning of 0th-frame (lasting 3 time-steps; CHECKED THIS)
                                    #This may not work for all DIN files; better to manually parse looking for FIRST ZERO/START TIME
    end_times = end_times[::3]      #Same check done

    print len(start_times), len(end_times)
    if len(end_times)<len(start_times):
        start_times=start_times[:len(end_times)]

    stim_times=[]
    for i in range(len(start_times)):
        stim_times.append([start_times[i], end_times[i]])
        
    stim_times = np.array(stim_times, dtype=np.float32)*1E-6 #Convert to seconds from usec
    return stim_times
    
def Load_din_driftbar(work_dir,file_name, experiment):

    din_data = np.fromfile(work_dir + file_name+'/'+file_name+'.din',dtype=np.int64)
    print din_data

    times = din_data[::2]       #Extract frame times
    frames = din_data[1::2]     #Extract frame #s ranging from 0 to 299 for MVI_1400_5s stimulus; each frame presented 3 times

    check=True
    start_index = []
    for i in range(len(frames)):
        if check and frames[i]==experiment:
            start_index.append(i)
            check=False
        if check==False and frames[i]!=experiment:
            check=True

    start_index = np.array(start_index)
    end_index = start_index + 801
    
    start_times = times[start_index]
    end_times = times[end_index]
    
    stim_times = np.column_stack((start_times, end_times))
    print stim_times*1E-6
    
    stim_times = np.array(stim_times, dtype=np.float32)*1E-6 #Convert to seconds from usec
    return stim_times




def Plot_firingrate(Sort1):

    #*********************** VLINE SETS OF PLOTS *****************************
    ax = plt.subplot(1, 1, 1)
    title("Original file: "+ Sort1.name, multialignment='center',fontsize=15)

    colors = ['blue','red', 'green', 'cyan', 'black']

    counter = 0 #These 2 variables used to track and color units that land on the same max channels
    depth = 0   #by cycling through the colors array above
    recording_length = 0
    for i in range(Sort1.n_units):
        if recording_length < max(Sort1.units[i]): recording_length = max(Sort1.units[i])

        #print "Plotting: ", i, "maxchan: ",Sort1.maxchan[i] #, Sort1.samplerate, Sort1.units[i][0:3],Sort1.units[i][-3:]
        x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5

        ymin=np.zeros(len(Sort1.units[i]))
        ymax=np.zeros(len(Sort1.units[i]))
        ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5 #Use x value to slightly offset units at same depth
        ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5-10.0 #lenght of tick line

        #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
        if depth == (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5): 
            counter+=1
        else:
            counter=0
        plt.vlines(x, ymin, ymax, linewidth=2, color=colors[mod(counter,5)])

        depth = (-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5)

    print recording_length #last spike in all rasters; gives sense of recording length;

    ##Compute state_vector
    #bin_width = 250. #bin width in ms; 4mil bins below should be enough for 1ms bins for 1 hr
    #n_bins = int(recording_length/Sort1.samplerate*1000/bin_width)+1
    #state_vector = np.zeros((Sort1.n_units, n_bins), dtype=np.int16) #Make array enough to fit 60mins 

    ##Compute and plot instantaneous firing rates: 
    #for i in range(Sort1.n_units):
        #for spike in Sort1.units[i]:
            #state_vector[i][int(spike/Sort1.samplerate*1000/bin_width)] +=1

    #state_vector = state_vector.T

    #inst_rate = np.zeros(n_bins,dtype=int16)
    #x = np.arange(0,n_bins)*bin_width/1000.
    #for i in range(n_bins):
        #inst_rate[i]=np.sum(state_vector[i])*10
    
    #plt.plot(x,inst_rate-1500, color='red')
    #for i in range(Sort1.n_units):
        #x = np.arange(0,int(recording_length/Sort1.samplerate*1000/bin_width)+1)*bin_width/1000.
        
        ##*bin_width/1000. #float(Sort1.samplerate)*2.5
        ##xmax = xmin+bin_width/1000.
        #print x
        #ymin=np.zeros(len(x))
        #ymax=np.zeros(len(x))
        #ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5 #Use x value to slightly offset units at same depth
        #ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5- \
        #10*np.array(state_vector[i],dtype=float32)#lenght of tick line
        
        ##r = matplotlib.patches.Rectangle((xmin, xmax), ymin, ymax, fill=True)
        ##ax.add_artist(r)
        ##ax.fill_between(x, 0, 1, where=y>theta, facecolor='green', alpha=0.5, transform=trans)
        #plt.vlines(x, ymin, ymax, linewidth=10, color='black', alpha=0.25)
    
    
    #for i in range(33000000/3):
        #plt.vlines(i*3/Sort1.samplerate, 0, -1200, linewidth=1, color='black')
    ax.xaxis.grid()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    #with open(Sort1.directory+ 'up_state_start_index', 'r') as f:
        #up_state_start_index = [line.rstrip('\n') for line in f]

    #with open(Sort1.directory+ 'up_state_end_index', 'r') as f:
        #up_state_end_index = [line.rstrip('\n') for line in f]

    #print up_state_start_index
    #print up_state_end_index

    #for i in range(len(up_state_start_index)):
        #p = ax.axvspan(up_state_start_index[i], up_state_end_index[i], facecolor='red', alpha=0.05)
        ## p = ax.axvspan(0.0, -2000.0, facecolor='1.0', alpha=0.0)

    plt.xlim(left=-0.03)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def Compute_upstate(Sort1):
    
    spike_matrix=[]
    tempint=0
    for i in range(Sort1.n_units):
        spike_matrix.append(Sort1.units[i])
        
        if tempint < max(Sort1.units[i]):
            tempint = max(Sort1.units[i])

    print tempint

    #Length of recording: tempint (in timesteps)
    time_index=5.0*Sort1.samplerate #Start at 5 seconds
    dt=int(0.030*Sort1.samplerate) #Use 30ms windows
    print "dt (in timesteps): ", dt
    n_spikes_t0=0
    n_spikes_t1=0
    n_spikes_t2=0

    up_state_flag = 0
    up_state_total_count=0

    up_state_start_index=[]
    up_state_end_index=[]
    up_state_count_index=[]

    #************************* LOOP OVER DATA *********************
    while True:
        n_spikes_t2=0
        for i in range(Sort1.n_units):
            a=np.array(Sort1.units[i])
            temp = np.where(np.logical_and(a>=time_index, a<=time_index+dt))[0]
            n_spikes_t2+=len(temp)


        #Up state start check
        if up_state_flag==1: 
            up_state_n_spikes+= n_spikes_t2
        else:
            if (n_spikes_t0 <2 and ((n_spikes_t1+n_spikes_t2)>10)):
                up_state_flag=1
                up_state_total_count+=1

                print "*****************"
                print "UP STATE #: ", up_state_total_count
                print "BEGIN: ", float(time_index-dt)/float(Sort1.samplerate), " sec."
                print n_spikes_t0
                print n_spikes_t1
                print n_spikes_t2
                up_state_time_index = time_index-dt
                up_state_n_spikes = n_spikes_t1+n_spikes_t2
                up_state_start_index.append(time_index-dt)

        #Up state ends
        if (up_state_flag==1  and n_spikes_t2<2):
            up_state_flag=0

            print "END: ", float(time_index)/float(Sort1.samplerate), " sec."
            print "LENGTH: ", float(time_index - up_state_time_index)/float(Sort1.samplerate)*1000, " ms."
            print "NO. SPIKES: ", up_state_n_spikes
            print n_spikes_t0
            print n_spikes_t1
            print n_spikes_t2
            print "*****************"
            print ""
            print ""
            up_state_end_index.append(time_index-dt)
            up_state_count_index.append(up_state_n_spikes-n_spikes_t2)

        time_index+=dt
        n_spikes_t0=n_spikes_t1
        n_spikes_t1=n_spikes_t2

        #print time_index, tempint
        if time_index>tempint:
            break

    print Sort1.directory

    with open(Sort1.directory+ 'up_state_start_index', 'w') as f:
        for s in up_state_start_index:
            f.write(str(float(s)/float(Sort1.samplerate)) + '\n')

    with open(Sort1.directory+'up_state_end_index', 'w') as f:
        for s in up_state_end_index:
            f.write(str(float(s)/float(Sort1.samplerate)) + '\n')

    with open(Sort1.directory+'up_state_count_index', 'w') as f:
        for s in up_state_count_index:
            f.write(str(s) + '\n')

def Compute_kuebler(Sort1, Sort2):
    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

    print "Computing Kuebler metric"

    bin_width = 0.1

    array=[]
    for i in range(Sort1.n_units):
        array.extend(Sort1.units[i])

    save_array = np.sort(np.array(array))/Sort1.samplerate

    colors = ['blue','red', 'black']
    yy = np.histogram(save_array, bins = np.arange(0,max(save_array),bin_width))

    Sort1.inst_firingrate=yy[0]/bin_width
    print "Length of firing rate vector: ", len(Sort1.inst_firingrate)
    return

    ax = plt.subplot(1, 2, 1)
    title("Original file: "+ Sort1.name + "  (bin_width = " + str(bin_width) + ")", multialignment='center',fontsize=15)
    plt.plot(yy[1][:-1],yy[0]/bin_width, linewidth=1, color='blue', alpha=0.750)


    max_x = max(save_array)
    plt.xlim(0, max_x)
    max_y = max(yy[0]/bin_width)
    plt.ylim(0, max_y)

    #plt.plot([0,max_x],[34, 34], linewidth=2, color='black', alpha=0.750)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Instantaneous Firing rate (#spikes/bin_width)',multialignment='center', fontsize=15)

    #Jittered spike trains
    jitter = int(max(save_array))

    ax = plt.subplot(1, 2, 2)
    title("Jittered spike trains ("+str(jitter) + " secs.)", multialignment='center',fontsize=15)

    for i in range(len(save_array)):
        save_array[i]=mod(save_array[i]+ jitter*random.random(), max_x)  

    colors = ['blue','red', 'black']
    yy = np.histogram(save_array, bins = np.arange(0,max(save_array),bin_width))
    plt.plot(yy[1][:-1],yy[0]/bin_width, linewidth=1, color='red', alpha=0.750)

    ave_rate = sum(yy[0]/bin_width)/len(yy[0]/bin_width)
 
    plt.plot([0,max_x],[ave_rate, ave_rate], linewidth=2, color='black', alpha=0.750)

    plt.xlim(0, max(save_array))
    plt.ylim(0, max_y)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Instantaneous Firing rate (#spikes/bin_width)',multialignment='center', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def Plot_upstate(Sort1, Sort2, lfp):
    
    #Compute SPECGRAM
    title(Sort1.directory + " (LFP events: " + str(len(Sort2.units)) + ")", multialignment='center',fontsize=15)
    if True:
        tsf_name = Sort1.directory + Sort1.filename+'.tsf'
        f = open(tsf_name, "rb")
        print "Loading "+ tsf_name + '.tsf'

        tsf = Tsf_file(f, Sort1.directory)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = Sort1.directory 
        tsf.tsf_name = tsf_name
        f.close()
        Sort1.tsf = tsf

    ax_image = plt.subplot(1,1,1) 
    img= lfp.specgram[::-1]
    rec_length = Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency
    im = ax_image.imshow(img, extent=[0,len(lfp.specgram.T)*lfp.tres,0,110], aspect='normal', 
    alpha=0.7, interpolation='none') #extent=lfp.extent, cmap=lfp.cm)

    #Compute ALL POWER from specgram
    allpower = []
    n_bands = 10
    power_specgram = np.zeros((n_bands,len(lfp.specgram.T)),dtype=np.float32)
    for i in range(len(lfp.specgram.T)):
        allpower.append(sum(lfp.specgram.T[i]))
        for j in range(n_bands):
            power_specgram[j][i]=sum(lfp.specgram.T[i][j*20:j*20+20])/60. #Have to multiply by 20 as 0.5Hz is the default stepsize
        
    colors=['blue','green','cyan','magenta','brown','pink','orange', 'black','red', 'yellow']
    
    #Compute MUA Histogram from all sorted units
    mua=[]
    for i in range(len(Sort1.units)):
        mua.append(Sort1.units[i])

    mua=np.hstack(mua)*25/(lfp.tres*Sort1.samplerate)

    y=np.histogram(mua, bins = np.arange(0,max(mua),50))
    plt.plot(y[1][0:-1]/100.,y[0]/3.5-15, linewidth=2, color='black')
    
    x = np.arange(0, len(lfp.specgram.T)*lfp.tres, lfp.tres)
    
    #mua_timeindices = muatimes.searchsorted(lfptimes)
    #muatimes[muatimeindices], lfptimes
    
    #for j in range(len(power_specgram)):
        #plt.plot(x, power_specgram[j]+j*10+20, color='red', linewidth=2)

    #Print Sort2 LFP rasters
    for i in range(len(Sort2.units)):
        print "LFP Unit: ", i, " # events: ", len(Sort2.units[i])
        ymin = -i*6-20
        ymax = ymin+5
        plt.vlines(np.array(Sort2.units[i])/lfp.SampleFrequency, ymin, ymax, linewidth=2, color=colors[i%10])
          
    plt.ylabel("LFP event rasters (multiple colours)   VS.  MUA (black)   VS.   Power in 10Hz increments (red/jet image)")
    plt.xlabel("Time (sec)")
    plt.show()
    
    
def Plot_LFP_triggered(Sort1, Sort2, lfp):

    
    #Load high pass filtered data if required
    if False:
        tsf_name = Sort1.directory + Sort1.filename+'.tsf'
        f = open(tsf_name, "rb")
        print "Loading "+ tsf_name + '.tsf'

        tsf = Tsf_file(f, Sort1.directory)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = Sort1.directory 
        tsf.tsf_name = tsf_name
        f.close()
        Sort1.tsf = tsf
    
    
    gs = gridspec.GridSpec(12,14)
    lfp_total=[]
    lfp_save =copy.deepcopy(lfp)
    for unit in range(len(Sort2.units)):  #Loop over population spikes
        print "Computing cumulative LFP: ", unit

        stim_array = np.array(Sort2.units[unit]) #Load spike times in ms
        
        window = 1000 #msec of window to average
        SampleFrequency=1000 #LFP Sample rate

        lfp_cumulative = np.zeros((lfp_save.n_electrodes,2*window), dtype=np.float32)
        
        #DO +/- window!!!
        for n in range(lfp_save.n_electrodes):
            counter=0
            for j in stim_array:
                start=max(0, int(j-window)) #Search backwards 1sec convert to ms
                end=min(int(start+2*window), lfp_save.n_vd_samples) #only search to end of recording; n_vd_samples already in ms (essentially timesteps)
                if (end-start)==2*window:
                    lfp_cumulative[n] += lfp_save.lfp_traces[n][start:end]
                    counter+=1
        lfp_cumulative /= counter
        lfp_total.append(lfp_cumulative)

        #Plot average LFP traces
        plt.suptitle(Sort1.directory + " (LFP events: " + str(len(Sort2.units)) + ") \n"
        + " LFP event: " + str(unit) + " # spikes: " + str(len(Sort2.units[unit])), multialignment='center',fontsize=15)

        ax_lfp = plt.subplot(gs[:10,7:14]) 
        for n in range(lfp.n_electrodes):
            ax_lfp.plot(lfp_total[unit][n]/10. - n*100, color='blue', linewidth=2)
        
        ax_lfp.set_ylim(-950., 50.)
        ax_lfp.plot([window,window], [0,-1000], 'r--', color='black')
        ax_lfp.xaxis.set_visible(False)
        ax_lfp.yaxis.set_visible(False)
    

        #Plot average specgram
        lfp.lfp_traces = lfp_total[unit]

        for i in range(10):
            ax_spec = plt.subplot(gs[i,:7]) 
            ax_spec.xaxis.set_visible(False)
            Plot_specgram(lfp, show_plot=False, spec_ch=i)

            img= lfp.specgram[::-1]
            im = ax_spec.imshow(img, extent=[0,2*window,0,60], aspect='normal')#,interpolation='none') #extent=lfp.extent, cmap=lfp.cm)
            if i==9: 
                ax_spec.xaxis.set_visible(True)
                ax_spec.set_xlabel("Time (ms)")
    
        #Plot average MUA
        mua_odd = []
        mua_even = []
        for n in range(len(Sort1.units)):
            for j in stim_array: #stim_array = Sort1.units[unit]
                
                temp2 = np.where(np.logical_and(np.array(Sort1.units[n])/Sort1.samplerate*1E3>=j-window, 
                np.array(Sort1.units[n])/Sort1.samplerate*1E3<=j+window))[0]
                for q in temp2 : 
                    if j%2 == 0:
                        mua_even.append(Sort1.units[n][q]/Sort1.samplerate*1E3-j+window)
                    else:
                        mua_odd.append(Sort1.units[n][q]/Sort1.samplerate*1E3-j+window)
                        
        ax_mua = plt.subplot(gs[10:12,7:14]) 

        yy_odd=np.histogram(np.array(mua_odd), bins = np.arange(0,2*window,20))
        plt.plot(yy_odd[1][0:-1],yy_odd[0], linewidth=1, color='blue', alpha=0.7)
        yy_even=np.histogram(np.array(mua_even), bins = np.arange(0,2*window,20))
        plt.plot(yy_even[1][0:-1],yy_even[0], linewidth=1, color='red', alpha=0.7)

        plt.plot(yy_even[1][0:-1],(yy_even[0]+yy_odd[0])/2., linewidth=4, color='black')

        ax_mua.set_xlabel("Time (ms)")

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()


    #Plot CSD triggered from different LFP events - Already have done this elsewhere though!
    #DO THIS +/- 1second though!!!


    #Plot synchrony index


def OLD_Plot_triggered_activity(Sort1,Sort2,lfp):

    #******** RASTER PLOTS *********
    if False:
        ax = plt.subplot(1, 1, 1)
        title("Original file: "+ Sort1.name +"\n Single unit spikes (blue)                 Population spikes (red) "
        + "\n Total units: " + str(Sort1.n_units) + "             Total pop spikes: " + str(len(Sort2.units[0])) , multialignment='center',fontsize=15)

        colors = ['blue','red', 'black']

        #Plot sorted spike rasters
        for i in range(Sort1.n_units):

            x = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate) #float(Sort1.samplerate)*2.5

            ymin=np.zeros(len(Sort1.units[i]))
            ymax=np.zeros(len(Sort1.units[i]))
            ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5
            ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]+Sort1.chanpos[Sort1.maxchan[i]][0]/56*5-8.0

            plt.vlines(x, ymin, ymax, linewidth=1.2, color='blue')

        #Plot pop spikes
        for i in range(Sort2.n_units):
            x = np.array(Sort2.units[i],dtype=float32)*20./float(Sort2.samplerate) #float(Sort1.samplerate)*2.5

            ymin=np.zeros(len(Sort2.units[i]))
            ymax=np.zeros(len(Sort2.units[i]))

            ymax+=-1400 #Sort2.chanpos[Sort1.maxchan[i]][1]+Sort2.chanpos[Sort1.maxchan[i]][0]/56*5-10.0

            plt.vlines(x, ymin, ymax, linewidth=0.5, color='red', alpha=0.5)

        ax.xaxis.grid()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

        plt.xlim(left=-0.03)

        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Depth from base of electrode (surface of cortex)',multialignment='center', fontsize=15)
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

    ##**********CROSS-CORRELOGRAM PLOTS ***************************************
    ##Plot cross-correlograms between pop spikes and each unit:
    ##Synchronized state from: 300sec -> 750sec
    #temp2 = np.array(Sort2.units[0],dtype=np.float32)*20./float(Sort2.samplerate)
    ##sync_array = temp2
    #sync_array = temp2[np.where(np.logical_and(temp2>=300., temp2<=750.))]

    #for i in range(Sort1.n_units):
        ##plt.title('Unit: ' + str(i))
        #print "plot: ", i
        #ax = plt.subplot(7, 10, i+1)
        ##ax =plt.subplot(1,1,1)
        #temp2 = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)
        #x1 = temp2[np.where(np.logical_and(temp2>=300., temp2<=750.))]
        #if len(x1)==0: continue
        #y = scipy.signal.correlate(x1,sync_array, mode='full')
        #y/=min(y)
        #x = np.arange(-len(y)/2,len(y)/2,1)
        #plt.bar(x,y)
        #plt.xlim(-500.,+500.)
    #figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    #plt.show()


    #*********** PLOT ALL UNIT RESPONSE TO EACH POP-SPIKE (all responses) *****************
    if False:
        for ss in range(len(Sort2.units)):
            
            #Make histograms -Xms .. +Xms from population spikes
            window=0.5        #Window to search in seconds
            bin_width = 0.020   #bin width for histograms in seconds
            min_spikes = 0      #min no. of spikes in window
            save_response_array=[]
            save_response_array1=[]
            save_response_array2=[]
            save_response_array3=[]
            save_response_array4=[]
            for i in range(Sort1.n_units):
                save_response_array.append([])
                save_response_array1.append([])
                save_response_array2.append([])
                save_response_array3.append([])
                save_response_array4.append([])

            #LOAD Pop Spike raster consisting of single-merged Unit
            pop_spikes = np.array(Sort2.units[ss],dtype=np.float32) 
            
            #Compute and plot pop-spike size distributions
            n_electrodes = 10
            ptp=[]
            maxchan=[]
            counter=0
            for s in pop_spikes:
                temp_ptp=[[0] for x in range(n_electrodes)]
                for i in range(n_electrodes):  
                    #Search 10 timesteps back and forth - NB: this is compressed LFP time
                    temp_ptp[i] = max(lfp.lfp_traces[i][s-10:s+10])-min(lfp.lfp_traces[i][s-10:s+10])

                ptp.append(max(temp_ptp))
                maxchan.append(np.argmax(temp_ptp))
                
                #print ptp[counter]
                #print maxchan[counter]
                counter+=1
            
            if False:
                yy=np.histogram(ptp, bins = np.arange(0,3000,50))

                plt.plot(yy[1][0:-1],yy[0], linewidth=2, color='blue')
                plt.show()

            temp_array = np.array(ptp, dtype=np.float32)
            temp2 = np.where(np.logical_and(temp_array>=150, temp_array<=1E6))[0]

            #Convert pop_spikes time back to same time scale as single units
            pop_spikes=pop_spikes[temp2]*100./float(Sort2.samplerate)  #This converts to Seconds

            rec_length = 0
            for i in range(len(Sort1.units)):
                if Sort1.units[i][-1]>rec_length: rec_length= Sort1.units[i][-1]
            rec_length/=float(Sort1.samplerate)
            print rec_length

            #Loop over all units in sorted data
            for i in range(Sort1.n_units):
                temp_array = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)  #This converts to seconds
                print "unit: ", i
                #Loop over all population spikes
                for j in range(len(pop_spikes)):
                    #Find all unit spikes that fall w/in +/- window of pop spike
                    temp2 = np.where(np.logical_and(temp_array>=pop_spikes[j]-window, temp_array<=pop_spikes[j]+window))[0]
                    if len(temp2)>min_spikes: #Min spikes in window 
                        save_response_array[i].extend(temp_array[temp2]-pop_spikes[j]) #Offset time of spike to t_0 
                        if j%2==0:  #Divide data sequentially
                            save_response_array1[i].extend(temp_array[temp2]-pop_spikes[j])
                        else:
                            save_response_array2[i].extend(temp_array[temp2]-pop_spikes[j])
                        if pop_spikes[j]< rec_length/2.:      #Divide data into 2 halves at half way point of pop_spikes
                            save_response_array3[i].extend(temp_array[temp2]-pop_spikes[j])
                        else:
                            save_response_array4[i].extend(temp_array[temp2]-pop_spikes[j])

            #def sinc_interp(x, s, u):
                #"""
                #Interpolates x, sampled at "s" instants
                #Output y is sampled at "u" instants ("u" for "upsampled")
                
                #from Matlab:
                #http://phaseportrait.blogspot.com/2008/06/sinc-interpolation-in-matlab.html        
                #"""
                
                #if len(x) != len(s):
                    #raise Exception, 'x and s must be the same length'
                
                ## Find the period    
                #T = s[1] - s[0]
                
                #sincM = tile(u, (len(s), 1)) - tile(s[:, newaxis], (1, len(u)))
                #y = dot(x, sinc(sincM/T))
                #return y

            #from numpy import exp
            #from scipy.optimize import curve_fit
            #from math import sin, pi

            #PLOT HISTOGRAMS - 3 UNITS
            if False:
                #aa=[4,11,25,15,21,24]
                aa=[17, 18, 23]
                scale=1
                plt.suptitle(Sort1.name+ "   Length of recording: " + str(rec_length) + " seconds. "+
                "\n plotted spikes: " + str(len(pop_spikes)) + "   of    Total pop spikes: " + str(len(Sort2.units[ss])), fontsize=20)

                for i in range(len(aa)):
                    ax = plt.subplot(math.ceil(len(aa)/3.), 3, i+1)
                    #ax = plt.subplot(1, 1, 1)
                    
                    plt.title("Unit: " + str(aa[i]))
                    plt.ylabel("# spikes bins ("+str(bin_width*1000.)+" ms)")
                    plt.xlabel("Seconds from pop event")

                    yy1=np.histogram(save_response_array1[aa[i]], bins = np.arange(-window,+window,bin_width))
                    yy2=np.histogram(save_response_array2[aa[i]], bins = np.arange(-window,+window,bin_width))
                    
                    plt.plot(yy1[1][0:-1],yy1[0], linewidth=1, color='blue')
                    plt.plot(yy2[1][0:-1],yy2[0], linewidth=1, color='red')

                    x=(yy1[0][0:-1]+yy2[0][0:-1])/2.
                    s=(yy1[1][0:-2]+yy2[1][0:-2])/2.
                    u=np.arange(-window,+window-0.020*3,.001)
                    fitdata=sinc_interp(x,s,u)
                    plt.plot(u,fitdata, linewidth=5, color='black')
                   
                    plt.plot([0,0],[0,1000], 'r--',color='black',linewidth=1)
                    plt.ylim(0,max(yy1[0])+max(yy1[0])*.2)

                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
                plt.show()
                #quit()

            #PLOT HISTOGRAMS - ALL UNITS
            counter=0
            scale=1
            plt.suptitle(Sort1.name+ "   Length of recording: " + str(rec_length) + " seconds. "+
            "\n plotted spikes: " + str(len(pop_spikes)) + "   of    Total pop spikes: " + str(len(Sort2.units[ss])), fontsize=20)
            peak=[]
            for i in range(Sort1.n_units):
                temp_ylim = 0

                ax = plt.subplot(int(sqrt(Sort1.n_units))+1, int(sqrt(Sort1.n_units))+1, i+1)
                print "Plotting Unit: ", i
                plt.ylabel("#"+str(i)+" : " + str(len(Sort1.units[i])), fontsize=14)

                if len(save_response_array[i])<min_spikes: continue
                yy1=np.histogram(save_response_array1[i], bins = np.arange(-window,+window,bin_width))
                plt.plot(yy1[1][0:-1],yy1[0], linewidth=1, color='blue')
                temp_ylim = max(max(yy1[0]),temp_ylim)
                
                yy2=np.histogram(save_response_array2[i], bins = np.arange(-window,+window,bin_width))
                plt.plot(yy2[1][0:-1],yy2[0], linewidth=1, color='red')
                temp_ylim = max(max(yy2[0]),temp_ylim)

                x=(yy1[0][0:-1]+yy2[0][0:-1])/2.
                s=(yy1[1][0:-2]+yy2[1][0:-2])/2.
                u=np.arange(-window,+window-0.020*3,.001)
                fitdata=sinc_interp(x,s,u)
                plt.plot(u,fitdata, linewidth=5, color='black')
                
                if abs(u[np.argmax(fitdata)])<0.1 :
                    temp_x = u[np.argmax(fitdata)]
                    temp_y = fitdata[np.argmax(fitdata)]
                    plt.plot([temp_x,temp_x],[0, temp_y+10], 'r--',color='green',linewidth=4)
                    peak.append(temp_x)
                else:
                    peak.append(0.0)
            
                #yy=np.histogram(save_response_array3[i], bins = np.arange(-window,+window,.020))
                #plt.plot(yy[1][0:-1],yy[0], linewidth=1, color='green')
                #temp_ylim = max(max(yy[0]),temp_ylim)

                #yy=np.histogram(save_response_array4[i], bins = np.arange(-window,+window,.020))
                #plt.plot(yy[1][0:-1],yy[0], linewidth=1, color='black')
                #temp_ylim = max(max(yy[0]),temp_ylim)

                plt.plot([0,0],[0,1000], 'r--',color='black',linewidth=1)
                plt.ylim(0,temp_ylim+2)

            mng = plt.get_current_fig_manager()
            mng.resize(*mng.window.maxsize())
            plt.show()


def Plot_templates(sim_dir, Sorts_lfp, lfp, recs, track_name):

    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    #Save template in list to make it easier to compute averages and std
    templates = []
    for unit in range(len(Sorts_lfp[0])):
        templates.append([])
        for ch in range(len(lfp.data[0])):        #Save multi-channel template array
            templates[unit].append([])
    
    
    #Filter lfp to remove 5Hz and lower freqs:
    if True:
        lowcut=5
        highcut=240
        fs=1000
        print "Filtering data"
        for rec in range(len(Sorts_lfp)):               #Loop over all recordings
            for ch in range(len(lfp.data[rec])):
                lfp.data[rec][ch]=butter_bandpass_filter(lfp.data[rec][ch], lowcut, highcut, fs, order = 2)

    print "Compiling averages"
    #Accumulate spike triggered traces
    for rec in range(len(Sorts_lfp)):               #Loop over all recordings
        for unit in range(len(Sorts_lfp[rec])):     #Loop over each unit in recording
            for spike in Sorts_lfp[rec][unit]:
                for ch in range(len(lfp.data[rec])):
                    temp_data = lfp.data[rec][ch][spike-100:spike+100]
                    if len(temp_data)==200:
                        templates[unit][ch].append(lfp.data[rec][ch][spike-100:spike+100])
    
    
    #Plot data
    xx = np.arange(0,200,1)
    scale_factor = .8
    for unit in range(len(Sorts_lfp[0])):
        ax=plt.subplot(1,len(Sorts_lfp[0]), unit+1)
        for ch in range(len(lfp.data[rec])):
            temp_ave = np.average(np.float32(templates[unit][ch]),axis=0)*scale_factor-ch*20000
            plt.plot(temp_ave, linewidth=5, color='black')

            y_vals = np.std(np.float32(templates[unit][ch]),axis=0)*scale_factor
            plt.fill_between(xx, temp_ave- y_vals,temp_ave, facecolor=colors[unit], alpha=0.4)
            plt.fill_between(xx, temp_ave+ y_vals,temp_ave, facecolor=colors[unit], alpha=0.4)         
        
        plt.title("Unit: "+str(unit+1)+"\n #: "+str(len(templates[unit][0])), fontsize=22)
        ax.yaxis.set_visible(False)
        
        x_original = np.arange(0,251,50)
        x_label = np.int32(np.arange(-100, 151,50))
        plt.xticks(x_original, x_label, fontsize=10)
        plt.plot([100,100],[-220000,30000*len(lfp.data[rec])],'r--', linewidth=3, color='black', alpha=0.4)
        
        plt.ylim(-220000,30000)
        plt.xlabel("Time (ms)", fontsize=22)
        ax.yaxis.set_visible(False)
    plt.suptitle(track_name, fontsize=20)
    plt.show()
    quit()

    print "# lfp cluters: ", end_lfp
    n_recs = len(recs)
    plotting=True
    start_window = -0.04
    end_window = 0.01
    
    si_limit = 0.7
    
    #Convert rec name from string to relative index in concatenated data; 
    #also compute # recordings w. minimum required lfp events
    rec_indexes=[]
    n_lfprecs = 0
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        if len(Sorts_lfp[lfp.rec.index(rec)][start_lfp])>5: n_lfprecs+=1
    if n_lfprecs==0:
        print "No lfp_event recordings!"
        return

def Plot_percentlock_synchronized_state_allrec(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):

    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    window=1
    start_lfp = 0; end_lfp = len(Sorts_lfp[0])  

    print "# lfp cluters: ", end_lfp
    n_recs = len(recs)
    plotting=True
    start_window = -0.03
    end_window = 0.01
    
    font_size=20    #Default fontsize value
    
    si_limit = 0.7
    
    #Convert rec name from string to relative index in concatenated data; 
    #also compute # recordings w. minimum required lfp events
    rec_indexes=[]
    n_lfprecs = 0
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        if len(Sorts_lfp[lfp.rec.index(rec)][start_lfp])>5: n_lfprecs+=1
    
    if n_lfprecs==0:
        print "No lfp_event recordings!"
        return

    #compute track recording length:
    track_length = 0
    for rec_index in rec_indexes:
        track_length+= lfp.t_end[rec_index] #Recording length in seconds

    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes
    
    #Make array to collect locked spikes for each unit across all recs
    total_locked_spikes_allrecs=np.zeros((end_lfp-start_lfp,total_units),dtype=np.float32)
    
    #Make array to collect all single unit spikes:
    sua_allspikes = np.zeros(total_units, dtype=np.float32)
                
    #Loop over all recordings and plot data where # lfp events > minimum (5)
    #gs = gridspec.GridSpec(end_lfp-start_lfp+2,n_lfprecs)       #Make rows = # of lfp clusters + 1 additional row for summary plots
    gs = gridspec.GridSpec(3,end_lfp-start_lfp)       #Make rows = # of lfp clusters + 1 additional row for summary plots

    #For each unit save list of non-locking spikes
    locked_spikes_times = []
    for i in range(total_units):
        locked_spikes_times.append([])

    #Pick a particular LFP event to lock to
    all_pop_spikes = 0 #add up all pop spikes over all recordings to use as control for final plot
    for ss in range(start_lfp, end_lfp, 1):
        print ""
        print ""
        print "LFP Cluster # : ", ss+1, " / ", len(Sorts_lfp[0])

        cluster_pop_spikes = 0  #Keep track of all pop spikes during sync periods for each LFP cluster
        
        #Make lists to hold # unique spikes that lock to lfp events
        locked_spikes = [[] for x in xrange(total_units)]
        
        #Make list to hold all spikes during sync periods;
        Sorts_sua_sync_spikes = np.zeros(total_units, dtype=np.float32)

        #Save total length of sync periods for each set of recordings; secs
        sync_period_total_time = 0 

        for rec_index in rec_indexes:
            print "...rec: ", rec_index, " / ", len(rec_indexes)
            if len(Sorts_lfp[rec_index][start_lfp])>5:  #If there are < 5 lfp events - skip recording;
                pass
            else:
                print "...too few lfp events... skipping."
                continue

            #Compute periods of synchrony from si index
            data_in = lfp.data[rec_index][9]
            si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)

            #Save total length of sync periods; it is the same for each lfp loop
            for k in range(len(sync_periods)):
                sync_period_total_time += sync_periods[k][1]-sync_periods[k][0]

            #Load pop events during synch periods (secs)
            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3  
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
                temp_list.extend(pop_spikes[indexes])
            pop_spikes=np.array(temp_list)

            #Remove any pop spikes that are w/in 50ms of another cluster pop spike:
            if True:
                remove_spikes = []
                for temp_ss in range(ss+1, end_lfp, 1):
                    if len(Sorts_lfp[rec_index][temp_ss])>0:
                        for k in range(len(pop_spikes)):
                            if np.min(np.abs(Sorts_lfp[rec_index][temp_ss]*1E-3-pop_spikes[k]))<abs(start_window):  #start_window is usually largest
                                remove_spikes.append(k)
                pop_spikes= np.delete(pop_spikes,remove_spikes)
            
            print "...no. of sync period pop events: ", len(pop_spikes)
            
            #Track all pop_spikes for individual clusters
            cluster_pop_spikes+= len(pop_spikes)
            
            #Track cumulative total of pop_spikes
            all_pop_spikes+= len(pop_spikes)
            
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                temp_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    temp_list.extend(spike_array[indexes])
                spike_array=np.array(temp_list)
                
                #Track cumulative total # spikes during synch period ACROSS ALL RECS; 
                if (ss-start_lfp)==0: sua_allspikes[unique_unit]+=len(spike_array) #Count only once over LFP cluster loops;
                
                #Save # of spikes during sync period; use sequential unit id - only used w/in a recording
                Sorts_sua_sync_spikes[unique_unit]+=len(spike_array)

                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if abs(pop_spikes[j+1]-pop_spikes[j])<abs(start_window): continue
                   
                    #Add # spikes occuring w/in window of lfp event
                    locked_indexes = np.where(np.logical_and(spike_array>=pop_spikes[j]+start_window, spike_array<=pop_spikes[j]+end_window))[0]
                    locked_spikes[unique_unit].extend(spike_array[locked_indexes])
                    
                    #Also save a master list of indexes of all locking spikes for each unit
                    locked_spikes_times[unique_unit].extend(spike_array[locked_indexes])

        #*****************************************************************************
        #************** SCATTER PLOT DISTRIBUTION OF LOCKED SPIKES TO LFP ************
        #*****************************************************************************
        #COMPUTE % LOCK AND FIRE RATE
        percent_lock = []
        fire_rate = []
        for unit in range(total_units):
            n_spikes_unique = len(np.unique(locked_spikes[unit]))
            if Sorts_sua_sync_spikes[unit]>0:
                percent_lock.append(float(n_spikes_unique)/Sorts_sua_sync_spikes[unit])
            else:
                percent_lock.append(0)
                
            fire_rate.append(float(n_spikes_unique)/sync_period_total_time)

            #Save number of unique spikes locked for each lfp cluster, each unit, and each recording
            total_locked_spikes_allrecs[ss-start_lfp][unit] = n_spikes_unique   
        
        #print percent_lock
        percent_lock = np.array(percent_lock)*1E2 #Convert to percent 
        fire_rate = np.array(fire_rate)

        ax1 = plt.subplot(gs[0,ss])
        plt.scatter(fire_rate, percent_lock, color=colors[ss%10])

        #Compute total % of recording that pop spikes occupy = n_pop_spikes x length_pop_spike/ recording_length
        exp_lock = cluster_pop_spikes*(end_window - start_window)/sync_period_total_time*100
        plt.plot([0,100],[exp_lock, exp_lock], 'r--', color='black',  linewidth=2., alpha=0.8)
        
        if (ss-start_lfp)==0: plt.title("rec: "+str(recs[rec_index]), fontsize=10)
        plt.plot([0,0], [-1,101], color='black')
        plt.plot([-1,100], [0,0], color='black')
                    
        ax1.set_yscale('symlog', linthreshy=0.1, basex=10)
        ax1.set_xscale('symlog', linthreshx=0.001, basex=10)
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        y_min = -0.0001
        y_max = 110
        x_min = -0.00001 #np.min(fire_rate[fire_rate!=0])*.8
        x_max = 100.
        plt.ylim(y_min,y_max)
        plt.xlim(x_min,x_max)
        plt.tick_params(axis='both', which='both', labelsize=font_size-10)


        plt.xlabel("Cell Firing Rate (Hz)", fontsize=font_size-3)
        plt.title("LFP Cluster: "+str(ss)+"\n# events: "+str(cluster_pop_spikes), fontsize=font_size-3)
        if (ss-start_lfp)==0: plt.ylabel("% Spikes Locked", fontsize=font_size)
        
        #***********************************************************
        #************** PLOT BAR HISTOGRAMS BY DEPTH ***************
        #***********************************************************
        #First, normalize the number of spikes for each unit:
        for unit in range(total_units):
            #print total_locked_spikes_allrecs[ss][unit], sua_allspikes[unit]
            if sua_allspikes[unit]>0: #make sure recording looped over contain spikes from unit - not all recs do
                total_locked_spikes_allrecs[ss][unit] = total_locked_spikes_allrecs[ss][unit] / sua_allspikes[unit]
            else:
                total_locked_spikes_allrecs[ss][unit] = 0
        
        #Plot bar graphs by depth of cell
        ax1 = plt.subplot(gs[1, :])
        indexes = np.arange(0,len(sua_allspikes),1) #Cells should have been exported by depth...
        
        cumulative=np.zeros(len(indexes),dtype=np.float32)
        x_count = 0
        for u in indexes: 
            for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                cumulative[u] += total_locked_spikes_allrecs[p][u]
            plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
            if (total_locked_spikes_allrecs[ss-start_lfp][u]*100. + cumulative[u]*100)>100:
                print sua_allspikes[u]
                print locked_spikes[u]
            x_count+=1
        #***********************************************************
        #************** PLOT BAR HISTOGRAMS BY FIR RATE ************
        #***********************************************************
        #Sort indexes in order of firing rate; use master list of sua_allspikes, otherwise incorrect;
        ax1 = plt.subplot(gs[2, :])
        indexes = sua_allspikes.argsort()
        cumulative=np.zeros(len(indexes),dtype=np.float32)
        x_count = 0
        for u in indexes: 
            for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                cumulative[u] += total_locked_spikes_allrecs[p][u]
            plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
           
            x_count+=1
                
    #Compute expected lock for allc ells also
    exp_lock = all_pop_spikes*(end_window-start_window)/sync_period_total_time*1E2
    print "total # pop spikes: ", all_pop_spikes
    print "track_length: ", track_length
    print "track_length - sync periods only: ", sync_period_total_time
    print "expected % lock: ", exp_lock

    #Put labels on depth histograms - FIRST ROW DOWN
    ax1 = plt.subplot(gs[1, :])
    plt.plot([0,len(sua_allspikes)],[exp_lock,exp_lock], 'r--', color='cyan', linewidth=2)

    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.tick_params(axis='x', which='both', labelsize=8)
    plt.ylabel("% Spikes Locked\n(LFP Cluster Contributions)", fontsize=font_size)
    plt.xlabel("Superficial <----- Depth of cell -----> Deeper ", fontsize=font_size)
    ax1.xaxis.set_ticks([])

    #Put labels on firing rate ordered histograms - SECOND ROW DOWN
    ax1 = plt.subplot(gs[2, :])
    plt.plot([0,len(sua_allspikes)],[exp_lock,exp_lock], 'r--', color='cyan', linewidth=2)

    #Compute ratio of expected lock to actual lock and plot on top of bar graphs
    cumulative=np.zeros(len(indexes),dtype=np.float32)
    for u in indexes:
        for p in range(0,ss-start_lfp+1):                         #Sum all locking percentages over all LFP spikes
            cumulative[u] += total_locked_spikes_allrecs[p][u]

    x_count=0
    for u in indexes:
        text_temp = str(int(round(cumulative[u]/exp_lock*100.)))
        plt.text(x_count, cumulative[u]*100.+2 , text_temp, fontsize=7)
        x_count+=1

    #Use firing rates for each cell to label axes; NB: use sua_allspikes which contains only spikes during sync periods
    x_label = []
    sua_allspikes = sua_allspikes[indexes] #Reorder all units by firing rate
    for k in range(len(sua_allspikes)):
        x_label.append(str(round(sua_allspikes[k]/track_length,3)))
    xx = np.linspace(0,len(sua_allspikes)-1,len(sua_allspikes))+0.5

    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.ylabel("%locked \n(%exp: "+str(round(exp_lock,2))+")", fontsize=font_size)

    plt.xticks(xx, x_label, fontsize=8)
    plt.xticks(rotation=90)
    plt.xlabel("Lower  <----- Cell firing rates (Hz) -----> Higher", fontsize=font_size)
    plt.suptitle(sim_dir+track_name + ", clusters: "+ str(len(Sorts_lfp[rec_index]))+", si_limit: "+str(si_limit)+ ", "+str(int(start_window*1E3))+
    "ms.."+str(int(end_window*1E3))+"ms.", fontsize=25)
            
    #SHOW PLOT
    plt.show()
    
    #******** SAVE NON-LOCKING SPIKES
    non_locked_spikes=[]
    locked_spikes=[]
    for k in range(len(Sorts_sua[rec_index].units)):
        unique_unit = Sorts_sua[rec_index].uid[k]
        print "Finding unique spikes unit: ", unique_unit

        locked_spikes.append(np.unique(np.array(locked_spikes_times[unique_unit])))                                     #times of all locked spikes; secs)
        all_spikes = np.float32(Sorts_sua[rec_index].units[k])/float(Sorts_sua[rec_index].samplerate)  #times of all spikes for unit; secs
        
        #find non_locking spikes; i.e. remove locked_spikes that are w/in 0.1 ms of all_spikes times:
        remove_spikes = []
        for s in range(len(locked_spikes[k])):
            closest_spike_index = np.argmin(np.abs(all_spikes-locked_spikes[k][s]))     #start_window is usually largest
            if abs(all_spikes[closest_spike_index]-locked_spikes[k][s])<0.0001:         #Need this to ensure float round off doesn't miss spike matches
                remove_spikes.append(closest_spike_index)
        non_locked_spikes.append(np.delete(all_spikes,remove_spikes))
    
    print Sorts_sua[rec_index].directory
    np.save(Sorts_sua[rec_index].directory+Sorts_sua[rec_index].filename+'_non_locked_spikes' , non_locked_spikes)
    np.save(Sorts_sua[rec_index].directory+Sorts_sua[rec_index].filename+'_locked_spikes' , locked_spikes)

    #********





def Plot_percentlock_synchronized_state(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):

    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    window=1
    small_window = 100 # ms of window for luczak plots
    large_window = window*1E3
    start_lfp = 0; end_lfp = len(Sorts_lfp[0])# 
    #start_lfp = 0; end_lfp = len(Sorts_lfp[0])# 

    print "# lfp cluters: ", end_lfp
    n_recs = len(recs)
    plotting=True
    start_window = -0.04
    end_window = 0.01
    
    si_limit = 0.7
    
    #Convert rec name from string to relative index in concatenated data; 
    #also compute # recordings w. minimum required lfp events
    rec_indexes=[]
    n_lfprecs = 0
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        if len(Sorts_lfp[lfp.rec.index(rec)][start_lfp])>5: n_lfprecs+=1
    if n_lfprecs==0:
        print "No lfp_event recordings!"
        return

    #compute track recording length:
    track_length = 0
    for rec_index in rec_indexes:
        track_length+= lfp.t_end[rec_index] #Recording length in seconds

    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes
    
    #Make array to collect locked spikes for each unit across all recs
    total_locked_spikes_allrecs=np.zeros((end_lfp-start_lfp, total_units),dtype=np.float32)
    
    #Make array to collect all single unit spikes:
    sua_allspikes = np.zeros(total_units, dtype=np.float32)
                
    #Loop over all recordings and plot data where # lfp events > minimum (5)
    gs = gridspec.GridSpec(end_lfp-start_lfp+2,n_lfprecs)       #Make rows = # of lfp clusters + 1 additional row for summary plots

    #Pick a particular LFP event to lock to
    all_pop_spikes = 0 #add up all pop spikes over all recordings to use as control for final plot
    sync_period_total_time = 0 #total period of synchronization across all recs (secs)
    for ss in range(start_lfp, end_lfp, 1):
        print ""
        print ""
        print "LFP Cluster # : ", ss+1, " / ", len(Sorts_lfp[0])

        cluster_pop_spikes = 0  #Keep track of all pop spikes during sync periods for each LFP cluster

        rec_counter=0
        for rec_index in rec_indexes:
            print "...rec: ", rec_index, " / ", len(rec_indexes)
            if len(Sorts_lfp[rec_index][start_lfp])>5:  #If there are < 5 lfp events - skip recording;
                pass
            else:
                print "...too few lfp events... skipping."
                continue

            #Compute periods of synchrony from si index
            data_in = lfp.data[rec_index][9]
            si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)

            #Save total length of sync periods; DO IT ONCE ONLY!
            if (ss-start_lfp)==0:
                for k in range(len(sync_periods)):
                    sync_period_total_time += sync_periods[k][1]-sync_periods[k][0]

            #Make lists to hold # unique spikes that lock to lfp events
            n_lock_spikes = [[] for x in xrange(total_units)]

            #Load pop events during synch periods (secs)
            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3  
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
                temp_list.extend(pop_spikes[indexes])
            pop_spikes=np.array(temp_list)
            
            print "...no. of sync period pop events: ", len(pop_spikes)
            
            #Track all pop_spikes for individual clusters
            cluster_pop_spikes+= len(pop_spikes)
            
            #Track cumulative total of pop_spikes
            all_pop_spikes+= len(pop_spikes)
            
            Sorts_sua_sync_spikes = np.zeros(total_units, dtype=np.float32)
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                temp_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    temp_list.extend(spike_array[indexes])
                spike_array=np.array(temp_list)
                
                #Track cumulative total # spikes during synch period across all recs; use unique unit as it's cumulative across all recs
                #NB: COUNT ONLY ONCE DURING MULTIPLE LFP LOOPS!!!
                if (ss-start_lfp)==0: sua_allspikes[unique_unit]+=len(spike_array)
                
                #Save # of spikes during sync period; use sequential unit id - only used w/in a recording
                Sorts_sua_sync_spikes[unit]=len(spike_array)
                
                #NB: May wish to remove LARGE SPIKE TIMES due to bug 1E+8 for secs should do it.
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                control_list=[]     #Collect all spikes from unit that locked to control spike
                control_n_spikes=0. #number of spikes locked to control lfp event
                
                #print cross-correlogram plots
                if False:
                    from brian import correlogram
                    cor_width = 30
                    cor_plot = correlogram(pop_spikes,spike_array,width=cor_width,bin=1,T=None)
                    
                    print len(pop_spikes)
                    print len(spike_array)
                    
                    #cor_plot = xcorr(pop_spikes, spike_array, [spike_array[0], spike_array[-1]])
                    #print cor_plot
                    plt.plot(cor_plot)
                    plt.plot([cor_width, cor_width],[0,max(cor_plot)*1.2], 'r--')
                    plt.ylim(bottom=0)
                    plt.title("LFP Cluster: "+str(ss)+"  Unit: "+str(unit)+ " / "+str(len(Sorts_sua[rec_index].units)))
                    plt.show()
                    #quit()
                
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)

                    #Add # spikes occuring w/in 25ms of lfp event
                    n_lock_spikes[unit].extend(spike_array[np.where(np.logical_and(spike_array>=pop_spikes[j]+start_window, spike_array<=pop_spikes[j]+end_window))[0]])

            #****************************************************************************************
            #************** SCATTER PLOT DISTRIBUTION OF LOCKED SPIKES TO LFP ***************************
            #****************************************************************************************

            #Plot blue-scatter plot: real fire rate vs. % of spikes 
            percent_lock = []
            for unit in range(len(Sorts_sua[rec_index].units)):
                n_spikes_unique = len(np.unique(n_lock_spikes[unit]))*100
                percent_lock.append(float(n_spikes_unique)/Sorts_sua_sync_spikes[unit])
                
                #Save number of unique spikes locked for each lfp cluster, each unit, and each recording
                total_locked_spikes_allrecs[ss-start_lfp][Sorts_sua[rec_index].uid[unit]] += n_spikes_unique*1E-2          #Add # locked spikes to total for each unique unit id
            
            percent_lock = np.array(percent_lock)
            
            #x-axis is firing rate of cell
            fire_rate = []
            for unit in range(len(Sorts_sua[rec_index].units)):
                fire_rate.append(float(len(Sorts_sua[rec_index].units[unit]))/lfp.t_end[rec_index])
            fire_rate = np.array(fire_rate)
            
            ax1 = plt.subplot(gs[ss-start_lfp, rec_counter])

            plt.scatter(fire_rate, percent_lock, color=colors[ss%10])

            #Compute total % of recording that pop spikes occupy = n_pop_spikes x length_pop_spike/ recording_length
            rec_length = lfp.t_end[rec_index] #Recording length in seconds
            exp_lock = len(pop_spikes)*(end_window - start_window)/rec_length*100
            plt.plot([0,100],[exp_lock, exp_lock], 'r--', color='black',  linewidth=2., alpha=0.8)
            
            if rec_counter ==0:
                pass
                #plt.ylabel("LFP Clstr: " + str(ss) + ", "+str(total_pop_spikes) + "\n%Spks Lock", fontsize=8)
            else:
                ax1.yaxis.set_visible(False)
            
            if (ss-start_lfp)==0: plt.title("rec: "+str(recs[rec_index]), fontsize=10)
            plt.plot([0,0], [-1,101], color='black')
            plt.plot([-1,100], [0,0], color='black')
            
            ax1.set_yscale('symlog', linthreshy=0.1, basex=10)
            ax1.set_xscale('symlog', linthreshx=0.001, basex=10)

            y_min = -0.01
            y_max = 110
            x_min = np.min(fire_rate[fire_rate!=0])*.8
            x_max = 100.
            plt.ylim(y_min,y_max)
            plt.xlim(x_min,x_max)
            plt.tick_params(axis='x', which='both', labelsize=8)

            rec_counter+=1
        
        #Go back and write the number of lfp events in ylabel for cluster
        ax1 = plt.subplot(gs[ss-start_lfp, 0])
        plt.ylabel(str(cluster_pop_spikes))


        #*******************************************************************************
        #Bottom Graphs: Bar graph of % of sua spikes for each lfp cluster **************
        #*******************************************************************************
        #First, normalize the number of spikes for each unit:
        for u in range(total_units):
            if sua_allspikes[u]>0: #make sure recording looped over contain spikes from unit - not all recs do
                total_locked_spikes_allrecs[ss-start_lfp][u] = total_locked_spikes_allrecs[ss-start_lfp][u] / sua_allspikes[u]
            else:
                total_locked_spikes_allrecs[ss-start_lfp][u] = 0
        
        #Plot bar graphs by depth of cell
        ax1 = plt.subplot(gs[end_lfp-start_lfp, :])
        indexes = np.arange(0,len(sua_allspikes),1)
        if (ss-start_lfp)==0: #Plot first bar graphs
            x_count = 0
            for u in indexes: 
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss])
                x_count+=1
        #Plot cumulative bar graphs for additional lfp clusters
        else:
            cumulative=np.zeros(len(indexes),dtype=np.float32)
            for u in indexes:
                for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                    cumulative[u] += total_locked_spikes_allrecs[p][u]
            x_count=0
            for u in indexes:
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
                x_count+=1

        #Sort indexes in order of firing rate; use master list of sua_allspikes, otherwise incorrect;
        ax1 = plt.subplot(gs[end_lfp-start_lfp+1, :])
        indexes = sua_allspikes.argsort()
        #Plot lfp cluster 0
        if (ss-start_lfp)==0: #Plot first bar graphs
            x_count = 0
            for u in indexes: 
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss])
                x_count+=1
        #Plot cumulative bar graphs for additional lfp clusters
        else:
            cumulative=np.zeros(len(indexes),dtype=np.float32)
            for u in indexes:
                for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                    cumulative[u] += total_locked_spikes_allrecs[p][u]
            x_count=0
            for u in indexes:
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
                x_count+=1
                

    #Put labels on depth histograms:
    ax1 = plt.subplot(gs[end_lfp-start_lfp, :])
    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.tick_params(axis='x', which='both', labelsize=8)
    plt.ylabel("%locked")
    plt.xlabel("Depth of cell (relative Units)", fontsize=15)
    

    #Put labels on firing rate ordered histograms
    ax1 = plt.subplot(gs[end_lfp-start_lfp+1, :])
    
    #Compute expected lock for allc ells also
    exp_lock = all_pop_spikes*(end_window-start_window)/sync_period_total_time*1E2
    print "total # pop spikes: ", all_pop_spikes
    print "track_length: ", track_length
    print "track_length - sync periods only: ", sync_period_total_time
    print "expected % lock: ", exp_lock
    plt.plot([0,len(sua_allspikes)],[exp_lock,exp_lock], 'r--', color='cyan', linewidth=2)

    #Compute ratio of expected lock to actual lock and plot on top of bar graphs
    cumulative=np.zeros(len(indexes),dtype=np.float32)
    for u in indexes:
        for p in range(0,ss-start_lfp+1):                         #Sum all locking percentages over all LFP spikes
            cumulative[u] += total_locked_spikes_allrecs[p][u]

    x_count=0
    for u in indexes:
        #text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)
        text_temp = str(int(round(cumulative[u]/exp_lock*100.)))
        plt.text(x_count, cumulative[u]*100.+2 , text_temp, fontsize=7)
        x_count+=1

    #Use firing rates for each cell to label axes; NB: use sua_allspikes which contains only spikes during sync periods
    x_label = []
    sua_allspikes = sua_allspikes[indexes] #Reorder all units by firing rate
    for k in range(len(sua_allspikes)):
        x_label.append(str(round(sua_allspikes[k]/track_length,3)))
    xx = np.linspace(0,len(sua_allspikes)-1,len(sua_allspikes))+0.5

    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.ylabel("%locked \n(%exp: "+str(round(exp_lock,2))+")")

    plt.xticks(xx, x_label, fontsize=8)
    plt.xticks(rotation=90)
    plt.xlabel("Firing rates - during sync periods (Hz)", fontsize=15)
    plt.suptitle(sim_dir+track_name + ", clusters: "+ str(len(Sorts_lfp[rec_index]))+", si_limit: "+str(si_limit)+ ", "+str(int(start_window*1E3))+
    "ms.."+str(int(end_window*1E3))+"ms.", fontsize=25)
    
        
    #SHOW PLOT
    plt.show()
    
    
def Plot_msl(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):

    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    start_lfp = 0
    end_lfp = 1 # len(Sorts_lfp[0])# 
    
    window=1.0  #Secs
    
    #Make list of rec_indexes
    rec_indexes=[]
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        
    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes

    #Pick a particular LFP event to lock to
    for ss in range(start_lfp, end_lfp, 1):
        print "LFP Cluster # : ", ss+1, " / ", len(Sorts_lfp[0])

        rec_counter=0
        fit_sum = np.zeros((len(Sorts_lfp[rec_index]),total_units,2000), dtype=np.float32)
        latencies=np.zeros((total_units, 2), dtype=np.float32)
        for rec_index in rec_indexes:
            print "...rec: ", rec_index+1, " / ", len(rec_indexes)

            #Compute periods of synchrony from si index
            data_in = lfp.data[rec_index][9]
            si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)

            #Load pop events during synch periods (secs)
            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3  
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
                temp_list.extend(pop_spikes[indexes])
            pop_spikes=np.array(temp_list)
            
            print " ... # pop events: ", len(pop_spikes)
            
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                temp_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    temp_list.extend(spike_array[indexes])
                spike_array=np.array(temp_list)
                
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)
                
                if len(xx_even)>0: xx_even = np.hstack(xx_even)
                if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)
                if len(xx1)>0: xx1 = np.hstack(xx1)

                #Convolve all spike train data w. 20ms gaussian  - Martin's suggestion; also in Luczak 2007, 2009;
                sig = 20
                fit_even = np.zeros(2000, dtype=np.float32)
                if len(xx_even)>0 : 
                    for g in range(len(xx_even)):
                        mu = np.array(xx_even[g]) 
                        fit_even += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)

                fit_odd = np.zeros(2000, dtype=np.float32)
                if len(xx_odd)>0 : 
                    for g in range(len(xx_odd)):
                        mu = xx_odd[g] 
                        fit_odd += gaussian(np.linspace(-1E3,1E3, 2000), mu, sig)

                #Plot control plots
                if False:
                    plt.plot(fit_even[960:1040], color='blue')
                    plt.plot(fit_odd[960:1040], color='red')
                    plt.plot([40,40],[0,max(np.max(fit_even),np.max(fit_odd))*1.2], 'r--', linewidth=3, color='black')
                    plt.title(track_name)
                    plt.show()

                fit_sum[ss][Sorts_sua[rec_index].uid[unit]] += (fit_even + fit_odd)/2.
        
            #Plot all units MSL for each reocrding
            n_lfp_clusters = end_lfp-start_lfp
            ax=plt.subplot(n_lfp_clusters,len(rec_indexes),(ss-start_lfp)*len(rec_indexes)+rec_index+1)
            img=[]
            lock_time=[]
            end_lock = 80
            start_lock = -end_lock
            firing_rate = []
            depth = []
            for unit in range(len(Sorts_sua[rec_index].units)):
                if np.max(fit_sum[ss][Sorts_sua[rec_index].uid[unit]][1000+start_lock:1000+end_lock])>0:
                    f_rate = len(Sorts_sua[rec_index].units[unit])/ lfp.t_end[rec_index]
                    if f_rate<1E-1: continue #Exclude cells firing slower than 0.1Hz
                    lock_time.append(np.argmax(fit_sum[ss][Sorts_sua[rec_index].uid[unit]][1000+start_lock:1000+end_lock]))
                    firing_rate.append(f_rate)
                    depth_temp = Sorts_sua[rec_index].chanpos[Sorts_sua[rec_index].maxchan[unit]][1]
                    print depth_temp
                    depth.append(depth_temp)
                    img.append(fit_sum[ss][Sorts_sua[rec_index].uid[unit]][1000+start_lock:1000+end_lock]/max(fit_sum[ss][Sorts_sua[rec_index].uid[unit]][1000+start_lock:1000+end_lock]))

            #inds = np.array(lock_time).argsort()
            #inds = np.array(firing_rate).argsort()
            inds = np.array(depth).argsort()
            img=np.array(img)[inds]
                
            for k in range(len(inds)-1):
                plt.plot([lock_time[inds[k]],lock_time[inds[k+1]]],[k,k+1], linewidth=10,color='black')

            print len(inds)
            im = ax.imshow(img, origin='upper', extent=[0,end_lock-start_lock, len(inds),0], aspect='auto', interpolation='none', alpha=0.5)
            #plt.plot([100,100],[1,total_units-1], 'r--', linewidth=2, color='white')
            xx = np.arange(0,end_lock-start_lock+1,end_lock)
            x_label = np.arange(-end_lock,end_lock-start_lock+1,end_lock)
            plt.ylim(len(inds),0)
            plt.xticks(xx,x_label, fontsize=25)
            plt.xlabel("Time (ms)", fontsize=30)
            #plt.ylabel("Neuron (lock order)", fontsize=30)
            #plt.ylabel("Neuron (firing rate order)", fontsize=25)
            plt.ylabel("Neuron (depth order)", fontsize=25)
            plt.plot([end_lock,end_lock],[0,total_units], 'r--', linewidth=4, color='white')
            ax.tick_params(axis='both', which='both', labelsize=25)
            ax.xaxis.labelpad = 0
            plt.title(track_name+ "  rec: "+str(recs[rec_index]), fontsize=25)
    plt.show()


def Plot_msl_compare_across_time(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):
    '''Align msl latencies by one LFP cluster and plot all others;
    '''
    
    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.0

    start_lfp = 3
    
    window=1.0  #Secs
    
    #Make list of rec_indexes
    rec_indexes=[]
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        
    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes

    #Fit all data across all lfp events; 
    fit_sum = np.zeros((1000,total_units,2000), dtype=np.float32)
    fit_even_array = np.zeros((1000,total_units,2000), dtype=np.float32)
    fit_odd_array = np.zeros((1000,total_units,2000), dtype=np.float32)

    rec_counter=0
    latencies=np.zeros((total_units, 2), dtype=np.float32)
    
    all_chunk_list = []
    all_chunk_events = []
    chunk_counter=0
    for rec_index in rec_indexes:
        print "...rec: ", rec_index+1, " / ", len(rec_indexes)

        #Compute periods of synchrony from si index
        data_in = lfp.data[rec_index][9]
        si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)
        print sync_periods

        #Load pop events during synch periods (secs)
        lfp_spikes = Sorts_lfp[rec_index][start_lfp]*1E-3  

        #Load data from sync periods into time chunks of X minutes
        chunk = 300 #seconds;
        chunk_list = []
        for p in range(len(sync_periods)):
            for c in range(1000):
                if (sync_periods[p][0]+(c+1)*chunk)>(sync_periods[p][1]): break  #Exit loop if chunked past the end point
                indexes = np.where(np.logical_and(lfp_spikes>=sync_periods[p][0]+c*chunk, lfp_spikes<=sync_periods[p][0]+(c+1)*chunk))[0]
                chunk_list.append(lfp_spikes[indexes])
                all_chunk_list.append(lfp_spikes[indexes])
        print "No. of ", chunk, "sec chunks = ", len(chunk_list)
        
        for t in range(len(chunk_list)):
            pop_spikes=np.array(chunk_list[t])
            
            print " ... chunk: ", t, " #events: ", len(pop_spikes)
            all_chunk_events.append(len(pop_spikes))
        
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                spike_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    spike_list.extend(spike_array[indexes])
                spike_array=np.array(spike_list)
                
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms

                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)
                
                if len(xx_even)>0: xx_even = np.hstack(xx_even)
                if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)

                #Convolve all spike train data w. 20ms gaussian  - Martin's suggestion; also in Luczak 2007, 2009;
                sig = 20
                fit_even = np.zeros(2000, dtype=np.float32)
                if len(xx_even)>0 : 
                    for g in range(len(xx_even)):
                        mu = np.array(xx_even[g]) 
                        fit_even += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)
                fit_even_array[chunk_counter][unique_unit] = fit_even
                
                fit_odd = np.zeros(2000, dtype=np.float32)
                if len(xx_odd)>0 : 
                    for g in range(len(xx_odd)):
                        mu = xx_odd[g] 
                        fit_odd += gaussian(np.linspace(-1E3,1E3, 2000), mu, sig)
                fit_odd_array[chunk_counter][unique_unit] = fit_odd

                fit_sum[chunk_counter][unique_unit] += (fit_even + fit_odd)/2.
    
            chunk_counter+=1
            
    #Plot msl using order from a particular cluster
    n_chunks = len(all_chunk_list)
    
    #Pick an chunk period to order rest of data: 
    ss = 7 
    img=[]
    lock_time=[]
    end_lock = 80
    start_lock = -end_lock
    for unit in range(total_units):
        if np.max(fit_sum[ss][unit][1000+start_lock:1000+end_lock])>0:
            lock_time.append(np.argmax(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
            img.append(fit_sum[ss][unit][1000+start_lock:1000+end_lock]/max(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
        else:
            lock_time.append(200)
            img.append(np.zeros(end_lock-start_lock, dtype=np.float32))
    inds = np.array(lock_time).argsort()

    for t in range(len(all_chunk_list)):
        ax=plt.subplot(2,n_chunks,t+1)
        img=[]
        lock_time=[]
        for unit in range(total_units):  #Looping over all units (even if 0 spikes occur)
            if np.max(fit_sum[t][unit][1000+start_lock:1000+end_lock])>0:
                lock_time.append(np.argmax(fit_sum[t][unit][1000+start_lock:1000+end_lock]))
                img.append(fit_sum[t][unit][1000+start_lock:1000+end_lock]/max(fit_sum[t][unit][1000+start_lock:1000+end_lock]))
            else:
                lock_time.append(200)
                img.append(np.zeros(end_lock-start_lock, dtype=np.float32))
        
        img=np.array(img)[inds]
        temp_img = []
        for i in range(len(img)):
            if np.max(img[i])!=0:
                temp_img.append(img[i])
        img=np.array(temp_img)

        print img.shape
        if img.shape[0]!=0: im = ax.imshow(img, origin='upper', extent=[0,end_lock-start_lock, len(img),0], aspect='auto', interpolation='none')
        #plt.plot([100,100],[1,total_units-1], 'r--', linewidth=2, color='white')

        xx = np.arange(0,end_lock-start_lock+1,end_lock)
        x_label = np.arange(-end_lock,end_lock-start_lock+1,end_lock)
        plt.xticks(xx,x_label, fontsize=25)

        yy = np.arange(0,len(img),10)
        y_label = np.arange(0,len(img),10)
        plt.yticks(yy, y_label, fontsize=25)
        
        plt.ylim(len(inds),0)
        
        plt.xlabel("Time (ms)", fontsize=18)
        plt.plot([end_lock,end_lock],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=18)
        ax.xaxis.labelpad = 0
        plt.title("#: "+str(all_chunk_events[t]), fontsize=25)


        #Plot latency 1st half vs 2nd:
        ax=plt.subplot(2,n_chunks,n_chunks+t+1)
        lock_width = 80
        for unit in range(len(Sorts_sua[rec_index].units)):
            if (np.max(fit_odd_array[t][unit])!=0) and (np.max(fit_even_array[t][unit])!=0):
                #if float(len(Sorts_sua[rec_index].units[unit]))/lfp.t_end[rec_index]>1E-3:      #Exclude cells firing <0.1Hz.
                plt.scatter(np.argmax(fit_odd_array[t][unit][1000-lock_width:1000+lock_width]), 
                np.argmax(fit_even_array[t][unit][1000-lock_width:1000+lock_width]), s=50, color=colors[t%10])
        plt.plot([0,lock_width*2],[0,lock_width*2], 'r--', linewidth=3, color='black', alpha=0.6)

        xx = np.arange(0,lock_width*2+1,lock_width)
        x_label = np.arange(-lock_width,lock_width*2+1,lock_width)
        plt.xticks(xx,x_label, fontsize=25)

        yy = np.arange(0,lock_width*2+1,lock_width)
        y_label = np.arange(-lock_width,lock_width*2+1,lock_width)
        plt.yticks(yy,y_label, fontsize=25)
        
        plt.xlim(0,lock_width*2)
        plt.ylim(0,lock_width*2)

        #plt.xlabel("Latency - 1st half (ms)", fontsize=30)
        #plt.ylabel("Latency - 2nd half (ms)", fontsize=30)
        plt.tick_params(axis='both', which='both', labelsize=18)
        #ax.xaxis.labelpad = 0
        #plt.title(track_name+ "  rec: "+str(recs[rec_index]+"\n"+"LFP "+str(ss)+ ", #: "+str(len(pop_spikes))), fontsize=25)

    plt.suptitle(track_name)
    plt.show()



def Plot_msl_compare_across_lfp_clusters(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):
    '''Align msl latencies by one LFP cluster and plot all others;
    '''
    
    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    start_lfp = 0
    end_lfp = len(Sorts_lfp[0])# 
    
    window=1.0  #Secs
    
    #Make list of rec_indexes
    rec_indexes=[]
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        
    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes

    #Pick a particular LFP event to lock to
    fit_sum = np.zeros((end_lfp-start_lfp,total_units,2000), dtype=np.float32)

    for ss in range(start_lfp, end_lfp, 1):
        print "LFP Cluster # : ", ss+1, " / ", len(Sorts_lfp[0])

        rec_counter=0
        latencies=np.zeros((total_units, 2), dtype=np.float32)
        for rec_index in rec_indexes:
            print "...rec: ", rec_index+1, " / ", len(rec_indexes)

            #Compute periods of synchrony from si index
            data_in = lfp.data[rec_index][9]
            si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)

            #Load pop events during synch periods (secs)
            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3  
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
                temp_list.extend(pop_spikes[indexes])
            pop_spikes=np.array(temp_list)
            
            print " ... # pop events: ", len(pop_spikes)
            
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                temp_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    temp_list.extend(spike_array[indexes])
                spike_array=np.array(temp_list)
                
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)
                
                if len(xx_even)>0: xx_even = np.hstack(xx_even)
                if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)
                if len(xx1)>0: xx1 = np.hstack(xx1)

                #Convolve all spike train data w. 20ms gaussian  - Martin's suggestion; also in Luczak 2007, 2009;
                sig = 20
                fit_even = np.zeros(2000, dtype=np.float32)
                if len(xx_even)>0 : 
                    for g in range(len(xx_even)):
                        mu = np.array(xx_even[g]) 
                        fit_even += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)

                fit_odd = np.zeros(2000, dtype=np.float32)
                if len(xx_odd)>0 : 
                    for g in range(len(xx_odd)):
                        mu = xx_odd[g] 
                        fit_odd += gaussian(np.linspace(-1E3,1E3, 2000), mu, sig)

                #Plot control plots
                if False:
                    plt.plot(fit_even[960:1040], color='blue')
                    plt.plot(fit_odd[960:1040], color='red')
                    plt.plot([40,40],[0,max(np.max(fit_even),np.max(fit_odd))*1.2], 'r--', linewidth=3, color='black')
                    plt.title(track_name)
                    plt.show()

                fit_sum[ss][unique_unit] += (fit_even + fit_odd)/2.
    
    #Plot msl using order from a particular cluster
    n_lfp_clusters = end_lfp-start_lfp
    
    #Pick an lfp cluster to order rest of data: 
    ss = 1 
    img=[]
    lock_time=[]
    end_lock = 80
    start_lock = -end_lock
    for unit in range(total_units):
        if np.max(fit_sum[ss][unit][1000+start_lock:1000+end_lock])>0:
            lock_time.append(np.argmax(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
            img.append(fit_sum[ss][unit][1000+start_lock:1000+end_lock]/max(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
        else:
            lock_time.append(200)
            img.append(np.zeros(end_lock-start_lock, dtype=np.float32))
    inds = np.array(lock_time).argsort()

    for ss in range(start_lfp, end_lfp, 1):
        print ss-start_lfp+1
        ax=plt.subplot(2,n_lfp_clusters,ss-start_lfp+1)
        img=[]
        lock_time=[]
        for unit in range(total_units):  #Looping over all units (even if 0 spikes occur)
            if np.max(fit_sum[ss][unit][1000+start_lock:1000+end_lock])>0:
                lock_time.append(np.argmax(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
                img.append(fit_sum[ss][unit][1000+start_lock:1000+end_lock]/max(fit_sum[ss][unit][1000+start_lock:1000+end_lock]))
            else:
                lock_time.append(200)
                img.append(np.zeros(end_lock-start_lock, dtype=np.float32))
        
        img=np.array(img)[inds]
        temp_img = []
        for i in range(len(img)):
            if np.max(img[i])!=0:
                temp_img.append(img[i])
        img=np.array(temp_img)

        im = ax.imshow(img, origin='upper', extent=[0,end_lock-start_lock, len(img),0], aspect='auto', interpolation='none')
        #plt.plot([100,100],[1,total_units-1], 'r--', linewidth=2, color='white')

        xx = np.arange(0,end_lock-start_lock+1,end_lock)
        x_label = np.arange(-end_lock,end_lock-start_lock+1,end_lock)
        plt.xticks(xx,x_label, fontsize=25)

        yy = np.arange(0,len(img),10)
        y_label = np.arange(0,len(img),10)
        plt.yticks(yy, y_label, fontsize=25)
        
        plt.ylim(len(inds),0)
        
        plt.xlabel("Time (ms)", fontsize=30)
        #plt.ylabel("Neuron (lock order)", fontsize=30)
        plt.plot([end_lock,end_lock],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.xaxis.labelpad = 0
        plt.title("lfp #: "+str(ss), fontsize=25)
    plt.show()


def Plot_msl_controls(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):
    '''Plot msl control comparisons; split data into odd/even or 1st and 2nd halves;
    nb: skipping cell firing <0.1Hz looks better; but not clear if useful;
    '''
    
    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    start_lfp = 0#
    end_lfp = len(Sorts_lfp[0])# 
    
    window=1.0  #Secs
    
    #Make list of rec_indexes
    rec_indexes=[]
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        
    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes

    #Pick a particular LFP event to lock to
    latencies_ave = np.zeros((end_lfp-start_lfp,total_units), dtype=np.float32)
    for ss in range(start_lfp, end_lfp, 1):
        print "LFP Cluster # : ", ss+1, " / ", len(Sorts_lfp[0])

        rec_counter=0
        fit_sum = np.zeros((len(Sorts_lfp[rec_index]),total_units,2000), dtype=np.float32)
        latencies=np.zeros((total_units, 2), dtype=np.float32)
        for rec_index in rec_indexes:
            print "...rec: ", rec_index+1, " / ", len(rec_indexes)

            #Compute periods of synchrony from si index
            data_in = lfp.data[rec_index][9]
            si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)

            #Load pop events during synch periods (secs)
            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3  
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
                temp_list.extend(pop_spikes[indexes])
            pop_spikes=np.array(temp_list)
            
            print " ... # pop events: ", len(pop_spikes)
            
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                #Load unique track-wide unit id 
                unique_unit = Sorts_sua[rec_index].uid[unit]

                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                temp_list = []
                for p in range(len(sync_periods)):
                    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                    temp_list.extend(spike_array[indexes])
                spike_array=np.array(temp_list)
                
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                for j in range(len(pop_spikes)):
                #for j in range(30):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    
                    #Split into even and odd LFP events; VS: Split into 1st half and 2nd half events:
                    if True:
                        if j % 2 == 0:
                            xx_even.append(x)
                        else:
                            xx_odd.append(x)
                        scatter_color='red'
                    else:
                        if j < (float(len(pop_spikes))/2):
                            xx_even.append(x)
                        else:
                            xx_odd.append(x)
                        scatter_color='blue'

                if len(xx_even)>0: xx_even = np.hstack(xx_even)
                if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)
                
                #Convolve all spike train data w. 20ms gaussian  - Martin's suggestion; also in Luczak 2007, 2009;
                sig = 20
                fit_even = np.zeros(2000, dtype=np.float32)
                if len(xx_even)>0 : 
                    for g in range(len(xx_even)):
                        mu = np.array(xx_even[g]) 
                        fit_even += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)

                fit_odd = np.zeros(2000, dtype=np.float32)
                if len(xx_odd)>0 : 
                    for g in range(len(xx_odd)):
                        mu = xx_odd[g] 
                        fit_odd += gaussian(np.linspace(-1E3,1E3, 2000), mu, sig)

                fit_sum[ss][unique_unit] += (fit_even + fit_odd)/2.
                #plt.plot(fit_sum[ss][Sorts_sua[rec_index].uid[unit]])
                #plt.show()
                
                end_lock = 80
                start_lock = -end_lock
                #Save latencies for controls:
                latencies[unit][0]=np.argmax(fit_even[1000+start_lock:1000+end_lock])
                latencies[unit][1]=np.argmax(fit_odd[1000+start_lock:1000+end_lock])
                latencies_ave[ss][unique_unit] = np.argmax(fit_sum[ss][unique_unit][1000+start_lock:1000+end_lock])
        
            #Plot latency 1st half vs 2nd:
            ax = plt.subplot(2,end_lfp-start_lfp,ss+1)
            for unit in range(len(Sorts_sua[rec_index].units)):
                if (latencies[unit][0]!=0) and (latencies[unit][1]!=0):
                    if float(len(Sorts_sua[rec_index].units[unit]))/lfp.t_end[rec_index]>1E-3:      #Exclude cells firing <0.1Hz.
                        plt.scatter(latencies[unit][0], latencies[unit][1], s=50, color=colors[ss])
            plt.plot([0,end_lock*2],[0,end_lock*2], 'r--', linewidth=3, color='black', alpha=0.6)

            xx = np.arange(0,end_lock-start_lock+1,end_lock)
            x_label = np.arange(-end_lock,end_lock-start_lock+1,end_lock)
            plt.xticks(xx,x_label, fontsize=25)

            yy = np.arange(0,end_lock-start_lock+1,end_lock)
            y_label = np.arange(-end_lock,end_lock-start_lock+1,end_lock)
            plt.yticks(yy,y_label, fontsize=25)
            
            plt.xlim(0,end_lock*2)
            plt.ylim(0,end_lock*2)

            #plt.xlabel("Latency - 1st half (ms)", fontsize=30)
            #plt.ylabel("Latency - 2nd half (ms)", fontsize=30)
            plt.tick_params(axis='both', which='both', labelsize=25)
            #ax.xaxis.labelpad = 0
            plt.title(track_name+ "  rec: "+str(recs[rec_index]+"\n"+"LFP "+str(ss)+ ", #: "+str(len(pop_spikes))), fontsize=25)
    plt.show()
    
    #Plot location of each cell across multiple LFP clusters
    counter=0
    for unit in range(total_units):
        line=[]
        for ss in range(start_lfp, end_lfp, 1):
            if latencies_ave[ss][unit]!=0.:
                #x_loc = unit*50 + latencies_ave[ss][unit]
                x_loc = counter*50 + latencies_ave[ss][unit]
                y_loc = 5*ss
                line.append([x_loc, y_loc])
        line=np.array(line)
        if len(line)==(end_lfp-start_lfp):
            for ss in range(start_lfp, end_lfp, 1):
                plt.scatter(line[ss][0],line[ss][1], color=colors[ss], s = 200)
            plt.plot(line[:,0], line[:,1], linewidth=2, color='black')
            plt.plot([counter*50+80,counter*50+80], [0,(ss)*5], 'r--', linewidth=2, color='black')
            counter+=1
        
    xx = np.arange(80,50*counter+51,50)
    x_label = np.arange(1,counter+1,1)
    plt.xticks(xx,x_label, fontsize=15)

    yy = np.arange(0,ss*5+5,5)
    y_label = np.arange(1,end_lfp-start_lfp+1,1)
    plt.yticks(yy,y_label, fontsize=25)
    plt.xlabel("Unit ID", fontsize=25)
    plt.ylabel("LFP Cluster", fontsize=25)
            
    plt.show()




def Plot_mls_lfp_trigger(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):

    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    window=1
    small_window = 100 # ms of window for luczak plots
    large_window = window*1E3
    start_lfp = 0
    end_lfp = len(Sorts_lfp[0])# 
    print start_lfp, end_lfp
    print "# lfp cluters: ", end_lfp
    n_recs = len(recs)
    plotting=True
    start_window = -0.030
    end_window = 0.010
    
    #Convert rec name from string to relative index in concatenated data; 
    #also compute # recordings w. minimum required lfp events
    rec_indexes=[]
    n_lfprecs = 0
    for rec in recs:
        rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        if len(Sorts_lfp[lfp.rec.index(rec)][start_lfp])>5: n_lfprecs+=1

    if n_lfprecs==0:
        print "No lfp_event recordings!"
        return

    #compute track recording length:
    track_length = 0
    for rec_index in rec_indexes:
        track_length+= lfp.t_end[rec_index] #Recording length in seconds

    #Search all SUA sorts to find max # units 
    total_units = 0
    for rec_index in rec_indexes:
        if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    total_units +=1 #Adjust for zero based indexes
    
    #Make array to collect locked spikes for each unit across all recs
    total_locked_spikes_allrecs=np.zeros((end_lfp-start_lfp, total_units),dtype=np.float32)
    
    #Make array to collect all single unit spikes:
    sua_allspikes = np.zeros(total_units, dtype=np.float32)
    for rec_index in rec_indexes:
        for unit in range(len(Sorts_sua[rec_index].units)):
            sua_allspikes[Sorts_sua[rec_index].uid[unit]]+=len(Sorts_sua[rec_index].units[unit])
                
    #Loop over all recordings and plot data where # lfp events > minimum (5)
    gs = gridspec.GridSpec(end_lfp-start_lfp+1,n_lfprecs)       #Make rows = # of lfp clusters + 1 additional row for summary plots

    #Pick a particular LFP event to lock to
    all_pop_spikes = 0 #add up all pop spikes over all recordings to use as control for final plot
    for ss in range(start_lfp, end_lfp, 1):
        print " pop event# : ", ss+1, " / ", len(Sorts_lfp[0])
        
        #Find total # of pop events for cluster ss across all recordings
        total_pop_spikes = 0
        for rec_index in rec_indexes:
            total_pop_spikes += len(Sorts_lfp[rec_index][ss])
        
        all_pop_spikes+=total_pop_spikes
        
        rec_counter=0
        for rec_index in rec_indexes:
            print "...rec: ", rec_index, " / ", len(rec_indexes)
            if len(Sorts_lfp[rec_index][start_lfp])>5:  #If there are < 5 lfp events - skip recording;
                pass
            else:
                continue

            peaks_depth=np.zeros((total_units,len(Sorts_lfp[rec_index])), dtype=np.float32)+100000. #Set all peaks to very large number
            peaks_order=np.zeros((total_units,len(Sorts_lfp[rec_index])), dtype=np.float32)+100000. #Set all peaks to very large number
            peaks_pval=np.zeros((total_units,len(Sorts_lfp[rec_index])), dtype=np.float32)+100000. #Set all peaks to very large number

            fit_sum = np.zeros((len(Sorts_lfp[rec_index]),total_units,1000), dtype=np.float32)
            p_val_array=np.zeros((len(Sorts_lfp[rec_index]),total_units), dtype=np.float64)+1. #Initialize p_val array

            #Make lists to hold unique spikes that lock to lfp events
            n_lock_spikes = [[] for x in xrange(total_units)]

            rec_length = lfp.t_end[rec_index] #Recording length in seconds

            pop_spikes = Sorts_lfp[rec_index][ss]*1E-3
            
            #Loop over all single units for each recording
            for unit in range(len(Sorts_sua[rec_index].units)):
                spike_array = np.array(Sorts_sua[rec_index].units[unit],dtype=float32)/float(Sorts_sua[rec_index].samplerate)  #This converts to seconds
                
                #NB: May wish to remove LARGE SPIKE TIMES due to bug 1E+8 for secs should do it.
                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                control_list=[]     #Collect all spikes from unit that locked to control spike
                control_n_spikes=0. #number of spikes locked to control lfp event
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue
                    
                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]+start_window, spike_array<=pop_spikes[j]+end_window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)

                    #Add # spikes occuring w/in 25ms of lfp event
                    #n_lock_spikes[unit].extend(spike_array[np.where(np.logical_and(spike_array>=pop_spikes[j]-lock_window, spike_array<=pop_spikes[j]+lock_window))[0]])
                    n_lock_spikes[unit].extend(spike_array[np.where(np.logical_and(spike_array>=pop_spikes[j]+start_window, spike_array<=pop_spikes[j]+end_window))[0]])


            #Plot Luczak plots - only locking units - FIRING Sequence Order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+1:int(sqrt(len(Sort1.units)))+5, 3:6])
                ax1.set_xlim(0,small_window*2)
                ax1.set_ylim(0,len(img))
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)
                                    
            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                    img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            for m in range(len(inds)):
                peaks_order[m][ss]=lock_time[inds[m]]

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (firing order w/in small window)")

            ##Plot Luczak plots - only locking units - Firing Sequence Order - PASSED PVAL TEST
            #if plotting:
                #ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+1:int(sqrt(len(Sort1.units)))+5, 6:9])
                #ax1.set_xlim(0,small_window*2)
                #ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            #img=[]
            #lock_time=[]
            #for unit in range(len(Sort1.units)):
                #if p_val_array[ss][unit]<= 0.01:
                    #lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                    #img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))

            #inds = np.array(lock_time).argsort()
            #img=np.array(img)[inds]

            #for m in range(len(inds)):
                #peaks_pval[m][ss]=lock_time[inds[m]]

            #if plotting and len(img)>0:
                #im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                #ax1.xaxis.labelpad = 0
                #plt.title("P_val < 0.01 units (firing order w/in small window)")
                #ax1.set_ylim(0,len(img))


            #Plot large windows *******************************************************
            #Depth
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 0:3])
                ax1.set_xlim(0,large_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (depth order; large window)")
                ax1.set_ylim(0,len(img))

            #Plot Luczak plots - only locking units - by firing order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 3:6])
                ax1.set_xlim(0,large_window*2)
                ax1.set_ylim(0,len(img))
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    lock_time.append(np.argmax(fit_sum[ss][unit]))
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            if plotting:
                if len(img)>0:
                    im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (firing order; large window)")

            #Plot Luczak plots - only locking units - by firing order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 6:9])
                ax1.set_xlim(0,large_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if p_val_array[ss][unit]<= 0.01:
                    lock_time.append(np.argmax(fit_sum[ss][unit]))
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')                
                ax1.xaxis.labelpad = 0
                plt.title("P_val < 0.01 units (depth order; large window)")
                ax1.set_ylim(0,len(img))

            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax1 = plt.subplot(gs[4:7, 11:14])
                xx1=[]
                for p in range(len(Sort1.units)):
                    ax1.scatter(peaks_depth[p][ss],-p*5, color='blue')
                    xx1.append(peaks_depth[p][ss])
                yy1 = np.arange(0,-len(Sort1.units)*5, -5)
                ax1.plot(xx1,yy1, color='blue')
                ax1.set_xlim(0, 2*small_window)
                ax1.set_ylabel("Depth Order")
                
            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax2 = plt.subplot(gs[7:10, 11:14])
                for p in range(len(Sort1.units)):
                    ax2.scatter(peaks_order[p][ss],-p*5, color='red')
                ax2.set_xlim(0, 2*small_window)                    
                ax2.set_ylabel("Firing Order")
                
            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax3 = plt.subplot(gs[10:13, 11:14])
                for p in range(len(Sort1.units)):
                    ax3.scatter(peaks_pval[p][ss],-p*5, color='black')
                ax3.set_xlim(0, 2*small_window)          
                ax3.set_ylabel("Firing Order (pval)")

            #SHOW PLOTS
            if plotting: 
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
                plt.show()    

    plt.show()

    if False:
        #PLotting cumulative locking times: *******************************************
        gs = gridspec.GridSpec(14,18)

        plt.suptitle(Sort1.directory[27:] + "  " + 
        str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
        + " mins) \n " + "window (x-axis): "+str(window*2*1E3) + " ms,   max spike count (y-axis): 100 spikes.", fontsize=20)

        ax1 = plt.subplot(gs[0:4, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_depth[p][ss]<1000:
                    ax1.scatter(peaks_depth[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_depth[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            yy1 = np.arange(0,-len(Sort1.units), -1)
            ax1.plot(xx1,yy1, color='blue')
            ax1.xaxis.set_visible(False)
            ax1.set_ylabel("Depth Order")
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)


        ax1 = plt.subplot(gs[4:8, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_order[p][ss]<1000:
                    ax1.scatter(peaks_order[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_order[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            yy1 = np.arange(0,-len(Sort1.units), -1)
            ax1.plot(xx1,yy1, color='red')
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)
            ax1.xaxis.set_visible(False)
            ax1.set_ylabel("Firing Order")

        ax1 = plt.subplot(gs[8:12, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_pval[p][ss]<1000:
                    ax1.scatter(peaks_pval[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_pval[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)
            yy1 = np.arange(0,-len(Sort1.units), -1)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            if all(x is None for x in xx1):
                plt.plot([100,100],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)
            else:
                ax1.plot(xx1,yy1, color='black')
            ax1.set_ylabel("Firing Order (pval)")
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)

        #SHOW PLOTS
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
        plt.show()   


def Plot_triggered_activity(Sort1,Sort2,lfp):

    #PLOT SINGLE UNIT RESPONSE TO EACH POP SPIKE
    unit_events = False     #Plot Single Unit response to all LFP events
    lfp_events = True       #Plot all unit responses to a single LFP event
    plotting = True     #Show plots
    start_unit = 0      #Starting single unit
    start_lfp = 10       #Starting lfp event
    window=1            #Window to search in seconds

    print "Loading TSF files"
    #Load required .tsf files for Sort1 and Sort2#Load high-pass filtered traces
    tsf_name = Sort1.name
    sim_dir = Sort1.directory
    if 'tim' in Sort1.directory:
        f = open(sim_dir + tsf_name + '_hp.tsf', "rb")
    else:
        f = open(sim_dir + tsf_name + '.tsf', "rb")

    Sort1.tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others

    #Load LFP traces
    tsf_name = Sort2.name
    sim_dir = Sort2.directory
    f = open(sim_dir + tsf_name + '.tsf', "rb")
    Sort2.tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others

    ss_lock_array=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float64)
    pie_array1=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float64)
    pie_array2=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float64)
    p_val_array=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float64)

    #Plot individual unit responses vs. looping over each 
    if unit_events:
        #Loop over all units - one plot for each unit
        for unit in range(start_unit, len(Sort1.units), 1):
            if plotting: 
                gs = gridspec.GridSpec(10,max(10,len(Sort2.units)+1))
                if len(Sort2.units)>10:
                    print "TOO MANY LFP UNITS FOR DISPLAY"
                    #quit()
                plt.suptitle(Sort1.directory[27:] + "      Recording Length: " + str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
                + " mins.", fontsize=20)
                
           #PLOT RASTERS + PSTHs + SPECGRAM 
            if plotting: 
                ax_image = plt.subplot(gs[6:,0:]) 
                
                
                img= np.rot90(lfp.specgram[::-1])
                rec_length = Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency
                ax_image.set_xlim(0,1000)
                ax_image.set_ylim(0,rec_length)
                im = ax_image.imshow(img, extent=[0,85,0,rec_length], aspect='normal') #extent=lfp.extent, cmap=lfp.cm)
                ax_image.yaxis.set_visible(False)

            spike_array = np.array(Sort1.units[unit],dtype=float32)/float(Sort1.samplerate)  #This converts to seconds

            #Loop over population spikes
            for ss in range(len(Sort2.units)):
                index_array=[] #Holds the indexes of the matching single units; important to remove duplicates;
                index_array_ttest=[]    
                                    
                #LOAD Pop Spike raster consisting of single-merged Unit
                if 'nick' in Sort2.directory:
                    pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)*100./float(Sort2.samplerate) #Convert to seconds
                else:
                    #Tim's data already converted to absolute timesteps
                    pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)/float(Sort2.samplerate)
            
                #Find max length of recording
                rec_length = max(pop_spikes) #Make sure to normalize depth plots also

                #temp_array = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)  #This converts to seconds
                print "unit: ", unit, " pop spike: ", ss
                #Loop over all population spikes
                xx1=[]
                xx_even=[]
                xx_odd=[]
                xx1_ttest=[]
                for j in range(len(pop_spikes)):

                    #PLOT RASTERS; find all unit spikes that fall w/in +/- window of pop spike
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]
                    
                    for i in temp2: #NB: Not excluding duplicate spikes from counts!!! make sure to take into account for analysis
                        index_array.append(i)

                        x=(spike_array[i]-pop_spikes[j])*1E3 #Offset time of spike to t_0 
                        xx1.append(x)
                        if j % 2 == 0:
                            xx_even.append(x)
                        else:
                            xx_odd.append(x)
                        ymin=pop_spikes[j]
                        ymax=ymin+10
                        if plotting:
                            x_temp = x/(float(window)*1.E3/36.)+40+102.*(ss+1)
                            ax_image.add_artist(Circle(xy=(x_temp,ymin), radius=0.5))

                    #Plot all LFP event spikes 
                    #if plotting: p = ax_image.plot([102.*(ss+1)-2,102*(ss+2)-18],[pop_spikes[j],pop_spikes[j]], linewidth=0.1,color='red', alpha=.25)
                    #p = ax_image.axhline(y=pop_spikes[j],xmin=102.*(ss+1)-2,xmax=102*(ss+2)-18,color='red', linewidth=100)

                ###Loop over ttest spikes
                #for j in range(len(ttest_spikes)):
                    #temp2_ttest = np.where(np.logical_and(spike_array>=ttest_spikes[j]-window, spike_array<=ttest_spikes[j]+window))[0]

                    #for i in temp2_ttest:
                        #if i in index_array_ttest:
                            ##print "Duplicate comparison spike"
                            #pass
                        #else:
                            #index_array_ttest.append(i)
                            #x=(spike_array[i]-ttest_spikes[j])*1E3 #Offset time of spike to t_0 
                            #xx1_ttest.append(x)
                
                #xx1=np.array(xx1)
                #xx1_ttest=np.array(xx1_ttest)
                
                if plotting: p = ax_image.axvspan(-18+102.*(ss+1), -2+102.*(ss+1), facecolor='black', alpha=0.25)

                #Plot PSTH on cumulative data from above: 3 methods available
                #1 Plot gaussian convolved PSTHs - Martin's suggestion
                def gaussian(x, mu, sig):
                    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
                sig = 20/(float(window)*1.E3/36.) #20 ms (scaled to fit window) - Martin's suggestion

                fitdata = np.zeros(1000, dtype=np.float32)
                for g in range(len(xx_even)):
                    mu = xx_even[g]/(float(window)*1.E3/36.) #Default of normalized height
                    fitdata += gaussian(np.linspace(-36.,36., 1000), mu, sig)
                u = np.linspace(-window*1E3, window*1E3, 1000)
                u_temp=u/(float(window)*1.E3/36.)+40+102.*(ss+1) 
                if plotting: ax_image.plot(u_temp,fitdata*5, linewidth=2, color='blue', alpha=0.7)

                fitdata = np.zeros(1000, dtype=np.float32)
                for g in range(len(xx_odd)):
                    mu = xx_odd[g]/(float(window)*1.E3/36.) #Default of normalized height
                    fitdata += gaussian(np.linspace(-36.,36., 1000), mu, sig)
                u = np.linspace(-window*1E3, window*1E3, 1000)
                u_temp=u/(float(window)*1.E3/36.)+40+102.*(ss+1) 
                if plotting: ax_image.plot(u_temp,fitdata*5, linewidth=2, color='red', alpha=0.7)


                #2 Plot bar graphs PSTH
                #bin_width = window/5. #Width of summation bin
                #yy1=np.histogram(np.array(xx1).flatten(), bins = np.arange(-window*1E3,+window*1E3,bin_width*1E3))
                #if plotting: ax_image.bar(yy1[1][0:-2]/(float(window)*1.E3/36.)+bin_width/2./(float(window)/36.)+40+102.*(ss+1),
                #yy1[0][0:-1]*10, bin_width/(float(window)/36.), color='blue', alpha=0.7)

                #3 Plot sinc interpolated PSTHs
                #x=yy1[0][0:-1]
                #s=yy1[1][0:-2]
                #u=np.arange(-window*1E3,+window*1E3,.001*1E3)
                #fitdata=sinc_interp(x,s,u)
                #u_temp = u/(float(window)*1.E3/36.)+bin_width*1E3/2./(float(window)*1.E3/36.)+40+102.*(ss+1) 
                #if plotting: ax_image.plot(u_temp,fitdata*5, linewidth=2, color='blue', alpha=0.7)

                #Plot max interploated value
                if fitdata[np.argmax(fitdata)]>4:  #Search for peaks with more than 5 spikes
                    temp_x = u_temp[np.argmax(fitdata)]
                    temp_y = max(fitdata)
                    #print temp_x, temp_y
                    if plotting: 
                        ax_image.plot([temp_x,temp_x],[0, temp_y+10000], 'r--',color='blue',linewidth=2, alpha=0.5)
                        ax_image.text(102.*(ss+1), -len(img)/13., "Peak: "+str(round(u[np.argmax(fitdata)],2)) + " ms. \n # spk: "+str(len(xx1))
                        + "\n duration: " + str(round(len(pop_spikes)*window*2/60.0,2)) + " mins." 
                        , color = 'black', fontsize = 10)
                    ss_lock_array[unit][ss]=u[np.argmax(fitdata)]
                else:
                    temp_x = u_temp[np.argmax(fitdata)]
                    temp_y = fitdata[np.argmax(fitdata)]
                    if plotting:
                        #ax_image.plot([temp_x,temp_x],[0, temp_y+10000], 'r--',color='blue',linewidth=2, alpha=0.5)
                        ax_image.text(102.*(ss+1), -len(img)/13., "\n # spk: "+str(len(xx1))
                        + "\n duration: " + str(round(len(pop_spikes)*window*2/60.0,2)) + " mins." 
                        , color = 'black', fontsize = 10)                    
                    ss_lock_array[unit][ss]=1000.0

                #Plot pie charts
                if plotting:
                    ax = plt.subplot(gs[5:6,ss+1]) 
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                    ax.set_xlim(-100, 100)
                    ax.set_ylim(-100, 100)
                    x=-65
                    y=40

                colors = ['green', 'black']
                a=min(float(len(xx1))/float(len(Sort1.units[unit])),1.)
                pie_array1[unit][ss]=round(a,4)
                b=1-a
                #if plotting: draw_pie(ax,    [a, b],     x,  y, size=100*6/len(Sort2.units), colors=colors)     

                colors = ['red', 'black']
                a=min(float(len(pop_spikes)*window*2.)/float(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency),1.)
                pie_array2[unit][ss]=round(a,4)
                b=1-a

                #Compute T-test values; #Limit analysis to at 0.1Hz firing cells
                if (len(Sort1.units[unit])/float(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency))> 0.05: 
                    temp_float=0.
                    for r in range(5):
                        xx1_ttest =  (np.random.random_sample((1000,))-0.5)*2*window*1E3
                        KS, p_val = sstats.ks_2samp(np.sort(xx1), np.sort(xx1_ttest)) #2 sample ks test
                        temp_float+=p_val
                    p_val=temp_float/5.
                    
                    p_val_array[unit][ss]=p_val
                    if plotting: 
                        plt.text(-95, -25, "P: %.1e" % p_val, fontsize = 10) #, ha='center', va='bottom')
                        plt.text(-95, 25, "KS: %.1e" % KS, fontsize = 10) #, ha='center', va='bottom')
                else:
                    p_val_array[unit][ss]=1.0
                        
                #if plotting: draw_pie(ax,    [a, b],     x,  y-70, size=100*6/len(Sort2.units), colors=colors)     

                #******PLOT LFP Spikes******* 
                if plotting: 
                    ax = plt.subplot(gs[0:5,ss+1])
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                    plt.title("LFP: " + str(ss)+"\n#: "+str(len(Sort2.units[ss]))
                    + " (" +str(round(float(len(Sort2.units[ss]))/(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency),2)) + "hz)"
                    ,weight = 'bold', color = 'black', fontsize = 11)

                n_points = 200
                x = np.zeros((Sort2.tsf.n_electrodes,n_points),dtype=np.float32)
                for i in range(Sort2.tsf.n_electrodes):
                    x[i]= Sort2.tsf.Siteloc[i*2] + np.array(arange(0,n_points,1))

                y1 = np.zeros((Sort2.tsf.n_electrodes,n_points), dtype=np.float32)
                counter=0
                spike_list = random.sample(range(1, len(Sort2.units[ss])), min(len(Sort2.units[ss]),50)-1)
                #spike_list = np.arange(0,min(len(Sort2.units[ss]),50),1) #Need this if some spikes too close to end of recording

                if plotting: 
                  for j in spike_list: #range(min(len(Sort1.units[unit]),50)): #Plot up to 100 spikes from unit
                    counter+=1
                    for i in range(Sort2.tsf.n_electrodes):
                        if Sort2.maxchan[ss]==i:
                            plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                            Sort2.units[ss][j]+n_points/2]/60.-Sort2.tsf.Siteloc[i*2+1]/5., color='blue', alpha=1.)
                        else:
                            plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                            Sort2.units[ss][j]+n_points/2]/60.-Sort2.tsf.Siteloc[i*2+1]/5., color='black', alpha=0.7)

                if plotting: 
                    ax.set_ylabel("Depth along probe (um)", weight = 'bold', color = 'black', fontsize = 12)
                    ax.set_xlim(0, 200)
                    ax.set_ylim(-20*Sort1.tsf.n_electrodes, 25)

            #Add one last grey bar to raster plots
            if plotting: p = ax_image.axvspan(-18+102.*(ss+2), -2+102.*(ss+8), facecolor='black', alpha=0.25)

            if (len(Sort1.units[unit])/float(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency))< 0.1: #Limit analysis to at 0.1Hz firing cells
                print "Skipping Unit: ", unit, " < 0.05Hz"

            #Plot custom xaxis labels
            x_origin=[]
            x_labels=[]
            for i in range(len(Sort2.units)):
                x_origin.append(40+102.*(i+1))
                x_labels.append(0.0)
                
                x_origin.append(-window*1E3/(float(window)*1.E3/36.)+40+102.*(i+1))
                x_labels.append(-window)

                x_origin.append(window*1E3/(float(window)*1.E3/36.)+40+102.*(i+1))
                x_labels.append(window)
                
            if plotting: 
                ax_image.set_xticks(x_origin)
                ax_image.set_xticklabels(x_labels)
            
            #PLOT SORTED UNIT
            points=40
            if plotting: ax = plt.subplot(gs[0:6,0]) 
            x = np.zeros((Sort1.tsf.n_electrodes,points),dtype=np.float32)
            for i in range(Sort1.tsf.n_electrodes):
                x[i]= Sort1.tsf.Siteloc[i*2] + np.array(arange(0,points,1))
                #print i, len(x[i])
                
            y1 = np.zeros((Sort1.tsf.n_electrodes,points), dtype=np.float32)
            counter=0
            spike_list = random.sample(range(0, len(Sort1.units[unit])), min(len(Sort1.units[unit]),20))

            if plotting: 
                for j in spike_list: #range(min(len(Sort1.units[unit]),50)): #Plot up to 100 spikes from unit
                    counter+=1
                    print "plotting spike: ", j, " time: ", Sort1.units[unit][j]
                    for i in range(Sort1.tsf.n_electrodes):
                        if Sort1.maxchan[unit]==i:
                            xx1 = x[i]
                            yy1 = Sort1.tsf.ec_traces[i][int(Sort1.units[unit][j])-points/2:int(Sort1.units[unit][j])+points/2]/10-Sort1.tsf.Siteloc[i*2+1]
                            plt.plot(xx1, yy1, color='blue', alpha=1.)
                        else:
                            xx1 = x[i]
                            yy1 = Sort1.tsf.ec_traces[i][int(Sort1.units[unit][j])-points/2:int(Sort1.units[unit][j])+points/2]/10-Sort1.tsf.Siteloc[i*2+1]
                            plt.plot(xx1,yy1, color='black', alpha=0.7)

                plt.title("Unit: " + str(unit) + " / "+str(len(Sort1.units))+"\n #spk: "+str(len(Sort1.units[unit])),weight = 'bold', color = 'black', fontsize = 14)
                ax.set_ylabel("Depth along probe (um)", weight = 'bold', color = 'black', fontsize = 12)
                ax.set_xlim(-65, 125)
                ax.set_ylim(top=50)
                plt.gca().xaxis.set_major_locator(plt.NullLocator())         

            #print p_val_array[unit]
            
            #SHOW PLOTS
            if plotting: 
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
                plt.show()
    
    #************************************* LFP BASED PLOTS ***********************************************
    #LFP vs. all unit plots
    window=1
    small_window = 100 # ms of window for luczak plots
    large_window = window*1E3
    start_lfp=0
    plotting=True

    if lfp_events:
        if plotting:
            print "size : ", int(sqrt(len(Sort1.units)))+1, int(sqrt(len(Sort1.units)))+2
            print "n units: ", len(Sort1.units)
            #gs = gridspec.GridSpec(int(sqrt(len(Sort1.units)))+10, int(sqrt(len(Sort1.units)))+3)
            gs = gridspec.GridSpec(25, 30)
                            
        peaks_depth=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float32)+100000. #Set all peaks to very large number
        peaks_order=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float32)+100000. #Set all peaks to very large number
        peaks_pval=np.zeros((len(Sort1.units),len(Sort2.units)), dtype=np.float32)+100000. #Set all peaks to very large number

        fit_sum = np.zeros((len(Sort2.units),len(Sort1.units),1000), dtype=np.float32)
        p_val_array=np.zeros((len(Sort2.units),len(Sort1.units)), dtype=np.float64)+1.

        for ss in range(start_lfp, len(Sort2.units), 1):
            print " pop spike: ", ss+1, " / ", len(Sort2.units), "  #spikes: ", len(Sort2.units[ss])

            rec_length = Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency #Recording length in seconds

            if plotting: 
                plt.suptitle("LFP #: " + str(ss) + " (" + Sort1.directory[27:] + "  " + 
                str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
                + " mins) \n " + "window (x-axis): "+str(window*2*1E3) + " ms,   max spike count (y-axis): 100 spikes.", fontsize=20)
                              
            #LOAD Pop Spike raster consisting of single-merged Unit

            #Cat data  use *100 conversion factor
            if 'nick' in Sort2.directory:
                pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)*100/float(Sort2.samplerate) #Convert to seconds

            #Mouse data use *20 conversion factor - NO! - CONVERT ALL RAW SPIKE STREAMS BACK TO ABSOLUTE TIME
            if 'tim' in Sort2.directory:
                #pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)*20/float(Sort2.samplerate) #Convert to seconds
                pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)/float(Sort2.samplerate) #Convert to seconds
            
            #Loop over all population spikes
            for unit in range(len(Sort1.units)):
                print "plotting unit: ", unit, " #spikes: ", len(Sort1.units[unit])
                
                if plotting: ax_image = plt.subplot(gs[unit/int(sqrt(len(Sort1.units))+1), unit % (int(sqrt(len(Sort1.units)))+1)+1]) 

                spike_array = np.array(Sort1.units[unit],dtype=float32)/float(Sort1.samplerate)  #This converts to seconds

                #print pop_spikes
                #print spike_array
                #quit()

                xx1=[]
                xx_even=[]
                xx_odd=[]
                xx1_ttest=[]
                for j in range(len(pop_spikes)):
                    #PLOT RASTERS; find all unit spikes that fall w/in +/- window of pop spike
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]
                    
                    #NB: Not excluding duplicate spikes from counts!!! make sure to take into account for analysis
                    #if len(temp2)>0:
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)
                            
                    #if plotting:
                    #    ax_image.vlines(x=x, ymin=pop_spikes[j],ymax=pop_spikes[j]+3, color='black',linewidth=2)
                
                xx1 = np.hstack(xx1)
                xx_even= np.hstack(xx_even)
                xx_odd = np.hstack(xx_odd)
                
                #Plot PSTH on cumulative data from above: 3 methods available
                #1 Plot gaussian convolved PSTHs - Martin's suggestion
                def gaussian(x, mu, sig):
                    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
                sig = 20 #20 ms (scaled to fit window) - Martin's suggestion

                fit_even = np.zeros(1000, dtype=np.float32)
                for g in range(len(xx_even)):
                    mu = np.array(xx_even[g])  #Default of normalized mean
                    fit_even += gaussian(np.linspace(-window*1E3, window*1E3, 1000), mu, sig)
                u = np.linspace(-window*1E3, window*1E3, 1000)
                if plotting: ax_image.plot(u,fit_even*10, linewidth=2, color='blue', alpha=0.7)

                fit_odd = np.zeros(1000, dtype=np.float32)
                for g in range(len(xx_odd)):
                    mu = xx_odd[g] #Default of normalized height
                    fit_odd += gaussian(np.linspace(-window*1E3,window*1E3, 1000), mu, sig)
                u = np.linspace(-window*1E3, window*1E3, 1000)
                if plotting: ax_image.plot(u,fit_odd*10, linewidth=2, color='red', alpha=0.7)
                
                fit_sum[ss][unit] = (fit_even + fit_odd)/2.

                #Compute T-test values; #Limit analysis to at 0.05Hz firing cells / try different firing rates and then show % locking
                p_val=1.0
                if (len(Sort1.units[unit])/float(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency))> 0.01: 
                    temp_float=0.
                    for r in range(5):
                        xx1_ttest =  (np.random.random_sample((1000,))-0.5)*2*window*1E3
                        KS, p_val = sstats.ks_2samp(np.sort(xx1), np.sort(xx1_ttest)) #2 sample ks test
                        temp_float+=p_val
                    p_val=temp_float/5.
                    p_val_array[ss][unit]=p_val
                    
                if plotting: 
                    plt.text(0, rec_length, "P: %.1e" % p_val, fontsize = 12)
                    plt.text(-window*1E3+window*1E3/20, rec_length, "U: "+str(unit), fontsize = 12, fontweight='bold') 

                    ax_image.set_xlim(-window*1E3, window*1E3)
                    ax_image.set_ylim(0,rec_length)
                    ax_image.xaxis.set_visible(False)
                    ax_image.yaxis.set_visible(False)

            #******PLOT LFP Spikes******* 
            if plotting: 
                ax1 = plt.subplot(gs[0: unit/int(sqrt(len(Sort1.units)))-1, 0]) 

            n_points = 200
            x = np.zeros((Sort2.tsf.n_electrodes,n_points),dtype=np.float32)
            for i in range(Sort2.tsf.n_electrodes):
                x[i]= Sort2.tsf.Siteloc[i*2] + np.array(arange(0,n_points,1))

            y1 = np.zeros((Sort2.tsf.n_electrodes,n_points), dtype=np.float32)

            spike_list = np.arange(0,min(len(Sort2.units[ss]),50),1) #Need this if some spikes too close to end of recording
            #if plotting: 
                #for j in spike_list: #range(min(len(Sort1.units[unit]),50)): #Plot up to 100 spikes from unit
                    #for i in range(Sort2.tsf.n_electrodes):
                        #if Sort2.maxchan[ss]==i:
                            #plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                            #Sort2.units[ss][j]+n_points/2]/600.-Sort2.tsf.Siteloc[i*2+1], color='blue', alpha=1.)
                        #else:
                            #plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                            #Sort2.units[ss][j]+n_points/2]/600.-Sort2.tsf.Siteloc[i*2+1], color='black', alpha=1.)

            #***********************************************************************************************
            #Make MLS Luczak plots
            #100ms wide windows
            #All units (exclude non-locking units)
            #Plot Luczak plots - only locking units - Original order - DEPTH
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+1:int(sqrt(len(Sort1.units)))+5, 0:3])
                ax1.set_xlim(0,small_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)
                
            img=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))
                
                    #Save peaks in depth order
                    peaks_depth[unit][ss]=(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))

               
            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (depth order)")
                ax1.set_ylim(0,len(img))

            #Plot Luczak plots - only locking units - FIRING Sequence Order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+1:int(sqrt(len(Sort1.units)))+5, 3:6])
                ax1.set_xlim(0,small_window*2)
                ax1.set_ylim(0,len(img))
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)
                                    
            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                    img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            for m in range(len(inds)):
                peaks_order[m][ss]=lock_time[inds[m]]

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (firing order w/in small window)")

            #Plot Luczak plots - only locking units - Firing Sequence Order - PASSED PVAL TEST
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+1:int(sqrt(len(Sort1.units)))+5, 6:9])
                ax1.set_xlim(0,small_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if p_val_array[ss][unit]<= 0.01:
                    lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                    img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            for m in range(len(inds)):
                peaks_pval[m][ss]=lock_time[inds[m]]

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("P_val < 0.01 units (firing order w/in small window)")
                ax1.set_ylim(0,len(img))


            #Plot large windows *******************************************************
            #Depth
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 0:3])
                ax1.set_xlim(0,large_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (depth order; large window)")
                ax1.set_ylim(0,len(img))

            #Plot Luczak plots - only locking units - by firing order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 3:6])
                ax1.set_xlim(0,large_window*2)
                ax1.set_ylim(0,len(img))
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if max(fit_sum[ss][unit])!= 0.:
                    lock_time.append(np.argmax(fit_sum[ss][unit]))
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            if plotting:
                if len(img)>0:
                    im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax1.xaxis.labelpad = 0
                plt.title("All units (firing order; large window)")

            #Plot Luczak plots - only locking units - by firing order
            if plotting:
                ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 6:9])
                ax1.set_xlim(0,large_window*2)
                ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

            img=[]
            lock_time=[]
            for unit in range(len(Sort1.units)):
                if p_val_array[ss][unit]<= 0.01:
                    lock_time.append(np.argmax(fit_sum[ss][unit]))
                    img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

            inds = np.array(lock_time).argsort()
            img=np.array(img)[inds]

            if plotting and len(img)>0:
                im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')                
                ax1.xaxis.labelpad = 0
                plt.title("P_val < 0.01 units (depth order; large window)")
                ax1.set_ylim(0,len(img))

            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax1 = plt.subplot(gs[4:7, 11:14])
                xx1=[]
                for p in range(len(Sort1.units)):
                    ax1.scatter(peaks_depth[p][ss],-p*5, color='blue')
                    xx1.append(peaks_depth[p][ss])
                yy1 = np.arange(0,-len(Sort1.units)*5, -5)
                ax1.plot(xx1,yy1, color='blue')
                ax1.set_xlim(0, 2*small_window)
                ax1.set_ylabel("Depth Order")
                
            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax2 = plt.subplot(gs[7:10, 11:14])
                for p in range(len(Sort1.units)):
                    ax2.scatter(peaks_order[p][ss],-p*5, color='red')
                ax2.set_xlim(0, 2*small_window)                    
                ax2.set_ylabel("Firing Order")
                
            #Plot all unit locking scattering times for each LFP event
            if plotting:
                ax3 = plt.subplot(gs[10:13, 11:14])
                for p in range(len(Sort1.units)):
                    ax3.scatter(peaks_pval[p][ss],-p*5, color='black')
                ax3.set_xlim(0, 2*small_window)          
                ax3.set_ylabel("Firing Order (pval)")

            #SHOW PLOTS
            if plotting: 
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
                plt.show()    

        #PLotting cumulative locking times: *******************************************
        gs = gridspec.GridSpec(14,18)

        plt.suptitle(Sort1.directory[27:] + "  " + 
        str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
        + " mins) \n " + "window (x-axis): "+str(window*2*1E3) + " ms,   max spike count (y-axis): 100 spikes.", fontsize=20)

        ax1 = plt.subplot(gs[0:4, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_depth[p][ss]<1000:
                    ax1.scatter(peaks_depth[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_depth[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            yy1 = np.arange(0,-len(Sort1.units), -1)
            ax1.plot(xx1,yy1, color='blue')
            ax1.xaxis.set_visible(False)
            ax1.set_ylabel("Depth Order")
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)


        ax1 = plt.subplot(gs[4:8, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_order[p][ss]<1000:
                    ax1.scatter(peaks_order[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_order[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            yy1 = np.arange(0,-len(Sort1.units), -1)
            ax1.plot(xx1,yy1, color='red')
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)
            ax1.xaxis.set_visible(False)
            ax1.set_ylabel("Firing Order")

        ax1 = plt.subplot(gs[8:12, 0:18])
        for ss in range(start_lfp, len(Sort2.units), 1):
            xx1=[]
            for p in range(len(Sort1.units)):
                if peaks_pval[p][ss]<1000:
                    ax1.scatter(peaks_pval[p][ss]+ss*small_window*2,-p, color='black')
                    xx1.append(peaks_pval[p][ss]+ss*small_window*2)
                else:
                    xx1.append(None)
            yy1 = np.arange(0,-len(Sort1.units), -1)

            ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

            if all(x is None for x in xx1):
                plt.plot([100,100],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)
            else:
                ax1.plot(xx1,yy1, color='black')
            ax1.set_ylabel("Firing Order (pval)")
            ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)

        #SHOW PLOTS
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
        plt.show()   
            
    quit()


    #********************************* PLOT Statistics on locking units from above ******************************
    
    fig, ax= plt.subplots()
    #plt.clf()
    #ax1 = plt.subplot(1,1,1)
    abs_log_array = np.absolute(np.log10(p_val_array))
    #for i in range(len(p_val_array)):
        #print "row: ", i
        #print p_val_array[i]
        #print abs_log_array[i]
    
    cax = ax.imshow(abs_log_array, cmap=plt.cm.Reds,   #TRY BLACK TO WHITE COLOUR CHART; OR BLUE TO WHITE
                    interpolation='None', vmin=1.3, vmax=6)
    #cbar = fig.colorbar(cax)
    cbar = fig.colorbar(cax, ticks=[2,3,4,5,6])

    cbar.ax.set_yticklabels(['1E-2', '1E-3', '1E-4','1E-5','1E-6'])

    cbar.set_label('KS (2 Sampel) P-Value')

    plt.title("P-Value Log Colour Map \n (i.e. abs(log(p-val)))" ,weight = 'bold', color = 'black', fontsize = 14)

    ax.set_xlabel("LFP Units", weight = 'bold', color = 'black', fontsize = 12)
    ax.set_ylabel("Single Units", weight = 'bold', color = 'black', fontsize = 12)

    #plt.suptitle(Sort1.directory[27:] + "      Recording Length: " + str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
    #+ " mins.", loc='left', fontsize=20)
    left, width = -.50, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    
    ax.text(left, 0.5*(bottom+top), Sort1.directory[27:] + "      Recording Length: " + str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
        + " mins    Window: " + str(window*2E3) + " ms.",
        horizontalalignment='right',
        verticalalignment='center',
        rotation='vertical',
        transform=ax.transAxes)
        
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
    plt.show()

    quit()

    if False:

        gs = gridspec.GridSpec(len(Sort2.units)+2,15)
        plt.suptitle(Sort1.directory[27:] + "      Recording Length: " + str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
        + " mins.", fontsize=20)

        #PLOT complied pie chart data from above plots (Shows % locked responses vs. time) 
        for i in range(len(Sort2.units)):
            ax = plt.subplot(gs[i:i+1,0:6]) 
            if i==0: plt.title("Overall locking percentages (green) and \n overall LFP event duration percentages (red)", weight = 'bold', color = 'black', fontsize = 14)
            ax.xaxis.set_visible(False)
            ax.set_ylabel("LFP #:" + str(i), weight = 'bold', color = 'black', fontsize = 12)

            max_yticks = 2
            yloc = plt.MaxNLocator(max_yticks)
            ax.yaxis.set_major_locator(yloc)
            plt.yticks(fontsize = 10)

            x = np.arange(0,len(Sort1.units),1)
            plt.bar(x,pie_array1.T[i], width=1, color='green')
            plt.bar(x,pie_array2.T[i], width=1, color='red', alpha=0.5)
            plt.ylim(0,1)
            plt.xlim(0,len(Sort1.units))

        #Plot total % response
        ax = plt.subplot(gs[len(Sort2.units):len(Sort2.units)+1,0:6]) 
        for i in range(len(Sort1.units)):
            plt.bar(x[i],sum(pie_array1[i]), width=1, color='green')
            plt.bar(x[i],sum(pie_array2[i]), width=1, color='red', alpha=0.5)
            plt.xlim(0,len(Sort1.units))

        plt.yticks(fontsize = 10)
        plt.xticks(fontsize=10)
        ax.set_ylabel("Totals", weight = 'bold', color = 'black', fontsize = 12)
        ax.xaxis.set_visible(False)

        #Plot individual LFP unit ratios
        ratios=[]
        for i in range(len(Sort2.units)):
            ax = plt.subplot(gs[i:i+1,6:13])
            if i==0: plt.title("Ratios of % of Total Spikes (green) / % of Total Time (red) \n (dashed line = 1, i.e. chance)", weight = 'bold', color = 'black', fontsize = 14)
            temp_ratio=[]
            for j in range(len(pie_array1)):
                a = pie_array1[j][i]
                b = pie_array2[j][i]
                if b ==0: 
                    temp_ratio.append(0)
                else:
                    temp_ratio.append(a/b)
                    
            print len(x), len(temp_ratio)
            plt.bar(x,temp_ratio, width=1, color='cyan')
            ratios.append(temp_ratio)
            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(False)
            plt.xlim(0,len(Sort1.units))
            plt.ylim(0,10)
            plt.plot([0,len(Sort1.units)],[1.0,1.0], 'r--', color = 'black', linewidth=4)
        
        plt.xlim(0,len(Sort1.units))
        plt.plot([0,len(Sort1.units)],[1.0,1.0], 'r--', color = 'black', linewidth=5)

        ax.set_xlabel("Single Unit IDs", weight = 'bold', color = 'black', fontsize = 12)

        #Plot Sum of ratios
        ax = plt.subplot(gs[len(Sort2.units):len(Sort2.units)+1,6:13]) 
        ratios=[]
        for i in range(len(Sort1.units)):
            a = float(sum(pie_array1[i]))
            b = float(sum(pie_array2[i]))
            if b ==0: 
                temp_ratio = 0
            else:
                temp_ratio = a/b
            plt.bar(x[i],temp_ratio, width=1, color='cyan')
            ratios.append(temp_ratio)

        ax.yaxis.set_visible(True)
        ax.set_ylabel("Averages", weight = 'bold', color = 'black', labelpad=-5,fontsize = 12)
        ax.yaxis.set_ticklabels([])    
        plt.xticks(fontsize=10)
        plt.ylim(0,10)
        plt.xlim(0,len(Sort1.units))
        plt.plot([0,len(Sort1.units)],[1.0,1.0], 'r--', color = 'black', linewidth=5)

        #ax.set_ylabel("Ratio", weight = 'bold', color = 'black', fontsize = 12)
        #ax.set_xlabel("Single Unit IDs", weight = 'bold', color = 'black', fontsize = 12)

        #Plot distribution of ratios:
        ax = plt.subplot(gs[len(Sort2.units)+1:len(Sort2.units)+2,6:13]) 
        bin_width=.2
        ax.set_xlabel("Distribution of ratios", weight = 'bold', color = 'black', fontsize = 14)
        yy=np.histogram(ratios, bins = np.arange(0,100,bin_width))
        plt.plot(yy[1][0:-1],yy[0], linewidth=4, color='cyan')
        plt.plot([1,1],[0,100.0], 'r--', color = 'black', linewidth=5)
        plt.ylim(0,max(yy[0]))
        plt.xlim(0,5)

        #Plot distribution of locking times ss_lock_array
        bin_width = 5
        for i in range(len(Sort2.units)):
            ax = plt.subplot(gs[i:i+1,13:15]) 
            if i ==0: plt.title("Dist. of peak \n locks (all units)", weight = 'bold', color = 'black', fontsize = 14)
            yy=np.histogram(ss_lock_array.T[i], bins = np.arange(-window*1E3,+window*1E3,bin_width))
            plt.plot(yy[1][0:-1],yy[0], linewidth=4, color='blue')

            index=[np.argmax(yy[0])]
            temp_x = yy[1][index]
            plt.plot([temp_x,temp_x],[0,yy[0][index]], 'r--', color='red', linewidth=2)
            plt.plot([0,0],[0,yy[0][index]], 'r--', color='black', linewidth=2)

            plt.ylim(bottom=0.0)
            plt.xlim(-50,50)

            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(False)

        ax.xaxis.set_visible(True)
        max_yticks = 3
        yloc = plt.MaxNLocator(max_yticks)
        ax.yaxis.set_major_locator(yloc)
        plt.yticks(fontsize = 10)
            
        ax.set_xlabel("Time of Peak (ms).", weight = 'bold', color = 'black', fontsize = 12)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
        plt.show()
        
    quit()
    
    plt.suptitle(Sort1.name+ "   Length of recording: " + str(rec_length) + " seconds. "+
    "\n plotted spikes: " + str(len(pop_spikes)) + "   of    Total pop spikes: " + str(len(Sort2.units[0])), fontsize=20)

    ax1 = plt.subplot(1, 2, 1)
    
    x = []
    y = []
    for i in range(len(Sort1.units)):
        x.append(peak[i])
        y.append(-Sort1.chanpos[Sort1.maxchan[i]][1])

    plt.plot([0,0],[0,-3000], 'r--',color='black',linewidth=1)

    ax1.plot(x,y,linewidth=3, color='blue')
    ax1.set_ylim(-2000,0)
    ax1.set_xlim(-0.075,0.075)
    
    plt.title("Peak of phase lock (xaxis) vs. depth of Unit (yaxis)" , fontsize=14)
    ax1.set_xlabel("Seconds from pop spike centre (secs)", fontsize=14)
    ax1.set_ylabel("Depth of unit from top of cortex/probe (um)", fontsize=14)
    
    #np.savetxt(Sort1.directory + Sort1.name+'_depthlock_x.csv', x, delimiter=",")
    #np.savetxt(Sort1.directory + Sort1.name+'_depthlock_y.csv', y, delimiter=",")

    ax2 = plt.subplot(1, 2, 2)
    yy1=np.histogram(peak, bins = np.arange(-0.150,+0.150,0.025))
    ax2.plot(yy1[1][0:-1],yy1[0], linewidth=3, color='red')
    plt.ylabel("#"+str(i)+" : " + str(len(Sort1.units[i])), fontsize=14)
    
    plt.title("Distribution of phase locking times (histogram across all units)" , fontsize=14)
    ax2.set_xlabel("Seconds from pop spike centre (secs)", fontsize=14)
    ax2.set_ylabel("Number of units in bin (25ms width)", fontsize=14)
    
    ax2.set_xlim(-0.075,0.075)
    plt.show()

    #sorted_file = ['61-tr5c-blankscreen', '63-tr5c-blankscreen', '67-tr5c-blankscreen', '68-tr5c-blankscreen_66Hz']
    
    #xd=[]
    #yd=[]
    #sim_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
    #for i in range(len(sorted_file)):
        #temp_x = np.loadtxt(sim_dir + sorted_file[i]+'/'+sorted_file[i]+'_depthlock_x.csv', dtype=np.float32, delimiter=",")
        #temp_y = np.loadtxt(sim_dir + sorted_file[i]+'/'+sorted_file[i]+'_depthlock_y.csv', dtype=np.float32, delimiter=",")
        #xd.append(temp_x)
        #yd.append(temp_y)

    #xd = np.hstack(np.array(xd))
    #yd = np.hstack(np.array(yd))

    #s=xd
    #xd=yd
    #u=np.arange(-0.75,+0.75,.001)
    #fitdata=sinc_interp(xd,s,u)
    #plt.plot(u,fitdata, linewidth=5, color='black')
    #plt.plot(x,y,linewidth=3, color='blue')

    #plt.show()
    quit()


    #PLOT RASTERS
    ax = plt.subplot(1, 2, 1)
    title("POPULATION SPIKE TRIGGERED ACTIVITY (ALL DATA)",fontsize=15)

    #Plot all spikes in response window
    counter=0
    scale=1
    for i in range(Sort1.n_units):
        print "Plotting Unit: ", counter
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1
        #yy = np.histogram(save_response_array[i], bins = np.arange(-window,+window,.010))
        #plt.vlines(yy[1], ymin, ymax, linewidth=3, color='blue', alpha=0.8)

        plt.vlines(save_response_array[i], ymin, ymax, linewidth=3, color='blue', alpha=0.8)

    #*************************SYNCHRONIZED STATE ONLY **************************
    ax = plt.subplot(1, 2, 2)
    title("POPULATION SPIKE TRIGGERED ACTIVITY (SYNCHRONIZED STATE ONLY)",fontsize=15)

    #Make histograms -Xms .. +Xms from population spikes
    window=0.500        #Window to search
    min_spikes = 0      #min no. of spikes in window
    save_response_array=[]
    save_response_array_individual=[]
    for i in range(Sort1.n_units):
        save_response_array.append([])
        save_response_array_individual.append([])

    temp2 = np.array(Sort2.units[0],dtype=np.float32)*20./float(Sort2.samplerate)
    #pop_spikes = temp2
    pop_spikes = temp2[np.where(np.logical_and(temp2>=300., temp2<=750.))]

    #Loop over all units in sorted data
    for i in range(Sort1.n_units):
        temp_array = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)

        #Loop over all population spikes
        for j in range(len(pop_spikes)):
            #Find all unit spikes that fall w/in +/- window of pop spike
            temp2 = np.where(np.logical_and(temp_array>=pop_spikes[j]-window, temp_array<=pop_spikes[j]+window))[0]
            if len(temp2)>min_spikes: #Min spikes in window 
                save_response_array[i].extend(temp_array[temp2]-pop_spikes[j]) #Offset time of spike to t_0 
                save_response_array_individual[i].append(temp_array[temp2]-pop_spikes[j])

    #Plot all spikes in response window
    counter=0
    scale=1
    for i in range(Sort1.n_units):
        print "Plotting Unit: ", counter
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1

        #yy = np.histogram(save_response_array[i], bins = np.arange(-window,+window,.010))
        #plt.vlines(yy[1], ymin, ymax, linewidth=3, color='blue', alpha=0.8)

        plt.vlines(save_response_array[i], ymin, ymax, linewidth=3, color='red', alpha=0.8)

    plt.show()


    fig, ax = plt.subplots()
    heatmap = ax.pcolor(np.array(save_response_array))
    plt.show()

    quit()
#*********************************************************************************************
    #Luczak sequence order plots

    ax = plt.subplot(1, 5, 1)

    bin_width=0.025

    save_array = np.sort(np.array(array))/Sort1.samplerate
    yy = np.histogram(save_array, bins = np.arange(0,max(save_array),bin_width))

    state=0
    state_start=[]
    state_end=[]
    for i in range(len(yy[1])-1):
        if (state==0 and ((yy[0][i]/bin_width)>ave_rate)):
            state=1
            state_start.append(yy[1][i])

        if (state==1 and ((yy[0][i]/bin_width)<ave_rate)):
            state=0
            state_end.append(yy[1][i])

    state_duration = np.array(state_end) - np.array(state_start[0:len(state_end)]) 

    title("Cumulative responses from -50ms to +50ms of up-phase start - # up-phases: "+str(len(state_start))
    +'\n'+Sort1.name, multialignment='center',fontsize=15)

    ###Plot state_durations
    #for i in range(500):
        #plt.fill_between([state_start[i],state_end[i]], 0, 500, color='black', alpha=0.25)
 
    #plt.plot([0,max_x],[ave_rate, ave_rate], linewidth=2, color='black', alpha=0.750)

    #plt.xlim(0, max(save_array))
    #plt.ylim(0, 1)

    #Plot firing rate histogram for each cell in data
    #Do it for 50ms before and after burst


    #Save arrays for data cut in 4 different ways
    save_response_array=[]
    save_response_array2=[]
    save_response_array_individual=[]
    save_response_array_individual1=[]
    save_response_array_individual2=[]
    save_response_array_individual3=[]
    save_response_array_individual4=[]
    for i in range(Sort1.n_units):
        save_response_array.append([])
        save_response_array2.append([])
        save_response_array_individual.append([])
        save_response_array_individual1.append([])
        save_response_array_individual2.append([])
        save_response_array_individual3.append([])
        save_response_array_individual4.append([])
    window=0.200


    ##Loop over all cells
    ##OLD WAY: Mean of means: Grand means
    #for i in range(Sort1.n_units):
        #temp_array = np.array(Sort1.units[i],dtype=np.float32)/Sort1.samplerate
        ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
        #for j in range(len(state_start)):
            #temp2 = np.where(np.logical_and(temp_array>=state_start[j], temp_array<=state_start[j]+window))[0]
            #if len(temp2)>0: #Min spikes in window after upphase start
                #save_response_array[i].extend(temp_array[temp2]-state_start[j]) #Offset time of spike to t_0 
                #save_response_array_individual[i].append(temp_array[temp2]-state_start[j])

    #Loop over all cells
    #New way: Sum spikes across all up-phase epochs

    #CHANGE THIS
    for i in range(Sort1.n_units):
        temp_array = np.array(Sort1.units[i],dtype=np.float32)/Sort1.samplerate
        #Loop over all bright flashes and grab response 3 seconds out from screen flash start
        for j in range(len(state_start)):
            temp2 = np.where(np.logical_and(temp_array>=state_start[j], temp_array<=state_start[j]+window))[0]
            if len(temp2)>0: #Min spikes in window after upphase start
                save_response_array[i].extend(temp_array[temp2]-state_start[j]) #Offset time of spike to t_0 
                save_response_array_individual[i].append(temp_array[temp2]-state_start[j])


    #Split response arrays into 2: alternating between even and odd phases counts
    #Split: May also wish to try splitting in halves
    temp=save_response_array_individual
    for i in range(Sort1.n_units):
        for j in range(len(temp[i])/2-1):
            if len(temp[i][j])>0:
                save_response_array_individual1[i].append(temp[i][2*j])  #Sequential split
                save_response_array_individual2[i].append(temp[i][2*j+1])
                save_response_array_individual3[i].append(temp[i][j])  #Half split
                save_response_array_individual4[i].append(temp[i][j+len(temp[i])/2])
    #    print "1: ", len(save_response_array_individual1[i])
    #    print "2: ", len(save_response_array_individual2[i])

    time.sleep(2)
    scale=1
    scale1=.1
    #mean_val_sort=[]
    #for i in range(Sort1.n_units):

        #x_list = [item for sublist in  save_response_array_individual1[i] for item in sublist]
        #x = np.array(x_list, dtype=float32)
        #x = np.sort(x)

        ##Plot mean of summed activity over all states
        #x = x[np.where(np.logical_and(x>=0.0, x<=0.1))[0]]
        #mean_val = np.mean(np.array(x, dtype=float64), axis=0)
        #print "Unit: ", i, "sum mean_val: ", mean_val
        
        #mean_val_sort.append(mean_val)

        #ymin=0
        #ymax=0
        #ymin+=i*scale
        #ymax+=i*scale+.7
   
    ##x = np.sort(mean_val_sort)
    #tracker = np.arange(0,len(mean_val_sort),1)
    #tracker = [y for (x,y) in sorted(zip(mean_val_sort,tracker))]

    #SORT DATA IN ORDER OF LATENCIES OF GRAND MEANS:
    mean_val_sort=[]
    counter=0
    for i in range(Sort1.n_units):  #plot grand meaned values
    #for i in tracker:  #plot sum meaned values;

        ymin=0
        ymax=0
        ymin+=counter*scale
        ymax+=counter*scale+.7
        counter+=1
        #plt.vlines(mean_val_sort[i], ymin, ymax, linewidth=3, color='red', alpha=0.8)

        #Plot mean of summed activity for each indivudal state transition
        x = save_response_array_individual1[i]
        mean_val_array=[]
        for j in range(len(x)): 
            x_individual = np.array(x[j],dtype=np.float64)
            #x_individual= x_individual[np.where(np.logical_and(x_individual>=0.0, x_individual<=0.1))[0]]
            #if len(x_individual)>0: #Check to see if any spikes occured in that upphase 

            mean_val = np.mean(np.array(x_individual, dtype=float64), axis=0)
            mean_val_array.append(mean_val)

            #plt.vlines(mean_val, ymin, ymax, linewidth=2, color='blue', alpha=0.1)
        
        if mean_val_array: 
            mean_val_array_value = np.mean(np.array(mean_val_array, dtype=float64), axis=0)
            #plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)
            mean_val_sort.append(mean_val_array_value)

            print "Mean of means: ", mean_val_array_value

    #x = np.sort(mean_val_sort)
    tracker = np.arange(0,len(mean_val_sort),1)
    tracker = [y for (x,y) in sorted(zip(mean_val_sort,tracker))]

    #PLOT DATA IN ORDER OF LATENCIES - SORTED BY GRAND MEAN:
    counter=0
    first_half_array=[]
    for i in range(len(tracker)):
        first_half_array.append([])

    for i in tracker:
        print "Plotting Unit: ", counter, " / ", len(tracker)
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1

        mean_val_array_value = mean_val_sort[i] #np.mean(np.array(mean_val_array[i], dtype=float64), axis=0)
        first_half_array[i]=mean_val_array_value
        plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)

        x = save_response_array_individual1[i]
        mean_val_array=[]
        for j in range(len(x)): 
            x_individual = np.array(x[j],dtype=np.float64)
            mean_val = np.mean(np.array(x_individual, dtype=float64), axis=0)
            plt.vlines(mean_val, ymin, ymax, linewidth=2, color='blue', alpha=0.1)

    x=[0,0]
    y=[0,Sort1.n_units]
    plt.plot(x,y,'r--', linewidth=1.5, color='red')

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Unit ID',multialignment='center', fontsize=15)

    plt.xlim(0, 0.200)
    plt.ylim(0, Sort1.n_units)

    #************************ PLOT 2nd HALF OF DATA ********************
    #Up state triggers
    ax = plt.subplot(1, 5, 2)
    print "Plotting 2nd HALF DATA"
    scale=1
    scale1=.1
    mean_val_sort2=[]
    ##Plot by summed mean
    #for i in tracker:

        #x_list = [item for sublist in  save_response_array_individual2[i] for item in sublist]
        #x = np.array(x_list, dtype=float32)
        #x = np.sort(x)

        ##Plot mean of summed activity over all states
        #x = x[np.where(np.logical_and(x>=0.0, x<=0.1))[0]]
        #mean_val = np.mean(np.array(x, dtype=float64), axis=0)
        #print "Unit: ", i, "mean_val: ", mean_val
        
        #mean_val_sort2.append(mean_val)

        #ymin=0
        #ymax=0
        #ymin+=i*scale
        #ymax+=i*scale+.7
   
    #DATA ALREADY SORTED ABOVE - USE SAME ORDER
    ##Sort by order of the 1st half of data
    #tracker = np.arange(0,len(mean_val_sort),1)
    #tracker = [y for (x,y) in sorted(zip(mean_val_sort,tracker))]

    #PLOT DATA IN ORDER OF LATENCIES:
    counter=0
    second_half_array=[]
    for i in range(len(tracker)):
        second_half_array.append([])
    for i in tracker:
        print "Plotting Unit: ", i, " / ", len(tracker)
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1
        #Plot mean of summed activity for each indivudal state transition
        x = save_response_array_individual2[i]
        mean_val_array=[]
        for j in range(len(x)):
            x_individual = np.array(x[j],dtype=np.float64)
            #x_individual= x_individual[np.where(np.logical_and(x_individual>=0.0, x_individual<=0.1))[0]]
            mean_val = np.mean(np.array(x_individual, dtype=float64), axis=0)
            plt.vlines(mean_val, ymin, ymax, linewidth=2, color='blue', alpha=0.1)
            
            mean_val_array.append(mean_val)

        if mean_val_array: 
            mean_val_array_value = np.mean(np.array(mean_val_array, dtype=float64), axis=0)
            plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)
            second_half_array[i]=mean_val_array_value
        else:
            second_half_array[i]=0

            #print "Mean of means: ", mean_val_array

        #mean_val_array_value = mean_val_sort2[i] #np.mean(np.array(mean_val_array[i], dtype=float64), axis=0)
        #plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)


    x=[0,0]
    y=[0,Sort1.n_units]
    plt.plot(x,y,'r--', linewidth=1.5, color='red')

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Unit ID',multialignment='center', fontsize=15)

    plt.xlim(0, 0.200)
    plt.ylim(0, Sort1.n_units)

    #************************ PLOT 2nd HALF OF DATA ********************
    #Up state triggers
    ax = plt.subplot(1, 5, 3)
    print "Plotting 3rd HALF DATA"
    scale=1
    scale1=.1
    mean_val_sort2=[]
    ##Plot by summed mean
    #for i in tracker:

        #x_list = [item for sublist in  save_response_array_individual2[i] for item in sublist]
        #x = np.array(x_list, dtype=float32)
        #x = np.sort(x)

        ##Plot mean of summed activity over all states
        #x = x[np.where(np.logical_and(x>=0.0, x<=0.1))[0]]
        #mean_val = np.mean(np.array(x, dtype=float64), axis=0)
        #print "Unit: ", i, "mean_val: ", mean_val
        
        #mean_val_sort2.append(mean_val)

        #ymin=0
        #ymax=0
        #ymin+=i*scale
        #ymax+=i*scale+.7
   
    #DATA ALREADY SORTED ABOVE - USE SAME ORDER
    ##Sort by order of the 1st half of data
    #tracker = np.arange(0,len(mean_val_sort),1)
    #tracker = [y for (x,y) in sorted(zip(mean_val_sort,tracker))]

    #PLOT DATA IN ORDER OF LATENCIES:
    counter=0
    third_half_array=[]
    for i in range(len(tracker)):
        third_half_array.append([])
    for i in tracker:
        print "Plotting Unit: ", i, " / ", len(tracker)
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1
        #Plot mean of summed activity for each indivudal state transition
        x = save_response_array_individual3[i]
        mean_val_array=[]
        for j in range(len(x)):
            x_individual = np.array(x[j],dtype=np.float64)
            #x_individual= x_individual[np.where(np.logical_and(x_individual>=0.0, x_individual<=0.1))[0]]
            mean_val = np.mean(np.array(x_individual, dtype=float64), axis=0)
            plt.vlines(mean_val, ymin, ymax, linewidth=2, color='blue', alpha=0.1)
            
            mean_val_array.append(mean_val)

        if mean_val_array: 
            mean_val_array_value = np.mean(np.array(mean_val_array, dtype=float64), axis=0)
            plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)
            third_half_array[i]=mean_val_array_value
        else:
            third_half_array[i]=0

            #print "Mean of means: ", mean_val_array

        #mean_val_array_value = mean_val_sort2[i] #np.mean(np.array(mean_val_array[i], dtype=float64), axis=0)
        #plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)


    x=[0,0]
    y=[0,Sort1.n_units]
    plt.plot(x,y,'r--', linewidth=1.5, color='red')

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Unit ID',multialignment='center', fontsize=15)

    plt.xlim(0, 0.200)
    plt.ylim(0, Sort1.n_units)

    #************************ PLOT 2nd HALF OF DATA ********************
    #Up state triggers
    ax = plt.subplot(1, 5, 4)
    print "Plotting 4th HALF DATA"
    scale=1
    scale1=.1
    mean_val_sort2=[]
    ##Plot by summed mean
    #for i in tracker:

        #x_list = [item for sublist in  save_response_array_individual2[i] for item in sublist]
        #x = np.array(x_list, dtype=float32)
        #x = np.sort(x)

        ##Plot mean of summed activity over all states
        #x = x[np.where(np.logical_and(x>=0.0, x<=0.1))[0]]
        #mean_val = np.mean(np.array(x, dtype=float64), axis=0)
        #print "Unit: ", i, "mean_val: ", mean_val
        
        #mean_val_sort2.append(mean_val)

        #ymin=0
        #ymax=0
        #ymin+=i*scale
        #ymax+=i*scale+.7
   
    #DATA ALREADY SORTED ABOVE - USE SAME ORDER
    ##Sort by order of the 1st half of data
    #tracker = np.arange(0,len(mean_val_sort),1)
    #tracker = [y for (x,y) in sorted(zip(mean_val_sort,tracker))]

    #PLOT DATA IN ORDER OF LATENCIES:
    counter=0
    fourth_half_array=[]
    for i in range(len(tracker)):
        fourth_half_array.append([])
    for i in tracker:
        print "Plotting Unit: ", i, " / ", len(tracker)
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1
        #Plot mean of summed activity for each indivudal state transition
        x = save_response_array_individual4[i]
        mean_val_array=[]
        for j in range(len(x)):
            x_individual = np.array(x[j],dtype=np.float64)
            #x_individual= x_individual[np.where(np.logical_and(x_individual>=0.0, x_individual<=0.1))[0]]
            mean_val = np.mean(np.array(x_individual, dtype=float64), axis=0)
            plt.vlines(mean_val, ymin, ymax, linewidth=2, color='blue', alpha=0.1)
            
            mean_val_array.append(mean_val)

        if mean_val_array: 
            mean_val_array_value = np.mean(np.array(mean_val_array, dtype=float64), axis=0)
            plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)
            fourth_half_array[i]=mean_val_array_value
        else:
            fourth_half_array[i]=0

            #print "Mean of means: ", mean_val_array

        #mean_val_array_value = mean_val_sort2[i] #np.mean(np.array(mean_val_array[i], dtype=float64), axis=0)
        #plt.vlines(mean_val_array_value, ymin, ymax, linewidth=3, color='green', alpha=0.8)


    x=[0,0]
    y=[0,Sort1.n_units]
    plt.plot(x,y,'r--', linewidth=1.5, color='red')

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Unit ID',multialignment='center', fontsize=15)

    plt.xlim(0, 0.200)
    plt.ylim(0, Sort1.n_units)

    #************************ PLOT DISTANCES BETWEEN 1st and 2nd HALVES ********************
    #Up state triggers
    ax = plt.subplot(1, 5, 5)
    print "Plotting 2nd HALF DATA"
    counter=0

    print first_half_array
    print second_half_array

    for i in tracker:
        print "Plotting Unit: ", i, " / ", len(tracker)
        ymin=counter*scale
        ymax=counter*scale+.7
        counter+=1

        #Plot mean of summed activity for each indivudal state transition

        plt.vlines(first_half_array[i], ymin, ymax, linewidth=2, color='red', alpha=0.9)
        plt.vlines(second_half_array[i], ymin, ymax, linewidth=2, color='blue', alpha=0.9)
        plt.vlines(third_half_array[i], ymin, ymax, linewidth=2, color='green', alpha=0.9)
        plt.vlines(fourth_half_array[i], ymin, ymax, linewidth=2, color='magenta', alpha=0.9)

        if abs(max(first_half_array[i],second_half_array[i],third_half_array[i],fourth_half_array[i])
        - min(first_half_array[i],second_half_array[i],third_half_array[i],fourth_half_array[i]))< 0.005:
            plt.axhspan(ymin, ymax, facecolor='0.5', alpha=0.5)

    plt.xlim(0, 0.200)
    plt.ylim(0, Sort1.n_units)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

    quit()

    #******************************PLOT DISTRIBUTION OF UPPHASE DURATION
    ax = plt.subplot(1, 1, 1)
    factor = 1e-1

    save_array = np.sort(state_duration)
    yy = np.histogram(save_array, bins = np.arange(0,max(save_array),bin_width))
    plt.plot(np.log(yy[1][1:-1]),np.log(yy[0][1:]+factor), linewidth=1, color='blue', alpha=0.750)

    #x = np.log(np.array(yy[1][:-1]).clip(1e-10),dtype=np.float32)
    #y = np.log(np.array(yy[0]).clip(1e-10),dtype=np.float32)
    x = np.log(yy[1][1:-1])
    y = np.log(yy[0][1:]+factor)

    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A,y)[0]
    #c+= -log(1e-10) #+2
    print m, c
    #x = yy[1][:-1]
    plt.plot(x, m*x + c, 'r', linewidth = 2, color = 'red', label='Fitted line')

    ax.set_yscale('linear')
    ax.set_xscale('linear')
    plt.xlim(left=x[1], right=4)
    plt.ylim(bottom=0, top=max(y))

    title("Distribution of length of Up Phases, slope = "+ str(m), multialignment='center', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def Plot_trackwide_LFP_triggered_activity(Sorts, LFP_sort, tsf):

    #************************************* LFP BASED PLOTS **************************************
    #Sorts contains a list of single unit recording objects
    #LFP_sort contains the track-wide LFP sort

    print "TRACKWIDE MSL PLOTS"

    sectioning = True
    sec_length = 30.*60.0  #length of partitioning of trackwide data in minutes converted to seconds

    window=1                        #Total window to compute MSL data
    small_window = 50              #Zoomed in window for MSL plots
    large_window = small_window*10  #Zoomed out windows for MLS plots
    start_lfp=0                    #The first LFP unit to start loop; default = 0
    end_lfp = len(LFP_sort.units)   #The last LFP unit to end loop; default = length of Sort1.units
    cutoff_freq = 0.001               #Minimum firing rate of cell to consider in analysis

    #Find total number of units in track wide sort:
    unit_ids = []
    for rec in range(len(Sorts)):
        for unit in range(len(Sorts[rec].units)):
            unit_ids.append([Sorts[rec].uid[unit]])

    unit_ids = np.unique(np.hstack(unit_ids))

    #Load this for display purposes
    LFP_sort.rec_start=[]
    LFP_sort.rec_end=[]
    LFP_sort.rec_start.append(0.0)
    LFP_sort.rec_end.append(float(Sorts[0].n_vd_samples)/float(Sorts[0].SampleFrequency))
    for i in range(len(Sorts)):
        LFP_sort.rec_start.append(LFP_sort.rec_end[i-1])
        LFP_sort.rec_end.append(LFP_sort.rec_end[i-1]+float(Sorts[i].n_vd_samples)/float(Sorts[i].SampleFrequency))

    #**** PARTITION DATA BY TIME CHUNKS ******
    if sectioning:
        #Partition analysis according to sec_length above
        #Find beginning of each recording for LFP_sort data from single unit Sorts indvidual lengths
        
        SUA_object = type('SUA_object', (object,), {})
        SUA_sort = SUA_object()

        SUA_sort.units = []
        SUA_sort.unit_ids = []         
        for i in range(max(unit_ids)+1):
            SUA_sort.units.append([])
            SUA_sort.unit_ids.append(i)

        recording_lengths = []
        #recording_lengths.append(0.0)
        for i in range(len(Sorts)):
            #for j in range(len(Sorts[i].units)): 
            for j in range(len(Sorts[i].units)): 
                #Add spike rasters only from unit spiking during particular recording
                SUA_sort.units[Sorts[i].uid[j]].extend(sum(recording_lengths)+ np.float32(Sorts[i].units[j])/float(Sorts[i].SampleFrequency))
            
            #Keep track of absolute .ptcs/.srf start time 
            recording_lengths.append(float(Sorts[i].n_vd_samples)/float(Sorts[i].SampleFrequency))
        
        #Save original recording lengths 
        SUA_sort.rec_len = recording_lengths

        #Calculate time sections starts/ends 
        temp_list = np.arange(0, sum(SUA_sort.rec_len), sec_length)
        temp_list = np.append(temp_list,sum(SUA_sort.rec_len)) #Add end of total recording also
        SUA_sort.sec_len=temp_list
        
    lfp_events=True
    plotting=True
    
    if lfp_events:
        
        #For every recording - must track MSL data for each unit
        peaks_depth=np.zeros((len(SUA_sort.sec_len),len(unit_ids),len(LFP_sort.units),small_window*2), dtype=np.float32) #Set all peaks to very large number
        peaks_depth_large=np.zeros((len(SUA_sort.sec_len),len(unit_ids),len(LFP_sort.units),small_window*10), dtype=np.float32) #Set all peaks to very large number

        peaks_order=np.zeros((len(SUA_sort.sec_len)-1,len(unit_ids),len(LFP_sort.units)), dtype=np.float32)+100000. #Set all peaks to very large number
        peaks_pval=np.zeros((len(SUA_sort.sec_len)-1,len(unit_ids),len(LFP_sort.units)), dtype=np.float32)+100000. #Set all peaks to very large number

        #Outer loop is over LFP units 
        for ss in range(start_lfp, end_lfp, 1):

            #Load pop_spikes raster
            pop_spikes_all = np.array(LFP_sort.units[ss],dtype=np.float32)*100/float(LFP_sort.samplerate) #Convert to seconds
            pop_spikes_timesteps = pop_spikes_all*1E3 #Convert back to 1Khz samplerate 

            if plotting:
                height = 25
                width_plots = int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.25)

                gs = gridspec.GridSpec(height, width_plots)

                #*************** SPECGRAM ************
                #Plot specgram
                ax_image = plt.subplot(gs[:5,int(width_plots*.05)+8:]) 
                #img= np.rot90(tsf.specgram[::-1])

                img= tsf.specgram[::-1]
                rec_length = tsf.n_vd_samples/tsf.SampleFrequency/60./60.
                #ax_image.set_xlim(0,1000)
                #ax_image.set_ylim(0,rec_length)
                im = ax_image.imshow(img, extent=[0,rec_length,0,110], aspect='normal') #, extent=tsf.extent, cmap=tsf.cm)
                #im = ax_image.imshow(img, aspect='normal') #extent=lfp.extent, cmap=lfp.cm)
                for i in range(len(Sorts)):
                    #plt.vlines(LFP_sort.rec_start[i]/60./60., 0, 110, linewidth=2, color='black', alpha=0.9)
                    plt.vlines(sum(SUA_sort.rec_len[0:i])/60./60., 0, 110, linewidth=2, color='black', alpha=0.9)

                ax_image.set_ylabel("Frequency (Hz)", weight='bold')
                ax_image.xaxis.set_visible(False)
                plt.yticks(fontsize=12)
            
                #*************** SAT RASTERS ***********
                #Plot track-wide pop rasters
                ax_image = plt.subplot(gs[5:7, int(width_plots*.05)+8:])

                #y = np.zeros(len(pop_spikes), dtype=np.float32) + 1.0
                ptp_lfpspike = []
                n_electrodes = 10
                width = 200
                for spike in pop_spikes_timesteps:
                    ptp_temp=[]
                    for i in range(n_electrodes):
                        ptp_temp.append(max(tsf.ec_traces[i][spike-width/2:spike+width/2])-min(tsf.ec_traces[i][spike-width/2:spike+width/2]))
                    ptp_lfpspike.append(max(ptp_temp))
                
                #for i in range(len(Sorts)):
                #    plt.vlines(sum(SUA_sort.rec_len[0:i])/60./60., 0, 5000, linewidth=2, color='black', alpha=0.9)

                if True:
                    for i in range(len(SUA_sort.sec_len)):
                        plt.vlines(sum(SUA_sort.sec_len[i])/60./60., 0, 5000, linewidth=2, color='red', alpha=0.7)
                    
                
                plt.scatter(pop_spikes_timesteps*1E-3/60./60., ptp_lfpspike, s=3, color='blue')
                plt.xlim(0,rec_length)
                ax_image.yaxis.set_ticks(np.arange(0,2000,500))
                ax_image.set_ylabel("PTP (uV)", weight='bold')
                ax_image.set_xlabel("Time (Hrs)", weight='bold')

                plt.ylim(0,2000)
                plt.yticks(fontsize=12)

                #*************** LFP TRACES **************
                #Plot LFP spikes
                ax_image = plt.subplot(gs[:7,0:int(width_plots*.05):])

                width = 200
                x = np.arange(0,width,1)
                for i in range(n_electrodes):
                    trace_sum=np.zeros(width, dtype=np.float32)
                    for spike in pop_spikes_timesteps[0:min(150,len(pop_spikes_timesteps))]:
                        trace_temp = tsf.ec_traces[i][spike-width/2:spike+width/2]*2
                        plt.plot(x, trace_temp - (i*2000), color='black', alpha=.5)
                        trace_sum+=trace_temp
                    plt.plot(x, trace_sum/min(150,len(pop_spikes_timesteps)) - (i*2000), color='red', linewidth=2, alpha=.8)

                plt.title("LFP Cluster: " + str(ss), fontsize=15)
                plt.ylim(-20000, 2000)
                ax_image.yaxis.set_visible(False)
                ax_image.xaxis.set_ticks(np.arange(0,201,100))
                ax_image.set_xlabel("Time (ms)", weight='bold')

            #*********** LOOP OVER SECTIONS OF RECORDING **************
            #Inner loop is over each recording in track-wide sorts
            #for rec in range(len(Sorts)): #Original recording lengths
            for sec in range(len(SUA_sort.sec_len)-1):

                fit_sum = np.zeros((len(LFP_sort.units),1000,1000), dtype=np.float32)
                p_val_array=np.zeros((len(LFP_sort.units),1000), dtype=np.float64)+1.

                #************** LOAD LFP AND SUA SPIKES DURING TIME CHUNK ************
                #Cat data  use *100 conversion factor
                if ('nick' in LFP_sort.directory):
                    temp = np.where(np.logical_and(pop_spikes_all>= SUA_sort.sec_len[sec], pop_spikes_all<=SUA_sort.sec_len[sec+1]))[0]

                    #pop_spikes contains lfp spikes that occur only during particular recording; NB: OFFSET TO 0 TO MATCH SINGLE UNIT RASTERS
                    #pop_spikes = pop_spikes_all[temp] - LFP_sort.rec_start[rec]  #USE THIS for .prf spikes which reset to 0.0 every epoch
                    pop_spikes = pop_spikes_all[temp]

                    if len(pop_spikes)==0: 
                        continue

                print "Section: ", sec, " pop spike: ", ss+1, " / ", len(LFP_sort.units), "  #spikes: ", len(pop_spikes), ' / ', len(LFP_sort.units[ss])                        

                ##Mouse data use *20 conversion factor - NO! - CONVERT ALL RAW SPIKE STREAMS BACK TO ABSOLUTE TIME
                #if 'tim' in LFP_sort.directory:
                #    pop_spikes = np.array(Sort2.units[ss],dtype=np.float32)/float(Sort2.samplerate) #Convert to seconds

                #Loop over particular pop spike times in particular recording 
                for unit in unit_ids:
                    temp = np.where(np.logical_and(SUA_sort.units[unit]>= SUA_sort.sec_len[sec], SUA_sort.units[unit]<=SUA_sort.sec_len[sec+1]))[0]
                    spike_array = np.array(SUA_sort.units[unit])[temp]

                    xx1=[]
                    xx_even=[]
                    xx_odd=[]
                    xx1_ttest=[]
                    for j in range(len(pop_spikes)):
                        #find all single unit spikes that fall w/in +/- window of pop spike
                        temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                        #NB: Not excluding duplicate spikes from counts!!! make sure to take into account for analysis
                        x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                        #print "locking spikes: ", x
                        xx1.append(x)
                        if j % 2 == 0:
                            xx_even.append(x)
                        else:
                            xx_odd.append(x)
                    
                    #Plot PSTH on cumulative data from above: 3 methods available
                    #1 Plot gaussian convolved PSTHs - Martin's suggestion
                    if len(xx_even)>0: xx_even = np.hstack(xx_even)
                    if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)
                    firing_rate = float(len(xx_even)+len(xx_odd))/float(SUA_sort.sec_len[sec+1] - SUA_sort.sec_len[sec])

                    if firing_rate > cutoff_freq:
                        #print "Unit ", unit, " firing rate: ", firing_rate
                        
                        sig = 20 #20 ms (scaled to fit window) - Martin's suggestion
                        fit_even = np.zeros(1000, dtype=np.float32)
                        if len(xx_even)>0 : 
                            for g in range(len(xx_even)):
                                mu = np.array(xx_even[g])  #Default of normalized mean
                                fit_even += gaussian(np.linspace(-window*1E3, window*1E3, 1000), mu, sig)
                            #u = np.linspace(-window*1E3, window*1E3, 1000)
                            #if plotting: ax_image.plot(u,fit_even*10, linewidth=2, color='blue', alpha=0.7)

                        fit_odd = np.zeros(1000, dtype=np.float32)
                        if len(xx_odd)>0 : 
                            for g in range(len(xx_odd)):
                                mu = xx_odd[g] #Default of normalized height
                                fit_odd += gaussian(np.linspace(-window*1E3,window*1E3, 1000), mu, sig)
                            #u = np.linspace(-window*1E3, window*1E3, 1000)
                            #if plotting: ax_image.plot(u,fit_odd*10, linewidth=2, color='red', alpha=0.7)

                        #fit_sum[ss][Sorts[rec].uid[unit]] = (fit_even + fit_odd)/2.
                        fit_sum[ss][SUA_sort.unit_ids[unit]] = (fit_even + fit_odd)/2.

                        if False:
                            plt.show()
                            
                            ax_image=plt.subplot(1,1,1)
                            
                            x = np.linspace(-window*1E3,window*1E3, 1000)
                            even_max = np.argmax(fit_even)
                            odd_max = np.argmax(fit_odd)
                            
                            plt.plot(x,fit_even, color='blue',linewidth=3, alpha=0.8)
                            plt.plot(x,fit_odd, color='red',linewidth=3, alpha=0.8)
                            
                            even_maxy = [0,fit_even[even_max]]
                            odd_maxy = [0,fit_odd[odd_max]]

                            even_max = [x[even_max], x[even_max]]
                            odd_max = [x[odd_max], x[odd_max]]
                            
                            print even_max, even_maxy
                            print odd_max, odd_maxy
                            
                            
                            
                            plt.plot(even_max,even_maxy,'r--', color = 'blue', linewidth=2, alpha=.9)
                            plt.plot(odd_max,odd_maxy,'r--', color = 'red', linewidth=2, alpha=.9)
                            
                            ax_image.tick_params(axis='both', which='major', labelsize=20)
                            
                            plt.xlabel("Time from LFP Event Peak (ms)", weight='bold', fontsize=30)
                            plt.xlim(-window*1E3/10, window*1E3/10)
                            ax_image.yaxis.set_visible(False)
                            plt.show()
                            #quit()
                        
                    else:
                        fit_sum[ss][SUA_sort.unit_ids[unit]] = np.zeros(1000, dtype=np.float32)
                    
                    #Compute T-test values; #Limit analysis to at 0.05Hz firing cells / try different firing rates and then show % locking
                    if False:
                        xx1 = np.hstack(xx1)
                        p_val=1.0
                        if (len(Sorts[rec].units[unit])/float(Sorts[rec].n_vd_samples/Sorts[rec].SampleFrequency))> 0.01: 
                            temp_float=0.
                            for r in range(5):
                                xx1_ttest =  (np.random.random_sample((1000,))-0.5)*2*window*1E3
                                KS, p_val = sstats.ks_2samp(np.sort(xx1), np.sort(xx1_ttest)) #2 sample ks test
                                temp_float+=p_val
                            p_val=temp_float/5.
                            p_val_array[ss][unit]=p_val

                    #if plotting: 
                        ##plt.text(0, rec_length, "P: %.1e" % p_val, fontsize = 12)
                        ##plt.text(-window*1E3+window*1E3/20, rec_length, "U: "+str(unit), fontsize = 12, fontweight='bold') 

                        #ax_image.set_xlim(-window*1E3, window*1E3)
                        ##ax_image.set_ylim(0,rec_length)
                        #ax_image.xaxis.set_visible(False)
                        #ax_image.yaxis.set_visible(False)

                ##******PLOT LFP Spikes******* 
                #if plotting: 
                    #ax1 = plt.subplot(gs[0: unit/int(sqrt(len(Sorts[rec].units)))-1, 0]) 

                #n_points = 200
                #x = np.zeros((LFP_sort.n_electrodes,n_points),dtype=np.float32)
                #for i in range(LFP_sort.n_electrodes):
                    #x[i]= LFP_sort.Siteloc[i*2] + np.array(arange(0,n_points,1))

                #y1 = np.zeros((LFP_sort.n_electrodes,n_points), dtype=np.float32)

                #spike_list = np.arange(0,min(len(LFP_sort.units[ss]),50),1) #Need this if some spikes too close to end of recording
                
                #if plotting: 
                    #for j in spike_list: #range(min(len(Sort1.units[unit]),50)): #Plot up to 100 spikes from unit
                        #for i in range(Sort2.tsf.n_electrodes):
                            #if Sort2.maxchan[ss]==i:
                                #plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                                #Sort2.units[ss][j]+n_points/2]/600.-Sort2.tsf.Siteloc[i*2+1], color='blue', alpha=1.)
                            #else:
                                #plt.plot(x[i], Sort2.tsf.ec_traces[i][Sort2.units[ss][j]-n_points/2:
                                #Sort2.units[ss][j]+n_points/2]/600.-Sort2.tsf.Siteloc[i*2+1], color='black', alpha=1.)

                #***********************************************************************************************
                #Make MLS plots - DEPTH ORDER
                #100ms wide windows
                #Plot Luczak plots - only locking units - Original order - DEPTH
                if plotting:
                    ax1 = plt.subplot(  gs[11:14,int(width_plots*.05)+8+sec*3:int(width_plots*.05)+11+sec*3])
                                        #gs[:5,int(width_plots*.05)+8:]
                    
                    ax1.set_xlim(0,small_window*2)

                img=[]
                for unit in range(len(unit_ids)):
                    #if unit in Sorts[rec].uid:
                    if unit in SUA_sort.unit_ids:
                        if max(fit_sum[ss][unit])!= 0.0:
                            img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))
                            #peaks_depth[sec][unit][ss]=(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                            peaks_depth[sec][unit][ss]= fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit])
                        else:
                            temp = np.zeros(small_window*2,dtype=np.float32)
                            peaks_depth[sec][unit][ss]= temp
                            img.append(temp)
                    else:
                        temp = np.zeros(small_window*2,dtype=np.float32)
                        img.append(temp)
                        peaks_depth[sec][unit][ss]= temp

                if plotting and len(img)>0:
                    im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                    ax1.xaxis.labelpad = 0
                    if sec>0: ax1.yaxis.set_visible(False)
                    ax1.xaxis.set_visible(False)
                    plt.title(str(sec), weight = 'bold', fontsize=10)
                    ax1.set_ylim(0,len(img))
                    #ax1.set_xlim(10,70)

                #***********************************************************************************************
                #Plot Luczak plots - FIRING Sequence Order
                if True:
                    if plotting:
                        ax1 = plt.subplot(gs[14:17,int(width_plots*.05)+sec*3+8:int(width_plots*.05)+11+sec*3])
                        ax1.set_xlim(0,small_window*2)
                        ax1.set_ylim(0,len(img))

                    img=[]
                    lock_time=[]
                    for unit in range(len(unit_ids)):
                        if unit in SUA_sort.unit_ids:
                            if max(fit_sum[ss][unit])!= 0.0:
                                lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                                img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))
                            else:
                                lock_time.append(1000)
                                temp = np.zeros(small_window*2,dtype=np.float32)
                                img.append(temp)
                        else: 
                            lock_time.append(1000)
                            temp = np.zeros(small_window*2,dtype=np.float32)
                            img.append(temp)

                    if True:
                    #if sec<1:
                        inds = np.array(lock_time).argsort()    #First section use just existing order
                    
                    img=np.array(img)[inds]
                    
                    #if sec>0:
                    #    img=np.array(img)[inds]
                    #    inds = np.array(lock_time).argsort()    #Other sections use previous order

                    #for m in range(len(inds)):
                    #    peaks_order[m][ss]=lock_time[inds[m]]

                    if plotting and len(img)>0:
                        im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                        ax1.xaxis.labelpad = 0
                        if sec>0: 
                            ax1.xaxis.set_visible(False)
                        ax1.yaxis.set_visible(False)
                        #ax1.set_xlim(10,70)

                #***********************************************************************************************
                #Plot Luczak plots - non-normalized
                if False:
                    if plotting:
                        ax1 = plt.subplot(gs[14:17,int(width_plots*.05)+sec*3+1:int(width_plots*.05)+4+sec*3])
                        ax1.set_xlim(0,small_window*2)

                    img=[]
                    for unit in range(len(unit_ids)):
                        if unit in Sorts[rec].uid:
                            if max(fit_sum[ss][unit])!= 0.0:
                                img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)])
                            else: 
                                temp = np.zeros(small_window*2,dtype=np.float32)
                                img.append(temp)
                        else: 
                            temp = np.zeros(small_window*2,dtype=np.float32)
                            img.append(temp)

                    #Clip all points to 1Hz max firing rate; 
                    img = np.clip(img,-1E7,1)
                    
                    if plotting and len(img)>0:
                        im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                        ax1.xaxis.labelpad = 0
                        ax1.xaxis.set_visible(False)
                        #if rec>0: ax1.yaxis.set_visible(False)
                        ax1.set_ylim(0,len(img))
                ##***********************************************************************************************
                ##Plot Luczak plots - only locking units - Firing Sequence Order - PASSED PVAL TEST
                #if plotting:
                    #ax1 = plt.subplot(gs[:5,7:10])
                    ##ax1 = plt.subplot(gs[int(sqrt(len(Sorts[rec].units)))+1:int(sqrt(len(Sorts[rec].units)))+5, 6:9])
                    #ax1.set_xlim(0,small_window*2)
                    #ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

                #img=[]
                #lock_time=[]
                #for unit in range(len(Sorts[rec].units)):
                    #if p_val_array[ss][unit]<= 0.01:
                        #lock_time.append(np.argmax(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]))
                        #img.append(fit_sum[ss][unit][500-int(small_window/window):500+int(small_window/window)]/max(fit_sum[ss][unit]))

                #inds = np.array(lock_time).argsort()
                #img=np.array(img)[inds]

                #for m in range(len(inds)):
                    #peaks_pval[m][ss]=lock_time[inds[m]]

                #if plotting and len(img)>0:
                    #im = ax1.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                    #ax1.xaxis.labelpad = 0
                    #plt.title("P_val < 0.01 units (firing order w/in small window)")
                    #ax1.set_ylim(0,len(img))

                ##***********************************************************************************************
                #Plot large windows *******************************************************
                #Depth
                if True:
                    if plotting:
                        #ax1 = plt.subplot(gs[18:21,7+rec*3:10+rec*3])
                        ax1 = plt.subplot(gs[19:22,int(width_plots*.05)+sec*3+8:int(width_plots*.05)+11+sec*3])

                        ax1.set_xlim(0,large_window*2)
                        if rec>0: ax1.yaxis.set_visible(False)

                    img=[]
                    for unit in range(len(unit_ids)):
                        if unit in SUA_sort.unit_ids:
                            if max(fit_sum[ss][unit])!= 0.0:
                                img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))
                                peaks_depth_large[sec][unit][ss]= fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit])
                            else: 
                                temp = np.zeros(large_window,dtype=np.float32)
                                img.append(temp)
                                peaks_depth_large[sec][unit][ss]= temp

                        else: 
                            temp = np.zeros(large_window,dtype=np.float32)
                            img.append(temp)
                            peaks_depth_large[sec][unit][ss]= temp

                    if plotting and len(img)>0:
                        im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                        ax1.xaxis.labelpad = 0
                        if sec>0: ax1.yaxis.set_visible(False)
                        ax1.set_ylim(0,len(img))
                        ax1.xaxis.set_visible(False)

                ##***********************************************************************************************
                #Plot Luczak plots - only locking units - by firing order
                if False:
                    if plotting:
                        ax1 = plt.subplot(gs[21:24,int(width_plots*.1)+sec*3+1:int(width_plots*.1)+4+sec*3])
                        ax1.set_xlim(0,large_window*2)
                        ax1.set_ylim(0,len(img))

                    img=[]
                    lock_time=[]
                    for unit in range(len(unit_ids)):
                        if unit in SUA_sort.unit_ids:
                            if max(fit_sum[ss][unit])!= 0.0:
                                lock_time.append(np.argmax(fit_sum[ss][unit]))
                                img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))
                            else: 
                                temp = np.zeros(large_window,dtype=np.float32)
                                img.append(temp)
                        else: 
                            temp = np.zeros(large_window,dtype=np.float32)
                            img.append(temp)
                            
                    inds = np.array(lock_time).argsort()
                    img=np.array(img)[inds]

                    if plotting:
                        if len(img)>0:
                            im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                        ax1.xaxis.labelpad = 0
                        if sec>0: 
                            ax1.yaxis.set_visible(False)
                            ax1.xaxis.set_visible(False)

                ##***********************************************************************************************
                ##Plot Luczak plots - only locking units - by firing order
                #if plotting:
                    #ax1 = plt.subplot(gs[int(sqrt(len(Sort1.units)))+6:int(sqrt(len(Sort1.units)))+10, 6:9])
                    #ax1.set_xlim(0,large_window*2)
                    #ax1.set_xlabel("Time (ms)", weight = 'bold', color = 'black', fontsize = 12)

                #img=[]
                #lock_time=[]
                #for unit in range(len(Sort1.units)):
                    #if p_val_array[ss][unit]<= 0.01:
                        #lock_time.append(np.argmax(fit_sum[ss][unit]))
                        #img.append(fit_sum[ss][unit][500-int(large_window/(2*window)):500+int(large_window/(2*window))]/max(fit_sum[ss][unit]))

                #inds = np.array(lock_time).argsort()
                #img=np.array(img)[inds]

                #if plotting and len(img)>0:
                    #im = ax1.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')                
                    #ax1.xaxis.labelpad = 0
                    #plt.title("P_val < 0.01 units (depth order; large window)")
                    #ax1.set_ylim(0,len(img))

                ##Plot all unit locking scattering times for each LFP event
                #if plotting:
                    #ax1 = plt.subplot(gs[4:7, 11:14])
                    #xx1=[]
                    #for p in range(len(Sort1.units)):
                        #ax1.scatter(peaks_depth[p][ss],-p*5, color='blue')
                        #xx1.append(peaks_depth[p][ss])
                    #yy1 = np.arange(0,-len(Sort1.units)*5, -5)
                    #ax1.plot(xx1,yy1, color='blue')
                    #ax1.set_xlim(0, 2*small_window)
                    #ax1.set_ylabel("Depth Order")
                    
                ##Plot all unit locking scattering times for each LFP event
                #if plotting:
                    #ax2 = plt.subplot(gs[7:10, 11:14])
                    #for p in range(len(Sort1.units)):
                        #ax2.scatter(peaks_order[p][ss],-p*5, color='red')
                    #ax2.set_xlim(0, 2*small_window)                    
                    #ax2.set_ylabel("Firing Order")
                    
                ##Plot all unit locking scattering times for each LFP event
                #if plotting:
                    #ax3 = plt.subplot(gs[10:13, 11:14])
                    #for p in range(len(Sort1.units)):
                        #ax3.scatter(peaks_pval[p][ss],-p*5, color='black')
                    #ax3.set_xlim(0, 2*small_window)          
                    #ax3.set_ylabel("Firing Order (pval)")

            #SHOW PLOTS
            if plotting: 
                #mng = plt.get_current_fig_manager()
                #mng.resize(*mng.window.maxsize())
                plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
                plt.suptitle(LFP_sort.directory+"\n Min firing rate: " + str(cutoff_freq) + " Hz   Time chunks: " + str(int(sec_length/60.))+" mins.", fontsize=20)
                plt.show()    

            #********************************** PLOT INDIVIDUAL UNIT RESPONSE PATTERNS ************************
            #gs = gridspec.GridSpec(25,40)

            counter=0
            for unit in unit_ids:
                ax = plt.subplot(int(sqrt(len(unit_ids)))+1, int(sqrt(len(unit_ids)))+1,counter+1)

                img=[]
                for sec in range(len(SUA_sort.sec_len)):
                    img.append(peaks_depth[sec][unit][ss])
                    #print peaks_depth[sec][unit][ss]
                    #time.sleep(1)
                print "LFP :" ,ss, " unit: ", unit

                counter+=1
                
                im = ax.imshow(img, extent=[0,small_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax.xaxis.labelpad = 0
                #ax.yaxis.set_visible(False)
                #ax.xaxis.set_visible(False)
                ax.set_ylim(0,len(img))
                
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
            plt.suptitle("Single unit responses ("+str(small_window*2)+"ms window) to LFP Cluster : "+str(ss))
            
            plt.show()

            counter=0
            for unit in unit_ids:
                ax = plt.subplot(int(sqrt(len(unit_ids)))+1, int(sqrt(len(unit_ids)))+1,counter+1)

                #temp_x = []
                img=[]
                for sec in range(len(SUA_sort.sec_len)):
                    img.append(peaks_depth_large[sec][unit][ss])
                print "LFP :" ,ss, " unit: ", unit

                counter+=1
                
                im = ax.imshow(img, extent=[0,large_window*2,0,len(img)], aspect='normal', interpolation='none')
                ax.xaxis.labelpad = 0
                ax.yaxis.set_visible(False)
                ax.xaxis.set_visible(False)
                ax.set_ylim(0,len(img))
                
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
            plt.suptitle("Single unit responses ("+str(large_window*2)+"ms window) to LFP Cluster : "+str(ss))
            
            plt.show()
            
            #quit()

            ##PLotting cumulative locking times: *******************************************

            #plt.suptitle(Sort1.directory[27:] + "  " + 
            #str(round(Sort1.tsf.n_vd_samples/Sort1.tsf.SampleFrequency/60.0,2))
            #+ " mins) \n " + "window (x-axis): "+str(window*2*1E3) + " ms,   max spike count (y-axis): 100 spikes.", fontsize=20)

            #ax1 = plt.subplot(gs[0:4, 0:18])
            #for ss in range(start_lfp, len(Sort2.units), 1):
                #xx1=[]
                #for p in range(len(Sort1.units)):
                    #if peaks_depth[p][ss]<1000:
                        #ax1.scatter(peaks_depth[p][ss]+ss*small_window*2,-p, color='black')
                        #xx1.append(peaks_depth[p][ss]+ss*small_window*2)
                    #else:
                        #xx1.append(None)

                #ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

                #yy1 = np.arange(0,-len(Sort1.units), -1)
                #ax1.plot(xx1,yy1, color='blue')
                #ax1.xaxis.set_visible(False)
                #ax1.set_ylabel("Depth Order")
                #ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)


            #ax1 = plt.subplot(gs[4:8, 0:18])
            #for ss in range(start_lfp, len(Sort2.units), 1):
                #xx1=[]
                #for p in range(len(Sort1.units)):
                    #if peaks_order[p][ss]<1000:
                        #ax1.scatter(peaks_order[p][ss]+ss*small_window*2,-p, color='black')
                        #xx1.append(peaks_order[p][ss]+ss*small_window*2)
                    #else:
                        #xx1.append(None)

                #ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

                #yy1 = np.arange(0,-len(Sort1.units), -1)
                #ax1.plot(xx1,yy1, color='red')
                #ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)
                #ax1.xaxis.set_visible(False)
                #ax1.set_ylabel("Firing Order")

            #ax1 = plt.subplot(gs[8:12, 0:18])
            #for ss in range(start_lfp, len(Sort2.units), 1):
                #xx1=[]
                #for p in range(len(Sort1.units)):
                    #if peaks_pval[p][ss]<1000:
                        #ax1.scatter(peaks_pval[p][ss]+ss*small_window*2,-p, color='black')
                        #xx1.append(peaks_pval[p][ss]+ss*small_window*2)
                    #else:
                        #xx1.append(None)
                #yy1 = np.arange(0,-len(Sort1.units), -1)

                #ax1.plot([100+ss*small_window*2,100+ss*small_window*2],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)

                #if all(x is None for x in xx1):
                    #plt.plot([100,100],[0,-len(Sort1.units)], 'r--', color='black', linewidth=2)
                #else:
                    #ax1.plot(xx1,yy1, color='black')
                #ax1.set_ylabel("Firing Order (pval)")
                #ax1.set_xlim(left=0, right=len(Sort2.units)*small_window*2)

            ##SHOW PLOTS
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
            #plt.show()   
        



def sinc_interp(x, s, u):
    """
    Interpolates x, sampled at "s" instants
    Output y is sampled at "u" instants ("u" for "upsampled")
    
    from Matlab:
    http://phaseportrait.blogspot.com/2008/06/sinc-interpolation-in-matlab.html        
    """
    
    if len(x) != len(s):
        raise Exception, 'x and s must be the same length'
    
    # Find the period    
    T = s[1] - s[0]
    
    sincM = tile(u, (len(s), 1)) - tile(s[:, newaxis], (1, len(u)))
    y = dot(x, sinc(sincM/T))
    return y


def Plot_SD_vs_specgram(Sort1, tsf, lfp):  #Nick - GRAY SCALE PLOTS!!!

        
    #Compute SD of each channel and plot
    bin_width = .00008 #bin width in seconds; NB: MAKE SURE INTEGER MULTIPLE OF TIME STEP
    step= bin_width*tsf.SampleFrequency
    #bin_average = 3 #steps to average at each point
    depth = []
    
    offset=5*tsf.SampleFrequency
    #offset_array=np.zeros(5*tsf.SampleFrequency,dtype=np.float32)
    print tsf.vscale_HP
    
    channel_sd_array=[]
    channel_val_array=[]
    channel_sd_mean_array=[]
    channel_sd_ave_array=[]
    for i in range(tsf.n_electrodes):
        print tsf.n_electrodes, i, " depth: ", tsf.Siteloc[i*2+1]
        channel_sd=[]
        channel_val=[]
        channel_sd_mean=[]
        channel_sd_ave=[]
        dd1=0
        for j in range(int(tsf.n_vd_samples/1000./(tsf.SampleFrequency*bin_width))):
            dd = tsf.ec_traces[i][offset+j*step:offset+(j+1)*step]
            channel_sd.append(np.std(dd))
            channel_val.append(dd)
            channel_sd_mean.append(dd-np.mean(dd))
            channel_sd_ave.append(dd-dd1)
            dd1=dd

        channel_sd_array.append(np.hstack(channel_sd))
        channel_val_array.append(np.hstack(channel_val))
        channel_sd_mean_array.append(np.hstack(channel_sd_mean))
        channel_sd_ave_array.append(np.hstack(channel_sd_ave))
        depth.append(tsf.Siteloc[i*2+1])

    inds=np.array(depth).argsort()

    temp_array=[]
    temp_array1=[]
    temp_array2=[]
    temp_array3=[]
    temp_array4=[]
    for i in range(len(inds)):
        temp_array.append(channel_sd_array[inds[i]])
        temp_array1.append(channel_val_array[inds[i]])
        temp_array3.append(channel_sd_mean_array[inds[i]])
        temp_array4.append(channel_sd_ave_array[inds[i]])

    img= lfp.specgram[::-1]
    im = plt.imshow(img, extent=[0,tsf.n_vd_samples/tsf.SampleFrequency,0,60], aspect='normal')#,interpolation='none') #extent=lfp.extent, cmap=lfp.cm)

    img0 = np.array(temp_array)
    im = plt.imshow(img0, extent=[0,tsf.n_vd_samples/1000./tsf.SampleFrequency,-100,0], aspect='normal', vmax=150,interpolation='none')
    
    img1 = np.array(temp_array1)
    im = plt.imshow(img1, extent=[0,tsf.n_vd_samples/1000./tsf.SampleFrequency,-200,-100], aspect='normal', vmax=500, cmap=cm.Greys,interpolation='none')

    img3 = np.array(temp_array3)
    im = plt.imshow(img3, extent=[0,tsf.n_vd_samples/1000./tsf.SampleFrequency,-300,-200], aspect='normal', vmax=100, cmap=cm.Greens,interpolation='none')
    
    img4 = np.array(temp_array4)
    im = plt.imshow(img4, extent=[0,tsf.n_vd_samples/1000./tsf.SampleFrequency,-400,-300], aspect='normal', vmax=100, cmap=cm.Reds,interpolation='none')

    plt.ylim(-300,60)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.title(tsf.tsf_name+"       depth: " + str(tsf.Siteloc[i*2+1])+ " um from top of probe")
    plt.show()
        
def Plot_MUA_vs_LFP_events(sim_dir, Sorts_sua, Sorts_lfp, lfp,  rec):

    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']
    
    #Load specgram

    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    rec_index = lfp.rec.index(rec)                      #compute relative index of recording in lists of concatenated recs
    print "rec_index: ", rec_index

    #Set window to plot data in seconds
    start_traces = 0                                    #Start of window; use 0 sec for start of rec
    end_traces = len(lfp.data[rec_index][0])*1E-3       #End of window; usually whole recording

    #Plot fft specgram + sync index
    if False: 
        #plt.close()
        font_size = 30
        ax = plt.subplot(1,1,1)
        channel = 9
        P, extent = Compute_specgram_signal(lfp.data[rec_index][channel][start_traces*1E6:end_traces*1E6],1000)
        plt.imshow(P, extent=extent, aspect='auto')

        data_in = lfp.data[rec_index][channel][start_traces*1E6:end_traces*1E6]
        si_limit=0.7
        si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)
        
        plt.plot(t, si*50-60, linewidth=5, color='black')
        plt.plot([0,max(t)],[-50*.3-10,-50*.3-10], 'r--', color='red', linewidth = 3, alpha=0.8)
        plt.plot([0,max(t)],[-10,-10], color='black', linewidth = 1, alpha=0.2)
        plt.plot([0,max(t)],[-60,-60], color='black', linewidth = 1, alpha=0.2)
        xx = np.linspace(0,max(t),int(max(t)/500.))
        x_label = np.round(np.linspace(0, max(t)/60.,int(max(t)/500.)))

        plt.xticks(xx, x_label, fontsize=20)
        ax.tick_params(axis='both', which='both', labelsize=font_size)
        #y_vals = si*0.
        #plt.fill_between(t, 50*si-50,-50+y_vals, facecolor='green', alpha=0.2)
        #plt.fill_between(t, 50*si-50,0, facecolor='green', alpha=0.2)                  
        plt.xlabel("Time (mins)", fontsize = font_size)
        plt.title(lfp.rec_name[rec_index], fontsize=font_size)
        plt.show()


    ax_image = plt.subplot(1,1,1)
    plt.title(lfp.rec_name[rec_index])
    
    plt.xlim(start_traces, end_traces)
    plt.ylim(-250, 110)
    plt.xlabel("Time (sec)")
    plt.ylabel("MUA Firing Rate*    LFP Cluster Raster   Frequency (Hz)", weight = 'bold')
    ax_image.get_xaxis().get_major_formatter().set_useOffset(False)
    
    #********** PLOT MTSpecgram ************
    #Find deepest channel by parsing sorted units:
    
    channel = 9  #max(Sorts_sua[rec_index].maxchan)
    print "deepest ch: ", channel
    file_name = glob.glob(sim_dir+rec+"*")[0]
    file_name = file_name.replace(sim_dir, '')
    
    fname = sim_dir+file_name+'/'+"specgram_ch_"+str(channel)+"_start_"+str(start_traces)+"_end_"+str(end_traces)
        
    if (os.path.exists(fname+".npy")==False):
        img, extent = Multitaper_specgram(lfp.data[rec_index][channel][start_traces*1E6:end_traces*1E6], 1E3)
        np.save(fname, img)
        np.save(fname+'_extent', extent)
    else: 
        img = np.load(fname+'.npy')
        extent = np.load(fname+'_extent.npy')

    extent = extent[0]+start_traces, extent[1]+start_traces, extent[2], extent[3]
    im = ax_image.imshow(img, extent=extent, aspect='auto') #, extent=tsf.extent, cmap=tsf.cm)
    
    #********** PLOT LFP Events ************
    print "No. LFP events: ",  len(Sorts_lfp[rec_index])
    for k in range(len(Sorts_lfp[rec_index])):
        spike_times = Sorts_lfp[rec_index][k]*1E-3  #Convert back to seconds from raw lfp event times

        print spike_times

        ymin = 0
        #ymax = -2*len(Sorts_sua[rec_index].units)
        ymax = -200 #*len(Sorts_sua[rec_index].units)
        
        plt.vlines(spike_times, ymin, ymax, linewidth=5, color=colors[k], alpha=0.2)

    #********* PLOT SUA spike locations ************
    #Order units by firing rate:
    if True:
        n_spikes = []
        for k in range(len(Sorts_sua[rec_index].units)):
            n_spikes.append(len(Sorts_sua[rec_index].units[k]))
        indexes = np.argsort(n_spikes)[::-1]

    counter=0
    for k in indexes:
        spike_times = np.float32(Sorts_sua[rec_index].units[k])/float(Sorts_sua[rec_index].samplerate) #convert to seconds
        temp_indexes = np.where(np.logical_and(spike_times>=float(start_traces), spike_times<=float(end_traces)))[0]

        yy = spike_times[temp_indexes]

        ymin = -counter
        ymax = ymin-1
        counter+=1
        
        plt.vlines(yy, ymin, ymax, linewidth=2, color=colors[k%10], alpha=1.0)

    #Order units by depth
    n_spikes = []
    for k in range(len(Sorts_sua[rec_index].units)):
        maxchan = Sorts_sua[rec_index].maxchan[0] #chanpos[0][1]
        n_spikes.append(Sorts_sua[rec_index].chanpos[maxchan][1])
    indexes = np.argsort(n_spikes)[::-1]

    plt.plot([0,len(lfp.data[rec_index][channel])],[-counter-2, -counter-2], color='black', linewidth = 3, alpha = 0.5)

    counter+=5
    for k in indexes:
        spike_times = np.float32(Sorts_sua[rec_index].units[k])/float(Sorts_sua[rec_index].samplerate) #convert to seconds
        temp_indexes = np.where(np.logical_and(spike_times>=float(start_traces), spike_times<=float(end_traces)))[0]

        yy = spike_times[temp_indexes]

        ymin = -counter
        ymax = ymin-1
        counter+=1
        
        plt.vlines(yy, ymin, ymax, linewidth=2, color=colors[k%10], alpha=1.0)
    
    #********** PLOT MUA HISTOGRAMS ********************
    if False:
        mua=[]
        for i in range(len(Sorts_sua[rec_index].units)):
            mua.extend(np.array(Sorts_sua[rec_index].units[i])/Sorts_sua[rec_index].SampleFrequency)

        mua_bin_width = 0.020
        mua = np.array(mua)

        mua_plot=np.histogram(mua, bins = np.arange(0,end_traces,mua_bin_width))
        plt.plot(mua_plot[1][0:-1],mua_plot[0]*2-counter-20, linewidth=3, color='black', alpha=0.3)

    ##******** Compute synchrony index ************
    data_in = lfp.data[rec_index][channel][start_traces*1E6:end_traces*1E6]
    si_limit=0.7
    si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)
    
    plt.plot(t, 50*si-counter-50, linewidth=5, color='green')
    y_vals = si*0.
    plt.fill_between(t, 50*si-counter-50,-counter-50+y_vals, facecolor='green', alpha=0.2)
    plt.fill_between(t, 50*si-counter-50,-counter, facecolor='green', alpha=0.2)
    plt.plot([0,max(t)],[-counter-50*.3,-counter-50*.3], 'r--', color='red', alpha=0.8)

        
    print indexes

    plt.show()

    quit()
    
    #************** Plot cross-correlograms between 
    #Find MUA over threhsold:
  
    
    x = np.where(y[0]*2> 20)[0]
    x = y[1][x]
    print "MUA times over threshold: ", x[0:10]
    print len(x)
    
    y = spike_times
    
    #from brian import *
    #xcor = correlogram(x*1E3,y*1E3,width=20,bin=1,T=None)
    
    yunbiased = y -np.mean(y)
    ynorm = np.sum(yunbiased**2)
    acor = np.correlate(x, yunbiased, "same")/ynorm
    # use only second half
    #acor = acor[len(acor)/2:]

    plt.plot(acor, color='black', linewidth = 3)
    plt.show()
    
    
def Multitaper_specgram(time_series, sampfreq):
    
    plotting=False
    if plotting:
        ax2 = plt.subplot(1,1,1)
        ax2.autoscale(enable=True, tight=True)
    
    #**************************************
    
    s = time_series
    print "Length of recording: ", len(s)

    #******************************************************
    #Plot multi taper specgram

    #General parametrs
    f0=0.1
    f1=110
    sampfreq=1000

    #time parameters
    t0 = None
    t1 = None
    ts = np.arange(0,len(s),1.0)/sampfreq

    if t0 == None:
        t0, t1 = ts[0], ts[-1] # full duration
    if t1 == None:
        t1 = t0 + 10 # 10 sec window
            
    #Multi taper function parameters
    N = 512     #Number of tapers 
    NW = 40     #
    step = 10   #shift
    k = 6       #
    tm = 6.0    #time support of the 
    Np = 201    #

    w = np.hamming(N)

    h,Dh,Th = hermf(N, k, tm)
    E,V   = dpss(N, NW, k)
    
    print "Computing trf_specgram..."
    spec = tfr_spec(s, N, step, Np, k, tm)
    print "..."
    #OTHER TYPES OF FUNCTIONS... Did not check in detail
    #mpsd  = mtm_psd(s, NW)
    #J     = mtfft(s, NW)
    #spec  = stft(s, w, step)
    #mspec = mtm_spec(s, N, step, NW)
    #tspec_zoom = tfr_spec(s, N, step, Np, k, tm, fgrid=np.linspace(0.1,0.475,512))
    #tspec_log = tfr_spec(s, N, step, Np, k, tm, fgrid=log_fgrid(0.1, 0.45, 256))
    

    extent = t0, t1, f0, f1

    lo = 0
    hi = int(float(N)/1000. * f1)
    print lo, hi

    spec = spec[lo:hi]

    zis = np.where(spec == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        spec[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float
        minnzval = spec.min() # get minimum nonzero value
        spec[zis] = minnzval # replace with min nonzero values
    spec = 10. * np.log10(spec) # convert power to dB wrt 1 mV^2?

    p0=-40
    p1=None
    if p0 != None:
        spec[spec < p0] = p0
    if p1 != None:
        spec[spec > p1] = p1

    if plotting:
        im = ax2.imshow(spec[::-1], extent=extent, aspect='auto', cmap=None) #, vmin=0, vmax=1e2)
        plt.show()

    return spec[::-1], extent

def Schneidman_vectors(sim_dir, file_name, bin_width):
    
    work_dir = sim_dir + file_name + "/"
    ptcs_flag = 1
    Sort = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
    Sort.name=file_name
    Sort.filename=file_name
    Sort.directory=work_dir
    
    #Find length of recording using largest spike time
    last_spikes = []
    for k in range(len(Sort.units)):
        temp_index = np.where(np.array(Sort.units[k])<1E10)[0] #Exclude very large spikes due to SS .ptcs spike time export error
        spikes = np.array(Sort.units[k])[temp_index]
        if len(spikes)>0: last_spikes.append(np.max(spikes))
    last_spike = np.max(last_spikes) #Don't use length of recording to make vector bins if .tsf is missing; use last spike

    #******** Order units by firing rate of firing rates
    n_spikes = []
    for f in range(len(Sort.units)): #Loop over all ptcs file-objects
        n_spikes.append(len(Sort.units[f]))
    indexes = np.argsort(n_spikes)[::-1]

    #******** Initialize and generate vectors
    n_units = len(Sort.units)
    n_bins = float(last_spike)/bin_width/Sort.samplerate +2 #Make sure use float otherwise rounding off
    print "n_units: ", n_units, "   n_bins: ", n_bins
    vectors = np.zeros((n_units,n_bins), dtype = np.int16)
    print "matrix shape: ", vectors.shape

    counter=0
    for k in indexes:
        bin_time = np.float32(Sort.units[k])/bin_width/Sort.samplerate
        for p in bin_time:
            if p <1E10:
                vectors[counter][p] += 1
        counter+=1

    return vectors, Sort, last_spike


def on_click(event):
    
    global coords, images_temp, ax, fig, cid, n_pix
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[200][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=128

        ax.imshow(images_temp[200])
        fig.canvas.draw()

    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
#def Define_generic_mask(images_processed, file_dir, file_name, window, img_rate, n_pixels):
def Define_generic_mask(images_processed, work_dir, file_name, img_rate, n_pixels):

    global coords, images_temp, ax, fig, cid, n_pix
    
    images_temp = images_processed.copy()
    fig, ax = plt.subplots()
    
    n_pix = n_pixels
    img_r = img_rate
    
    coords=[]

    ax.imshow(images_processed[200])#, vmin=0.0, vmax=0.02)
    ax.set_title("Compute generic (outside the brain) mask")
    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()

    #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
    #Search points outside and black them out:
    all_points = []
    for i in range(len(images_processed[0][0])):
        for j in range(len(images_processed[0][0])):
            all_points.append([i,j])

    all_points = np.array(all_points)
    vertixes = np.array(coords) 

    from matplotlib.path import Path
    vertixes_path = Path(vertixes)
    
    mask = vertixes_path.contains_points(all_points)
    counter=0
    coords_save=[]
    for i in range(len(images_processed[0][0])):
        for j in range(len(images_processed[0][0])):
            if mask[counter] == False:
                images_processed[200][i][j]=0
                coords_save.append([i,j])
            counter+=1

    fig, ax = plt.subplots()
    ax.imshow(images_processed[200])
    plt.show()
   
    genericmask_file = work_dir + 'genericmask.txt'
    np.savetxt(genericmask_file, coords_save)

    print "Finished Making General Mask"
    
    return coords_save


