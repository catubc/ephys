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

#************************** LOAD PTCS FILES ****************************

#np.random.seed(12345678)  #fix random seed to get the same result
#np.random.seed(seed=None)

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
    

    
#***************************************************************************************
#************************************TIM DATA******************************************
#***************************************************************************************

ptcs_flag = 1

#TIM'S DATA
#sim_dirs = [
#'/media/cat/4TB/in_vivo/tim/2015-5-6/',
#'/media/cat/4TB/in_vivo/tim/2015-5-12/',
#'/media/cat/4TB/in_vivo/tim/2015-5-13/',
#'/media/cat/4TB/in_vivo/tim/2015-7-8/',
#'/media/cat/4TB/in_vivo/tim/2015-7-14/',
#'/media/cat/4TB/in_vivo/tim/2015-7-22/',
#'/media/cat/4TB/in_vivo/tim/2015-7-23/',
#'/media/cat/4TB/in_vivo/tim/2015-7-30/'
#]

#sim_dir = '/media/cat/4TB/in_vivo/tim/2015-7-30/2015-7-30-1/'

#file_name = '2015-11-18-18-deep-iso0_hp'
#file_name = '2015-11-18-14-allelectrodein-iso0_hp'

#file_name = '2015-11-18-8-9electrodein-iso0_hp'

#sim_dir = '/media/cat/12TB/in_vivo/tim/2015-11-18/' +file_name[0:-3] +'/'

#file_name = '2015-11-27-13-deep-iso1_hp'
#file_name = '2015-12-2-14-allelectrodeinthalamus-is0_hp'
#sim_dir = '/media/cat/12TB/in_vivo/tim/2015-11-27/' +file_name+'/'



file_dir = '/media/cat/12TB/in_vivo/tim/2016-1-14/' #/media/cat/12TB/in_vivo/tim/2015-12-11/'  #*************************************
file_names = [
'2016-1-14-13-allelectrodeinthalamus-iso1'

]

#Save single units into individual files (seconds)
if True:

    for file_name in file_names:
        
        sim_dir = file_dir +file_name +'/'

        work_dir = sim_dir
        Sort1 = Loadptcs(file_name+'_hp', work_dir, ptcs_flag, save_timestamps=False)
        Sort1.name=file_name
        Sort1.filename=file_name
        Sort1.directory=work_dir
                            
        for j in range(len(Sort1.units)):
            with open(sim_dir +"/unit_"+str(j).zfill(2)+"_channel_"+ str(Sort1.maxchan[j]+1).zfill(2)+"_ptp_"+str(int(Sort1.ptp[j])).zfill(3)+".csv", "w") as f:
                writer = csv.writer(f)
                for k in range(len(Sort1.units[j])):
                    writer.writerow([round(Sort1.units[j][k]/(1E+4*2),6)])     

quit()



min_ptp = 0
min_spikes = 0

if False:        #Firing rate analysis

    counter=0
    burstiness=[]
    burstiness_ptps=[]
    burstiness_x=[]
    bin_width = .001
    length = 1.

    isi=np.zeros((length/bin_width-1),dtype=np.float32)
    isi_log=np.zeros((length/bin_width-1),dtype=np.float32)

    plotting = True
    recording_location = "thalamus"

    cell_counter=0
    fire_rate=[]
    for sim_dir in sim_dirs:
        dir_counter=0

        for sim_dir, subdirs, files in os.walk(sim_dir):
            files = os.listdir(sim_dir)
            sim_dir = sim_dir + '/'

            for i in range(len(files)):
                cell_type=[]
                print sim_dir + files[i]
                if "ptcs" in files[i]:

                    with open(sim_dir +"recording_loc", "r") as f:
                        recording_type = f.read().strip()
                    
                    #print recording_type
                    #print recording_location
                    #quit()
                    if recording_type in recording_location:
                        file_name = files[i][:-5]
                        work_dir = sim_dir
                        Sort1 = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
                        Sort1.name=file_name
                        Sort1.filename=file_name
                        Sort1.directory=work_dir
                        
                        #Load TSF file and plot specgrams for all channels
                        if True:
                            tsf_name = sim_dir + files[i][:-5]+'_hp.tsf'
                            f = open(tsf_name, "rb")
                            print "Loading "+ sim_dir + files[i][:-5] + '_hp.tsf'
                            print "Recording: ", counter
                            counter+=1
                            tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
                            tsf.sim_dir = sim_dir
                            tsf.tsf_name = tsf_name
                            f.close()
                            Sort1.tsf = tsf
                                
                        for unit in range(len(Sort1.units)):
                            if len(Sort1.units[unit])< min_spikes: continue
                            if Sort1.ptp[unit]< min_ptp: continue
                            
                            fire_rate.append(float(len(Sort1.units[unit]))/(float(Sort1.tsf.n_vd_samples)/float(Sort1.tsf.SampleFrequency)))

    print "# cells: ", len(fire_rate)
    bin_width2 = .5
    ax = plt.subplot(1,1,1)
    plt.title("Firing Rate Distributions   #" + str(len(fire_rate)) + " cells")
    ax.set_xlabel("Firing Rate (Hz)")
    
    y=np.histogram(np.float32(fire_rate), bins = np.arange(0.005,50,bin_width2))
    #y=np.histogram(np.float32(fire_rate), bins=np.logspace(1E-2,1E+2,2000))
    plt.bar(y[1][:-1],y[0], bin_width2, color='blue')
    #ax.set_xscale('log')
    ax.set_xlim(0,40)
    plt.show()
    quit()

min_ptp = 100
min_spikes = 100

colors = np.array(['blue', 'red', 'green', 'magenta', 'cyan', 'black', 'yellow'])

if True:        #Do ISI analysis

    counter=0
    burstiness=[]
    burstiness_ptps=[]
    burstiness_x=[]
    bin_width = .001
    length = 1.

    isi=np.zeros((length/bin_width-1),dtype=np.float32)
    isi_log=np.zeros((length/bin_width-1),dtype=np.float32)

    plotting = True

    cell_counter=0
    for sim_dir in sim_dirs:
        dir_counter=0

        for sim_dir, subdirs, files in os.walk(sim_dir):
            files = os.listdir(sim_dir)
            sim_dir = sim_dir + '/'

            for i in range(len(files)):
                cell_type=[]
                print sim_dir + files[i]
                if "ptcs" in files[i]:

                    if "18" in files[i]: continue
                    file_name = files[i][:-5]
                    work_dir = sim_dir
                    Sort1 = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
                    Sort1.name=file_name
                    Sort1.filename=file_name
                    Sort1.directory=work_dir

                    for unit in range(len(Sort1.units)):

                        #SELECT PARTICULAR UNIT IN SORT
                        #if unit ==29: 
                            #pass
                        #else:
                            #continue

                        if Sort1.ptp[unit]<min_ptp: continue        #Exclude units w. < min_ptp amplitude
                        if len(Sort1.units[unit])<min_spikes: continue     #Exclude units w. less 100 spikes
                        cell_counter+=1
                        temp_isi = []
                        temp_isi2 = []
                        isi_index = 0
                        for q in range(1,len(Sort1.units[unit])-2,1):
                            pre_spike = (Sort1.units[unit][q]-Sort1.units[unit][q-1])/Sort1.samplerate*1E3 #Convert to ms time
                            post_spike = (Sort1.units[unit][q+1]-Sort1.units[unit][q])/Sort1.samplerate*1E3
                            if pre_spike == 0: continue     #Duplicates in data
                            if post_spike == 0: continue    #Dupliates in data
                            temp_isi.append(pre_spike)
                            temp_isi2.append(post_spike)
                            
                            if pre_spike<0.01 or post_spike<0.01:
                                isi_index+=1
                                
                        y_log = np.histogram(temp_isi, bins=np.logspace(0, length, length/bin_width)) # 
                        y = np.histogram(temp_isi, bins = np.arange(0,length,bin_width)) # 

                        isi+= np.hstack(y[0])
                        isi_log+=np.hstack(y_log[0])

                        plt.title(file_name+ " unit: " + str(unit) + " no. of spikes: " + str(len(Sort1.units[unit]))
                        + "  % burstiness: " + str(float(isi_index)/float(len(Sort1.units[unit]))) +
                        "\n color order" + str(colors))

                        burstiness.append(sum(y[0][0:5])/len(Sort1.units[unit])) # Percent spikes in first 5ms; 

                        #Cluster the log-log ISI plots
                        #K Means
                        if False:
                            n_clusters = 7
                            from sklearn import cluster, datasets
                            
                            clusters = cluster.KMeans(n_clusters)
                            
                            temp_data=np.zeros((len(temp_isi),2), dtype=np.float64)
                            for p in range(len(temp_isi)):
                                temp_data[p]=[log(temp_isi[p]),log(temp_isi2[p])]

                            clusters.fit(temp_data)

                        #Mean shift
                        if True:
                            from sklearn.cluster import MeanShift, estimate_bandwidth
                            from sklearn.datasets.samples_generator import make_blobs

                            ###############################################################################
                            # Generate sample data
                            #centers = [[1, 1], [-1, -1], [1, -1]]
                            #X, _ = make_blobs(n_samples=10000, centers=centers, cluster_std=0.6)

                            X=np.zeros((len(temp_isi),2), dtype=np.float64)
                            for p in range(len(temp_isi)):
                                X[p]=[log(temp_isi[p]),log(temp_isi2[p])]
                                
                            ###############################################################################
                            # Compute clustering with MeanShift

                            # The following bandwidth can be automatically detected using
                            bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=len(X))

                            ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
                            ms.fit(X)
                            labels = ms.labels_
                            cluster_centers = ms.cluster_centers_

                            labels_unique = np.unique(labels)
                            n_clusters = len(labels_unique)
                            #print labels_unique
                            #print n_clusters_
                            #print labels

                            print("number of estimated clusters : %d" % n_clusters)
                            #quit()

                        ax = plt.subplot(1,1,1)

                        #plt.scatter(np.array(temp_isi), np.array(temp_isi2), s=5, color=colors[clusters.labels_])
                        plt.scatter(np.array(temp_isi), np.array(temp_isi2), s=5, color=colors[labels])
                        mng = plt.get_current_fig_manager()
                        mng.resize(*mng.window.maxsize())
                        ax.set_xscale('log')
                        ax.set_yscale('log')
                        ax.set_xlim(1E0,1E4)
                        ax.set_ylim(1E0,1E4)

                        for k in range(n_clusters):
                            #cluster_indexes = np.hstack(np.where(clusters.labels_==k)[0])
                            cluster_indexes = np.hstack(np.where(labels==k)[0])
                            cluster_spikes = []
                            for s in range(len(cluster_indexes)):
                                cluster_spikes.append(Sort1.units[unit][cluster_indexes[s]])
                            print cluster_spikes
                            
                            #Save single units into individual files (seconds)
                            if False:
                                np.savetxt(sim_dir + 'unit_'+str(unit)+'_isi_cluster_'+colors[k], np.array(cluster_spikes)/(1.E+4*2.), delimiter="/n")

                        plt.show()
                        quit()



                        #CROSS CORRELOGRAM ANALYSIS
                        if False:
                            temp_array = np.int64(np.array(Sort1.units[unit])*40)
                            print temp_array
                            time_array = np.int64([-1E6,1E6])
                            print time_array
                            
                            acorr = xcorrm(temp_array,temp_array, time_array)*1E-3 #Convert back to ms

                            print "ACORR: ", acorr
                            acorr = np.histogram(acorr, bins = np.arange(-1E3,1E3,1)) # 
                            #print acorr.shape
                            acorr[0][1000]=0
                            acorr[0][999]=0
                            plt.bar(acorr[1][:-1],acorr[0],1, color='red')
                            plt.show()
                        
                        #plt.plot(y_log[0]/float(max(y_log[0])), color='blue', linewidth=2)
                        #CDF data analysis
                        cdf_sum = np.cumsum(y[0]*bin_width)
                        cdf_sum = cdf_sum/float(max(cdf_sum))
                        inds = np.argmax(cdf_sum>0.50)
                        burstiness_x.append(bin_width*inds)

                        if False: 
                            ax = plt.subplot(1,1,1)
                            plt.title("Unit: " + str(unit) + " no. of spikes: " + str(len(Sort1.units[unit])))
                            #plt.bar(y[1][:-1],y[0]/max(y[0]), bin_width, color='red', linewidth=2)
                            #plt.bar(y_log[1][:-1],y_log[0]/max(y_log[0]), bin_width, color='blue', linewidth=2)
                            plt.plot(y[1][:-1],y[0]/float(max(y[0])),  color='red', linewidth=2)

                            plt.plot([bin_width*inds, bin_width*inds], [0,1], 'r--', color='black', linewidth=2)
                            plt.plot(y[1][:-1],cdf_sum, color='green', linewidth=2)

                            mng = plt.get_current_fig_manager()
                            mng.resize(*mng.window.maxsize())
                            #ax.set_xscale('log')
                            plt.show()

                        #burstiness.append(float(sum(y[0][0:10]))/float(len(Sort1.units[unit]))) #Store % spikes w/in 1 second of each other.
                        #a = float(sum(y[0][1])+1)
                        #b = float(sum(y[0][0])+1)
                        #burstiness.append(10**(b/a)) #percent spike sin first bin vs first second
                        #burstiness_ptps.append(Sort1.ptp[unit])

                        #print isi

    #print isi
    print cell_counter
    bin_width2 = 0.001
    y=np.histogram(burstiness, bins = np.arange(0,1,bin_width2))
    plt.bar(y[1][:-1], y[0], bin_width2, color='black')

    #y=np.histogram(burstiness_x, bins = np.arange(0,1,bin_width2))
    #plt.bar(y[1][:-1], y[0], bin_width2, color='black')
    
    #ax = plt.subplot(1,1,1)
    #plt.plot(isi/max(isi),color='red',linewidth=2)
    #plt.plot(isi_log/max(isi_log), color='blue', linewidth=2)
    ##plt.xlim(0,1)
    #ax.set_xscale('log')
    plt.show()

    quit()
    #plt.scatter(burstiness, burstiness_ptps, color='blue')
    #plt.show()


peakwidth_data=[]
troughwidth_data=[]
prepeakheight_data=[]

min_ptp = 0
min_spikes = 20

#Compute spike shapes and plot cell-type distributions
if True:

    counter=0
    burstiness=[]
    burstiness_ptps=[]
    burstiness_x=[]
    bin_width = .001
    length = 1.
    cell_counter=0
    
    recording_location = "cortex"
    
    for sim_dir in sim_dirs:
        dir_counter=0

        for sim_dir, subdirs, files in os.walk(sim_dir):
            files = os.listdir(sim_dir)
            sim_dir = sim_dir + '/'

            for i in range(len(files)):
                cell_type=[]
                print sim_dir + files[i]
                if "ptcs" in files[i]:
                    
                    with open(sim_dir +"recording_loc", "r") as f:
                        recording_type = f.read()

                    #print recording_type                   
                    
                    if recording_location in recording_type:
                        print recording_location, recording_type
                        #quit()
                        file_name = files[i][:-5]
                        work_dir = sim_dir
                        Sort1 = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
                        Sort1.name=file_name
                        Sort1.filename=file_name
                        Sort1.directory=work_dir
                            
                        #Save single units into individual files (seconds)
                        if False:
                            for j in range(len(Sort1.units)):
                                print sim_dir
                                print files[i]
                                with open(sim_dir +"/unit_"+str(j)+"_channel_"+ str(Sort1.maxchan[j]+1)+".csv", "w") as f:
                                    writer = csv.writer(f)
                                    for k in range(len(Sort1.units[j])):
                                        writer.writerow([round(Sort1.units[j][k]/(1E+4*2),6)])        
                        #quit()
                        #Load TSF file and plot specgrams for all channels
                        if True:
                            tsf_name = sim_dir + files[i][:-5]+'_hp.tsf'
                            f = open(tsf_name, "rb")
                            print "Loading "+ sim_dir + files[i][:-5] + '_hp.tsf'
                            print "Recording: ", counter
                            counter+=1
                            tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
                            tsf.sim_dir = sim_dir
                            tsf.tsf_name = tsf_name
                            f.close()
                            Sort1.tsf = tsf

                        if False:
                            spec_ch = []
                            for chs in range(len(tsf.ec_traces)):
                                spec_ch.append(chs)
                            Plot_quickspecgram(tsf, show_plot=True, spec_ch=spec_ch)
                        
                        #PLOT SORTED UNIT
                        for unit in range(len(Sort1.units)):

                            if Sort1.ptp[unit]<min_ptp: continue
                            points=99

                            x = np.zeros((Sort1.tsf.n_electrodes,points),dtype=np.float32)
                            for n in range(Sort1.tsf.n_electrodes):
                                x[n]= Sort1.tsf.Siteloc[n*2] + np.array(arange(0,points,1))
                                
                            y1 = np.zeros((Sort1.tsf.n_electrodes,points), dtype=np.float32)
                            y1_norm = np.zeros((Sort1.tsf.n_electrodes,points), dtype=np.float32)

                            spike_list = np.arange(0,min(500,len(Sort1.units[unit])-1)) #DON"T GO TO END OF RECORD, use -1

                            #print spike_list
                            for j in spike_list: #range(min(len(Sort1.units[unit]),50)): #Plot up to 100 spikes from unit
                                #Exclude spikes too close to beginning or end;
                                if (int(Sort1.units[unit][j])+points*2/3) > Sort1.tsf.n_vd_samples: continue
                                if (int(Sort1.units[unit][j])-points/3) < 1: continue
                                for n in range(Sort1.tsf.n_electrodes):
                                        yy1 = Sort1.tsf.ec_traces[n][int(Sort1.units[unit][j])-points/3:int(Sort1.units[unit][j])+points*2/3]
                                        y1[n]+= yy1 #/10-Sort1.tsf.Siteloc[n*2+1]
                                        
                                        #if plotting: 
                                            #xx1 = x[n]
                                            #if Sort1.maxchan[unit]==n-1:
                                                #plt.plot(xx1, yy1, color='blue', alpha=1.) #/10-Sort1.tsf.Siteloc[n*2+1], color='blue', alpha=1.)
                                            #else:
                                                #plt.plot(xx1, yy1, color='black', alpha=0.7) #/10-Sort1.tsf.Siteloc[n*2+1], color='black', alpha=0.7)
                            y1 /= len(spike_list)
                            
                            plotting = False
                            if plotting: 
                                plt.suptitle(Sort1.name + " # units: " + str(len(Sort1.units)), fontsize=20)
                                ax = plt.subplot(int(sqrt(len(Sort1.units))), int(sqrt(len(Sort1.units)))+2, unit+1)

                            #for n in range(Sort1.tsf.n_electrodes):
                            n = Sort1.maxchan[unit]
                            y1[n]=-y1[n]/min(y1[n])*10

                            #Interpolate spike using sinc function
                            resolution = .05
                            u=np.arange(0,points,resolution)
                            y_sinc=sinc_interp(np.array(y1[n]),np.array(x[n]),u)
                            
                            if plotting: 
                                plt.plot(u,y_sinc, color='black', linewidth=2, alpha=0.8)

                            #Trough minimum
                            trough_index = np.argmin(y_sinc)

                            #Peak following trough
                            peak_index = np.argmax(y_sinc[trough_index:])+trough_index
                            
                            #End of peak (zero crossing)
                            endpeak_index = np.where(y_sinc[peak_index:]<=0.0)+peak_index
                            if len(endpeak_index[0])>0: 
                                endpeak_index=endpeak_index[0][0] #Take first index of value < 0.0 after positive peak
                            else:
                                endpeak_index = int(points/resolution)-1
                            
                            #Begnning of peak (zero crossing)
                            startpeak_index = np.where(y_sinc[trough_index:peak_index]>=0.0)+trough_index
                            if len(startpeak_index[0])>0: 
                                startpeak_index=startpeak_index[0][0] #Take first index of value < 0.0 after positive peak
                            else:
                                startpeak_index = 1

                            #Old "cap" phase, or pre-Na peak - if it exists 
                            prepeak_index = np.argmax(y_sinc[trough_index-10/resolution:trough_index])+trough_index-10/resolution #Search only 5 timesteps backwards = 500ms @ 20Khz
                            if y_sinc[prepeak_index]>0.5 : 
                                prepeak_index=prepeak_index #Take first index of value < 0.0 after positive peak
                            else:
                                prepeak_index= 1
                            
                            #Beginning of trough
                            starttrough_index = np.where(y_sinc[:trough_index]>=0.0)
                            if len(starttrough_index[0])>0: 
                                if (trough_index - starttrough_index[0][-1]) < 6/resolution: #Must be located w/in 250ms - maybe even lesss
                                    starttrough_index=starttrough_index[0][-1] #Take the last index of value 
                                else:
                                    starttrough_index = trough_index-(startpeak_index-trough_index)
                            else:
                                starttrough_index = trough_index-(startpeak_index-trough_index)

                            ##Check to see if positive peak exists
                            for m in range(startpeak_index,endpeak_index,1):
                                if y_sinc[m]>(y_sinc[peak_index]/4.): 
                                    temp1=m
                                    break
                            for m in range(endpeak_index,startpeak_index,-1):
                                if y_sinc[m]>(y_sinc[peak_index]/4.): 
                                    temp2=m
                                    break
                            peakwidth_data.append(temp2-temp1)
                                
                            if(temp2-temp1)>475:
                                cell_type.append("pyramid")
                            else:
                                cell_type.append("non-pyramid")
                                
                            #Check trough width
                            if starttrough_index!=1 and startpeak_index!=1: 
                                for m in range(starttrough_index,startpeak_index,1):
                                    if y_sinc[m]<(y_sinc[trough_index]/4.): 
                                        temp3=m
                                        break
                                for m in range(startpeak_index,starttrough_index,-1):
                                    if y_sinc[m]<(y_sinc[trough_index]/4.): 
                                        temp4=m
                                        break
                                troughwidth_data.append(temp4-temp3)
                            else:
                                troughwidth_data.append(0)                    
                            
                            ##Check to see if prepeak exists
                            prepeakheight_data.append(y_sinc[startpeak_index])
                            #prepeakheight_data.append(y_sinc[endpeak_index]-y_sinc[startpeak_index])
                            #prepeakheight_data.append(abs(y_sinc[prepeak_index]-y_sinc[peak_index])/(y_sinc[prepeak_index]+y_sinc[peak_index]))
                            #prepeakheight_data.append(y_sinc[peak_index])

                            if plotting:
                                #ax.set_xlim(0, points)
                                ax.set_title(str(unit+1) + " ("+str(len(Sort1.units[unit]))+")", fontsize = 12, fontweight='bold') 
                                ax.set_ylim(-12, 12) #top=-Sort1.tsf.Siteloc[(Sort1.maxchan[unit]+1)*2+1]+150, bottom=-Sort1.tsf.Siteloc[(Sort1.maxchan[unit]+1)*2+1]-150)
                        
                            #Plotting waveshape features
                            if plotting:
                               
                                plt.plot([trough_index*resolution,trough_index*resolution],[-10000,10000], color='red')
                                plt.plot([peak_index*resolution,peak_index*resolution],[-10000,10000], color='blue')
                                plt.plot([endpeak_index*resolution,endpeak_index*resolution],[-10000,10000], color='cyan')
                                plt.plot([startpeak_index*resolution,startpeak_index*resolution],[-10000,10000], color='cyan')
                                plt.plot([prepeak_index*resolution,prepeak_index*resolution],[-10000,10000], color='magenta')
                                plt.plot([starttrough_index*resolution,starttrough_index*resolution],[-10000,10000], color='green')
                                plt.plot([temp2*resolution,temp1*resolution],[y_sinc[temp2],y_sinc[temp1]], color = 'black', linewidth=2)
                                plt.plot([temp4*resolution,temp3*resolution],[y_sinc[temp4],y_sinc[temp3]], color = 'black', linewidth=2)
                                ax.xaxis.set_visible(False)


                        if plotting: 
                            mng = plt.get_current_fig_manager()
                            mng.resize(*mng.window.maxsize())
                            plt.show()

                        #np.savetxt(sim_dir + files[i][:-5] + '_type.csv', np.array(cell_type), delimiter="/n")
                        with open(sim_dir + files[i][:-5] + '_type.csv', "w") as f:
                            writer = csv.writer(f)
                            for j in range(len(cell_type)):
                                writer.writerow([cell_type[j]])
                                    


x = peakwidth_data
y = troughwidth_data
z = prepeakheight_data
colors=[]

#plt.scatter(peakwidth_data, burstiness, color='black')
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(peakwidth_data)):
    colors.append('b')
    #if x[i]<160:
    #    colors.append('r')
    #else:
    #    colors.append('b')
        
print "TOTAL # CELLS: ", len(peakwidth_data)

ax.scatter(x, y, z, c=colors, marker='o',s=15)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
plt.show()

        #quit()
            
        ##SHOW PLOTS
        #if plotting:
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.subplots_adjust(left=0.07, right=0.93, top=0.91, bottom=0.06)
            #plt.show()
        
            
quit()

#***************************************************************************************
#******************************** LOAD LFP & MAKE COMPRESSED LFP ***********************
#***************************************************************************************

if True:
    fname = sim_dir + sorted_file+'/' + sorted_file + '.lfp'
    lfp = Load_lfp(fname, Sort1)

    #Plot_PSD(lfp,'PSD')
    Plot_specgram(lfp, show_plot=True)
    

#***************************************************************************************
#*********************************** LOAD 2ND PTCS SORT*********************************
#***************************************************************************************


# **** Load LFP Spikes ****
if True:
    lfp_file=sorted_file+'_scaled_notch_multiunit'
    work_dir = sim_dir + sorted_file + "/"

    print "Loading: ", lfp_file

    Sort2 = Loadptcs(lfp_file, work_dir, 1, save_timestamps=True) #Auto load flag for Nick's data
    Sort2.name=lfp_file
    Sort2.filename=lfp_file
    Sort2.directory=work_dir




#***************************************************************************************
#********************************** ANALYSIS METHODS ************************************
#***************************************************************************************

#Plot_rasters(Sort1)
#Plot_firingrate(Sort1)

#Sort2=False
#Plot_LFP(lfp) #,Sort2)

#Plot_CSD(lfp)

#Plot_specgram(lfp)

#Plot_PSD(lfp,'PSD')

Plot_triggered_activity(Sort1,Sort2,lfp)

#Compute_kuebler(Sort1, Sort2)

#Compute_upstate(Sort1)
