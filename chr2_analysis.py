'''Compute SUA and MUA histograms as well as LFP and CSD for ChR2 experiments
'''

from tsf_ptcs_classes import *
import matplotlib.pyplot as plt
import numpy as np
from math import *
import imp

#areas = ['HL_L','V1_R','V1_R','PT_L','M1_R','M1_R','BC_L','V2_R',
#'V2_R','FL_L','BC_R','BC_R','RS_L','FL_R','FL_R','V1_L','HL_R','HL_R',
#'V2_L','PT_R','PT_R','M1_L','RS_R','RS_R']

#******** FILE NAMES *******
#2015_05_03
#dir_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/'
#file_name = '2016_05_03_tr2_deep_laser535_sorted_160503_171740'
#file_name = '2016_05_03_tr2_deep_laser535_sorted_mua_160503_171740'
#file_name = '2016_05_03_tr2_deep_laser535_2_sorted_160503_172757'
#file_name = '2016_05_03_tr2_deep_laser535_3_sorted_160503_173756'
#file_name = '2016_05_03_tr2_deep_laser535_4_sorted_160503_174635'

#dir_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/'
#file_name = '2016_05_03_tr3_deep_laser535_sorted_160503_185155'
#file_name = '2016_05_03_tr3_deep_laser535_2_160503_192518'

#2015_05_25
#dir_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25/'
#file_name = '2016_05_25_deep_laser535_tr1_160525_174613'
#dir_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/tsf_files/'
#file_name = '2016_05_25_tr1_deep_laser535_10ms_5sec_5rep_160525_174613.tsf'


#*****2016-05-25
dir_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/'
#file_name = '2016_05_25_tr1_deep_laser535_10ms_5sec_5rep_160525_174613'
#file_name = '2016_05_25_tr2_deep_laser535_10ms_5sec_5rep_160525_193958'
file_name = '2016_05_25_tr31_deep_laser535_7ms_5sec_5rep_160525_214419'
#file_name = '2016_05_25_tr32_deep_laser535_7ms_1sec_15rep_160525_224253'


#******* READ TRACK NUMBER ***** (Essentially recording number)
start = file_name.find('tr')
end = file_name[start:].find('_')
track = file_name[start:start+end]
print track


#******************************************************************************
#******************************READ LASER DATA ********************************
#******************************************************************************
laser_file = dir_name+ '/laser_files/'+track + '_laser_codes.txt'
laser_times = np.load(dir_name+ '/laser_files/'+track + '_laser_times.npy')
meta_data = np.load(dir_name+ '/laser_files/'+track + '_meta_data.npy')

#Load laser times == 1 
indexes = np.where(laser_times!=0)[0]   #Remove all OFF values

laser_t = []
laser_t.append(indexes[0])
counter = 0
for k in range(len(indexes)):
    if indexes[k]>(laser_t[counter]+10000):
        laser_t.append(indexes[k])
        counter+=1
laser_t = np.array(laser_t)
print "Length of laser_t arraY: ", len(laser_t)

#Load laser locations
laser_data = np.loadtxt(laser_file, dtype=np.str)
laser_data = np.delete(laser_data,0,1)
laser_data = np.delete(laser_data,1,1)

#Load locations into grid; 132 locations in this case;
grid_laser_times = []
for k in range(132):
    grid_laser_times.append(np.where(laser_data[:,1]==str(k))[0])

print "Examples of laser times for particular area:..."
print laser_t[grid_laser_times[0]]
print laser_t[grid_laser_times[131]]
#quit()


#******************************************************************************
#******************************SUA AND MUA PLOTS ******************************
#******************************************************************************
if True: 
    work_dir = dir_name+'ptcs_files/'
    sort_file = file_name+"_MUA"
    Sort = Loadptcs(sort_file, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
    Sort.name=file_name
    Sort.filename=sort_file
    Sort.directory=work_dir

    print "# units: ", len(Sort.units)
    print np.array(Sort.units[0])/Sort.samplerate

    counter=0
    stim_areas_list = np.arange(0,len(grid_laser_times),1)
    for i in range(len(stim_areas_list)):
        #loc = 121 - (i%12)*11+int(i/12+1)   #Need to flip this to go upside down;
        loc = (i%12)*11+int(i/12+1)   #Need to flip this to go upside down;
        ax=plt.subplot(12,11,loc)
        print "...plot loc: ", i, loc
        
        area_histogram= []
        for unit in range(len(Sort.units)):
            spike_times=[]
            for k in range(len(laser_t[grid_laser_times[i]])):
                spike_array = np.array(Sort.units[unit])/Sort.samplerate
                laser_flash = float(laser_t[grid_laser_times[i]][k])/Sort.samplerate
                temp2 = np.where(np.logical_and(spike_array>=laser_flash-1.0, spike_array<=laser_flash+1.0))[0]
                spike_times.extend(spike_array[temp2]-laser_flash)
            #print np.sort(np.array(spike_times))
            #quit()
            #x = np.linspace(-1,1,50)
            #y = np.histogram(spike_times, bins = x)
            #plt.plot(x[:-1], -(y[0]-np.average(y[0]))/10.+unit, linewidth=1,alpha=1)
            area_histogram.extend(spike_times)

        #plt.close()
        #x = np.linspace(-1,1,50)
        #plt.tick_params(axis='both', which='major', labelsize=10)
        #y = np.histogram(area_histogram, bins = x)
        #plt.plot(x[:-1], y[0], color='blue', linewidth=4)
        #plt.show()
        
        
        x = np.linspace(-1,1,50)
        y = np.histogram(area_histogram, bins = x)
        #plt.plot(x[:-1], -(y[0]-np.average(y[0])), linewidth=1,alpha=1)
        plt.plot(x[:-1], y[0], color='blue', linewidth=4)

        plt.ylim(len(Sort.units),0)
        plt.ylim(0,200)
        plt.xlim(-0.100, 0.100)

        #plt.title("Area: "+str(i), fontsize=15)
        #if counter==0: plt.ylabel("Unit # and Relative Firing Rate\n(deeper<----- unit depth ------>shallower)", fontsize=15)
        p = plt.axvspan(0.0, 0.010, facecolor='0.5', alpha=0.3)
        plt.tick_params(axis='both', which='major', labelsize=8)

        #if a>=12: plt.xlabel("Time from pulse (sec)", fontsize=8)
        #plt.plot([0,0],[0,1.43],'r--', linewidth=3, color='black', alpha=0.8)
        #plt.ylim(1.44,0)
        #plt.tick_params(axis='both', which='major', labelsize=8)
        ax.get_yaxis().set_visible(False)
        #ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        #plt.title(areas[stim_areas_list[i]],fontsize=15)
        #if (counter==0) or (counter==12):
        #    plt.ylabel("Depth along probe (mm)", fontsize=20)
        #    ax.get_yaxis().set_visible(True)
        #else:
        #    ax.get_yaxis().set_visible(False)
        if i%12==11: 
            x_labels = [-0.1,0.0,0.1]
            plt.xticks(x_labels, x_labels)

        counter+=1



    plt.suptitle(file_name, fontsize=30)
    plt.show()


    #All area plots
    area_histogram=[]
    for a in range(24):
        area_histogram.append([])
        for unit in range(len(Sort.units)):
            spike_times=[]
            for k in range(len(area_repeats[a])):
                spike_array = np.array(Sort.units[unit])/Sort.samplerate
                temp2 = np.where(np.logical_and(spike_array>=area_repeats[a][k]-1.0, spike_array<=area_repeats[a][k]+1.0))[0]
                spike_times.extend(spike_array[temp2]-area_repeats[a][k])
            area_histogram[a].extend(spike_times)

        ax= plt.subplot(4,6,a+1)

        x = np.linspace(-1,1,50)
        plt.title("Area: "+areas[a])
        plt.tick_params(axis='both', which='major', labelsize=10)
        y = np.histogram(area_histogram[a], bins = x)
        plt.plot(x[:-1], y[0], color='blue', linewidth=4)
        plt.ylim(0,600)

        p = plt.axvspan(0.0, 0.190, facecolor='0.5', alpha=0.5)

    plt.suptitle(file_name, fontsize=30)
    plt.show()

#******************************************************************************
#******************************LFP / CSD ANALYSIS ******************************
#******************************************************************************
if True:
    
    #***************Load original (unsorted/unfiltered) raw .tsf file
    file2_name = file_name.replace('sorted_','')
    work_dir=dir_name
    tsf_name = work_dir + 'tsf_files/'+file2_name+'.tsf'
    f = open(tsf_name, "rb")
    tsf = Tsf_file(f, work_dir+file2_name+'/')  
    f.close()


    #PLOT HISTOGRAMS BY AREA: 
    stim_areas_list = np.arange(0,len(grid_laser_times),1)
    
    lfp_ave = np.zeros((len(stim_areas_list),len(tsf.ec_traces),2*tsf.SampleFrequency),dtype=np.float32)

    plotting = False
    counter=0
    for a in stim_areas_list:
        ax=plt.subplot(111)
        print "...computing grid area: ", a
        #ax=plt.subplot(1,3,counter+1)
        #ax= plt.subplot(2,12,a+1)
        #plt.title("Area: "+areas[a], fontsize=15)
        for ch in range(len(tsf.ec_traces)):
            for t in grid_laser_times[a]:   #Times are already in sampletime steps
                t0 = laser_t[t]
                lfp_ave[counter][ch]+= tsf.ec_traces[ch][t0-tsf.SampleFrequency:t0+tsf.SampleFrequency][0:50000]
                
            x = np.linspace(-tsf.SampleFrequency,tsf.SampleFrequency,tsf.SampleFrequency*2)/tsf.SampleFrequency
            lfp_ave[counter][ch]=lfp_ave[counter][ch]/len(grid_laser_times[a])
        
            if plotting: plt.plot(x, -lfp_ave[counter][ch]/5.+tsf.Siteloc[ch*2+1])

        #if counter==0: plt.ylabel("Unit # and Relative Firing Rate\n(deeper<----- unit depth ------>shallower)", fontsize=15)
        if plotting:
            p = plt.axvspan(0.0, 0.190, facecolor='0.5', alpha=0.2)
            plt.tick_params(axis='both', which='major', labelsize=10)
            plt.ylim(1500,-100)
            #if a>=12: plt.xlabel("Time from pulse (sec)", fontsize=8)
            if counter==0: plt.ylabel("Electrode position (um)", fontsize=20)
            plt.show()

        counter+=1
        
    #********Compute CSD
    probe_layout = tsf.Siteloc[1::2]

    z = probe_layout*1E-3 #Convert to mm size
    csdmod = imp.load_source('csd_est_funds','csd_est_funcs.py')
    
    sigma = 0.3  # extracellular conductivity (mS/mm)
    b = 1.0   # assumed radius of the column (mm)
    SNR = 6.0    # average ratio of signal to noise on each channel

    [A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
    [W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
    [W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check

    #[fig100,fig101] = csdmod.show_operators(A,A0,W,W0,R,R0)    # display forward and inverse operators 
    #plt.show()

    no_lfps=1
    counter=0
    for i in range(len(stim_areas_list)):
    #for i in range(24):
        #loc = 121 - (i%12)*11+int(i/12+1)   #Need to flip this to go upside down;
        loc = (i%12)*11+int(i/12+1)   #Need to flip this to go upside down;
        ax=plt.subplot(12,11,loc)
        print "...plot loc: ", i, loc
        #ax=plt.subplot(111)

        CSD=np.dot(W,lfp_ave[counter])   # units: mS/mm^3*mV = uA/mm^3
        t=np.arange(-.1,0.1,1./(tsf.SampleFrequency*2))
        
        #x_labels = [-0.1,0.0,0.3]
        #plt.xticks(x_labels, x_labels)
        
        #print CSD.shape
        CSD = CSD[:,22500:27500]    #This sets data from -100ms to +300ms
        plt.imshow(CSD, vmin = -250, vmax = 250,extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto',
                   interpolation='sinc'); 
        plt.plot([0,0],[0,1.43],'r--', linewidth=3, color='black', alpha=0.8)
        plt.ylim(1.44,0)
        #plt.tick_params(axis='both', which='major', labelsize=8)
        ax.get_yaxis().set_visible(False)
        #ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        #plt.title(areas[stim_areas_list[i]],fontsize=15)
        #if (counter==0) or (counter==12):
        #    plt.ylabel("Depth along probe (mm)", fontsize=20)
        #    ax.get_yaxis().set_visible(True)
        #else:
        #    ax.get_yaxis().set_visible(False)
        if i%12==11: 
            x_labels = [-0.1,0.0,0.1]
            plt.xticks(x_labels, x_labels)
        
        counter+=1
        
    plt.suptitle(file_name, fontsize=25)
    plt.show()
    
    for i in [4,16,22]:
        CSD=np.dot(W,lfp_ave[i])   # units: mS/mm^3*mV = uA/mm^3
        t=np.arange(-.1,0.3,1./(tsf.SampleFrequency*2))
        
        x_labels = [-0.1,0.0,0.1,0.2,0.3]
        plt.xticks(x_labels, x_labels)
        
        print CSD.shape
        CSD = CSD[:,22500:32500]
        plt.imshow(CSD, vmin = -250, vmax = 250,extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto',
                   interpolation='sinc'); 
        plt.plot([0,0],[0,1.43],'r--', linewidth=3, color='black', alpha=0.6)
        plt.plot([.1,.1],[0,1.43],'r--', linewidth=3, color='black', alpha=0.25)
        plt.plot([.2,.2],[0,1.43],'r--', linewidth=3, color='black', alpha=0.2)

        plt.ylim(1.44,0)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.title(areas[stim_areas_list[i]],fontsize=15)
        plt.ylabel("Depth along probe (mm)", fontsize=20)
        plt.xlabel("Time from laser pulse (sec)",fontsize=20)

        plt.suptitle(file_name, fontsize=25)
        plt.show()

