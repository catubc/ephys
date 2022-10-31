from tsf_ptcs_classes import *
from distributions import *

from scipy import stats
from scipy.interpolate import interp1d
import scipy.optimize
import scipy.io
import numpy as np
import time, math
import sys
import os.path
from pylab import *
import struct, array, csv
import pandas as pd
import matplotlib.mlab as mlab
import itertools
import cPickle as pkl
#import pickle

def Plot_fullraster(LGN_array, stim_array1, stim_array2):

    #*********************** VLINE SETS OF PLOTS *****************************
    ax = plt.subplot(1, 1, 1)
    title("DAN - LGN Rasters & Stimulus + PETH",fontsize=15)

    colors = ['blue','red', 'black']

    for i in range(len(LGN_array)):

        x = np.array(LGN_array[i],dtype=float32)

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=i
        ymax+=i-.8

        #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
        plt.vlines(x, ymin, ymax, linewidth=1, color='black')

    xx = list(itertools.chain(*LGN_array))
    
    yy = np.histogram(xx, bins = np.arange(0,260,.010))
    plt.plot(yy[1][:-1], yy[0])

    #Plot bright stimulus: stim_array1
    for i in range(min(1000, len(stim_array1))): #len(stim_array)):
        plt.fill_between([stim_array1[i][0],stim_array1[i][1]], 0, 180, color='red', alpha=0.25)

    #Plot dark stimulus: stim_array2
    #for i in range(min(1000, len(stim_array2))): #len(stim_array)):
        #plt.fill_between([stim_array2[i][0],stim_array2[i][1]], 0, 180, color='black', alpha=0.25)


    plt.ylim(0,len(LGN_array)+2)
    #plt.xlim(0,260)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Cell Index (not ordered)',multialignment='center', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def Plot_3secraster(save_array_bright, save_array_dark, file_name):

    #*********************** PLOT BRIGHT FLASH RESPONSE *****************************
    ax = plt.subplot(1, 1, 1)
    title("1.0sec Responses to Bright Flashes (left) and Dark Flashes (right) (40 trials each summed) & PSTH\n"
    +"Recording:  "+file_name,fontsize=15)

    save_array = save_array_bright
    colors = ['blue','red', 'black']

    scale=1
    scale1=0.025

    for i in range(len(save_array)):

        x = np.array(save_array[i],dtype=float32)
        x = np.sort(x)

        #Take only up to 1.0sec data;
        if len(x)>0:
            for j in range(len(x)):
                if x[j]>1.0:
                    break
            x = x[0:j-1]
        
        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=i*scale
        ymax+=i*scale+.1

        #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
        plt.vlines(x, ymin, ymax, linewidth=1, color='black')

        if max(save_array[i])>0:
            print max(save_array[i])
            yy = np.histogram(save_array[i], bins = np.arange(0,1.0,.010))
            plt.plot(yy[1][:-1], yy[0]*scale1+i*scale, linewidth=1, color='black', alpha=0.750)

    #xx = list(itertools.chain(*LGN_array))
    
    #Plot bright stimulus: stim_array1
    p = ax.axvspan(0.0, 0.05, facecolor='red', alpha=0.25)
    p = ax.axvspan(0.05, 0.10, facecolor='blue', alpha=0.25)
    p = ax.axvspan(0.10, 0.15, facecolor='green', alpha=0.25)
    p = ax.axvspan(0.15, 0.20, facecolor='magenta', alpha=0.25)
    #p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)


    plt.plot([1.,1.],[0,len(save_array)*scale+2],linewidth=5,color='black',alpha=.9)
    ##Plot dark stimulus: stim_array2
    #for i in range(min(1000, len(stim_array2))): #len(stim_array)):
        #plt.fill_between([stim_array2[i][0],stim_array2[i][1]], 0, 180, color='black', alpha=0.25)

    plt.ylim(0,len(save_array)*scale)
    plt.xlim(0,2.0)


    #*********************** PLOT DARk FLASH RESPONSE *****************************
    #ax = plt.subplot(1, 2, 2)
    #title("3 Seconds rasters + PETH",fontsize=15)

    save_array = save_array_dark

    scale=1
    xshift = 1.0
    for i in range(len(save_array)):

        x = np.array(save_array[i],dtype=float32)
        x = np.sort(x)

        #Take only up to 1.0sec data;
        if len(x)>0:
            for j in range(len(x)):
                if x[j]>1.0:
                    break
            x = x[0:j-1]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=i*scale
        ymax+=i*scale+.1

        #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
        plt.vlines(x+xshift, ymin, ymax, linewidth=1, color='black')

        if max(save_array[i])>0:
            print max(save_array[i])
            yy = np.histogram(save_array[i], bins = np.arange(0,1.0,.010))
            plt.plot(yy[1][:-1]+xshift, yy[0]*scale1+i*scale,linewidth=1, color='black', alpha=0.750)
    #xx = list(itertools.chain(*LGN_array))
    
    #Plot bright stimulus: stim_array1
    p = ax.axvspan(xshift+0.0,xshift+ 0.05, facecolor='black', alpha=0.25)
    p = ax.axvspan(xshift+0.05, xshift+0.10, facecolor='blue', alpha=0.25)
    p = ax.axvspan(xshift+0.10, xshift+0.15, facecolor='green', alpha=0.25)
    p = ax.axvspan(xshift+0.15, xshift+0.20, facecolor='magenta', alpha=0.25)
    #p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)

    ###Plot dark stimulus: stim_array2
    ##for i in range(min(1000, len(stim_array2))): #len(stim_array)):
        ##plt.fill_between([stim_array2[i][0],stim_array2[i][1]], 0, 180, color='black', alpha=0.25)

    #plt.ylim(0,len(save_array)*scale+2)
    #plt.xlim(0,1.0)

    #plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Cell Index (no particular order)',multialignment='center', fontsize=15)

    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    plt.show()

def Plot_3secraster_brightonly(save_array_bright, save_array_dark, file_name, no_flashes):

    #*********************** PLOT BRIGHT FLASH RESPONSE *****************************
    ax = plt.subplot(1, 2, 1)
    title("1.0sec Responses to Bright Flashes (left) and Dark Flashes (right) (40 trials each summed) & PSTH\n"
    +"Recording:  "+file_name,fontsize=15)

    #save_array_bright contains the pooled responses for each cell across all stimulus trials
    #sev uses 3 second stimulus windows (?); Dan uses 1 sec but alternates bright to dark flashes
    save_array = save_array_bright
    colors = ['blue','red', 'black']

    scale=1
    bin_width = 0.010 #10ms bin width
    scale1=.1/float(no_flashes)/float(bin_width)

    
    #Looping over the # of cells 
    for i in range(len(save_array)):
        x = np.array(save_array[i],dtype=float32)
        x = np.sort(x)

        #Take only up to 1.0sec data;
        if len(x)>0:
            for j in range(len(x)):
                if x[j]>1.0:
                    break
            x = x[0:j-1]
        
        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=i*scale
        ymax+=i*scale+.1

        #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
        plt.vlines(x, ymin, ymax, linewidth=1, color='black')

        if save_array[i]:           #Checks if list is empty
            print max(save_array[i])
            yy = np.histogram(save_array[i], bins = np.arange(0,1.0,bin_width))
            plt.plot(yy[1][:-1], yy[0]*scale1+i*scale, linewidth=1, color='black', alpha=0.750)

    #xx = list(itertools.chain(*LGN_array))
    
    #Plot bright stimulus: stim_array1
    p = ax.axvspan(0.0, 0.05, facecolor='red', alpha=0.25)
    p = ax.axvspan(0.05, 0.10, facecolor='blue', alpha=0.25)
    p = ax.axvspan(0.10, 0.15, facecolor='green', alpha=0.25)
    p = ax.axvspan(0.15, 0.20, facecolor='magenta', alpha=0.25)
    #p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)

    plt.plot([1.,1.],[0,len(save_array)*scale+2],linewidth=5,color='black',alpha=.9)

    plt.ylim(0,len(save_array)*scale)
    plt.xlim(0,0.5)

    plt.ylabel('Cell Index (no particular order)',multialignment='center', fontsize=15)
    plt.xlabel('Time post stimulus: seconds')
    #ax2 = ax.twinx()

    #label_y = np.arange(0,100*len(save_array),50)
    
    #my_yticks1 = ['-', '-']
    #my_yticks= my_yticks1
    #for i in range(len(save_array)):
        #my_yticks = np.concatenate((my_yticks,my_yticks1))

    #plt.yticks(label_y, my_yticks)

    #ax2.set_ylabel('Average Firing Rate (PSTH)', fontsize=16)
    #ax2.set_ylim(0, 100*len(save_array))

    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    plt.show()

#************************************** READ DAN DATA ******************************************
def Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array):
    f = open(sim_dir+file_name,'r')
    object_file=pkl.load(f)

    array_keys= object_file.keys()
    #print array_keys

    for i in range(len(array_keys)):
        if ((array_keys[i]<> 'bright') or (array_keys[i]<>'dark')):
            LGN_array.append(np.array(object_file[array_keys[i]])/25000.)

    #print LGN_array[0]

    temp_array =  np.array(object_file['bright'])/25000.
    for i in range(len(temp_array)):
        bright_array.append([temp_array[i],temp_array[i]+0.050])

    temp_array = np.array(object_file['dark'])/25000.
    for i in range(len(temp_array)):
        dark_array.append([temp_array[i],temp_array[i]+0.050])

    #print bright_array
    #print dark_array

#***************************LOAD PTCS FILE  **************************

np.random.seed(12345678)  #fix random seed to get the same result
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

##Load spike times from .ptcs file
#sorted_file = '1CDS1.dat_LGN'
#sim_dir = '/home/cat/neuron/in_vivo/dan/M150566/2015-04-01_16-41-14_flashes_vis/'
sim_dir = '/media/cat/4TB/in_vivo/sev/M74/'

ptcs_file = 'CSD1.dat_LGN'
sorted_file = ptcs_file
print "Loading: ", ptcs_file
fname = sim_dir + ptcs_file + '.ptcs'
f = open(fname, "rb")
work_dir = sim_dir+sorted_file+'/'
Sort1 = Loadptcs(sorted_file, sim_dir)
Sort1.name=sorted_file
Sort1.filename=sorted_file
Sort1.directory=work_dir
if sorted_file == 'truth':
    Sort1.flag=0 #This flag indicates the sort is blind and there is no ground truth
else:
    Sort1.flag=1
f.close()

LGN_array=[]
for i in range(Sort1.n_units):
    LGN_array.append(np.array(Sort1.units[i])/Sort1.samplerate)

#***********  SEV'S STIMULUS TIMES ************
#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M72/'
#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M73/'
##sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M74/'

####READ SPIKE TIME .mat DATA

#file_name = 'SpT_in_sec_M72_CSD1.mat'
#file_name = 'SpT_in_sec_CSD_M73.mat'
##file_name = 'SpT_in_sec_M74_CSD1.mat'

##Load Sev's sorted data spike times
##mat = scipy.io.loadmat(sim_dir+file_name)
##indexes = mat.keys()
##array = np.array(mat[indexes[1]])
##LGN_array=[]
##for i in range(array.shape[1]):
    ##LGN_array.append(np.hstack(array[0][i]))

#********SEV'S STIMULUS FILES FROM .MAT FORMAT ************
file_name = 'CSD1.SweepTime.mat'
mat = scipy.io.loadmat(sim_dir+file_name)

indexes = mat.keys()

array = np.array(mat[indexes[1]][0][0][4][0][1]) #Cat: many wrappers on these data structures 
bright_array= array/20000. #Convert to seconds; DON"T HARDWIRE THIS


no_flashes = len(bright_array)
file_name_sum = file_name

save_array=[]
test_array=[]
lgn_array=[]
for j in range(6):
    test_array.append([])
save_array1={} #Make lists to capture 1second of expeirmental data across all units
save_array2={} #same as above

for i in range(len(LGN_array)):
    save_array.append([])
    for j in range(6):
        test_array[j].append([])
#    save_array1.append([])
#    save_array2.append([])

#Loop over all LGN cells; LGN_array contains complete recording 
for i in range(len(LGN_array)):
    temp_array = np.array(LGN_array[i],dtype=np.float32)
    #Loop over all bright flashes and grab response 3 seconds out from screen flash start
    for j in range(len(bright_array)):
        #Search for spikes within 1.0 sec of stimulus start time.
        temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+1.0))[0]
        save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

        #Save raster to a total LGN raster list
        lgn_array.append(temp_array[temp2]-bright_array[j][0])

        #This extracts the LGN responses to the 2nd and 3rd stimulus;
        if j == 1:
            save_array1[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[0][i].extend(temp_array[temp2]-bright_array[j][0])
        if j == 2:
            save_array2[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[1][i].extend(temp_array[temp2]-bright_array[j][0])
        if j == 3:
            #save_array2[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[2][i].extend(temp_array[temp2]-bright_array[j][0])
        if j == 55:
            #save_array2[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[3][i].extend(temp_array[temp2]-bright_array[j][0])
        if j == 95:
            #save_array2[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[4][i].extend(temp_array[temp2]-bright_array[j][0])
        if j == 130:
            #save_array2[str(i)]=(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 
            test_array[5][i].extend(temp_array[temp2]-bright_array[j][0])

temp_array = save_array #Save copy of array as the plotting routines below overwrite it

print len(lgn_array)
print len(bright_array)
lgn_array= np.array(lgn_array)
np.savez(sim_dir+'lgn_array_M74',lgn_array)

quit()

#Plot responses to individual stimuli [names above];
if False:
    names = [1,2,3,55,95,130]
    #Looping over the # of cells
    scale=1
    scale1=1
    bin_width = 0.010 #10ms bin width
    
    for j in range(6):
        ax = plt.subplot(2, 3, j+1)
        save_array = test_array[j]
        plt.title('Stimulus #: ' + str(names[j]))
        #yy_sum saves cumulative rasters
        yy_sum = np.arange(0,1.0-bin_width,bin_width)*0.
        for i in range(len(save_array)):
            x = np.array(save_array[i],dtype=float32)
            x = np.sort(x)
    
            #Take only up to 1.0sec data;
            if len(x)>0:
                for j in range(len(x)):
                    if x[j]>1.0:
                        break
                x = x[0:j-1]
            
            ymin=np.zeros(len(x))
            ymax=np.zeros(len(x))
            ymin+=i*scale
            ymax+=i*scale+.1
    
            #if Sort1.chanpos[Sort1.maxchan[i]][0]==0:
            plt.vlines(x, ymin, ymax, linewidth=1, color='black')
    
            if save_array[i]:           #Checks if list is empty
                print max(save_array[i])
                yy = np.histogram(save_array[i], bins = np.arange(0,1.0,bin_width))
                plt.plot(yy[1][:-1], yy[0]*scale1+i*scale, linewidth=1, color='black', alpha=0.750)
                #print yy[0]
                yy_sum += yy[0]*scale1
                xx = yy[1][:-1]
                #print xx
                
        plt.plot(xx, yy_sum, linewidth='2', color='red')
        
        plt.ylim(0,150)
    plt.show()

#np.savez(sim_dir+file_name+'_first.csv',**save_array1)
#np.savez(sim_dir+file_name+'_second.csv',**save_array2)

#quit()
save_array_bright = temp_array

##*************** DAN'S STIMULUS TIMES ***************

#file_name = '/home/cat/neuron/in_vivo/dan/M150566/2015-04-01_16-41-14_flashes_vis/stim_times.txt'

#bright_array = np.genfromtxt(file_name, dtype='float32') #, delimiter=",")
#bright_array = bright_array[::2]

#no_flashes = len(bright_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+1.0))[0] #Search forward + 1.0sec
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#save_array_bright = save_array



#**************************************************************************************


Plot_3secraster_brightonly(save_array_bright, save_array_bright, file_name, no_flashes)

quit()

#**************************************************

#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M72/'
sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M73/'
#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M74/'

###READ SPIKE TIME .mat DATA

#file_name = 'SpT_in_sec_M72_CSD1.mat'
file_name = 'SpT_in_sec_CSD_M73.mat'
#file_name = 'SpT_in_sec_M74_CSD1.mat'

mat = scipy.io.loadmat(sim_dir+file_name)
indexes = mat.keys()
array = np.array(mat[indexes[1]])

LGN_array=[]
for i in range(array.shape[1]):
    LGN_array.append(np.hstack(array[0][i]))

#READ SWEEP TIME .mat DATA
file_name = '1CSD1.SweepTime.mat'
mat = scipy.io.loadmat(sim_dir+file_name)

indexes = mat.keys()

array = np.array(mat[indexes[1]][0][0][4][0][1]) #Cat: many wrappers on these data structures 
bright_array= array/20000. #Convert to seconds

file_name_sum += file_name

save_array=[]
for i in range(len(LGN_array)):
    save_array.append([])

#Loop over all LGN cells
for i in range(len(LGN_array)):
    temp_array = np.array(LGN_array[i],dtype=np.float32)
    #Loop over all bright flashes and grab response 3 seconds out from screen flash start
    for j in range(len(bright_array)):
        temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

for i in range(len(save_array)):
    save_array_bright.append(save_array[i])

file_name_sum += ', '+ file_name


#**************************************************

sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M72/'
#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M73/'
#sim_dir = '/media/cat/Data1/in_vivo/sev/LGN/M74/'

###READ SPIKE TIME .mat DATA

file_name = 'SpT_in_sec_M72_CSD1.mat'
#file_name = 'SpT_in_sec_CSD_M73.mat'
#file_name = 'SpT_in_sec_M74_CSD1.mat'

mat = scipy.io.loadmat(sim_dir+file_name)
indexes = mat.keys()
array = np.array(mat[indexes[1]])

LGN_array=[]
for i in range(array.shape[1]):
    LGN_array.append(np.hstack(array[0][i]))

#READ SWEEP TIME .mat DATA
file_name = 'CSD1.SweepTime.mat'
mat = scipy.io.loadmat(sim_dir+file_name)

indexes = mat.keys()

array = np.array(mat[indexes[1]][0][0][4][0][1]) #Cat: many wrappers on these data structures 
bright_array= array/20000. #Convert to seconds

file_name_sum += file_name

save_array=[]
for i in range(len(LGN_array)):
    save_array.append([])

#Loop over all LGN cells
for i in range(len(LGN_array)):
    temp_array = np.array(LGN_array[i],dtype=np.float32)
    #Loop over all bright flashes and grab response 3 seconds out from screen flash start
    for j in range(len(bright_array)):
        temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

for i in range(len(save_array)):
    save_array_bright.append(save_array[i])

file_name_sum += ', '+ file_name

Plot_3secraster_brightonly(save_array_bright, save_array_bright, file_name_sum, no_flashes)


quit()

##*********************** READ DAN DATA *****************
#sim_dir = '/media/cat/Data1/in_vivo/dan/LGN/'
#file_name = '149873_lgnSpikes.pkl'
#file_name_sum = file_name
#LGN_array=[]
#bright_array = []
#dark_array=[]

#Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#save_array_bright = save_array

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all dark flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(dark_array)):
        #temp2 = np.where(np.logical_and(temp_array>=dark_array[j][0], temp_array<=dark_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-dark_array[j][0]) #Offset time of spike to t_0 

#save_array_dark = save_array

##******************************** LOAD ANOHTER FILE *************

#sim_dir = '/media/cat/Data1/in_vivo/dan/LGN/'
#file_name = '177366_lgnSpikes.pkl'
#LGN_array=[]
#bright_array = []
#dark_array=[]

#Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_bright.append(save_array[i])

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all dark flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(dark_array)):
        #temp2 = np.where(np.logical_and(temp_array>=dark_array[j][0], temp_array<=dark_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-dark_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_dark.append(save_array[i])
#file_name_sum += ', '+ file_name

##******************************** LOAD ANOHTER FILE *************

#sim_dir = '/media/cat/Data1/in_vivo/dan/LGN/'
#file_name = '179401_lgnSpikes.pkl'
#LGN_array=[]
#bright_array = []
#dark_array=[]

#Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_bright.append(save_array[i])

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all dark flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(dark_array)):
        #temp2 = np.where(np.logical_and(temp_array>=dark_array[j][0], temp_array<=dark_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-dark_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_dark.append(save_array[i])
#file_name_sum += ', '+ file_name

##******************************** LOAD ANOHTER FILE *************

#sim_dir = '/media/cat/Data1/in_vivo/dan/LGN/'
#file_name = '179402_lgnSpikes.pkl'
#LGN_array=[]
#bright_array = []
#dark_array=[]

#Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_bright.append(save_array[i])

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all dark flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(dark_array)):
        #temp2 = np.where(np.logical_and(temp_array>=dark_array[j][0], temp_array<=dark_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-dark_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_dark.append(save_array[i])
#file_name_sum += ', '+ file_name

##******************************** LOAD ANOHTER FILE *************

#sim_dir = '/media/cat/Data1/in_vivo/dan/LGN/'
#file_name = '181420_lgnSpikes.pkl'
#LGN_array=[]
#bright_array = []
#dark_array=[]

#Read_dan(sim_dir, file_name, LGN_array, bright_array, dark_array)

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all bright flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(bright_array)):
        #temp2 = np.where(np.logical_and(temp_array>=bright_array[j][0], temp_array<=bright_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-bright_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_bright.append(save_array[i])

#save_array=[]
#for i in range(len(LGN_array)):
    #save_array.append([])

##Loop over all LGN cells
#for i in range(len(LGN_array)):
    #temp_array = np.array(LGN_array[i],dtype=np.float32)
    ##Loop over all dark flashes and grab response 3 seconds out from screen flash start
    #for j in range(len(dark_array)):
        #temp2 = np.where(np.logical_and(temp_array>=dark_array[j][0], temp_array<=dark_array[j][0]+3.0))[0]
        #save_array[i].extend(temp_array[temp2]-dark_array[j][0]) #Offset time of spike to t_0 

#for i in range(len(save_array)):
    #save_array_dark.append(save_array[i])

#file_name_sum += ', '+ file_name

#Plot_3secraster(save_array_bright, save_array_dark, file_name_sum)


##***************************** PLOTTING SUBROUTINES ******************************

#Plot_fullraster(LGN_array, bright_array, dark_array)

#Plot_3secraster(save_array_bright, save_array_dark)
