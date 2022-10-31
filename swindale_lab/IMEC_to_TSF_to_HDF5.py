import struct
import numpy as np
import matplotlib.pyplot as plt
from array import array
import time
import sys, csv
import cPickle as pkl
import marshal
import json
import OpenEphys as oe
import h5py
import scipy.io

#************ SET FILE PARAMETERS ***************
header = 'Test spike file '
iformat = 1002
n_electrodes = 128
SampleFrequency = 25000
vscale_HP = .2
SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array

for i in range (n_electrodes):
    SiteLoc[i][0]=i%2 * 22
    SiteLoc[i][1]=i/2*22

print "No. channels: ", n_electrodes
#print SiteLoc

#IMEC PROBE LAYOUT:
left_shank=[ 64,  62,  60,  58,  56,  54,  52,  50,  48,  46,  44,  42,  40,
         38,  36,  34,  32,  30,  28,  26,  24,  22,  20,  18,  16,  14,
         12,  10,   8,   6,   4,   2, 128, 126, 124, 122, 120, 118, 116,
        114, 112, 110, 108, 106, 104, 102, 100,  98,  96,  94,  92,  90,
         88,  86,  84,  82,  80,  78,  76,  74,  72,  70,  68,  66]

right_shank=[ 65,  67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89,
         91,  93,  95,  97,  99, 101, 103, 105, 107, 109, 111, 113, 115,
        117, 119, 123, 121, 127, 125,   1,   3,   5,   7,   9,  11,  13,
         15,  19,  17,  23,  21,  27,  25,  31,  29,  35,  33,  39,  37,
         43,  41,  47,  45,  51,  49,  55,  53,  59,  57,  63,  61]

#interweave the channel locations;
probe128 = np.vstack((left_shank,right_shank)).ravel([-1])
#probe128 = np.arange(1,129)
print probe128

#LOAD directories and files
file_dir = '/media/cat/4TB/in_vivo/dan/'

mouse_dir = 'M150566'
exp_dir = '2015-04-01_16-41-14_flashes_vis'

#file_dir = 'M149873/2015-02-06_13-02-32_flash_CSD/'

#file_dir = 'M150566/2015-04-01_16-41-14_flashes_vis/'

#file_dir = 'M150566/2015-04-01_20-55-32_flashes_vis/'

#file_dir = 'M181420/181420_2015-04-07_13-25-07_flashes_vis/'

#file_dir = 'M181420/181420_2015-04-07_18-10-46_flash_vis/'


##********************** LOAD STIMULUS TIMES AND SAVE THEM *****************************
if False:
    file_name = '103_PAI7.StimTimeshere.continuous'
    data=[]

    data = oe.loadContinuous(file_dir+mouse_dir+'/'+exp+dir+'/'+file_name)
    dt = np.array(data['data'])

    stim_array=[]
    state=0
    for i in range(len(dt)):
        if state==0:
            if dt[i]>1:
                print "ON"
                temp=i
                state=1
        if state==1:
            if dt[i]<1:
                print "OFF"
                stim_array.append([temp,i])
                state=0

    del stim_array[0]
    st= np.array(stim_array)/float(SampleFrequency)

    file_name = file_dir+mouse_dir+'/'+exp+dir+'/'+'stim_times.txt'

    np.savetxt(file_name,st, fmt='%8.3f')

#******************************* LOAD RAW DATA AND PREP FOR .TSF FILE ***************************
if True:
    data=[]
    if 'M149873' in mouse_dir: 
        for i in range(128):
            file_name = file_dir+mouse_dir+'/'+exp_dir+'/' + '100_CH'+ str(i+1) +'.continuous'
            data.append(oe.loadContinuous(file_name)['data'])
    else:
        for i in range(128):
            file_name = file_dir+mouse_dir+'/'+exp_dir+'/' + '100_HS_CH'+ str(i) +'.continuous'
            data.append(oe.loadContinuous(file_name)['data'])

    print "Convert to numpy"
    dt = np.array(data, dtype=np.float32)
    print "Convert to numpy - int16"
    dt = np.array(dt/vscale_HP, dtype=np.int16)

    print dt.shape
    ####Plot first 10000 steps from channel 0
    #x = np.arange(0,len(dt1),1.)/SampleFrequency
    #plt.plot(x,dt1, 'r-', color='black',linewidth=1)
    #plt.show()

    n_electrodes = dt.shape[0]
    n_vd_samples= dt.shape[1]
    print "n_electrodes: ", n_electrodes 
    print "n_vd_samples: ", n_vd_samples

    ##***************************************** SAVE .TSF FILE *************************************
    file_name = file_dir+mouse_dir+'/'+exp_dir+'/' + 'V1.tsf'
    dat_name = file_dir+mouse_dir+'/'+exp_dir+'/' + 'V1.dat'

    f = open(file_name, 'wb')
    f_dat = open(dat_name, 'wb')

    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP))
    for i in range (n_electrodes):
        #Siteloc[i][0].tofile(f)
        #Siteloc[i][1].tofile(f)
        f.write(struct.pack('h', SiteLoc[i][0]))
        f.write(struct.pack('h', SiteLoc[i][1]))
        f.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

    #IMEC data
    print "Writing data"
    for i in probe128:
        print "Writing ch: ", i, "Probe ch: ", i
        dt[i-1].tofile(f)
        dt[i-1].tofile(f_dat)
        
    #for i in range(n_electrodes):
        #print "Writing ch: ", i, "Probe ch: ", probe128[i]
        #dt[probe128[i]-1].tofile(f)
        
    tempint = 0
    f.write(struct.pack('i', tempint)) #Write # of fake spikes
    f.close()
    f_dat.close()

quit()
#************************** SAVE HDF5 FILE FORMAT *******************************

#Create animal specific hdf5 file
exp_name = 'V1'  #Change this value to correct animal
ofname = file_dir+mouse_dir+'/'+exp_dir+'/' #'/media/cat/4TB/in_vivo/dan/M149873/2015-02-06_13-02-32_flash_CSD/' 

hdf5file = h5py.File(ofname+mouse_dir+'_'+exp_dir+'.hdf5','w',libver='latest')

#SampleFrequency=25000.

print "writing hdf5 file: ", ofname
    
#Load stimulus times
file_name = ofname+'stim_times.txt' #exp_name+'/'+str(i)+'/'+str(i)+'.SweepTime.mat'
temp_array=[]
with open(file_name, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        temp_array.append([float(str.split(row[0])[0]),float(str.split(row[0])[1])])

#Split stimulus times into bright and dark arrays
bright_array=[]
dark_array=[]
for i in range(len(temp_array)/2):
    bright_array.append(temp_array[i*2])
    dark_array.append(temp_array[i*2+1])
    
bright_array=np.array(bright_array,dtype=np.float32)
dark_array=np.array(dark_array,dtype=np.float32)

#Load raw .dat files #CHANGE THIS
file_name = ofname+'V1.dat'
print "processing: ", file_name
ecp = np.fromfile(file_name,dtype=np.int16)
ecp=ecp.reshape(128,ecp.shape[0]/128)
print ecp.shape


#Create Recording group
grp=hdf5file.create_group(mouse_dir+'_'+exp_dir)
grp.attrs['recording name']= 'Dan Denman CSD experiment ' + ofname+mouse_dir+'_'+exp_dir
grp.attrs['date of recording']= 'Feb 6 2015'
grp.attrs['length of recording (sec)']= str(ecp.shape[0]/SampleFrequency)

#Save Raw ecp
#x = np.arange(0,len(ecp.T[0][0:100000]),1.)/SampleFrequency
#plt.plot(x,ecp.T[0][100000:200000]*vscale_HP, 'r-', color='red',linewidth=1)
#plt.plot(x,ecp.T[1][100000:200000]*vscale_HP+1000, 'r-', color='blue',linewidth=1)
#plt.plot(x,ecp.T[100][100000:200000]*vscale_HP+2000, 'r-', color='green',linewidth=1)
#plt.show()
#quit()

grp.create_dataset('raw_ecp',data=ecp.T)
grp.attrs['raw_ecp']=  'raw extracellular voltages (dtype=np.int16) \
in array of size = [no. of channels, length of recording], (NB: the values need to be multiplied by vscale).'

#Save time steps for Sergey
x = np.arange(0.,1.,1./25000.)
grp.create_dataset('time_points',data=x)
grp.attrs['time_points']=  'corresponding time steps in seconds (dtype=np.float32) for 1sec (total) saved traces'

#Save time steps for Sergey
x = np.zeros(128,dtype=np.int32)+1
temp_array=np.array([16,31,46,63,80,95,110,126,127],dtype=np.int32)
x[temp_array-1]=0
print x
grp.create_dataset('channel_quality',data=x)
grp.attrs['channel_quality']=  'channel quality: 1=good channel, 0=bad channel (dtype=np.int32)'


#vscale
grp.create_dataset('vscale',data=vscale_HP)
grp.attrs['vscale']= 'scaling factor needed to convert (i.e. multiply) raw voltages (int16) to uV units'

#probe layout
grp.create_dataset('probe_layout',data=SiteLoc)
grp.attrs['probe_layout']= 'array of x, y locations for each electrode site from top of probe down (in um)'

#stim times - bright flashes
grp.create_dataset('stimulus_times_bright',data=bright_array)
grp.attrs['stimulus_times_bright']= 'array containing begining and end of full-screen-flash (in seconds) for \
bright flashes'

#stim times - bright flashes
grp.create_dataset('stimulus_times_dark',data=dark_array)
grp.attrs['stimulus_times_dark']= 'array containing begining and end of full-screen-flash (in seconds) for \
dark flashes'


#stim triggered ave ecp - BRIGHT ARRAY
window = 1 #secs of window to average
print "Computing cumulative LFP for: ", window, " sec."
lfp_cumulative = np.zeros((n_electrodes,window*SampleFrequency), dtype=np.float32)
for i in range(n_electrodes):
    for j in bright_array:
        start=int(j[0]*SampleFrequency)
        end=start+window*SampleFrequency
        #print lfp_cumulative[i].shape
        #print ecp.T[Save_probe[i]][start:end].shape
        lfp_cumulative[i] += ecp[i][start:end]
lfp_cumulative /= len(bright_array)

grp.create_dataset('stim_triggered_ecp_bright',data=lfp_cumulative)
grp.attrs['stim_triggered_ecp_bright']= 'Stimulus triggered average ecp for 1 second after onset of \
bright stimulus (summed over no. of stimuli and divided by no. of stimuli; unfiltered)'

#stim triggered ave ecp - DARK ARRAY
window = 1 #secs of window to average
print "Computing cumulative LFP for: ", window, " sec."
lfp_cumulative = np.zeros((n_electrodes,window*SampleFrequency), dtype=np.float32)
for i in range(n_electrodes):
    for j in dark_array:
        start=int(j[0]*SampleFrequency)
        end=start+window*SampleFrequency
        #print lfp_cumulative[i].shape
        #print ecp.T[Save_probe[i]][start:end].shape
        lfp_cumulative[i] += ecp[i][start:end]
lfp_cumulative /= len(dark_array)

grp.create_dataset('stim_triggered_ecp_dark',data=lfp_cumulative)
grp.attrs['stim_triggered_ecp_dark']= 'Stimulus triggered average ecp for 1 second after onset of \
bright stimulus (summed over no. of stimuli and divided by no. of stimuli; unfiltered)'


hdf5file.close()



