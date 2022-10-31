import struct
import numpy as np
import matplotlib.pyplot as plt
from array import array
import time
import sys
import h5py
import scipy.io



#************* Convert Sev's H32 (32 channel) or A64 + Edge32 (99 channel) recordings to .tsf *******

# LOAD .DAT FILE 
print "Loading Traces"

#file_name = '/media/cat/4TB/in_vivo/sev/M71/CSD3/CSD3.dat'
file_name = '/media/cat/4TB/in_vivo/sev/M71/DeepCSD4/DeepCSD4.dat'

dt = np.fromfile(file_name,dtype=np.int16)

dt=dt.reshape(dt.shape[0]/99,99)



#Find length of recording in # of samples
n_vd_samples= dt.shape[0]

###Plot first 10000 steps from channel 0
#plt.plot(dt[0:100000,0], 'r-', color='black',linewidth=1)
#plt.show()


dt=np.array(dt, dtype=np.int16)

print dt[0][0:100]
print dt[1][0:100]
print dt[2][0:100]

print dt[0:100][0]
print dt[0:100][1]
print dt[0:100][2]

quit()

#x = np.arange(0,5000, .05)
#len(x)
#len(dt[0:100000,0])
#plt.plot(x, dt[0:100000,0], 'r-', linewidth=1)
#plt.plot(x, dt[0:100000,1], 'r-', linewidth=1)
#plt.plot(x, dt[0:100000,2], 'r-', linewidth=1)

#plt.show()


#******************************* SELECT A PROBE *********************************
#A64 - 6 shank (10chs each) LGN; Edge32 - linear 32 channels in V1; H32 - older recordings (e.g. M34)

A64 = False
Edge32 = True
H32 = False

#* V1 H32 PROBES **
#Sev's V1 probe data: M34_20132808;
if H32: 
    probe32 = np.array([16,6,5,15,4,7,3,8,2,9,1,10,14,13,12,11,22,21,20,19,23,25,24,18,26,17,27,29,28,31,30,32])
    adapter32 = np.array([23,24,25,26,27,28,29,30,31,32,17,18,21,22,19,20,13,14,11,12,15,16,1,2,3,4,5,6,7,8,9,10])

    for i in range(4):
        SiteLoc[0][i*8]=-24+i*200
        SiteLoc[1][i*8]=0

        SiteLoc[0][i*8+1]=-18+i*200
        SiteLoc[1][i*8+1]=60

        SiteLoc[0][i*8+2]=-12+i*200
        SiteLoc[1][i*8+2]=120

        SiteLoc[0][i*8+3]=-6+i*200
        SiteLoc[1][i*8+3]=180

        SiteLoc[0][i*8+4]=0+i*200
        SiteLoc[1][i*8+4]=210

        SiteLoc[0][i*8+5]=6+i*200
        SiteLoc[1][i*8+5]=150

        SiteLoc[0][i*8+6]=12+i*200
        SiteLoc[1][i*8+6]=90

        SiteLoc[0][i*8+7]=18+i*200
        SiteLoc[1][i*8+7]=30

#******************************** A64 PROBE LAYOUT #******************************** 
if A64: 
    n_electrodes = 6 * 10
    SiteLoc = np.zeros((2,n_electrodes), dtype=np.int16) #Read as 1D array

    for i in range(6):
        SiteLoc[0][i*10]=-54+i*200
        SiteLoc[1][i*10]=0

        SiteLoc[0][i*10+1]=+48+i*200
        SiteLoc[1][i*10+1]=30

        SiteLoc[0][i*10+2]=-42+i*200
        SiteLoc[1][i*10+2]=60

        SiteLoc[0][i*10+3]=+36+i*200
        SiteLoc[1][i*10+3]=90

        SiteLoc[0][i*10+4]=-30+i*200
        SiteLoc[1][i*10+4]=120

        SiteLoc[0][i*10+5]=+24+i*200
        SiteLoc[1][i*10+5]=150

        SiteLoc[0][i*10+6]=-18+i*200
        SiteLoc[1][i*10+6]=180

        SiteLoc[0][i*10+7]=+12+i*200
        SiteLoc[1][i*10+7]=210

        SiteLoc[0][i*10+8]=-6+i*200
        SiteLoc[1][i*10+8]=240

        SiteLoc[0][i*10+9]=0+i*200
        SiteLoc[1][i*10+9]=270

#******************************** EDGE 32 PROBE LAYOUT #******************************** 
if Edge32: 
    n_electrodes = 32
    #Make Edge 32 probe
    SiteLoc = np.zeros((2,n_electrodes), dtype=np.int16) #Read as 1D array
    for i in range(n_electrodes):
        SiteLoc[0][i]=0
        SiteLoc[1][i]=20*i

print "No. channels: ", len(SiteLoc)
print SiteLoc

#******************************** LGN A64 PROBE MAPPING ************************************
#Sev's LGN probe data: M72, M73, M74

probe1 = [11,13,27,15,5,21,19,23,1,31]
probe2 = [22,3,28,26,18,7,24,30,14,9]
probe3 = [8,12,2,0,6,16,10,4,25,20]
probe4 = [50,39,62,35,46,47,58,17,42,29]
probe5 = [61,40,36,34,57,44,32,38,55,48]
probe6 = [51,53,49,37,43,59,41,45,33,63]

LGNprobe = []
LGNprobe.extend(probe1)
LGNprobe.extend(probe2)
LGNprobe.extend(probe3)
LGNprobe.extend(probe4)
LGNprobe.extend(probe5)
LGNprobe.extend(probe6)

print "LGNprobe: ", LGNprobe
#***************************  CRAP CHS ************
probe7 = [52,54,60,56]  #Caca, i.e. poopy

#******************************* V1 EDGE PROBE MAPPING *****************************
probe8 = [83,91,90,82] #V1 probes on Edge 32 probe; 
probe9 = [89,92,88,93] 
probe10 = [87,94,86,95]
probe11 = [85,84,81,80]
probe12 = [79,78,75,74]
probe13 = [64,66,65,77]
probe14 = [67,76,68,70]
probe15 = [69,72,71,73]
V1probe = []
V1probe.extend(probe8)
V1probe.extend(probe9)
V1probe.extend(probe10)
V1probe.extend(probe11)
V1probe.extend(probe12)
V1probe.extend(probe13)
V1probe.extend(probe14)
V1probe.extend(probe15)

V1probe=V1probe[::-1] #Invert probe - i.e. it's upside down
print "V1probe: ", V1probe

#quit()

#************************* STIMULUS CHS, SCREEN REFRESH ETC **********************
probe16 = [96,97,98] #This contains stimulus info, monitor refresh rate and X?
#**************************************************************

#**************SAVE STRAIGHT .DAT FILE FOR COSTAS ****************

#file_name = file_name+'_V1.dat'

#f = open(file_name, 'wb')

#for i in range(n_electrodes):
    #print ""
    #print ""
    #print "Writing ch: ", i+1
    #print "Probe ch: ", V1probe[i]
    #dt[0:n_vd_samples,V1probe[i]].tofile(f)

#quit()

#*********** SAVE .TSF FILES ***********
if Edge32:
    file_name = file_name[:-4]+'_V1.tsf'
    Save_probe = V1probe

if A64: 
    file_name = file_name[:-4]+'_LGN.tsf'
    Save_probe = LGNprobe

header = 'Test spike file '
iformat = 1002
SampleFrequency = 20000
vscale_HP = .25

if Edge32:
    f = open(file_name, 'wb')
    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP))
    for i in range(len(Save_probe)):
        f.write(struct.pack('h', SiteLoc[0][i]))
        f.write(struct.pack('h', SiteLoc[1][i]))
        f.write(struct.pack('i', i+1)) #Use '1' based indices or SS will crash
    print "Writing data"

    for i in range(n_electrodes):
        print ""
        print ""
        print "Writing ch: ", i+1
        print "Probe ch: ", V1probe[i]
        dt[0:n_vd_samples,V1probe[i]].tofile(f)

    f.write(struct.pack('i', 0)) #Write # of fake spikes
    f.close()



#Save A64 probe
if A64:

    file_name = file_name+'_A64_LGN_singleshank.tsf'
    #Save_probe = probe1 

    f = open(file_name, 'wb')
    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP))
    for i in range(len(Save_probe)):
        f.write(struct.pack('h', SiteLoc[0][i]))
        f.write(struct.pack('h', SiteLoc[1][i]))
        f.write(struct.pack('i', i+1)) #Use '1' based indices or SS will crash
    print "Writing data"

    #Sev's V1 data from 32 ch Edge probe; LGN Files: M72-74
    for i in Save_probe:
        print ""
        print ""
        print "Writing Probe ch: ", i
        (dt[:,i]).tofile(f)

    tempint = 0
    f.write(struct.pack('i', tempint)) #Write # of fake spikes
    f.close()

if H32:
    f = open(file_name, 'wb')
    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP))
    for i in range(len(Save_probe)):
        f.write(struct.pack('h', SiteLoc[0][i]))
        f.write(struct.pack('h', SiteLoc[1][i]))
        f.write(struct.pack('i', i+1)) #Use '1' based indices or SS will crash
    print "Writing data"

    #Sev's V1 data on H64 probe
    for i in range (16,n_electrodes+16):
        print ""
        print ""
        print "Writing ch: ", i+1
        print "Probe ch: ", probe32[i]
        ch_index = adapter32[probe32[i]-1]  #np.where(adapter32==(probe32[i]-1))[0][0]
        print "Adapter ch: ", ch_index
        dt[0:n_vd_samples,ch_index-1].tofile(f)

    #Sev's V1 data from 32 ch Edge probe; LGN Files: M72-74
    for i in Save_probe:
        print ""
        print ""
        print "Writing Probe ch: ", i
        dt[:,i].tofile(f)

    tempint = 0
    f.write(struct.pack('i', tempint)) #Write # of fake spikes
    f.close()

#************************** SAVE HDF5 FILE FORMAT *******************************

#quit()

#Create animal specific hdf5 file
exp_name = 'M71'  #Change this value to correct animal
ofname = '/media/cat/4TB/in_vivo/sev/' 
hdf5file = h5py.File(ofname+exp_name+'/'+exp_name+'.hdf5','w',libver='latest')

SampleFrequency=20000.

print "writing hdf5 file: ", ofname
for i in ['CSD2','CSD3','DeepCSD4','DeepCSD5']: #Change list to required list
    
    ##Load stimulus times - SEV DATA
    #file_name = ofname+exp_name+'/'+str(i)+'/'+str(i)+'.SweepTime.mat'
    #mat = scipy.io.loadmat(file_name)
    #bright_array = np.array(mat[mat.keys()[1]][0][0][4][0][1])/SampleFrequency  

    #Load raw .dat files
    file_name = ofname+exp_name+'/'+str(i)+'/'+str(i)+'.dat'
    print "processing: ", file_name
    ecp = np.fromfile(file_name,dtype=np.int16)
    ecp=ecp.reshape(ecp.shape[0]/99,99)
    print ecp[:,Save_probe].T[0][0:100]
    quit()
    
    #Create Recording group
    grp=hdf5file.create_group(str(i))
    grp.attrs['recording name']= 'Severine Durand CSD experiment'
    grp.attrs['date of recording']= 'Sep 26 2014'
    grp.attrs['length of recording (sec)']= str(ecp.shape[0]/SampleFrequency)

    #Raw ecp
    grp.create_dataset('raw_ecp',data=ecp[:,Save_probe].T)
    grp.attrs['raw_ecp']=  'raw extracellular voltages in dtype=np.int16 \
    (scaled by vscale) in array of size = [no. of channels, length of recording]'

    #vscale
    grp.create_dataset('vscale',data=vscale_HP)
    grp.attrs['vscale']= 'scaling factor used to convert raw data into int16 format'

    #probe layout
    grp.create_dataset('probe_layout',data=SiteLoc.T)
    grp.attrs['probe_layout']= 'array of x, y locations for each electrode site from top of probe down (in um)'

    #stim times
    grp.create_dataset('stimulus_times',data=bright_array)
    grp.attrs['stimulus_times']= 'array containing beginnig and end of full-screen-flash (in seconds)'

    #stim triggered ave ecp
    window = 1 #secs of window to average
    print "Computing cumulative LFP for: ", window, " sec."
    lfp_cumulative = np.zeros((SiteLoc.shape[1],window*SampleFrequency), dtype=np.float32)
    for i in range(len(Save_probe)):
        for j in bright_array:
            start=int(j[0]*SampleFrequency)
            end=start+window*SampleFrequency
            #print lfp_cumulative[i].shape
            #print ecp.T[Save_probe[i]][start:end].shape
            lfp_cumulative[i] += ecp.T[Save_probe[i]][start:end]
    lfp_cumulative /= len(bright_array)

    grp.create_dataset('stim_triggered_ecp',data=lfp_cumulative)
    grp.attrs['stim_triggered_ecp']= 'Stimulus triggered average ecp for 1 second after onset of \
    stimulus (summed over no. of stimuli and divided by no. of stimuli; unfiltered)'


hdf5file.close()

