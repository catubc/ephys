import struct
import numpy as np
import matplotlib.pyplot as plt
from array import array
import time
import sys
import cPickle as pkl
import marshal
import json
import OpenEphys as oe
import h5py

#************ SET FILE PARAMETERS ***************
header = 'Test spike file '
iformat = 1002
n_electrodes = 32
SampleFrequency = 32000
vscale_HP = 1.
SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array


##IMEC PROBE LAYOUT:
#for i in range (n_electrodes):
#    SiteLoc[i][0]=i%2 * 22
#    SiteLoc[i][1]=i/2*22
#left_shank=[ 64,  62,  60,  58,  56,  54,  52,  50,  48,  46,  44,  42,  40,
         #38,  36,  34,  32,  30,  28,  26,  24,  22,  20,  18,  16,  14,
         #12,  10,   8,   6,   4,   2, 128, 126, 124, 122, 120, 118, 116,
        #114, 112, 110, 108, 106, 104, 102, 100,  98,  96,  94,  92,  90,
         #88,  86,  84,  82,  80,  78,  76,  74,  72,  70,  68,  66]
#right_shank=[ 65,  67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89,
         #91,  93,  95,  97,  99, 101, 103, 105, 107, 109, 111, 113, 115,
        #117, 119, 123, 121, 127, 125,   1,   3,   5,   7,   9,  11,  13,
         #15,  19,  17,  23,  21,  27,  25,  31,  29,  35,  33,  39,  37,
         #43,  41,  47,  45,  51,  49,  55,  53,  59,  57,  63,  61]
##interweave the channel locations;
#probe128 = np.vstack((left_shank,right_shank)).ravel([-1])

probe32 = np.linspace(0,32,32)

#**************************** LOAD RAW DATA AND PREP FOR .TSF FILE ***************************

#file_name = '/media/cat/4TB/in_vivo/peter/lp.h5/lp.h5'

file_name = '/media/cat/4TB/in_silico/einevoll/ViSAPy_nonfiltered.h5'
f=h5py.File(file_name, 'r')
print f.keys()


#dt=np.array(f[f.keys()[0]]).T.copy()*f[f.keys()[2]][0]
dt=np.array(f[f.keys()[0]]).T.copy()

####Plot first 10000 steps from channel 0
#dt = dt[0][1000000:1300000]
#x = np.arange(0,len(dt),1.)/SampleFrequency
#print len(x)
#plt.plot(x,dt)
#plt.show()

dt = np.array(dt*1000, dtype=np.int16)
n_electrodes = dt.shape[0]
n_vd_samples= dt.shape[1]
print "n_electrodes: ", n_electrodes 
print "n_vd_samples: ", n_vd_samples

for i in range (n_electrodes):
    SiteLoc[i][0]=0
    SiteLoc[i][1]=i*22

#***************************************** SAVE .TSF FILE *************************************
file_name = file_name+'.tsf'

f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples))
f.write(struct.pack('f', vscale_HP))
for i in range (n_electrodes):
    f.write(struct.pack('h', SiteLoc[i][0]))
    f.write(struct.pack('h', SiteLoc[i][1]))
    f.write(struct.pack('i', i+1))

#IMEC data
print "Writing data"
for i in range(n_electrodes):
    print "Writing ch: ", i, "Probe ch: ", probe32[i]
    dt[i].tofile(f)

#dt.tofile(f)


tempint = 0
f.write(struct.pack('i', tempint)) #Write # of fake spikes
f.close()
