import numpy as np
import matplotlib.pyplot as plt
from array import array
import time
import sys
import random as rnd

from scipy import stats
from scipy.interpolate import interp1d
import scipy.optimize
import scipy.io
import numpy as np
import time, math
import os.path
from pylab import *
import struct, array, csv

file_name = '/home/cat/neuron/src_save_haycells_work/modelsdir/input/lam_electrode_rat.npy'

#IMEC PROBE: 22um spacing 
n_electrodes = 30
electrode_x=[]
electrode_y=[]
electrode_z=[]

for i in range(n_electrodes):
    electrode_x.append(0.)

for i in range(n_electrodes):
    electrode_y.append((i%2)*22.)
    electrode_z.append(1050+int(i/2)*22.)

temp=n_electrodes

###H32 PROBE LAYOUT - 8 Channels
n_electrodes = 8
for i in range(n_electrodes):
    electrode_x.append(0.)

electrode_y.append(-24.)     #Start ~L4 where pyr4s are
electrode_z.append(1100.)
electrode_y.append(18.)
electrode_z.append(1130.)
electrode_y.append(-18.)
electrode_z.append(1160.)
electrode_y.append(12.)
electrode_z.append(1190.)
electrode_y.append(-12.)
electrode_z.append(1220.)
electrode_y.append(6.)
electrode_z.append(1250.)
electrode_y.append(-6.)
electrode_z.append(1280.)
electrode_y.append(0.)
electrode_z.append(1310.)

temp+=n_electrodes

#NEURO-NEXUS ~65um spacing probe
n_electrodes = 5
for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(0.)
    electrode_z.append(1050+i*65.)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(65.)
    electrode_z.append(1050 + i*65.+32.5)
temp+=n_electrodes

#New IMEC probe 4 column layout: 
n_electrodes = 8
for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(0.)
    electrode_z.append(1050 + i*44)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(22.)
    electrode_z.append(1050 + 22 + i*44)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(44.)
    electrode_z.append(1050 + i*44)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(0.)
    electrode_y.append(66.)
    electrode_z.append(1050 + 22+ i*44)
temp+=n_electrodes


#Theoretical 3 shank 22um vertical spacing x 2 rows:
n_electrodes = 15
for i in range(n_electrodes):
    electrode_x.append(-50.)
    electrode_y.append(-100.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(-50.)
    electrode_y.append(0.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(-50.)
    electrode_y.append(+100.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(+50.)
    electrode_y.append(-100.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(+50.)
    electrode_y.append(0.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

for i in range(n_electrodes):
    electrode_x.append(+50.)
    electrode_y.append(+100.)
    electrode_z.append(1050 + i*22)
temp+=n_electrodes

    
e = np.array([electrode_x, electrode_y, electrode_z])

np.save(file_name, e)

print e
print temp
