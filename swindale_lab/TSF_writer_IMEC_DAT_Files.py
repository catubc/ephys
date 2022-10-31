import struct, h5py
import numpy as np
import matplotlib.pyplot as plt
from array import array
import time

print "Loading Traces"

#file_name = '/home/cat/neuron/in_vivo/dan/M149873/2015-02-06_13-02-32_flash_CSD/data/openephysLGN.dat'
file_name = '/home/cat/neuron/in_vivo/dan/M150566/2015-04-01_16-41-14_flashes_vis/openephysLGN.dat'

##Read h5py files
#f=h5py.File(file_name, 'r')
#print f.keys()[0]
#ec_traces=f[f.keys()[0]]
#print ec_traces.shape
##ec_traces= np.array(ec_traces).reshape(1386872704,0)
##ec_traces.reshape(128,10834943)

#with open(file_name, 'rb') as f:
    #data = f.read()

##print data
#quit()

ec_traces = np.fromfile(file_name, dtype=np.int16)
ec_traces = ec_traces.reshape(len(ec_traces)/128, 128)
ec_traces = ec_traces.T
print ec_traces
print ec_traces.shape

##Plot first 10000 steps from channel 0
#x = np.arange(0,25000,1.)/25000.
#for i in range(128):
    #plt.plot(x,ec_traces[i][0:25000]+100*i, 'r-', color='black',linewidth=1)

#plt.show()

n_electrodes = ec_traces.shape[0]
n_vd_samples = ec_traces.shape[1]


#**************** Generate Header Info ****************************
header = 'Test spike file '
iformat = 1002
SampleFrequency = 25000
vscale_HP = .15
Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array

for i in range (n_electrodes-1,0,-1):
    Siteloc[i][0]=i%2 * 22
    Siteloc[i][1]=i/2*22

print Siteloc
#*********** SAVE .TSF FILE *
#read(14)SampleFrequency,n_channels,n_vd_samples,vscale_HP
#read(14)(SiteLoc(1,n),SiteLoc(2,n),electrode_sort_order(n),n=1,n_channels)

file_name = file_name+'.tsf'

f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples))
f.write(struct.pack('f', vscale_HP))
for i in range (n_electrodes):
    #Siteloc[i][0].tofile(f)
    #Siteloc[i][1].tofile(f)
    f.write(struct.pack('h', Siteloc[i][0]))
    f.write(struct.pack('h', Siteloc[i][1]))
    f.write(struct.pack('i', n_electrodes-1-i))
print "Writing data"
for i in range (n_electrodes):
    print "Writing ch: ", i
    ec_traces[i].tofile(f)
tempint = 0
f.write(struct.pack('i', tempint)) #Write # of fake spikes
f.close()
