import struct
import numpy as np
import matplotlib.pyplot as plt
from array import array

print "Loading Traces"

file_name = '/media/cat/Data1/in_vitro/all_cell_random_depth/Traces/James_in-vitro_extracellular_traces.tsf'
#file_name = '/media/cat/Data1/in_vivo/openephysLGN.tsf'

#******************************* LOAD TSF FILE *****************************
#Cat: this code should be simplified/made more pythonic; 
with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]  
    print iformat
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency, currently 10KHz
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples (varries)
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1
    if iformat==1001:
        Siteloc = np.zeros((2*56), dtype=np.int16)
        Siteloc = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
    if iformat==1002:
        Siteloc = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)
    ec_traces.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0]  #Read # of ground truth spikes
    print n_cell_spikes

    if n_cell_spikes>0:
        if (iformat==1001):
            vertical_site_spacing = struct.unpack('i',fin.read(4))[0] 
            n_cell_spikes = struct.unpack('i',fin.read(4))[0] 

        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

##Plot first 10000 steps from channel 0
#plt.plot(ec_traces[0][0:10000], 'r-', color='black',linewidth=1)
#plt.show()

#**************** Generate Header Info If Required **************************
#header = 'Test spike file '
#iformat = 1001
#SampleFrequency = SampleFrequency*2
#vscale_HP = 1.0 
#Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array

#for i in range (n_electrodes):
    #Siteloc[i][0]=i%2 * 22
    #Siteloc[i][1]=i/2*22

print "Saving Traces"
print Siteloc
print iformat

print ec_traces.shape

#************ SAVE .BIN FILE
file_name = file_name+'.bin'

fout = open(file_name, 'wb')
ec_traces.tofile(fout)

quit()
#*********** SAVE .TSF FILE *
SampleFrequency=20000
file_name = file_name+'_20Khz.tsf'

fout = open(file_name, 'wb')
fout.write(header)
fout.write(struct.pack('i', iformat))
fout.write(struct.pack('i', SampleFrequency))
fout.write(struct.pack('i', n_electrodes))
fout.write(struct.pack('i', n_vd_samples))
fout.write(struct.pack('f', vscale_HP))
if iformat==1002:
    for i in range (n_electrodes):
        fout.write(struct.pack('h', Siteloc[i*2]))
        fout.write(struct.pack('h', Siteloc[i*2+1]))
        fout.write(struct.pack('i', i))
else:
    fout.write(struct.pack(str(2*56)+'h', *Siteloc))

ec_traces.tofile(fout)

fout.write(struct.pack('i', n_cell_spikes))
print "No. cell spikes: ", n_cell_spikes
if (n_cell_spikes>0):
    if (iformat==1001):
        fout.write(struct.pack('i', vertical_site_spacing)) # = struct.unpack('i',fin.read(4))[0] 
        fout.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    fake_spike_times.tofile(fout)
    fake_spike_assignment.tofile(fout) 
    fake_spike_channels.tofile(fout) 
fout.close()
