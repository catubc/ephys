import struct
import numpy as np
import matplotlib.pyplot as plt

print "Loading Traces"

file_name = '/media/cat/Data1/in_vivo/sev/M34_20132808/2SPARSE1.dat_shank3.tsf'

#Read in-vivo data first
#Cat: this code should be simplified/made more pythonic; 
with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency, currently 10KHz
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples (varries)
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1

    if iformat==1001:
        Siteloc1 = np.zeros((2*56), dtype=np.int16)
        Siteloc1 = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
    if iformat==1002:
        Siteloc1 = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc1[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc1[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ec_traces1 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)
    ec_traces1.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        if (iformat==1001):
            vertical_site_spacing = struct.unpack('i',fin.read(4))[0] 
            n_cell_spikes = struct.unpack('i',fin.read(4))[0] 

        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

print "Sample freq - in vivo - : ", SampleFrequency

print Siteloc1

#********************************************************************************************

file_name = '/media/cat/Data1/in_vitro/20Khz_10cells/20Khz_10cells_truth.tsf'

#Read in-vitro data
with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency, currently 10KHz
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples (varries)
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1

    if iformat==1001:
        Siteloc2 = np.zeros((2*56), dtype=np.int16)
        Siteloc2 = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
    if iformat==1002:
        Siteloc2 = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc2[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc2[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ec_traces2 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)
    ec_traces2.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        if (iformat==1001):
            vertical_site_spacing = struct.unpack('i',fin.read(4))[0] 
            n_cell_spikes = struct.unpack('i',fin.read(4))[0] 

        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

print "Sample freq - in vitro - : ", SampleFrequency


print Siteloc2

mapping =([0,2,4,6,7,5,3,1])

#[-24   0 -18  60 -12 120  -6 180   0 210   6 150  12  90  18  30]

#-24, 0, 
#18, 30, 
#-18, 60, 
#12, 90, 
#-12, 120, 
#6, 150, 
#-6, 180, 
#0, 210

#Plot first 10000 steps from channel 0
#x = np.arange(0,5000, .05)
#plt.plot(x, ec_traces[0][0:100000], 'r-', color='black',linewidth=1)
#plt.show()


n_electrodes=8

#*********** SAVE .TSF FILE *
file_name = '/media/cat/Data1/in_vivo/hybrid/2SPARSE1.dat_shank3_20Khz_10cells_truth.tsf'
iformat = 1002
f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples))
f.write(struct.pack('f', vscale_HP))
for i in range (n_electrodes):
    f.write(struct.pack('h', Siteloc1[i*2]))
    f.write(struct.pack('h', Siteloc1[i*2+1]))
    f.write(struct.pack('i', Readloc[i]))
print "Writing data"
for i in range (n_electrodes):
    print "Writing ch: ", i
    (ec_traces1[i][0:n_vd_samples]+ec_traces2[mapping[i]][0:n_vd_samples]).tofile(f)

#n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes

print "No. cell spikes: ", n_cell_spikes
if (n_cell_spikes>0):
    if (iformat==1001):
        f.write(struct.pack('i', vertical_site_spacing)) # = struct.unpack('i',fin.read(4))[0] 
        f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    fake_spike_times.tofile(f)
    fake_spike_assignment.tofile(f) 
    fake_spike_channels.tofile(f) 
f.close()
