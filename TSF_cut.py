import struct
import numpy as np
import matplotlib.pyplot as plt


file_name = '/media/cat/12TB/in_silico/july_23_7200chs/ECP_0.tsf'

#f = open(file_name[:-4]+'.bin', 'rb')
#bin_data = np.fromfile(f, dtype=np.int16).reshape(7500,-1)
#plt.plot(bin_data[0][:10000])
#plt.plot(bin_data[1][:10000]+1000)
#plt.show()
#quit()


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
    
    print Siteloc1
    np.savetxt(file_name[:-4]+'_layout.txt', Siteloc1)
    print "...loading traces..."
    ec_traces1 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)

    #IN SILICO DATA ONLY (NB: Add uncorrelated noise before separating data for .tsf and .dat file; so that traces are identical)
    if 'silico' in file_name:  
        print "...generating noise..."
        noise = np.array(np.random.randn(len(ec_traces1))*5, dtype=np.int16) #insert normalized noise to 5uV SD; also multiply by scaling factor of 10
        print "...adding noise..."
        ec_traces1 = ec_traces1 + noise
    
    ec_traces1.shape = n_electrodes, n_vd_samples
    print ec_traces1.shape
                
    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        if (iformat==1001):
            vertical_site_spacing = struct.unpack('i',fin.read(4))[0] 
            n_cell_spikes = struct.unpack('i',fin.read(4))[0] 

        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

print "...saving .bin file..."
f = open(file_name[:-4]+'.bin', 'wb')
bin_data.tofile(f)

quit()
#*********** SAVE .TSF FILE *
file_name = file_name+'_noise.tsf'
iformat = 1002
print vscale_HP
print np.random.rand(n_vd_samples)
n_electrodes = 100
electrode_list = np.arange(0,n_electrodes, 2)
#quit()
Siteloc1*=10


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
    print ec_traces1[i][0:10000]
    np.array((ec_traces1[i]+(np.random.standard_normal(n_vd_samples)*5)),dtype=np.int16).tofile(f)

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
