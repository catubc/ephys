import struct
import numpy as np
import matplotlib.pyplot as plt

print "Loading Traces"

#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr3/2016_05_03_tr3_spont_160503_182605_lp_compressed_30chs.tsf'
file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr1_spontaneous_cortex_01_160527_150948.tsf'
camera_onoff = '/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr1_spontaneous_cortex_01_160527_150948_camera_onoff.npy'

#Cat: this code should be simplified/made more pythonic; 
with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
    vscale_HP1 = struct.unpack('f',fin.read(4))[0]          #Scaling of int2 values below to save space, currently 0.1
    print vscale_HP1
    
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

    ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples) 

    ec_datfile = ec_traces.copy()

    #IN SILICO DATA ONLY (NB: Add uncorrelated noise before separating data for .tsf and .dat file; so that traces are identical)
    if False:  
        noise = np.array(np.random.randn(len(ec_traces))*5, dtype=np.int16) #insert normalized noise to 5uV SD; also multiply by scaling factor of 10
        ec_traces = ec_traces + noise
    
    #Reshape traces
    ec_traces.shape = n_electrodes, n_vd_samples
    
    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

    #n_bytes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    #text = fin.read(n_bytes)
    fin.close()

camera_data = np.load(camera_onoff)

plt.plot(ec_traces[0])
plt.plot(camera_data*10, linewidth=2)
plt.show()

print n_bytes
print text
            
#Plot first 10000 steps from channel 0
plt.plot(ec_traces[0][0:10000], 'r-', color='black',linewidth=1)
plt.show()
