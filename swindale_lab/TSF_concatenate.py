import struct
import numpy as np
import matplotlib.pyplot as plt
from tsf_ptcs_classes import *
print "Loading Traces"


#concanteate 17th channel from multi .mcd files
if True:
    file_name = '/media/cat/12TB/in_vivo/tim/2016-1-14/2016-1-14-3-10electrodein cortex-iso1/2016-1-14-3-10electrodein cortex-iso1_1.mcd'
    
    data = MCD_read_imagingtimes_old(file_name)

    print "Finding beginning and end of imaging trigger on channel 17th .mcd file"
    temp_data = []
    for i in range(data['extra'].item_count):
        temp_data.append(data['extra'].get_data(i)) 
    temp_data = np.array(temp_data)[:,0]    #Select time column only from 17th channel

    start_array = []
    end_array = []
    start_array.append(temp_data[0])
    for i in range(1,len(temp_data)-1,1):
        if temp_data[i+1]-temp_data[i]>1.0:
            end_array.append(temp_data[i])
            start_array.append(temp_data[i+1])
    end_array.append(temp_data[-1])
    
    ##Visualize different epoch starts/ends - FOR SINGLE .mcd file
    #data_array = np.vstack((start_array,end_array))
    #data_array = data_array.T
    #print data_array
    #for i in range(len(start_array)):
        #plt.plot([start_array[i],start_array[i]],[0,len(temp_data)])
        #plt.plot([end_array[i],end_array[i]],[0,len(temp_data)])
        #plt.axvspan(start_array[i], end_array[i], color='red', alpha=0.35)
    #plt.show()
    #quit()
    
    n_epochs = len(start_array)-1

    file_name = '/media/cat/12TB/in_vivo/tim/2016-1-14/2016-1-14-3-10electrodein cortex-iso1/2016-1-14-3-10electrodein cortex-iso1_2.mcd'
    data = MCD_read_imagingtimes_old(file_name)

    temp_data = []
    for i in range(data['extra'].item_count):
        temp_data.append(data['extra'].get_data(i)) 
    temp_data = np.array(temp_data)[:,0]    #Select time column only from 17th channel

    start_array.append(temp_data[0]+end_array[n_epochs])
    for i in range(1,len(temp_data)-1,1):
        if temp_data[i+1]-temp_data[i]>1.0:
            end_array.append(temp_data[i]+end_array[n_epochs])
            start_array.append(temp_data[i+1]+end_array[n_epochs])
    end_array.append(temp_data[-1]+end_array[n_epochs])

    data_array = np.vstack((start_array,end_array))
    data_array = data_array.T
    #data_array = np.delete(data_array,2, axis=0)
    
    
    print data_array
    #quit()

    #Visualize different epoch starts/ends
    #plt.plot(temp_data)
    for i in range(len(start_array)):
        plt.plot([start_array[i],start_array[i]],[0,len(temp_data)])
        plt.plot([end_array[i],end_array[i]],[0,len(temp_data)])
        plt.axvspan(start_array[i], end_array[i], color='red', alpha=0.35)
    plt.show()
    quit()
    

#****************Read first dataset****************
file_name = '/media/cat/12TB/in_vivo/tim/2016-1-14/2016-1-14-3-10electrodein cortex-iso1/2016-1-14-3-10electrodein cortex-iso1_raw_1.tsf'

with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples1 = struct.unpack('i',fin.read(4))[0]        #No. of samples 
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

    ec_traces1 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples1)
    ec_traces1.shape = n_electrodes, n_vd_samples1

    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

print "Length of rec1: ", n_vd_samples1/SampleFrequency

#****************Read second dataset****************
file_name = '/media/cat/12TB/in_vivo/tim/2016-1-14/2016-1-14-3-10electrodein cortex-iso1/2016-1-14-3-10electrodein cortex-iso1_raw_2.tsf'


with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    n_vd_samples2 = struct.unpack('i',fin.read(4))[0]       #No. of samples (varries)
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Assume same scaling - but could be df

    if iformat==1001:                                       #Assuming same siteloc mapping
        Siteloc1 = np.zeros((2*56), dtype=np.int16)
        Siteloc1 = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
    if iformat==1002:
        Siteloc1 = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc1[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc1[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ec_traces2 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples2)
    ec_traces2.shape = n_electrodes, n_vd_samples2

    n_cell_spikes2 = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes2
    if (n_cell_spikes2>0):
        fake_spike_times2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)
        fake_spike_assignment2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)
        fake_spike_channels2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)

print "Length of rec2: ", n_vd_samples2/SampleFrequency


#*********** SAVE .TSF FILE *
file_name = file_name+'_combined.tsf'
iformat = 1002
print vscale_HP

f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples1+n_vd_samples2))
f.write(struct.pack('f', vscale_HP))

for i in range (n_electrodes):
    f.write(struct.pack('h', Siteloc1[i*2]))
    f.write(struct.pack('h', Siteloc1[i*2+1]))
    f.write(struct.pack('i', Readloc[i]))

print "Writing data"

for i in range (n_electrodes):
    print "Writing channel: ", i
    ec_traces1[i].tofile(f)
    ec_traces2[i].tofile(f)

f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes

print "No. cell spikes: ", n_cell_spikes
if (n_cell_spikes>0):
    if (iformat==1001):
        f.write(struct.pack('i', vertical_site_spacing)) # = struct.unpack('i',fin.read(4))[0] 
        f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    fake_spike_times.tofile(f)
    (fake_spike_assignment+1).tofile(f) 
    fake_spike_channels.tofile(f) 
f.close()
