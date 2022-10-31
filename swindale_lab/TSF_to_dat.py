import struct
import numpy as np
import matplotlib.pyplot as plt
from tsf_ptcs_classes import *

print "Loading Traces"

###****************Read raw .dat file****************
#file_name = '/media/cat/4TB/in_silico/ucsd_Sep6_rat_3k_20Khz/ECP_0.dat'

#with open(file_name, 'rb') as f:
    #data = np.fromfile(file_name, dtype=np.int16)

#data=data.reshape(data.shape[0]/30,30)
#data = data.T
#print data.shape
##print data[0][0:1000]
#for i in range(len(data)):
    #plt.plot(data[i][0:10000]+100*i)

#plt.show()
#quit()


#****************Read first dataset****************
#file_name = '/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/ECP_1.tsf'
#file_name = '/media/cat/4TB/in_vitro/20Khz_10cells/20Khz_10cells.tsf'
#file_name = '/media/cat/4TB/in_vivo/sev/M34_20132808/1SFTF1.dat_chs_2ndShank_4min.tsf'

#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/track5c_59-76_notch_1_240Hz_5subsample/track5c_59-76_notch_1_240Hz_5subsample_uncompressed.tsf'
#file_name = '/media/cat/4TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/2015-11-27-16-deep-iso0_lp.tsf'
#file_name = '/media/cat/4TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-2-10electrodein-iso1/2015-11-27-2-10electrodein-iso1_lp.tsf'

#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/57-tr5c-mseq32_40ms/57-tr5c-mseq32_40ms.tsf'
#file_name = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-4-10electrodein-iso0/2015-11-27-4-10electrodein-iso0_raw.tsf'
file_name = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-2-10electrodein-iso1/2015-11-27-2-10electrodein-iso1_raw.tsf'

if True:
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

#**** Save .lfp.zip file
if True:
   
    data = ec_traces
    t0 = 0
    t1 = len(data[0])*1E6/SampleFrequency  #time of end of file in usec 
    tres = 1000     #1Khz sample rate
    uVperAD = 1.0
    chans = np.arange(0, n_electrodes, 1)
    #print Siteloc1
    chanpos = Siteloc1
    #quit()
    
    out_file = file_name+'.lfp.zip'
    
    #If file already low pass filtered:
    if False:
        np.savez(out_file, data=data, t0=t0, t1=t1, tres=tres)
    else:
        temp=[]
        for k in range(len(data)):      #Subsample to 1Khz and notch filter
            temp.append(np.array(filter.notch(data[k][::SampleFrequency/1000])[0], dtype=np.int16)) # remove 60 Hz mains noise
        
        #Band pass below 240Hz
        lowcut = 0.1; highcut=240; fs=1000
        for c in range(len(temp)):      #Subsample down to 1Khz for the filter
            temp[c]= np.array(butter_bandpass_filter(temp[c], lowcut, highcut, fs, order = 2), dtype=np.int16)
                
        
        np.savez(out_file, chans=chans, chanpos=chanpos, data=temp, t0=t0, t1=t1, tres=tres, uVperAD=uVperAD)
    
    #np.savez(out_file, data=data[::20], t0=t0, t1=t1, tres=tres) #Use to convert from raw data;
    quit()

file_name = '/media/cat/4TB/in_vivo/nick/ptc21/57-tr5c-mseq32_40ms/57-tr5c-mseq32_40ms.lfp.zip'
data_in = np.load(file_name)
ec_traces = data_in['data']
print data_in.keys()

temp=[]
for k in range(len(ec_traces)):      #No need to subsample Martin's LFP records before notch filter
    temp.append(np.array(filter.notch(ec_traces[k])[0], dtype=np.int16)) # remove 60 Hz mains noise
ec_traces = np.array(temp)

#Pad beginning of record to match high-pass record
start_offset = np.zeros((len(data_in['data']), data_in['t0']*1E-3), dtype=np.int16)
ec_traces = np.concatenate((start_offset, ec_traces), axis=1)

#Set Siteloc values
n_electrodes=10
Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
for i in range(n_electrodes):
    Siteloc[i][0]=0
    Siteloc[i][1]=i*100

#Subsample datarecord for LFP sort compression
temp = []
for k in range(len(ec_traces)):
    temp.append(ec_traces[k][::5])
ec_traces=np.int16(temp)

#*********** SAVE TIME COMPRESSED .TSF FILE *
if True:  #save new .tsf file with additional uncorrelated noise
    
    file_name = file_name+'_all_notch_5subsample.tsf'

    header = 'Test spike file '
    iformat = 1002
    SampleFrequency = 20000
    n_electrodes = 10
    n_vd_samples = len(ec_traces[0])
    vscale_HP = data_in['uVperAD']
    n_cell_spikes = 0

    f = open(file_name, 'wb')
    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP))

    for i in range (n_electrodes):
        f.write(struct.pack('h', Siteloc[i][0]))
        f.write(struct.pack('h', Siteloc[i][1]))
        f.write(struct.pack('i', i+1))

    print "Writing data"

    for i in range (n_electrodes):
        print "Writing channel: ", i
        ec_traces[i].tofile(f)

    #f.write(struct.pack('i', 0)) #Write zero spikes
    print "No. cell spikes: ", n_cell_spikes

    f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    if (n_cell_spikes>0):
        fake_spike_times.tofile(f)
        (fake_spike_assignment+1).tofile(f) 
        fake_spike_channels.tofile(f) 
        
    f.close()
