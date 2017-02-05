import struct
import numpy as np
import matplotlib.pyplot as plt

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
#file_name = '/media/cat/4TB/in_silico/ucsd_Sep13_rat_3k_20khz/ECP_2.tsf'
file_name = '/media/cat/4TB/in_silico/ucsd_Sep6_rat_3k_20Khz/ECP_0.tsf'

with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
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
    ec_temp = ec_traces1.copy()
    ec_traces1.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)

rec1_length = n_vd_samples

print Siteloc1
quit()

print "Length of rec1: ", rec1_length
print "vscale_HP: ", vscale_HP
print "n_electrodes: ", n_electrodes 

print "spike times: ", fake_spike_times
print "# total spikes: ", len(fake_spike_times)
print "spike ids: ", fake_spike_assignment
print "# units: ", len(np.unique(np.array(fake_spike_assignment)))
print "spike channels: ", fake_spike_channels

#quit()

#Save .dat file for KlustaKwik KK
if True: 

    #Add uncorrelated noise
    noise = np.array(np.random.randn(len(ec_temp))*5, dtype=np.int16) #insert normalized noise to 5uV SD; also multiply by scaling factor of 10
    ec_temp = ec_temp + noise

    ec_temp.shape = n_electrodes, n_vd_samples
    kk_traces = ec_temp.T.copy()
    
    #for i in range(len(kk_traces)):
    f = open(file_name[:-4]+'.dat', 'wb')
    print "Saving .dat file"
    (kk_traces).tofile(f)
    f.close()
    print "Done writing .dat"
    quit()

##****************Read second dataset****************
#file_name = '/media/cat/4TB/in_vivo/peter/lp.h5/lp.h5_truth.tsf'

#with open(file_name,'rb') as fin:
    #header = fin.read(16)                                   #Header info; not used
    #iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    #SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    #n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes, currently 8, Buzsaki H32 single shank
    #n_vd_samples2 = struct.unpack('i',fin.read(4))[0]       #No. of samples (varries)
    #vscale_HP = struct.unpack('f',fin.read(4))[0]           #Assume same scaling - but could be df

    #if iformat==1001:                                       #Assuming same siteloc mapping
        #Siteloc1 = np.zeros((2*56), dtype=np.int16)
        #Siteloc1 = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
    #if iformat==1002:
        #Siteloc1 = np.zeros((2*n_electrodes), dtype=np.int16)
        #Readloc = np.zeros((n_electrodes), dtype=np.int32)
        #for i in range(n_electrodes):
            #Siteloc1[i*2] = struct.unpack('h', fin.read(2))[0]
            #Siteloc1[i*2+1] = struct.unpack('h', fin.read(2))[0]
            #Readloc[i] = struct.unpack('i', fin.read(4))[0]

    #ec_traces2 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples2)
    #ec_traces2.shape = n_electrodes, n_vd_samples2

    #n_cell_spikes2 = struct.unpack('i',fin.read(4))[0] 
    #print "No. cell spikes: ", n_cell_spikes2
    #if (n_cell_spikes2>0):
        #fake_spike_times2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)
        #fake_spike_assignment2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)
        #fake_spike_channels2 =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes2)

#rec2_length = n_vd_samples2
#print "Length of rec1: ", rec2_length

#*********** SAVE .TSF FILE *
file_name = file_name[:-4]+'_noise.tsf'
iformat = 1002
print vscale_HP

##Save only single column of IMEC probe:
#n_vd_samples2=0
#n_cell_spikes2=0
#fake_spike_channels*=1

f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples))
f.write(struct.pack('f', vscale_HP))

for i in range (n_electrodes):
    print "Writing electrode: ", i, Siteloc1[i*2], Siteloc1[i*2+1], Readloc[i]+1
    f.write(struct.pack('h', Siteloc1[i*2]))
    f.write(struct.pack('h', Siteloc1[i*2+1]))
    f.write(struct.pack('i', Readloc[i]+1))

print "Writing data"

for i in range (n_electrodes):
    print "Writing channel: ", i
    np.array((ec_traces1[i]+np.array(np.random.randn(len(ec_traces1[i])))*5),dtype=np.int16).tofile(f)

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
