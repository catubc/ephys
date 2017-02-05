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
#file_name = '/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/ECP_1.tsf'
#file_name = '/media/cat/4TB/in_vitro/20Khz_10cells/20Khz_10cells.tsf'
#file_name = '/media/cat/4TB/in_vivo/sev/M34_20132808/1SFTF1.dat_chs_2ndShank_4min.tsf'

file_name = '/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/vivo_1.tsf'


with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
    vscale_HP1 = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1
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

    ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples) #Increase in vivo data by factor of 2.

    #IN SILICO DATA ONLY (NB: Add uncorrelated noise before separating data for .tsf and .dat file; so that traces are identical)
    if False:  
        noise = np.array(np.random.randn(len(ec_traces))*5, dtype=np.int16) #insert normalized noise to 5uV SD; also multiply by scaling factor of 10
        ec_traces = ec_traces + noise
    
    
    #Reshape traces
    ec_traces.shape = n_electrodes, n_vd_samples
    
    #Reorganize traces to match Recording #2
    temp_array = [0,7,1,6,2,5,3,4]
    ec_traces = ec_traces[temp_array]
    #for i in range(len(temp_array)):
    #    temp_traces = ec_traces[
    
    #Cut off part of traces if required:
    #print "Length of recording (mins): ", n_vd_samples / SampleFrequency / 60
    #n_vd_samples = 4*60*SampleFrequency 
    #print "Choped down to 4 mins =  ", n_vd_samples , " timesteps"
    #ec_traces = ec_traces[:,0:n_vd_samples]

    


    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)


#************ Load 2nd .tsf file **************

#****************Read first dataset****************
#file_name = '/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/ECP_1.tsf'
#file_name = '/media/cat/4TB/in_vitro/20Khz_10cells/20Khz_10cells.tsf'

file_name = '/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/ECP_1.tsf'

with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1
    print vscale_HP

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
    print Readloc
    print Siteloc1

    ec_traces2 =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)*5 #x 10 / 2 to decrease amplitude of cells to be closer to in vivo
    
    #Reshape traces
    ec_traces2.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)


#****** Save .dat file for KlustaKwik KK ********
if False: 
    ec_datfile.shape = n_electrodes, n_vd_samples
    ec_datfile = ec_datfile.T.copy()
    
    #for i in range(len(kk_traces)):
    f = open(file_name[:-4]+'_noise.dat', 'wb')
    print "Saving .dat file"
    (ec_datfile).tofile(f)
    f.close()
    print "Done writing .dat"
    #quit()

#*********** SAVE .TSF FILE *
if True:  #save new .tsf file with additional uncorrelated noise
    file_name = file_name[:-4]+'_4min_hybrid.tsf'
    iformat = 1002
    print vscale_HP1

    f = open(file_name, 'wb')
    f.write(header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', SampleFrequency))
    f.write(struct.pack('i', n_electrodes))
    f.write(struct.pack('i', n_vd_samples))
    f.write(struct.pack('f', vscale_HP1))

    for i in range (n_electrodes):
        print "Writing electrode: ", i, Siteloc1[i*2], Siteloc1[i*2+1], Readloc[i]+1
        f.write(struct.pack('h', Siteloc1[i*2]))
        f.write(struct.pack('h', Siteloc1[i*2+1]))
        f.write(struct.pack('i', Readloc[i]))

    print "Writing data"

    for i in range (n_electrodes):
        print "Writing channel: ", i
        print ec_traces[i][0:100]
        print ec_traces2[i][0:100]
        (ec_traces[i]+ec_traces2[i]).tofile(f)


    #f.write(struct.pack('i', 0)) #Write zero spikes
    print "No. cell spikes: ", n_cell_spikes

    f.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    if (n_cell_spikes>0):
        fake_spike_times.tofile(f)
        (fake_spike_assignment+1).tofile(f) 
        fake_spike_channels.tofile(f) 
        
    f.close()
