import struct
import numpy as np
import matplotlib.pyplot as plt
from array import array
from scipy.signal import butter, lfilter, filtfilt
import filter
from scipy.signal import hilbert

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

print "Loading Traces"

#file_name = '/media/cat/12TB/in_vitro/all_cell_random_depth/Traces/James_in-vitro_extracellular_traces.tsf'
#file_name = '/media/cat/Data1/in_vivo/openephysLGN.tsf'
#file_name = '/media/cat/12TB/in_vivo/nick/ptc17/57-tr2b-MVI_1400_5s/57-tr2b-MVI_1400_5s.tsf'

#PTC17 - tr2b
#file_name = '/media/cat/12TB/in_vivo/nick/ptc17/57-tr2b-MVI_1400_5s/57-tr2b-MVI_1400_5s.lfp.zip'
file_name = '/media/cat/12TB/in_vivo/nick/ptc17/44-tr2b-blankscreen/44-tr2b-blankscreen.lfp.zip'

#PTC21 - tr5c
#file_name = '/media/cat/12TB/in_vivo/nick/ptc21/61-tr5c-blankscreen/61-tr5c-blankscreen.lfp.zip'


#Dongsheng - 2015-11-27
#file_name = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-2-10electrodein-iso1/2015-11-27-2-10electrodein-iso1_raw.tsf'

if 'tim' in file_name:
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

else:
    lfp_data = np.load(file_name)
    print lfp_data.keys()
    print lfp_data['uVperAD']
    
    ec_traces = lfp_data['data'] #*1.0 #float(lfp_data['uVperAD'])
    header = 'Test spike file '
    SampleFrequency = 1000         #Martin's lfp data
    
    n_vd_samples = len(ec_traces[0])
    n_electrodes = 10
    vscale_HP = 1.0
    n_cell_spikes = 0

    Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
    for i in range (n_electrodes):
        Siteloc[i][0]=0
        Siteloc[i][1]=100*i

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

print ec_traces.shape

#*********** SAVE FILTERED .TSF FILE *

SampleFrequency=SampleFrequency
fs = SampleFrequency

for cut in range(10):
    lowcut = cut*10
    highcut = (cut+1)*10
    #file_name = file_name+'_'+str(lowcut)+'_'+str(highcut)+'Hz_bandpass.tsf'

    range_x1 = 0 
    range_x2 = range_x1 + 100 * SampleFrequency     #Number of seconds of data to display
    y_offset = -1500

    t = np.linspace(0, (range_x2-range_x1),(range_x2-range_x1))/SampleFrequency

    save_array=[]
    for k in range(n_electrodes):
        print k
        temp = filter.notch(ec_traces[k][range_x1:range_x2])[0]
        temp = np.array(butter_bandpass_filter(temp, lowcut, highcut, fs, order = 2), dtype=np.int16)

        analytic_signal = hilbert(temp).copy()
        
        save_array.append(analytic_signal[::1000])
        
        #plt.plot(t, temp+k*y_offset, color='black')
        #plt.plot(t, np.abs(analytic_signal)+k*y_offset, color='blue', linewidth=2)
    #plt.show()

    save_array = np.array(save_array)

    np.save(file_name+"_envelope_"+str(highcut)+"hz", save_array)

quit()

fout = open(file_name, 'wb')
fout.write(header)
fout.write(struct.pack('i', 1002))
fout.write(struct.pack('i', SampleFrequency))
fout.write(struct.pack('i', n_electrodes))
fout.write(struct.pack('i', n_vd_samples))
fout.write(struct.pack('f', vscale_HP))
for i in range (n_electrodes):
    fout.write(struct.pack('h', Siteloc[i][0]))
    fout.write(struct.pack('h', Siteloc[i][1]))
    fout.write(struct.pack('i', i+1))

    
#Apply butterworth filters - if required
for i in range(n_electrodes):
    ec_traces[i] = filter.notch(ec_traces[i])[0]
    ec_traces[i] = np.array(butter_bandpass_filter(ec_traces[i], lowcut, highcut, fs, order = 2), dtype=np.int16)
                        
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
