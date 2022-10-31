from neo import io
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter
import filter
from tsf_ptcs_classes import *


file_dir = '/media/cat/12TB/in_vivo/tim/2015-11-18/'
file_name = '2015-11-18-2-9electrodein'

r = io.AxonIO(filename=file_dir+file_name+'/'+file_name+'.abf')
bl = r.read_block(lazy=False, cascade=True)

signal = bl.segments[0].analogsignals[1]

print bl.segments[0].analogsignals
print bl.segments[0].eventarrays
quit()
plt.plot(signal)
plt.show()

sample_rate = 9803.92156863

times = np.arange(0, len(signal), 1.)/sample_rate

#start = 0.024
#end = 0.030
#step = end-start
#step_freq = 1./6.1

#print start, end

#check_signal = np.int16(times)*0
#for i in range(30):
    #print [int((start + step_freq*i)*sample_rate), int((start + step_freq*i + step)*sample_rate)]
    #signal[int((start + step_freq*i)*sample_rate) : int((start + step_freq*i + step)*sample_rate)]=0
    #check_signal[int((start + step_freq*i)*sample_rate) : int((start + step_freq*i + step)*sample_rate)]=1


header = 'Test spike file '
iformat = 1002
n_electrodes = 1
n_vd_samples= len(signal)
vscale_HP = 1.0
SampleFrequency = sample_rate
SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array

for i in range (n_electrodes):
    SiteLoc[i][0]=i%2 * 22
    SiteLoc[i][1]=i/2*22

##***************************************** SAVE .TSF FILE *************************************
file_name = file_dir+file_name+'/'+file_name+'_eeg.tsf'

f = open(file_name, 'wb')
f.write(header)
f.write(struct.pack('i', iformat))
f.write(struct.pack('i', SampleFrequency))
f.write(struct.pack('i', n_electrodes))
f.write(struct.pack('i', n_vd_samples))
f.write(struct.pack('f', vscale_HP))

for i in range (n_electrodes):
    f.write(struct.pack('h', SiteLoc[i][0]))
    f.write(struct.pack('h', SiteLoc[i][1]))
    f.write(struct.pack('i', i))

#IMEC data
print "Writing data"
for i in range(n_electrodes):
    np.array(signal*1E4, dtype=np.int16).tofile(f)

tempint = 0
f.write(struct.pack('i', tempint)) #Write # of fake spikes
f.close()
quit()

signal = signal [0:30000]
check_signal = check_signal[0:30000]
times = times[0:30000]
plt.plot(times, signal)
plt.plot(times, check_signal, color='red')
plt.show()
quit()

Plot_PSD_signal(signal, sample_rate)
quit()      
Plot_specgram_signal(signal, sample_rate)
quit()



fs = sample_rate
lowcut = 0.1
highcut = 20.0
signal = butter_bandpass_filter(signal, lowcut, highcut, fs, order = 2)


plt.plot(times, signal)
plt.show()
quit()
print bl.segments[0].analogsignals
print bl.segments[0].eventarrays
