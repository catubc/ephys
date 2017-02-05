#Function to generate compressed LFP data for spike sorting from raw (i.e. full bandwidth and samplerate) data
#TODO: add functionality for Nick's lfp data which is 1Khz and has particular format including offsets.
#TODO: add functionality for Martin's LFP data which is already available

import numpy as np
import matplotlib.pyplot as plt
import struct
from tsf_ptcs_classes import *

#************DONGSHENG RECS *******************
#data_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/'
#file_name = '2015-11-27-2-10electrodein-iso1'
#file_name = '2015-11-27-4-10electrodein-iso0'

#data_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-1/'
#file_name = '2015-12-1-1-10electrodeiniso1'

#data_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-2/'
#file_name = '2015-12-2-2-10electrodeincortex-iso1'

#************CAT RECORDINGS **************

#data_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr1/'
#file_name = '2016_05_03_tr1_spont_afterinsertion_noisyamp_160503_145951'

#data_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr1/'
#file_name = '2016_05_03_tr1_spont_deep_noisyamp_160503_152854'

data_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr2/'
file_name = '2016_05_03_tr2_spont_afterinsertion_160503_162252'

#data_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr3/'
#file_name = '2016_05_03_tr3_spont_160503_182605'


if True:
    if 'dongsheng' in data_dir: 
        file_raw = data_dir + file_name+'/'+file_name+'_raw.tsf'
        out_file = file_raw[:-8]+'.lfp.zip'
    else: 
        file_raw = data_dir + file_name+'.tsf'
        out_file = file_raw[:-4]+'.lfp.zip'

    #*********LOAD TSF FILE********
    with open(file_raw,'rb') as fin:
        header = fin.read(16)                                   #Header info; not used
        iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
        SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
        n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
        n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
        vscale_HP = struct.unpack('f',fin.read(4))[0]          #Scaling of int2 values below to save space, currently 0.1
        print vscale_HP

        Siteloc = np.zeros((2*n_electrodes), dtype=np.int16)   #NB iformat 1001 is slightly different
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

        ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples) 
        ec_traces.shape = n_electrodes, n_vd_samples
        
    #**** Save .lfp.zip file
    data = ec_traces
    t0 = 0
    t1 = len(data[0])*1E6/SampleFrequency  #time of end of file in usec 
    tres = 1000     #1Khz sample rate
    uVperAD = 1.0
    chans = np.arange(0, n_electrodes, 1)
    chanpos = Siteloc
    
    
    #************BUTTERWORTH FILTER DATA***********
    if True:
        #Notch filter
        temp=[]
        for k in range(len(data)):      #Subsample to 1Khz and notch filter
            temp.append(np.array(filter.notch(data[k][::SampleFrequency/1000])[0], dtype=np.int16)) # remove 60 Hz mains noise

        #Band pass below 240Hz
        lowcut = 0.1; highcut=240; fs=1000
        lowcut = 5; highcut=240; fs=1000    #Use 5Hz low cutoff to reduce slower oscillations for spike sorting;
        for c in range(len(temp)):      #Subsample down to 1Khz for the filter
            temp[c]= np.array(butter_bandpass_filter(temp[c], lowcut, highcut, fs, order = 2), dtype=np.int16)
        
    print "Saving LFP : ", out_file
    np.savez(out_file, chans=chans, chanpos=chanpos, data=temp, t0=t0, t1=t1, tres=tres, uVperAD=uVperAD)
    os.rename(out_file+'.npz', out_file)
    
    #****************COMPRESS LFP *********
    if 'dongsheng' in data_dir: file_lfp = data_dir + file_name+'/'+file_name
    else: file_lfp = data_dir + file_name

    print "Loading LFP (.lfp.zip) file: ", file_lfp
    
    data_in  = np.load(file_lfp+'.lfp.zip')
    t_start = data_in['t0']*1E-6      #Convert to seconds
    t_end = data_in['t1']*1E-6        #Convert to seconds
    start_offset = np.zeros((len(data_in['data']), int(t_start*1E3)), dtype=np.int16)  #make padding at front of data
    
    #Notch filter 60Hz frequency
    temp_data = data_in['data']
    for k in range(len(temp_data)):
        temp_data[k] = np.array(filter.notch(temp_data[k])[0], dtype=np.int16) # remove 60 Hz mains noise
    data_out = temp_data

    #Combine filtered data and front end padding
    data_out = np.concatenate((start_offset, data_out), axis=1)

    #*********** SAVE COMPRESSED LOW PASS .TSF FILE *********
    header = 'Test spike file '
    iformat = 1002
    n_electrodes = len(data_in['data'])
    Compress_factor = 20
    Subsample = .1
    SampleFrequency = 1000 * Compress_factor
    vscale_HP = vscale_HP #use the same as above
    overall_sample = Compress_factor*Subsample
    overall_compression = Compress_factor/Subsample
    #header = 'Test spike file '
    #header = str('compress:  '+str(overall_compression)).zfill(16)
    print header, len(header)

    n_vd_samples = len(data_out[0][::overall_sample])

    SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
    for i in range (n_electrodes):
        SiteLoc[i][0]=data_in['chanpos'][i*2]
        SiteLoc[i][1]=data_in['chanpos'][i*2+1]
        
    file_name = file_lfp + '_lp_compressed.tsf'
    f1 = open(file_name, 'wb')

    f1.write(header)
    f1.write(struct.pack('i', iformat))
    f1.write(struct.pack('i', SampleFrequency))
    f1.write(struct.pack('i', n_electrodes))
    f1.write(struct.pack('i', n_vd_samples))
    f1.write(struct.pack('f', vscale_HP))
    for i in range (n_electrodes):
        f1.write(struct.pack('h', SiteLoc[i][0]))
        f1.write(struct.pack('h', SiteLoc[i][1]))
        f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

    print "Writing data"
    for i in range(n_electrodes):
        data_out[i][::overall_sample].tofile(f1)

    f1.write(struct.pack('i', 0)) #Write # of fake spikes
    text = "Compression:" + str(overall_compression)
    n_bytes = len(text)
    f1.write(struct.pack('i', n_bytes))
    f1.write(text)
    f1.close()


    #*********** SAVE TOP CHANNELS ONLY LOW PASS .TSF FILE *********
    n_electrodes = 30
       
    file_name = file_lfp + '_lp_compressed_30chs.tsf'
    f1 = open(file_name, 'wb')

    f1.write(header)
    f1.write(struct.pack('i', iformat))
    f1.write(struct.pack('i', SampleFrequency))
    f1.write(struct.pack('i', n_electrodes))
    f1.write(struct.pack('i', n_vd_samples))
    f1.write(struct.pack('f', vscale_HP))
    for i in range (n_electrodes):
        f1.write(struct.pack('h', SiteLoc[i][0]))
        f1.write(struct.pack('h', SiteLoc[i][1]))
        f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

    print "Writing data"
    for i in range(n_electrodes):
        data_out[i][::overall_sample].tofile(f1)

    f1.write(struct.pack('i', 0)) #Write # of fake spikes
    text = "Compression:" + str(overall_compression)
    n_bytes = len(text)
    f1.write(struct.pack('i', n_bytes))
    f1.write(text)
    f1.close()
        
    print "DONE"
    #quit()
