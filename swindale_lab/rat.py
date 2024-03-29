#Animal class and subclasses 
#Code: Cat 

from animal_analysis import *

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
import cPickle as pickle
import gc
from skimage.measure import block_reduce
import shutil
from load_intan_rhd_format import *
 

class Rat(object):      
    
    def __init__(self, animal_name, home_dir):
        
        self.name = animal_name
        self.home_dir = home_dir

        self.load_filenames()
        
        #self.probe = Probe()   #Load intan probe map;
    
    def load_filenames(self):
        print "...loading filenames..."
        self.filenames = glob.glob(self.home_dir+self.name+'*.tsf')   #use .tif extension otherwise will load .npy files
        for f in range(len(self.filenames)):
            self.filenames[f] = self.filenames[f][:-4]
        
        #self.tsf_filenames = glob.glob( self.home_dir+self.name+'/tsf_files/*.tsf')   #use .tif extension otherwise will load .npy files
    
    
    def load_tsf_header(self, file_name):
        self.tsf = Tsf_file(file_name)


    def load_tsf(self, file_name):
        self.tsf = Tsf_file(file_name)
        self.tsf.read_ec_traces()
       
        
    def load_channel(self, file_name, channel):
        
        self.tsf = Tsf_file(file_name)
        self.tsf.file_name = file_name
        self.tsf.read_trace(channel)
    

    def load_lfp_channel(self, file_name, channel):     #Nick/Martin data has different LFP structure to their data.
        
        class Object_empty(object):
            def __init__(self):
                pass
        
        self.tsf = Object_empty()
        self.tsf.file_name = file_name
        
        data_in = np.load(file_name)
        self.tsf.SampleFrequency = data_in['tres']
        print "freq: ", self.tsf.SampleFrequency
        self.tsf.ec_traces = data_in['data']
        self.tsf.ec_traces = self.tsf.ec_traces * data_in['uVperAD']
        print self.tsf.ec_traces.shape
        print self.tsf.ec_traces[0]
        
        self.tsf.ec_traces = self.tsf.ec_traces[channel]

        #Home made notch filter; filter.notch doesn't seem to work...
        notched = np.int16(butter_bandpass_filter(self.tsf.ec_traces, 59.75, 60.25, self.tsf.SampleFrequency, order = 2))

        self.tsf.ec_traces = self.tsf.ec_traces-notched

        #compute LFP offset relative high-pass and pad LFP record.
        print "offseting lfp data to highpass"
        self.tsf.ec_traces = np.append(np.zeros(int(data_in['t0']*1E-3), dtype=np.int16), self.tsf.ec_traces)
        #self.tsf.ec_traces = np.concatenate((start_offset, self.tsf.ec_traces), axis=0)
               
    def load_lfp_all(self, file_name):     #Nick/Martin data has different LFP structure to their data.
        
        class Object_empty(object):
            def __init__(self):
                pass
        
        self.tsf = Object_empty()
        self.tsf.file_name = file_name
        
        data_in = np.load(file_name)
        self.tsf.SampleFrequency = data_in['tres']
        print "freq: ", self.tsf.SampleFrequency
        self.tsf.ec_traces = data_in['data']
        self.tsf.ec_traces = self.tsf.ec_traces * data_in['uVperAD']
        print self.tsf.ec_traces.shape
        print self.tsf.ec_traces[0]
        
        #self.tsf.ec_traces = self.tsf.ec_traces[channel]

        #Home made notch filter; filter.notch doesn't seem to work...
        offset = np.zeros(int(data_in['t0']*1E-3), dtype=np.int16)
        for k in range(len()):
            notched = np.int16(butter_bandpass_filter(self.tsf.ec_traces[k], 59.75, 60.25, self.tsf.SampleFrequency, order = 2))

            self.tsf.ec_traces[k] = self.tsf.ec_traces[k]-notched

            #compute LFP offset relative high-pass and pad LFP record.
            print "offseting lfp data to highpass"
            self.tsf.ec_traces[k] = np.append(offset, self.tsf.ec_traces[k])
            #self.tsf.ec_traces = np.concatenate((start_offset, self.tsf.ec_traces), axis=0)
        
        


    def tsf_to_lfp(self):
        '''Read .tsf files - subsample to 1Khz, save as *_lp.tsf
        '''
        
        print "...making low-pass tsf files (1Khz sample rates)..."

        for file_name in self.filenames:
            
            file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp.tsf'
            if os.path.exists(file_out)==True: continue

            file_in = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_raw.tsf'

            print "Processing: \n", file_in

            self.load_tsf(file_in)
            #self.load_channel(file_in, 0)
            print self.tsf.Siteloc.shape
            
            n_vd_samples = len(self.tsf.ec_traces[0]); print "Number of samples: ", n_vd_samples
            
            print "...converting raw to .lfp (1Khz) sample rate tsf files ..."
            temp_array=[]
            lowcut = 0.1; highcut=110; fs=1000
            #import matplotlib.pyplot as plt
            #lowcut = 5; highcut=240; fs=1000    #Use 5Hz low cutoff to reduce slower oscillations for spike sorting;
            for k in range(len(self.tsf.ec_traces)):      #Subsample to 1Khz and notch filter
                print "ch: ", k
                temp = np.array(butter_bandpass_filter(self.tsf.ec_traces[k][::int(self.tsf.SampleFrequency/1000)], lowcut, highcut, fs, order = 2), dtype=np.int16)
                #Home made notch filter; filter.notch doesn't seem to work...
                notch = np.array(butter_bandpass_filter(temp, 59.9, 60.1, fs, order = 2), dtype=np.int16)
                temp = temp-notch
                temp_array.append(temp)
            
            ec_traces = np.int16(temp_array) 
            print ec_traces.shape

            #SAVE RAW DATA ******NB:  SHOULD CLEAN THIS UP: the write function should be shared by all, just data is changing so no need to repeat;
            print "Writing raw data ..."

            fout = open(file_out, 'wb')
            fout.write(self.tsf.header)
            fout.write(struct.pack('i', 1002))
            fout.write(struct.pack('i', 1000))
            fout.write(struct.pack('i', self.tsf.n_electrodes))
            fout.write(struct.pack('i', len(ec_traces[0])))
            fout.write(struct.pack('f', self.tsf.vscale_HP))
            
            for i in range (self.tsf.n_electrodes):
                fout.write(struct.pack('h', self.tsf.Siteloc[i*2]))
                fout.write(struct.pack('h', self.tsf.Siteloc[i*2+1]))
                fout.write(struct.pack('i', i+1))

            for i in range(self.tsf.n_electrodes):
                print i,
                ec_traces[i].tofile(fout)  #Frontside

            fout.write(struct.pack('i', self.tsf.n_cell_spikes))
            fout.close()
            #quit()

    def lfp_compress(self):
        
        print "...making compressed lfp files ..."

        for file_name in self.filenames:
            
            file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp_compressed.tsf'
            if os.path.exists(file_out)==True: continue

            #Load low-pass .tsf file
            file_in = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp.tsf'
            print "Processing: \n", file_in
            self.load_tsf(file_in)
            print self.tsf.Siteloc.shape

            #Save .lfp.zip file; Martin format; still used by some routines (may wish to remove eventually)
            out_file = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'.lfp.zip'
            print "Saving LFP : ", out_file
            t0 = 0
            t1 = len(self.tsf.ec_traces[0])*1E6/self.tsf.SampleFrequency  #time of end of file in usec 
            tres = 1000     #1Khz sample rate
            uVperAD = 1.0
            chans = np.arange(0, self.tsf.n_electrodes, 1)
            chanpos = self.tsf.Siteloc
    
            np.savez(out_file, chans=chans, chanpos=chanpos, data=self.tsf.ec_traces, t0=t0, t1=t1, tres=tres, uVperAD=uVperAD)
            os.rename(out_file+'.npz', out_file)
            
            #COMPRESS LFP 
            #PAD DATA - Required for Nick/Martin recs as LFP starts bit after highpass data;
            #print "Loading LFP (.lfp.zip) file: ", file_lfp
            #data_in  = np.load(file_lfp+'.lfp.zip')
            #t_start = data_in['t0']*1E-6      #Convert to seconds
            #t_end = data_in['t1']*1E-6        #Convert to seconds
            #start_offset = np.zeros((len(data_in['data']), int(t_start*1E3)), dtype=np.int16)  #make padding at front of data
            #data_out = np.concatenate((start_offset, data_out), axis=1)

            #*********** SAVE COMPRESSED LOW PASS .TSF FILE *********
            header = 'Test spike file '
            iformat = 1002
            n_electrodes = self.tsf.n_electrodes
            SampleFrequency = 50000
            vscale_HP = self.tsf.vscale_HP #use the same as above
            #n_vd_samples = len(self.tsf.ec_traces[0][::Compress_factor])
            n_vd_samples = len(self.tsf.ec_traces[0])

            f1 = open(file_out, 'wb')
            f1.write(header)
            f1.write(struct.pack('i', iformat))
            f1.write(struct.pack('i', SampleFrequency))
            f1.write(struct.pack('i', n_electrodes))
            f1.write(struct.pack('i', n_vd_samples))
            f1.write(struct.pack('f', vscale_HP))
            for i in range (n_electrodes):
                f1.write(struct.pack('h', self.tsf.Siteloc[i*2]))
                f1.write(struct.pack('h', self.tsf.Siteloc[i*2+1]))
                f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

            print "Writing data"
            for i in range(n_electrodes):
                self.tsf.ec_traces[i].tofile(f1)

            f1.write(struct.pack('i', 0)) #Write # of fake spikes
            #text = "Compression:" + str(overall_compression)
            #n_bytes = len(text)
            #f1.write(struct.pack('i', n_bytes))
            #f1.write(text)
            f1.close()

    def write_tsf(self):
        pass



class Probe(object):      

    def __init__(self):

        print "...loading probe..."

        self.name = "NeuroNexus 64Ch probe"         #Hardwired, but should add options here...
    
        self.load_layout()
    
    def load_layout(self):
        ''' Load intan probe map layout 
        '''

        self.n_electrodes = 64

        #Fixed location array for NeurNexus probe layotus
        self.Siteloc = np.zeros((self.n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (self.n_electrodes):
            self.Siteloc[i][0]=30*(i%2)
            self.Siteloc[i][1]=i*23


        #A64 Omnetics adaptor
        adaptor_map = []
        adaptor_map.append([34,35,62,33,60,54,57,55,10,8,11,5,32,3,30,31])
        adaptor_map.append([64,58,63,56,61,59,52,50,15,13,6,4,9,2,7,1])
        adaptor_map.append([53,51,49,47,45,36,37,38,27,28,29,20,18,16,14,12])
        adaptor_map.append([48,46,44,42,40,39,43,41,24,22,26,25,23,21,19,17])

        adaptor_layout1=[]      #Concatenated rows
        for maps in adaptor_map:
            adaptor_layout1.extend(maps)

        #Intan adapter - if inserted right-side up
        intan_map = []
        intan_map.append(list(reversed([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])))     #NB: need to reverse these arrays:  list(reversed(...))
        intan_map.append(list(reversed([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])))
        intan_map.append(list(reversed([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])))
        intan_map.append(list(reversed([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])))

        intan_layout1=[]
        for maps in intan_map:
            intan_layout1.extend(maps)


        #Intan adapter - if inserted upside-down; no need to reverse
        intan_map = []
        intan_map.append([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])
        intan_map.append([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])
        intan_map.append([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])
        intan_map.append([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])

        intan_layout2=[]
        for maps in intan_map:
            intan_layout2.extend(maps)


        #A1x64 probe layout
        a = [27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,28,29,30,31,32]
        b = [37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,36,35,34,33]
        probe_map = a+b
        probe_map[::2] = a
        probe_map[1::2] = b
        
        self.layout = []
        for i in range(len(probe_map)):
            self.layout.append(intan_layout1[adaptor_layout1.index(probe_map[i])])



def find_previous(self, array, value):
    temp = (np.abs(array-value)).argmin()
    if array[temp]>value: return temp-1
    else: return temp

      
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y
       
