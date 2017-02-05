import numpy as np
import os
import glob
import struct

class Tsf_file(object):

    def __init__(self, file_name):
        
        self.read_header(file_name)
        
    def read_header(self, file_name):
        
        self.fin = open(file_name, "rb")
        
        self.header = self.fin.read(16)
        self.iformat = struct.unpack('i',self.fin.read(4))[0] 
        self.SampleFrequency = struct.unpack('i',self.fin.read(4))[0] 
        self.n_electrodes = struct.unpack('i',self.fin.read(4))[0] 
        self.n_vd_samples = struct.unpack('i',self.fin.read(4))[0] 
        self.vscale_HP = struct.unpack('f',self.fin.read(4))[0] 

        if self.iformat==1001:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Siteloc = struct.unpack(str(2*self.n_electrodes)+'h', self.fin.read(2*self.n_electrodes*2))
        if self.iformat==1002:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
            for i in range(self.n_electrodes):
                self.Siteloc[i*2] = struct.unpack('h', self.fin.read(2))[0]
                self.Siteloc[i*2+1] = struct.unpack('h', self.fin.read(2))[0]
                self.Readloc[i] = struct.unpack('i', self.fin.read(4))[0]

    def read_ec_traces(self):
        print " ... reading data, #chs: ", self.n_electrodes, " nsamples: ", self.n_vd_samples, " len: ", float(self.n_vd_samples)/float(self.SampleFrequency), " sec."
        self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        self.ec_traces.shape = self.n_electrodes, self.n_vd_samples

        self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 
        
        #print "No. ground truth cell spikes: ", self.n_cell_spikes
        if (self.n_cell_spikes>0):
            if (self.iformat==1001):
                self.vertical_site_spacing = struct.unpack('i',self.fin.read(4))[0] 
                self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 

            self.fake_spike_times =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_assignment =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_channels =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
        
        #self.fin.close()

    def read_trace(self, channel):
        #Load single channel 

        indent = 16+20+self.n_electrodes*8

        self.fin.seek(indent+channel*2*self.n_vd_samples, os.SEEK_SET)         #Not 100% sure this indent is correct.
        self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_vd_samples)
        
        #self.fin.close()
    
    def save_tsf(self, file_name):
        
        fout = open(file_name, 'wb')
        print "...saving: ",  file_name
        fout.write(self.header)
        fout.write(struct.pack('i', self.iformat))
        fout.write(struct.pack('i', self.SampleFrequency))
        fout.write(struct.pack('i', self.n_electrodes))
        fout.write(struct.pack('i', self.n_vd_samples))
        fout.write(struct.pack('f', self.vscale_HP))
        
        for i in range (self.n_electrodes):
            fout.write(struct.pack('h', self.Siteloc[i*2]))
            fout.write(struct.pack('h', self.Siteloc[i*2+1]))
            fout.write(struct.pack('i', i+1))                 #CAREFUL, SOME FILES MAY USE ReadLoc values..

        self.ec_traces.tofile(fout)

        fout.write(struct.pack('i', self.n_cell_spikes))

        #try:
            #self.subsample
        #except NameError:
            #self.subsample = 1.0

        #fout.write(struct.pack('i', self.subsample))

        fout.close()
        
#DUPLICATE FUNCTION WITH TSF CLASS FUNCTION; May still need it for stand alone functions; but LIKELY OBSOLETE... ERASE!!!!!!!!!!!!
def save_tsf_single(tsf, file_name):
    
    fout = open(file_name, 'wb')
    print file_name
    fout.write(tsf.header)
    fout.write(struct.pack('i', tsf.iformat))
    fout.write(struct.pack('i', tsf.SampleFrequency))
    fout.write(struct.pack('i', tsf.n_electrodes))
    fout.write(struct.pack('i', tsf.n_vd_samples))
    fout.write(struct.pack('f', tsf.vscale_HP))
    for i in range (tsf.n_electrodes):
        fout.write(struct.pack('h', tsf.Siteloc[i*2]))
        fout.write(struct.pack('h', tsf.Siteloc[i*2+1]))
        fout.write(struct.pack('i', i+1))                 #CAREFUL, SOME FILES MAY USE ReadLoc values..

    tsf.ec_traces.tofile(fout)

    fout.write(struct.pack('i', tsf.n_cell_spikes))
    fout.close()

        
#**********************************************************
filename = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/ephys/hp_wavelet/track_1_hp_wavelet_alltrack.tsf'

tsf = Tsf_file(filename)

tsf.temp_traces = []
for channel in range(29,64):
    print "...reading ch: ", channel
    tsf.read_trace(channel)
    tsf.temp_traces.append(tsf.ec_traces)

tsf.ec_traces = np.array(tsf.temp_traces)
tsf.n_electrodes = 10
tsf.n_cell_spikes = 0

save_tsf_single(tsf, filename[:-4]+"_bottom_35ch.tsf")





