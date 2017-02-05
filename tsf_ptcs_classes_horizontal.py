import struct, array, csv
import numpy as np
import os
import math, operator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
import time
from matplotlib_venn import * #venn2, venn2_circles, venn3, venn3_circles

colors = ['blue','red', 'black']

class Tsf_file(object):

    def __init__(self, fin, sim_dir):

        self.read_tsf(fin,sim_dir)

    def read_tsf(self,fin,sim_dir):

        self.header = fin.read(16)
        self.iformat = struct.unpack('i',fin.read(4))[0] 
        self.SampleFrequency = struct.unpack('i',fin.read(4))[0] 
        self.n_electrodes = struct.unpack('i',fin.read(4))[0] 
        self.n_vd_samples = struct.unpack('i',fin.read(4))[0] 
        self.vscale_HP = struct.unpack('f',fin.read(4))[0] 
        print "iformat: ", self.iformat

        if self.iformat==1001:
            self.Siteloc = np.zeros((2*56), dtype=np.int16)
            self.Siteloc = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
        if self.iformat==1002:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
            for i in range(self.n_electrodes):
                self.Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
                self.Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
                self.Readloc[i] = struct.unpack('i', fin.read(4))[0]

        self.ec_traces =  np.fromfile(fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        self.ec_traces.shape = self.n_electrodes, self.n_vd_samples
        self.ec_distribution = self.ec_traces

        self.n_cell_spikes = struct.unpack('i',fin.read(4))[0] 
        print "No. cell spikes: ", self.n_cell_spikes
        if (self.n_cell_spikes>0):
            if (self.iformat==1001):
                self.vertical_site_spacing = struct.unpack('i',fin.read(4))[0] 
                self.n_cell_spikes = struct.unpack('i',fin.read(4))[0] 

            self.fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)

        #print self.n_cell_spikes

        if self.n_cell_spikes>0:
            #print self.fake_spike_assignment
            self.n_cells=max(self.fake_spike_assignment)
            print "No. of ground truth cell: ",self.n_cells
            self.cell_spikes=[[] for x in range(self.n_cells)]  #make lists of cells; .tsf files don't store #cells; could just read file twice;
            self.cell_ptp=[[] for x in range(self.n_cells)] 
            self.cell_size=[[] for x in range(self.n_cells)] 
            self.cell_maxchan=[[] for x in range(self.n_cells)] 
            self.cell_ptpsd=[[] for x in range(self.n_cells)] 

            for k in range(len(self.fake_spike_times)):
                self.cell_spikes[self.fake_spike_assignment[k]-1].append(self.fake_spike_times[k])
           
            ##Compute PTP for ground truth data
            if (os.path.exists(sim_dir+'/cell_ptp.csv')==False):
                for i in range(self.n_cells):
                    #print "Computing PTP cell: ", i
                    self.compute_ptp(self.cell_spikes[i], self.n_electrodes, self.ec_traces, self.SampleFrequency, self.vscale_HP)
                    self.cell_ptp[i]=self.ptp
                    self.cell_size[i]=len(self.cell_spikes[i])
                    self.cell_maxchan[i]=self.maxchan
                    self.cell_ptpsd[i]=self.ptpsd
                    print "Cell: ", i , " ptp: ", self.cell_ptp[i], "ptp_sd: ", self.cell_ptpsd[i], " size: ", self.cell_size[i], " max ch: ", self.cell_maxchan[i]

                np.savetxt(sim_dir+'/cell_ptp.csv', self.cell_ptp, delimiter=",")
                np.savetxt(sim_dir+'/cell_ptpsd.csv', self.cell_ptpsd, delimiter=",")

    def compute_ptp(self, unit, n_electrodes, ec_traces, SampleFrequency, vscale_HP):
        #print unit[0:10]
        #print ec_traces[0][unit[0:10]]
        ptp=[[0] for x in range(n_electrodes)]
        for i in range(n_electrodes):  
            for k in unit:
                #print "Looking up spike: ", k
                ptp[i] += vscale_HP*(max(ec_traces[i][k-10:k+10])-min(ec_traces[i][k-10:k+10]))

        self.ptp = max(ptp)/float(len(unit))
        self.maxchan = np.argmax(ptp)
        
        #Compute SD of ptp on maxchan:
        ptp=np.zeros(len(unit), dtype=np.float32) 
        tempint=0
        for i in [self.maxchan-1]:
            for k in unit:
                ptp[tempint] = vscale_HP*(max(ec_traces[i][k-10:k+10])-min(ec_traces[i][k-10:k+10]))            
                tempint+=1

        self.ptpsd = np.std(ptp)

    def make_DCoffset_tsf(self):

        self.ec_traces
        self.fake_spike_times

        self.ec_traces_DCoffset =  np.zeros((self.n_electrodes, self.n_vd_samples), dtype=np.int16)

        for p in range(len(self.fake_spike_times)):
            print "Inserting spike: ", p, " of total: ", len(self.fake_spike_times)
            k = self.fake_spike_times[p]
            ws=5  #window start
            we=5  #window end
            for j in range(min(100, len(self.fake_spike_times)-p-1)):
                if (self.fake_spike_times[p+j]-k)>ws:
                    #ws+=20  #window start
                    we+=5  #window end

            for i in range(self.n_electrodes):  
                self.ec_traces_DCoffset[i][k-ws:k+we] = self.ec_traces[i][k-ws:k+we] - self.ec_traces[i][k-ws]

        #*********** SAVE .TSF FILE *
        file_name = '/media/cat/Data1/in_silico/12.0sec_poly_3_sum.tsf_truth_old_DCOffset.tsf'
        fout = open(file_name, 'wb')
        fout.write(self.header)
        fout.write(struct.pack('i', self.iformat))
        fout.write(struct.pack('i', self.SampleFrequency))
        fout.write(struct.pack('i', self.n_electrodes))
        fout.write(struct.pack('i', self.n_vd_samples))
        fout.write(struct.pack('f', self.vscale_HP))
        fout.write(struct.pack(str(2*56)+'h', *self.Siteloc))
        self.ec_traces_DCoffset.tofile(fout)
        fout.write(struct.pack('i', self.n_cell_spikes))
        if self.n_cell_spikes>0:
            fout.write(struct.pack('i', self.vertical_site_spacing))
            fout.write(struct.pack('i', self.n_cell_spikes))
            fout.write(struct.pack(str(self.n_cell_spikes)+'i', *self.fake_spike_times))
            fout.write(struct.pack(str(self.n_cell_spikes)+'i', *self.fake_spike_assignment))
            fout.write(struct.pack(str(self.n_cell_spikes)+'i', *self.fake_spike_channels))
        fout.close()


class Loadptcs(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, f, work_dir):
        # call the appropriate method:
        self.VER2FUNC = {1: self.readHeader, 2: self.readHeader, 3: self.readHeader}

        self.readHeader(f)
        self.loadData(self.nsamplebytes, f, work_dir)

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def readHeader(self, f):
        """Read in neuron record of .ptcs file version 3. 'zpos' field was replaced
        by 'sigma' field.
        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
        (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        sigma: float64 (um) (Gaussian spatial sigma)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
        (template waveform data, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
        (template waveform standard deviation, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """

        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.nneurons = int(np.fromfile(f, dtype=np.uint64, count=1)) # nneurons
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        self.nsamplebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsamplebytes
        self.samplerate = int(np.fromfile(f, dtype=np.uint64, count=1)) # samplerate
        self.npttypebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # npttypebytes

        print "self.samplerate: ", self.samplerate

        self.pttype = f.read(self.npttypebytes).rstrip('\0 ') # pttype

        self.nptchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nptchans
        self.chanpos = np.fromfile(f, dtype=np.float64, count=self.nptchans*2) # chanpos
        self.chanpos.shape = self.nptchans, 2 # reshape into rows of (x, y) coords
        self.nsrcfnamebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsrcfnamebytes
        self.srcfname = f.read(self.nsrcfnamebytes).rstrip('\0 ') # srcfname
        # maybe convert this to a proper Python datetime object in the Neuron:
        self.datetime = float(np.fromfile(f, dtype=np.float64, count=1)) # datetime (days)
        self.ndatetimestrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndatetimestrbytes
        self.datetimestr = f.read(self.ndatetimestrbytes).rstrip('\0 ') # datetimestr

    def loadData(self, n_bytes, f, work_dir):
        #pass
        ## call the appropriate method:
        #self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        #self.header = header
        self.nsamplebytes = n_bytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.nsamplebytes]

        self.n_units=self.nneurons
        print self.nneurons
        self.units=[None]*self.n_units
        self.n_sorted_spikes = [None]*self.n_units
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size = []
        self.maxchan = []
        
        print "No. of units sorted: ", self.nneurons

        for k in range(self.n_units):
            print "Loading unit: ", k
            self.readUnit(f)
            self.units[k]= self.spikes
            self.n_sorted_spikes[k] = len(self.units[k])
            self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes to usec
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
            self.ptp[k]=max(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]])- \
            min(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) #compute PTP of template;
        f.close()

        if (os.path.exists(work_dir+'ptps.csv')==False):
            np.savetxt(work_dir+'ptps.csv', self.ptp, delimiter=",")

    def readUnit(self,f):
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr
        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zpos = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) # chanids
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid
        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wavedata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        # spike timestamps (us):
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes)

        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.int64)

    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        # nwavedata/nwavestd bytes, padded:
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1))
        fp = f.tell()
        count = nbytes // self.nsamplebytes # trunc to ignore any pad bytes
        X = np.fromfile(f, dtype=self.wavedtype, count=count) # wavedata/wavestd (uV)
        if nbytes != 0:
            X.shape = self.nchans, self.nt # reshape
        f.seek(fp + nbytes) # skip any pad bytes
        return nbytes, X

    def rstrip(s, strip):
        """What I think str.rstrip should really do"""
        if s.endswith(strip):
            return s[:-len(strip)] # strip it
        else:
            return s

    def read(self):
        self.nid = self.parse_id()
        with open(self.fname, 'rb') as f:
            self.spikes = np.fromfile(f, dtype=np.int64) # spike timestamps (us)
        self.nspikes = len(self.spikes)

    def read_tsf(self,f):
        pass

class Loadcsv(object):

    def __init__(self, fname, tsf, work_dir, csv_version):
        
        self.loadSpikes(fname, tsf, work_dir, csv_version)

    def loadSpikes(self, fname, tsf, work_dir, csv_version):

        if csv_version==2:
            lines=csv.reader(open(fname))
            self.units=[[] for x in range(sum(1 for row in lines))]
            lines=csv.reader(open(fname))

            counter=0
            for row in lines:
                #print row[0]
                #print len(row[1:])
                self.units[counter]=(row[1:])
                self.units[counter]=np.array([float(i) for i in self.units[counter]])
                print type(self.units[counter][0])
                print self.units[counter][0:10]
                counter+=1

            self.n_units=counter
            print "No. of units in .csv: ", self.n_units

            self.n_sorted_spikes = [None]*self.n_units
            for k in range(self.n_units):
                self.n_sorted_spikes[k] = len(self.units[k])

        elif csv_version==1:
            f = open(fname, "r")
            spiketimes = np.genfromtxt(f, delimiter=',', dtype=float32, usecols=(0,))
            f.close()
            f = open(fname, "r")
            spikeids =  np.genfromtxt(f, delimiter=',', dtype=int32, usecols=(1,))

            self.n_sorted_spikes = len(spiketimes)
            self.n_units = max(spikeids)
            self.units=[[] for x in range(self.n_units)]

            for k in range(len(spiketimes)):
                self.units[spikeids[k]-1].extend([spiketimes[k]])

            #print type(self.units[0][0])
            #time.sleep(3)
            self.n_sorted_spikes = [None]*self.n_units
            for k in range(self.n_units):
                self.n_sorted_spikes[k] = len(self.units[k])

            #James spiketimes need to be scaled as below
            #Scale spiketimes; DO THIS IN NUMPY?! couldn't get it to work earlier.
            for k in range(self.n_units):
                self.units[k]=[i * tsf.SampleFrequency for i in self.units[k]] 

        elif csv_version==3:
            f = open(fname, "r")
            spiketimes = np.genfromtxt(f, delimiter=',', dtype=float32, usecols=(0,))
            f.close()
            f = open(fname, "r")
            spikeids =  np.genfromtxt(f, delimiter=',', dtype=int32, usecols=(1,))

            self.n_sorted_spikes = len(spiketimes)
            self.n_units = max(spikeids)
            self.units=[[] for x in range(self.n_units)]

            for k in range(len(spiketimes)):
                self.units[spikeids[k]-1].extend([spiketimes[k]])

            self.n_sorted_spikes = [None]*self.n_units
            for k in range(self.n_units):
                self.n_sorted_spikes[k] = len(self.units[k])

            #James spiketimes need to be scaled as below
            #Scale spiketimes; DO THIS IN NUMPY?! couldn't get it to work earlier.
            for k in range(self.n_units):
                self.units[k]=[i * tsf.SampleFrequency*1E-9 for i in self.units[k]] 

        print self.units[0]

        #Compute PTP values from original .tsf file; needed for .csv sorted files
        print "Computing PTP" #unit is an array of spikes
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size=np.zeros((self.n_units), dtype=np.float32)
        for i in range(self.n_units):
            self.size[i]=len(self.units[i])

        self.maxchan=np.zeros((self.n_units), dtype=np.float32)

        if (os.path.exists(work_dir+'ptps.csv')==False):
            for i in range(self.n_units):
                print "Computing ptp for unit: ", i
                tsf.compute_ptp(self.units[i], tsf.n_electrodes, tsf.ec_traces, tsf.SampleFrequency, tsf.vscale_HP)
                self.ptp[i]=tsf.ptp
                self.maxchan[i]=tsf.maxchan
                print "Unit: ", i , " ptp: ", self.ptp[i], " size: ", self.size[i], " max ch: ", tsf.maxchan

            np.savetxt(work_dir+'ptps.csv', self.ptp, delimiter=",")

def Check_sort(tsf, Sort1, sim_dir):

    if (os.path.exists(sim_dir+'/purity.csv')==False):
        #Look for best match w/in 10 timesteps; then re-run matcher just for particular cell and unit w. 2-3ms window in case alignemnt is wrong
        dist=2.0 * tsf.SampleFrequency*1.0E-3 #time difference in 100 usec blocks allowable between units and cells

        #NB: dist is a tricky parameter; for in-vitro datasets where at least 5-10ms between spikes it's ok to use 3ms, but for other data WILL BE PROBLEM!

        n_cells=tsf.n_cells
        sort = np.zeros((Sort1.n_units, tsf.n_cells), dtype=np.int32)
        duplicate=np.zeros((Sort1.n_units, tsf.n_cells), dtype=np.int32)
        duplicate_mainmatch=np.zeros(Sort1.n_units, dtype=np.int32)

        sort_flags=[]
        for i in range(tsf.n_cells):
            sort_flags.append([0]*len(tsf.cell_spikes[i]))

        spike_membership=[]
        for i in range(Sort1.n_units):
            spike_membership.append([999]*len(Sort1.units[i]))   #Keep track of the membership of every single spike;
                                                            #999 = noise
                                                            #other values = cell no. to which spike belongs

        purity=[]
        completeness=[]
        maxchan=[]
        sortedspikes_fromcells=[]   #Keeps track of spikes in units that are from other cells 
                                    #(i.e not from main unit and not noise
        bestmatch=[]

        for k in range(Sort1.n_units): #ptcs_file.nneurons):
            print ""
            print "**************************************************************"
            print "Checking Unit: ", k, " # spikes: ", Sort1.n_sorted_spikes[k]
            for j in range(tsf.n_cells):
                for p in range(Sort1.n_sorted_spikes[k]):
                    #print "Looking for spike: ", Sort1.units[k][p]
                    #print tsf.cell_spikes[j][0:50]
                    idx=find_nearest(tsf.cell_spikes[j], Sort1.units[k][p], dist)
                    if (idx<>None):
                        if (sort_flags[j][idx]==0):     #Ensure spike hasn't already been detected
                            sort_flags[j][idx]=1        #Set flag to 0 to ensure only single detection
                            sort[k][j]+=1
                            spike_membership[k][p]=j
                        else:
                            duplicate[k][j]+=1
                    else:
                        pass
                    #time.sleep(1)

                print "Unit: ", k, "(",Sort1.n_sorted_spikes[k],")  vs. Cell: ", j, "(", len(tsf.cell_spikes[j]), ") matches : ", sort[k][j], " duplicates: ", duplicate[k][j]
            print "Matched spikes: ", sort[k]
            print "Duplicates    : ", duplicate[k]

            #Rerun comparison routine just for highest match; use 2 x size of window in case overlapping units... (Q: in-vitro data no temporal overlap)
            #This is important for spikes that will be temporally overlapping as might be case in silico; 
            ind = np.argmax(np.array(sort[k],dtype='int')) #Get index of cell with best match to unit
            print "Rechecking unit: ", k, " size: ", Sort1.n_sorted_spikes[k], " with max match cell: ", ind
            bestmatch.append(ind)
            print "1st iteration # spikes match: ", max(sort[k])
            sort[k][ind]=0
            sort_flag2=np.zeros(len(tsf.cell_spikes[ind]), dtype=int32)  #Set flags back to zero

            for p in range(Sort1.n_sorted_spikes[k]):
                idx=find_nearest(tsf.cell_spikes[ind], Sort1.units[k][p], dist*2)
                if (idx<>None):
                    #print idx, Sort2.n_sorted_spikes[k]
                    if (sort_flag2[idx]==0): #Ensure spike hasn't already been detected
                        sort_flag2[idx]=1   #Set flag to ensure single detection
                        sort[k][ind]+=1
                    else:
                        #NOT COMPLETELY CORRECT HERE; ON RECHECK CAN FIND ADDITIONAL DUPLICATES... MUST KEEP TRACK
                        duplicate_mainmatch[k]+=1
                else:
                    pass
            print "2nd iteration # spikes match: ", sort[k][ind] , " duplicate_mainmatch: ", duplicate_mainmatch[k]
            duplicate[k][ind]=duplicate_mainmatch[k]
            purity.append(float(sort[k][ind]-duplicate[k][ind])/float(Sort1.n_sorted_spikes[k])*100) #max # spikes matched/Total # spikes in unit
            completeness.append(float(sort[k][ind])/len(tsf.cell_spikes[ind])*100) #max # spikes matched/Total # spikes in closest cell
            sortedspikes_fromcells.append(sum(sort[k]))

            print "Total spikes coming from all cells: ", sortedspikes_fromcells[k]
            print "Unit %d, PTP: %.1f uV, Cell_match: %d, purity: %.1f %%, comp.: %.1f %%" \
            %(k, Sort1.ptp[k], np.argmax(sort[k])+1, purity[k], completeness[k])
            #time.sleep(5)

        print "Bestmatch array: ", bestmatch
        print sort
        sort-=duplicate
        print sort

        #np.savetxt(sim_dir+'ptps.csv', ptps, delimiter=",")
        np.savetxt(sim_dir+'purity.csv', purity, delimiter=",")
        np.savetxt(sim_dir+'completeness.csv', completeness, delimiter=",")
        np.savetxt(sim_dir+'size.csv', Sort1.n_sorted_spikes, delimiter=",")
        np.savetxt(sim_dir+'maxchan.csv', Sort1.maxchan, delimiter=",")
        np.savetxt(sim_dir+'sortedspikes_fromcells.csv', sortedspikes_fromcells, delimiter=",")
        np.savetxt(sim_dir+'bestmatch.csv', bestmatch, delimiter=",")
        np.savetxt(sim_dir+'full_sort.csv', sort, delimiter=",")
        np.savetxt(sim_dir+'duplicate.csv', duplicate, delimiter=",")
        np.savetxt(sim_dir+'duplicate_mainmatch.csv', duplicate_mainmatch, delimiter=",")

        with open(sim_dir + 'spike_membership.csv', "w") as f:
            writer = csv.writer(f)
            writer.writerows(spike_membership)

    else:
        print "Skipping detection (loading previous ground truth check data)"


def Compare_sorts(Sort1, Sort2, tsf):

    if (os.path.exists(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv')==False):
        #Look for best match w/in 10 timesteps; then re-run matcher just for particular cell and unit w. 2-3ms window in case alignemnt is wrong
        dt=2.0*float(tsf.SampleFrequency)*1.0E-3   #time difference in ms allowable between units and cells
        distance=70.0   #distance in um allowable between units:

        sort = np.zeros((Sort1.n_units, Sort2.n_units), dtype=np.int32)
        sort_flags = []
        for i in range(Sort1.n_units):
            sort_flags.append([0]*Sort1.n_sorted_spikes[i])

        sortedspikes_fromcells=[None]*Sort1.n_units     #Keeps track of spikes in units that are from other cells 
                                                        #(i.e not from main unit and not noise
        comparematch = [] #Should keep flags '999' if there is no closest match - very rare occassions;
        duplicate = np.zeros(Sort1.n_units, dtype=np.int32)

        for k in range(Sort1.n_units): 
            print ""
            print "**************************************************************"
            print "Checking Unit: ", k, " # spikes: ", Sort1.n_sorted_spikes[k]
            for j in range(Sort2.n_units):
                #Cat: use distance between max channel for compared units; DOESN"T WORK FOR .CSV FILES; ASK FOR CHANNEL LOCATION IN FUTURE
                #if (sqrt((Sort1.chanpos[Sort1.maxchan[k]][0]-Sort2.chanpos[Sort2.maxchan[j]][0])**2+
                    #(Sort1.chanpos[Sort1.maxchan[k]][1]-Sort2.chanpos[Sort2.maxchan[j]][1])**2))<distance:
                    for p in range(Sort1.n_sorted_spikes[k]):
                        idx=find_nearest(Sort2.units[j], Sort1.units[k][p], dt)
                        if (idx<>None):
                            if (sort_flags[k][p]==0): #Ensure spike hasn't already been detected
                                sort_flags[k][p]=1   #Set flag to 0 to ensure only single detection
                                sort[k][j]+=1
                            else:
                                duplicate[k]+=1
                                #print "****************DUPLICATE SPIKE******************"
                        else:
                            pass
                    print "Vs. Unit: ", j, " # spikes: ", len(Sort2.units[j]), " matches : ", sort[k][j], " duplicates: ", duplicate[k]
            print sort[k]
            sortedspikes_fromcells.append(sum(sort[k]))
            print "Total spikes matched : ", sortedspikes_fromcells[k]
            
            #Rerun comparison routine just for highest match; use 2 x size of window in case overlapping units... (Q: in-vitro data no temporal overlap)
            #This is important for spikes that will be temporally overlapping as might be case in silico; 
            ind = np.argmax(np.array(sort[k],dtype='int')) #Get index of cell with best match to unit
            print "Rechecking unit: ", k, " size: ", Sort1.n_sorted_spikes[k], " with max match unit: ", ind
            comparematch.append(ind)
            print "1st iteration # spikes match: ", max(sort[k])
            sort[k][ind]=0 
            sort_flag2=np.zeros(Sort2.n_sorted_spikes[ind], dtype=int32)  #Set flags back to zero

            for p in range(Sort1.n_sorted_spikes[k]):
                idx=find_nearest(Sort2.units[ind], Sort1.units[k][p], dt*2)
                if (idx<>None):
                    if (sort_flag2[idx]==0): #Ensure spike hasn't already been detected
                        sort_flag2[idx]=1   #Set flag to ensure single detection
                        sort[k][ind]+=1
                    #else:
                        #NOT COMPLETELY CORRECT HERE; ON RECHECK CAN FIND ADDITIONAL DUPLICATES... MUST KEEP TRACK
                else:
                    pass

            print "2nd iteration # spikes match: ", sort[k][ind]

        print "Bestmatch array: ", comparematch
        print "Compare array: ", sort

        np.savetxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', comparematch, delimiter=",")
        np.savetxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', sort, delimiter=",")
        np.savetxt(Sort1.directory+'/duplicate_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', duplicate, delimiter=",")
        with open(Sort1.directory+'/sort_flags_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', "w") as f:
            writer = csv.writer(f)
            writer.writerows(sort_flags)        
    else:

        print "Skipping detection (loading previous unit compare data)"

def Compare_sorts_plot(Sort1, Sort2, tsf):

    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = np.array(Sort1.size) #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = np.array(Sort2.size) #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.suptitle('Compare Analysis for: \n"' + Sort1.filename[0:8] + '"  AND  "' + 
    Sort2.filename[0:8] + '" (in vitro)', fontsize = 16)

    scale=1

    #*********************** FIRST PLOT *****************************
    ax = plt.subplot(1, 1, 1)
    title(Sort1.filename[0:8] + ' vs. ' + Sort2.filename[0:8] , fontsize=15)

    #plt.ylabel('% of Total Spikes in Recording', fontsize=17)

    colors = ['blue','red', 'black']
    plt.xlabel('Average PTP (uV)',fontsize=15)
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = sorted(x)
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))] #Keeps track of ids sorted by PTP value; use this in place of index

    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=i*scale+scale/10. #comparesort[i][comparematch[i]]/size[i]*100

    a=np.array(a)*scale
    #y=a
    plt.scatter(a,y, s=size, alpha=.35, color='blue')

    ss = []
    for i in range(len(ptps)):
        ss.append(size2[int(comparematch[i])])
    plt.scatter(a,y, s=ss, alpha=.35, color='red')

    plt.xlim(0,max(a)+scale/10.)
    plt.ylim(0,max(y)+scale/10.)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()


    #***********************************************************************************************************
    #***************************************** SECOND PLOT - ORDERED PLOTS**************************************
    #***********************************************************************************************************

    ax = plt.subplot(1, 2, 1)
    title(Sort1.filename[0:8] + ' vs. ' + Sort2.filename[0:8] , fontsize=15)

    #plt.ylabel('% of Total Spikes in Recording', fontsize=17)

    colors = ['blue','red', 'black']
    plt.xlabel('Average PTP (uV)',fontsize=15)
      
    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort1.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort2.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    title(Sort1.filename[0:8] + ' vs. ' + Sort2.filename[0:8] , fontsize=15)

    #plt.ylabel('% of Total Spikes in Recording', fontsize=17)

    colors = ['blue','red', 'black']
    plt.xlabel('Average PTP (uV)',fontsize=15)
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size)
    
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=a[i] #tracker.index(i)*scale*10+scale/10.
    #y=a

    a=np.array(a)*scale
    plt.scatter(a,y, s=size/2, alpha=.35, color='blue')
    
    sort1_x=a #Identical to y
    sort1_y=y
    sort1_comparematch=comparematch
    sort1_comparesort=comparesort
    sort1_size = size
    scale2 = max(ptps)

    #************************* SECOND PLOT NO SCALING *****************************
    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort2.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort2.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort1.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.xlabel('Average PTP (uV)',fontsize=15)
    plt.ylabel('Average PTP (uV)',fontsize=15)
    #ax.set_yticks([])    
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size) 
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=a[i] #tracker.index(i)*scale*10+scale/10. 

    a=np.array(a)*scale
    plt.scatter(a+int(math.ceil(scale2/50.0))*50,y, s=size/2, alpha=.35, color='red')

    for i in range(len(sort1_comparesort)): #Loop over all units
        for j in range(len(sort1_comparesort[i])): #Loop over each match in each unit
            if sort1_comparesort[i][j]>0:
                xx = (sort1_x[i], a[j]+int(math.ceil(scale2/50.0))*50)
                yy = (sort1_y[i], y[j])
                alpha_val = (min(sort1_comparesort[i][j]/sort1_size[i]+.05,1.0))
                plt.plot(xx,yy, 'r-', color='black',linewidth=2, alpha=alpha_val)

    #Labels and limits
    plt.xlim(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50.0))*50)
    label_x = np.arange(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50))*50,50)
    my_xticks1 = np.arange(0,int(math.ceil(scale2/50))*50,50).astype('str')
    my_xticks2 = np.arange(0,int(math.ceil(max(ptps)/50))*50, 50).astype('str')
    my_xticks = np.concatenate((my_xticks1,my_xticks2))

    plt.xticks(label_x, my_xticks)
    plt.ylim(bottom=0)

    for i in range(len(label_x)):
        xx=[i*50,i*50]
        yy=[0,10000]
        plt.plot(xx,yy, 'r--', color='black',linewidth=2, alpha=0.5)

    xx=[int(math.ceil(scale2/50))*50,int(math.ceil(scale2/50))*50]
    yy=[0,1000]
    plt.plot(xx,yy,'r--', color='black',linewidth=3)

    #***********************************************************************************************************
    #****************************** NODE Comparisons - PTP BASED ************************************
    #***********************************************************************************************************

    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort1.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort2.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    ax = plt.subplot(1, 2, 2)
    title(Sort1.filename[0:8] + ' vs. ' + Sort2.filename[0:8] , fontsize=15)
    ax.set_yticklabels([])
    
    #plt.ylabel('% of Total Spikes in Recording', fontsize=17)

    colors = ['blue','red', 'black']
    plt.xlabel('Ordered (PTP, no units)',fontsize=15)
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size)
    
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=tracker.index(i)*scale*50+100
    #y=a

    a=np.array(a)*scale
    plt.scatter(a,y, s=size/2, alpha=.35, color='blue')
    
    sort1_x=a #Identical to y
    sort1_y=y
    sort1_comparematch=comparematch
    sort1_comparesort=comparesort
    sort1_size = size
    scale2 = max(ptps)

    #************************* SWITCH PLOTS *****************************
    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort2.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort2.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort1.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.xlabel('Average PTP (uV)',fontsize=15)
    plt.ylabel('Ordered by PTP (no units)',fontsize=15)
    #ax.set_yticks([])    
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size) 
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=tracker.index(i)*scale*50+100

    a=np.array(a)*scale
    plt.scatter(a+int(math.ceil(scale2/50.0))*50,y, s=size/2, alpha=.35, color='red')

    for i in range(len(sort1_comparesort)): #Loop over all units
        for j in range(len(sort1_comparesort[i])): #Loop over each match in each unit
            if sort1_comparesort[i][j]>0:
                xx = (sort1_x[i], a[j]+int(math.ceil(scale2/50.0))*50)
                yy = (sort1_y[i], y[j])
                alpha_val = (min(sort1_comparesort[i][j]/sort1_size[i]+.05,1.0))
                plt.plot(xx,yy, 'r-', color='black',linewidth=2, alpha=alpha_val)

    #Labels and limits
    plt.xlim(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50))*50)
    label_x = np.arange(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50))*50,50)
    my_xticks1 = np.arange(0,int(math.ceil(scale2/50))*50,50).astype('str')
    my_xticks2 = np.arange(0,int(math.ceil(max(ptps)/50))*50, 50).astype('str')
    my_xticks = np.concatenate((my_xticks1,my_xticks2))

    plt.xticks(label_x, my_xticks)
    plt.ylim(bottom=0)

    for i in range(len(label_x)):
        xx=[i*50,i*50]
        yy=[0,10000]
        plt.plot(xx,yy, 'r--', color='black',linewidth=2, alpha=0.5)

    xx=[int(math.ceil(scale2/50))*50,int(math.ceil(scale2/50))*50]
    yy=[0,1000]
    plt.plot(xx,yy,'r--', color='black',linewidth=3)

    #************************************************

    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = np.array(Sort1.size) #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = np.array(Sort2.size) #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.suptitle('Compare Analysis for: \n"' + Sort1.filename[0:8] + '"  AND  "' + 
    Sort2.filename[0:8] + '" : '+tsf.tsf_name, fontsize = 16)

    #VEN-LIKE DIAGRAMS:
    #Need 2 of them still; x-axis = PTP of Sort1, y-axis = % matching of best unit i.e. comparematch.csv file
    #Find radical line and centres of circles to be plotted; 
    #Combine this for up to 3-4 different matching units: i.e. circles should have bubles on them;

    scale=1

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()


    #***********************************************************************************************************
    #********************************** NODE Comparisons - PIE CHARTS ******************************************
    #***********************************************************************************************************

    comparematch = np.genfromtxt(Sort1.directory+'/comparematch_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort1.directory+'/comparesort_vs_'+Sort2.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort1.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort2.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    ax = plt.subplot(1, 1, 1)
    plt.suptitle('Compare Analysis for: \n"' + Sort1.filename[0:8] + '"  &  "' + 
    Sort2.filename[0:8] + '" : '+tsf.tsf_name, fontsize = 16)

    if Sort1.flag == 1:
        title(Sort1.filename[0:8] + ' vs. ' + Sort2.filename[0:8] + ' - Blind Sort ', fontsize=15)
        #ax.set_yticklabels([])
    else:
        title('Cell Completeness (Blue) & Unit Purity (Green) ', fontsize=15)
        #ax.set_yticklabels([])
    
    colors = ['blue','red', 'black']
    plt.xlabel('Ordered (PTP, no units)',fontsize=15)
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size)
    
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=a[i] #tracker.index(i)

    a=sorted(a)*scale

    a = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        a[i]=tracker.index(i)*scale*5+10 #tracker[i]*10 #comparesort[i][comparematch[i]]/size[i]*100

    #Cat: must always use "truth" in Sort1 - otherwise must make the loop below conditional also (as the one in the next block)
    colors = ['blue','red', 'black']
    if Sort1.flag == 1:
        colors = ['cyan','red', 'black']

    y_scaling=1
    if Sort2.n_units>Sort1.n_units:
        y_scaling = Sort2.n_units/Sort1.n_units
    for i in range(len(a)):
        #Cat: Not enough to find maximum match for each unit/cell;
        #but need to ensure that the matching unit/cell represents the current unit
        #and doesn't have other - larger unit/cells - in it

        itemindex = int(comparematch[i])
        tempint=0
        for j in range(len(comparematch)):
            if tempint < comparesort[j][itemindex]:
                tempint = comparesort[j][itemindex]

        if tempint<=comparesort[i][itemindex]:
            a1=float(comparesort[i][itemindex]/size[i])
        else: 
            a1=0
        print ""

        b1=float(sum(comparesort[i])/size[i]-a1)
        c1=1-a1-b1

        #            pie slices,   x loc       y loc   
        draw_pie(ax,[a1, b1, c1], a[i], y[i], size=size[i]/8, colors=colors)

    #plt.scatter(a,y, s=size/2, alpha=.35, color='blue')

    sort1_x=a #Identical to y
    sort1_y=y
    sort1_comparematch=comparematch
    sort1_comparesort=comparesort
    sort1_size = size
    scale2 = max(ptps)

    #************************* SWITCH PLOTS *****************************
    comparematch = np.genfromtxt(Sort2.directory+'/comparematch_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort2.directory+'/comparesort_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort2.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort2.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort1.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.xlabel('Average PTP (uV)',fontsize=15)
    plt.ylabel('Ordered by PTP (no units)',fontsize=15)
    #ax.set_yticks([])    
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size) 

    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=a[i] #tracker.index(i)*50+100
    #y=a

    a = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        a[i]=tracker.index(i)*scale*5+10 #tracker[i]*10 #comparesort[i][comparematch[i]]/size[i]*100

    colors = ['green','red', 'black']

    y_scaling=1
    if Sort1.n_units>Sort2.n_units:
        y_scaling = float(Sort1.n_units)/float(Sort2.n_units)

    if Sort1.flag==1: 
        colors = ['cyan','red', 'black'] #Blind COMPARISON; DO NOT SUBTRACT AB

        for i in range(len(a)):
            #Cat: Not enough to find maximum match for each unit/cell;
            #but need to ensure that the matching unit/cell represents the current unit
            #and doesn't have other - larger unit/cells - in it

            itemindex = int(comparematch[i])
            tempint=0
            for j in range(len(comparematch)):
                if tempint < comparesort[j][itemindex]:
                    tempint = comparesort[j][itemindex]

            if tempint<=comparesort[i][itemindex]:
                a1=float(comparesort[i][itemindex]/size[i])
            else: 
                a1=0
            print ""

            b1=float(sum(comparesort[i])/size[i]-a1)
            c1=1-a1-b1

            #            pie slices,   x loc       y loc   
            draw_pie(ax,[a1, b1, c1], a[i]+int(math.ceil(scale2/50.0))*50, y[i], size=size[i]/8, colors=colors)

    else:
        for i in range(len(a)):
            a1=float(max(comparesort[i])/size[i])
            b1=float(sum(comparesort[i])-max(comparesort[i]))/size[i]
            c1=1-a1-b1
            #            pie slices,   x loc       y loc   
            draw_pie(ax,[a1, b1, c1], a[i]+int(math.ceil(scale2/50.0))*50, y[i], size=size[i]/8, colors=colors)


    #plt.scatter(a+int(math.ceil(scale2/50.0))*50,y, s=size/2, alpha=.35, color='red')

    for i in range(len(sort1_comparesort)): #Loop over all units
        for j in range(len(sort1_comparesort[i])): #Loop over each match in each unit
            if sort1_comparesort[i][j]>0:
                xx = (sort1_x[i], a[j]+int(math.ceil(scale2/50.0))*50)
                yy = (sort1_y[i], y[j])
                alpha_val = (min(sort1_comparesort[i][j]/sort1_size[i]+.05,1.0))
                plt.plot(xx,yy, 'r-', color='black',linewidth=2, alpha=alpha_val)

    #Labels and limits
    plt.xlim(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50))*50)
    label_x = np.arange(0,int(math.ceil(scale2/50))*50+int(math.ceil(max(ptps)/50))*50,50)
    my_xticks1 = np.arange(0,int(math.ceil(scale2/50))*50,50).astype('str')
    my_xticks2 = np.arange(0,int(math.ceil(max(ptps)/50))*50, 50).astype('str')
    my_xticks = np.concatenate((my_xticks1,my_xticks2))

    plt.xticks(label_x, my_xticks)
    plt.ylim(bottom=0)

    for i in range(len(label_x)):
        xx=[i*50,i*50]
        yy=[0,10000]
        plt.plot(xx,yy, 'r--', color='black',linewidth=2, alpha=0.5)

    xx=[int(math.ceil(scale2/50.0))*50,int(math.ceil(scale2/50.0))*50]
    yy=[0,1000]
    plt.plot(xx,yy,'r--', color='black',linewidth=3)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

    quit()
    #******************RIGHT PLOT******************************

    comparematch = np.genfromtxt(Sort2.directory+'/comparematch_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort2.directory+'/comparesort_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort2.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort2.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort1.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")

    ax = plt.subplot(1, 2, 2)
    if Sort1.flag == 1:
        title(Sort2.filename[0:8] + ' vs. ' + Sort1.filename[0:8] + ' - Blind Sort ', fontsize=15)
        ax.set_yticklabels([])
    else:
        title(Sort2.filename[0:8] + ' vs. ' + Sort1.filename[0:8] + ' - Completeness ', fontsize=15)
        ax.set_yticklabels([])
    
    colors = ['blue','red', 'black']
    plt.xlabel('Ordered (PTP, no units)',fontsize=15)
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size)
    
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=tracker.index(i)*scale*50+100
    #y=a

    a=np.array(a)*scale
    colors = ['green','red', 'black']

    if Sort1.flag==1: 
        colors = ['cyan','red', 'black'] #Blind COMPARISON; DO NOT SUBTRACT AB

        for i in range(len(a)):
            #Cat: Not enough to find maximum match for each unit/cell;
            #but need to ensure that the matching unit/cell represents the current unit
            #and doesn't have other - larger unit/cells - in it

            itemindex = int(comparematch[i])
            tempint=0
            for j in range(len(comparematch)):
                if tempint < comparesort[j][itemindex]:
                    tempint = comparesort[j][itemindex]

            if tempint<=comparesort[i][itemindex]:
                a1=float(comparesort[i][itemindex]/size[i])
            else: 
                a1=0
            print ""

            b1=float(sum(comparesort[i])/size[i]-a1)
            c1=1-a1-b1

            #            pie slices,   x loc       y loc   
            draw_pie(ax,[a1, b1, c1], a[i], y[i], size=size[i]/8, colors=colors)

    else:
        for i in range(len(a)):
            a1=float(max(comparesort[i])/size[i])
            b1=float(sum(comparesort[i])-max(comparesort[i]))/size[i]
            c1=1-a1-b1
            #            pie slices,   x loc       y loc   
            draw_pie(ax,[a1, b1, c1], a[i], y[i], size=size[i]/8, colors=colors)

    #plt.scatter(a,y, s=size/2, alpha=.35, color='blue')
    
    sort1_x=a #Identical to y
    sort1_y=y
    sort1_comparematch=comparematch
    sort1_comparesort=comparesort
    sort1_size = size
    scale2 = max(ptps)

    #************************* SWITCH PLOTS *****************************
    comparematch = np.genfromtxt(Sort2.directory+'/comparematch_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    comparesort = np.genfromtxt(Sort2.directory+'/comparesort_vs_'+Sort1.name+'_'+tsf.tsf_name+'.csv', dtype='float32', delimiter=",")
    ptps = np.genfromtxt(Sort1.directory+'/ptps.csv', dtype='float32', delimiter=",")
    size = Sort1.size #np.genfromtxt(Sort2.directory+'/size.csv', dtype='float32', delimiter=",")
    size2 = Sort2.size #np.genfromtxt(Sort1.directory+'/size.csv', dtype='float32', delimiter=",")

    plt.xlabel('Average PTP (uV)',fontsize=15)
    plt.ylabel('Ordered by PTP (no units)',fontsize=15)
    #ax.set_yticks([])    
    
    #X-axis; SORT ALL DATA BY PTPS ORDER
    x = ptps
    X = x #Use this array to sort all other data by PTP 
    a = x
    tracker = np.arange(0,len(ptps),1)
    tracker = [y for (x,y) in sorted(zip(X,tracker))]

    size = np.array(size) 
    #Y-axis
    y = np.zeros(len(ptps), dtype=np.float32)
    for i in range(len(ptps)):
        y[i]=tracker.index(i)*scale*50+100

    a=np.array(a)*scale


    plt.scatter(a+int(math.ceil(scale2/50.0))*50,y, s=size/2, alpha=.35, color='red')


    for i in range(len(sort1_comparesort)): #Loop over all units
        for j in range(len(sort1_comparesort[i])): #Loop over each match in each unit
            if sort1_comparesort[i][j]>0:
                xx = (sort1_x[i], a[j]+int(math.ceil(scale2/50.0))*50)
                yy = (sort1_y[i], y[j])
                alpha_val = (min(sort1_comparesort[i][j]/sort1_size[i]+.05,1.0))
                plt.plot(xx,yy, 'r-', color='black',linewidth=2, alpha=alpha_val)

    #Labels and limits
    plt.xlim(0,int(math.ceil(scale2/50.0))*50+int(math.ceil(max(ptps)/50.0))*50.0)
    label_x = np.arange(0,int(math.ceil(scale2/50.0))*50+int(math.ceil(max(ptps)/50.0))*50,50)
    my_xticks1 = np.arange(0,int(math.ceil(scale2/50.0))*50,50).astype('str')
    my_xticks2 = np.arange(0,int(math.ceil(max(ptps)/50.0))*50, 50).astype('str')
    my_xticks = np.concatenate((my_xticks1,my_xticks2))

    plt.xticks(label_x, my_xticks)
    plt.ylim(bottom=0)

    for i in range(len(label_x)):
        xx=[i*50,i*50]
        yy=[0,10000]
        plt.plot(xx,yy, 'r--', color='black',linewidth=2, alpha=0.5)

    xx=[int(math.ceil(scale2/50.0))*50,int(math.ceil(scale2/50.0))*50]
    yy=[0,1000]
    plt.plot(xx,yy,'r--', color='black',linewidth=3)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()


    quit()

    #*********************** VLINE SETS OF PLOTS *****************************
    #ax = plt.subplot(1, 2, 1)
    title("Original file: "+ Sort1.rawfile+ '\n'+Sort1.name+ '(Red)           '+Sort2.name+' (Blue)' , 
    multialignment='center',fontsize=15)
    comparematch = np.genfromtxt(Sort1.directory+'/comparematch.csv', dtype='float32', delimiter=",")

    #with open(Sort1.directory+'/sort_flags.csv', "w") as f:

    colors = ['blue','red', 'black']
    
    for i in range(Sort1.n_units):
        print "Plotting: ", i, "maxchan: ",Sort1.maxchan[i] #, Sort1.samplerate, Sort1.units[i][0:3],Sort1.units[i][-3:]

        x = np.array(Sort1.units[i],dtype=float32)*1E-4 #float(Sort1.samplerate)*2.5
        
        ymin=np.zeros(len(Sort1.units[i]))
        ymax=np.zeros(len(Sort1.units[i]))
        ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]
        ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]-5.0

        plt.vlines(x, ymin, ymax, linewidth=.5, color='red')

    for i in range(Sort2.n_units):
        print "Plotting: ", i, "maxchan: ",Sort2.maxchan[i] #, Sort1.samplerate, Sort1.units[i][0:3],Sort1.units[i][-3:]

        x = np.array(Sort2.units[i],dtype=float32)*1E-4 #/float(Sort1.samplerate)*2.5

        ymin=np.zeros(len(Sort2.units[i]))
        ymax=np.zeros(len(Sort2.units[i]))
        ymin+=-Sort2.chanpos[Sort2.maxchan[i]][1]-5.0
        ymax+=-Sort2.chanpos[Sort2.maxchan[i]][1]-10.0

        plt.vlines(x, ymin, ymax, linewidth=0.5, color='blue')

    plt.xlim(left=-0.03)

    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Depth from top of electrode',multialignment='center', fontsize=15)
    #plt.savefig(Sort1.directory+'comparison.pdf', format='pdf', dpi=6000)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

    #******************** PLOT MATCHING UNIT RASTERS *********************

    title("Original file: "+ Sort1.rawfile+ '\n'+Sort1.name+ '(Red)           '+Sort2.name+' (Blue)' , 
    multialignment='center',fontsize=15)

    #comparematch = np.genfromtxt(Sort2.directory+'/comparematch.csv', dtype='float32', delimiter=",")
    #comparesort = np.genfromtxt(Sort2.directory+'/comparesort.csv', dtype='float32', delimiter=",")


    for i in range(Sort1.n_units):
        print "Plotting: ", i, "maxchan: ", Sort1.maxchan[i] 

        x = np.array(Sort1.units[i],dtype=float32)*1E-4 #float(Sort1.samplerate)*2.5
        
        ymin=np.zeros(len(Sort1.units[i]))
        ymax=np.zeros(len(Sort1.units[i]))
        ymin+=-Sort1.chanpos[Sort1.maxchan[i]][1]
        ymax+=-Sort1.chanpos[Sort1.maxchan[i]][1]-5.0

        plt.vlines(x, ymin, ymax, linewidth=.5, color='red')

        #PLOT 2nd UNIT BEST MATCH
        print int(comparematch[i])
        print "Plotting: ", i, "maxchan: ",Sort2.maxchan[int(comparematch[i])] #, Sort1.samplerate, Sort1.units[i][0:3],Sort1.units[i][-3:]

        x = np.array(Sort2.units[int(comparematch[i])],dtype=float32)*1E-4 #/float(Sort1.samplerate)*2.5

        ymin=np.zeros(len(Sort2.units[int(comparematch[i])]))
        ymax=np.zeros(len(Sort2.units[int(comparematch[i])]))
        ymin+=-Sort2.chanpos[Sort2.maxchan[int(comparematch[i])]][1]-5.0
        ymax+=-Sort2.chanpos[Sort2.maxchan[int(comparematch[i])]][1]-10.0

        plt.vlines(x, ymin, ymax, linewidth=0.5, color='blue')
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()


def Check_sort_plot(tsf, Sort1, sim_dir, ptcs_name):

    ptps = np.genfromtxt(sim_dir+ptcs_name+'/ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+ptcs_name+'/purity.csv', dtype='float32', delimiter=",")
    completeness = np.genfromtxt(sim_dir+ptcs_name+'/completeness.csv', dtype='float32', delimiter=",")
    size = np.genfromtxt(sim_dir+ptcs_name+'/size.csv', dtype='float32', delimiter=",")
    maxchan = np.genfromtxt(sim_dir+ptcs_name+'/maxchan.csv', dtype='float32', delimiter=",")
    sd = np.genfromtxt(sim_dir+ptcs_name+'/sd.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+ptcs_name+'/sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    bestmatch = np.genfromtxt(sim_dir+ptcs_name+'/bestmatch.csv', dtype='float32', delimiter=",")
    cell_ptp = np.genfromtxt(sim_dir+'/cell_ptp.csv', dtype='float32', delimiter=",")
    cell_ptpsd = np.genfromtxt(sim_dir+'/cell_ptpsd.csv', dtype='float32', delimiter=",")
    #cell_ptp = np.genfromtxt(sim_dir+'Traces/'+'/cell_ptp.csv', dtype='float32', delimiter=",")
    #cell_ptpsd = np.genfromtxt(sim_dir+'Traces/'+'/cell_ptpsd.csv', dtype='float32', delimiter=",")
    full_sort = np.genfromtxt(sim_dir+ptcs_name+'/full_sort.csv', dtype='float32', delimiter=",")
    duplicate = np.genfromtxt(sim_dir+ptcs_name+'/duplicate.csv', dtype='float32', delimiter=",")

    #*************** Compute some stats ***********************************
    print "Total ground truth spikes: ", tsf.n_cell_spikes
    n_sorted_spikes=np.array(Sort1.n_sorted_spikes)
    purity=np.array(purity)
    print "Total sorted spikes: ", sum(Sort1.n_sorted_spikes)
    print "Total correct sorted spikes: ", int(sum(Sort1.n_sorted_spikes*purity/100))
    tempreal = sum(Sort1.n_sorted_spikes*purity/100)/sum(Sort1.n_sorted_spikes)*100
    print 'Percent sorted correct: %.2f' % tempreal
    tempreal = sum(purity)/len(purity)
    print 'Unit based percent correct: %.2f' % tempreal

    n_cells = tsf.n_cells
    n_units = Sort1.n_units
    cells = tsf.cell_spikes

    #****************  PLOT ROUTINES **************************************

    plt.suptitle('Sort Analysis for: "' + ptcs_name + '" (In-Vitro Dataset; Rodrigo and Costas)', fontsize = 20)

    #*********************** FIRST PLOT *****************************
    ax = plt.subplot(1, 3, 1)
    title('Cells - Spike Detection Rates' , fontsize=20)

    plt.ylim((0,100))
    plt.xlim((0,250))

    colors = ['blue','red', 'black']

    plt.xlabel('Average PTP (uV)',fontsize=17)
    plt.ylabel('% Completeness\n (Total Spikes Detected From a Cell)',multialignment='center', fontsize=17)

    #****************************** COMPUTE OVERPLIT INDEX ************************************
    a_sum=0
    b_sum=0
    c_sum=0
    oversplitx=[]
    oversplity=[]
    for i in range(n_cells): #Loop over cells
        #print "********************************************************************"
        #print "Searching for cell: ", i
        a = 0.0
        b = 0.0
        c = 1.0
        tempint=0 #cumulative spike count for detected spikes NOT assigned to main unit
        tempint3=0 #Keeps track of # spikes of largest matching unit to cell
        tempint2=1000 #Keeps track of index of max spike unit; 1000 = flag for no match;
        for j in range(n_units): #Loop over units
            if(full_sort[j][i]>0): #Check if any detected spikes in unit from cell 'i'
                tempint += full_sort[j][i] #Accumulate all sorted spikes belonging to cell 'i'
                #if full_sort[j][i]>tempint3: #Look for largest unit among all units
                if bestmatch[j]==i:
                    if full_sort[j][i]>tempint3: #Look for best match unit w. largest # cell spikes
                        tempint2 = j #Keep track of index of max spike unit
                        tempint3=full_sort[j][i] #keep track of # cell spikes in largest, best match unit 

        #print "Total spikes in cell: ", len(cells[i])
        #print "Best matching unit: ", tempint2, " # spikes: ", full_sort[tempint2][i]
        #Cell spikes assigned to best matching unit
        if tempint2<1000:
            a = float(full_sort[tempint2][i])/float(len(cells[i]))
        else:
            a = 0
        a_sum+=a
        
        #tempint-=duplicate[i]

        #print "Other sorted spikes (Overplit Index): ", tempint - full_sort[tempint2][i]
        #Cell spikes detected but assigned to non-best matching unit
        b = float(tempint)/float(len(cells[i]))-a
        b_sum+=b

        #print "Remaining Spikes: ", len(cells[i]) - tempint
        c = 1-a-b
        c_sum+=c

        draw_pie(ax,[a, b, c], cell_ptp[i], float(a)*100, size=125, colors=colors)

        oversplity.append(b)
        oversplitx.append(cell_ptp[i])

    #Compute oversplit sums to show dashed line on graph
    oversplit_sums = np.array([oversplitx, oversplity]).T  #Load oversplit data into 2 columns
    oversplit_sums = oversplit_sums[oversplit_sums[:,0].argsort()]  #Sort data by PTP amplitude (1st column)
    #print "Sorted: ", oversplit_sums

    oversplity_sum=[] #Cumulatively add all oversplit values backwards;
    for i in range(len(oversplity)):
        tempreal=0.0
        for j in range(i,len(oversplity)):
            tempreal+= oversplit_sums[j][1]
        oversplity_sum.append(tempreal)

    oversplity_sum = np.array(oversplity_sum)/float(n_cells)*100
    oversplitx=oversplit_sums[:,0] #Select sorted 1st column of PTP values;
    dash_patch, = ax.plot(oversplitx, oversplity_sum, '--k', color='red', linewidth=5)

    temp_array = [oversplitx, oversplity_sum]
    np.savetxt(sim_dir+ptcs_name+'/oversplit.csv', temp_array, delimiter=",")

    #Plot pie-charts
    blue_patch = mpatches.Patch(color='blue')
    red_patch = mpatches.Patch(color='red')
    black_patch = mpatches.Patch(color='black')

    labels = ['% Correctly assigned', '% Oversplit', '% Missed', '% Oversplit vs PTP ']

    ax.legend([blue_patch, red_patch, black_patch, dash_patch], labels, fontsize=12, loc=0, 
    title="Cell Spikes - Detected")

    #Plot Large pie-chart
    pie_data = a_sum/n_cells, b_sum/n_cells, c_sum/n_cells
    colors = ['blue', 'red', 'black']

    draw_pie(ax,pie_data, 200, 65, size=2500, colors=colors)

    p = ax.axhspan(50.0, 100, facecolor='blue', alpha=0.05)
    p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)

    x = (50,50)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (100,100)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (70,70)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (80,80)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (90,90)
    plt.plot(x,y, 'r--', color='black',linewidth=1)


    #*********************************************** SECOND PLOT *****************************************
    ax = plt.subplot(1, 3, 2)
    title('Units - Contents & Purity' , fontsize=20)
    colors = ['green','red','black','magenta','purple']

    size_scale = 0

    #Plot scatter plots by # of spikes in unit (size)
    if size_scale==0: 
        s = [int(float(x)/5+1) for x in size]
    else:
        #Plot scatter plots by SD of max channel of unit
        s = [int((float(x)-min(sd))*60) for x in sd]

    #ax.scatter(ptps, purity,s=s, alpha=0.4, color='black')

    undersplity=[]
    undersplitx=[]
    noise=[]
    for i in range (len(ptps)):
        a = purity[i]/100.0 #Purity
        b = (sortedspikes_fromcells[i] - n_sorted_spikes[i]*purity[i]/100.0)/ n_sorted_spikes[i] #undersplit
        c = 1-a-b
        draw_pie(ax,[a,b,c], ptps[i], purity[i], size=125, colors=colors)

        undersplity.append(b+c) #*** MUST ADD UNDERSPLIT ERRORS AND NOISE TOGETHER FOR COMPOSITE ERROR METRIC *** 
        undersplitx.append(ptps[i])
        noise.append(c)

    #Compute oversplit sums to show dashed line on graph
    undersplit_sums = np.array([undersplitx, undersplity]).T  #Load oversplit data into 2 columns
    undersplit_sums = undersplit_sums[undersplit_sums[:,0].argsort()]  #Sort data by PTP amplitude (1st column)
    #print "Sorted: ", undersplit_sums

    undersplity_sum=[] #Cumulatively add all oversplit values backwards;
    for i in range(len(undersplity)):
        tempreal=0.0
        for j in range(i,len(undersplity)):
            tempreal+= undersplit_sums[j][1]
        undersplity_sum.append(tempreal)

    undersplity_sum = np.array(undersplity_sum)/float(n_units)*100
    undersplitx=undersplit_sums[:,0] #Select sorted 1st column of PTP values;
    dash_patch, = ax.plot(undersplitx, undersplity_sum, '--k', color='red', linewidth=5)

    temp_array = [undersplitx, undersplity_sum]
    np.savetxt(sim_dir+ptcs_name+'/undersplit.csv', temp_array, delimiter=",")

    green_patch = mpatches.Patch(color='green', label='Single Source Spikes')
    red_patch = mpatches.Patch(color='red', label='Multi Source Spikes')
    black_patch = mpatches.Patch(color='black', label='Noise')

    labels = ['% Correct Spikes', '% Undersplit', '% Noise', '% Undersplit vs. PTP']

    green_patch = mpatches.Patch(color='green', label='% Spikes from Single-Source')
    yellow_patch = mpatches.Patch(color='yellow', label='% Spikes from Multi-Source (Undersplitting)')
    red_patch = mpatches.Patch(color='red', label='Noise')
    ax.legend([green_patch, red_patch, black_patch, dash_patch], labels, loc=4, prop={'size':12}, 
    title='Sorted Spikes - Source Cells')

    #Plot Large pie charts spikes
    a = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))/100.
    b = float(sum(sortedspikes_fromcells))/float(sum(n_sorted_spikes))-a
    c = 1-a-b

    pie_data = float(a), float(b),float(c)
    colors = ['green', 'red', 'black']

    draw_pie(ax,pie_data, 200, 65, size=2500, colors=colors)

    p = ax.axhspan(90.0, 100, facecolor='green', alpha=0.05)
    p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)

    plt.ylim((0,100))
    plt.xlim((0,250))
    plt.xlabel('Average PTP (uV)',fontsize=17)
    plt.ylabel('% Purity\n (Spikes in Unit from Unique Cell)',multialignment='center', fontsize=17)

    x = (50,50)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (100,100)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (70,70)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (80,80)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (90,90)
    plt.plot(x,y, 'r--', color='black',linewidth=1)


    #********************************* SIZE PIECHARTS *******************************

    ax = plt.subplot(1, 3, 3)
    plt.ylim((0,100))
    plt.xlim((0,250))

    title('Units - Purity vs. Size of Unit' , fontsize=20)
    colors = ['green','red','black','magenta','purple']

    size_scale = 0

    #Plot scatter plots by # of spikes in unit (size)
    if size_scale==0: 
        s = [int(float(x)/5+1) for x in size]
    else:
        #Plot scatter plots by SD of max channel of unit
        s = [int((float(x)-min(sd))*60) for x in sd]

    #ax.scatter(ptps, purity,s=s, alpha=0.4, color='black')
    size = np.genfromtxt(sim_dir+ptcs_name+'/size.csv', dtype='float32', delimiter=",")

    undersplity=[]
    undersplitx=[]
    noise=[]
    for i in range (len(ptps)):
        a = purity[i]/100.0 #Purity
        b = (sortedspikes_fromcells[i] - n_sorted_spikes[i]*purity[i]/100.0)/ n_sorted_spikes[i] #undersplit
        c = 1-a-b
        #draw_pie(ax,[a,b,c], ptps[i], purity[i], size=size[i]/10., colors=colors)
        plt.scatter(ptps[i],purity[i], s=size[i]/5., color='blue', alpha=0.5)
        undersplity.append(b+c) #*** MUST ADD UNDERSPLIT ERRORS AND NOISE TOGETHER FOR COMPOSITE ERROR METRIC *** 
        undersplitx.append(ptps[i])
        noise.append(c)

    x = (50,50)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (100,100)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (70,70)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (80,80)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (90,90)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

    plt.show()


#************************************ PURITY VS. COMPLETENESS PLOTS ***********************************************

    plt.suptitle('Sort Analysis for: "' + ptcs_name + '" (In-Vitro Dataset; Rodrigo and Costas)', fontsize = 20)

    #************************************************************* THIRD PLOT *****************************
    ax = plt.subplot(1, 2, 1)
    plt.ylim((0,100))
    plt.xlim((0,250))

    colors = ['cyan','red']

    title('Units - Reliability Metric' , fontsize=20)

    a_sum = 0
    for i in range (len(completeness)):
        a = completeness[i]/100*purity[i]/100.0
        a_sum+= a
        b = 1-a #
        #print b
        draw_pie(ax,[a, b], ptps[i], a*100, size=125)

    blue_patch = mpatches.Patch(color='cyan', label='% Detected spikes')
    white_patch = mpatches.Patch(color='red', label='% Missed spikes')
    #yellow_patch = mpatches.Patch(color='yellow', label='Multi Source Spikes')
    #red_patch = mpatches.Patch(color='red', label='Noise')

    ax.legend([blue_patch, white_patch], ['% Cell Spikes Sorted into Unique Unit', '% Other'], loc=4,
    fontsize=12)

    pie_data = float(a_sum)/len(completeness), 1-float(a_sum)/len(completeness)
    colors = ['cyan', 'red']

    draw_pie(ax,pie_data, 200, 65, size=2500)
    #ax.pie(frac,colors=colors ,labels=labels, autopct='%1.1f%%')

    p = ax.axhspan(50.0, 100, facecolor='cyan', alpha=0.05)
    p = ax.axvspan(0.0, 200.0, facecolor='1.0', alpha=0.0)

    plt.xlabel('Average PTP (uV)',fontsize=17)
    plt.ylabel('% Purity*Completeness Composite\n (Completeness of Cell Activity in Unit)',multialignment='center', fontsize=17)

    x = (50,50)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (100,100)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (70,70)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (80,80)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (90,90)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    ##*********************** FOURTH PLOT ***************************
    ax = plt.subplot(1, 2, 2)

    plt.ylim((0,100))
    plt.xlim((0,100))
    plt.xlabel('% Purity',fontsize=17)
    plt.ylabel('% Completeness',multialignment='center', fontsize=17)

    colors = ['green','blue', 'black']

    a_sum = 0
    print purity
    print completeness

    for i in range (len(completeness)):
        a = purity[i]/200.
        b = completeness[i]/200.
        c = 1 - a - b
        draw_pie(ax,[a, b, c], purity[i], completeness[i], size=125)


    x = (80,80)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (90,90)
    y = (0,100)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (80,80)
    plt.plot(x,y, 'r--', color='black',linewidth=1)

    x = (0,250)
    y = (90,90)
    plt.plot(x,y, 'r--', color='black',linewidth=1)



    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

#*******************************************************************************************************
#***************************************** MULTISORT PLOTS *****************************
#*******************************************************************************************************

    plt.suptitle('MultiSort Comparison (In-Vitro Dataset; Rodrigo and Costas)', fontsize = 20)
    sim_dir = '/media/cat/Data1/in_vitro/all_cell_random_depth/'

    ax = plt.subplot(1, 2, 1)
    title('Oversplit Errors (Cells)' , fontsize=20)

    ax.set_xlim(0,150)
    ax.set_ylim((0,100))
    ax.set_xlabel('Average Peak-to-Peak (PTP) of Spikes in Cell (uV)',fontsize=15)
    ax.set_ylabel('% Oversplit - Cell spikes split across multiple units', fontsize=17)

    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")

    all_cell='james9/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    james_patch, = ax.plot(oversplit[0], oversplit[1], color='green', linewidth=3)

    all_cell='nick/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    nick_patch, = ax.plot(oversplit[0], oversplit[1], color='blue', linewidth=3)

    all_cell='dan/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    dan_patch, = ax.plot(oversplit[0], oversplit[1], color='orange', linewidth=3)

    all_cell='cat/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    cat_patch, = ax.plot(oversplit[0], oversplit[1], color='red', linewidth=3)

    all_cell='dirty_auto/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    dirty_patch, = ax.plot(oversplit[0], oversplit[1], color='cyan', linewidth=3)

    all_cell='clean_auto/'
    oversplit = np.genfromtxt(sim_dir+all_cell+'oversplit.csv', dtype='float32', delimiter=",")
    clean_patch, = ax.plot(oversplit[0], oversplit[1], color='magenta', linewidth=3)


    #************************************ 2ND OVERLAPPING FIGURES ****************************
    ax2 = ax.twinx()

    ax2.set_ylabel('Number of Cells from electrode with spikes > PTP (Normalized)', fontsize=16)
    ax2.set_xlim(0,150)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel('PTP (uV)',fontsize=20)

    plt.yticks([0,43,50,86,100])

    x = np.arange(7.,150.,.01)
    y = (3500/(x-7)+5).tolist()

    dash_patch, = ax2.plot(x,y, '--k', color='black',  linewidth=2)

    x = (150, x[min(range(len(y)), key=lambda i: abs(y[i]-43))])
    a = (43,43)
    ax2.plot(x,a, 'r--', color='black',linewidth=1)

    x = (100,100)
    a = (0,43)
    ax2.plot(x,a, 'r--', color='black', linewidth=1)

    x = np.arange(7.,150.,.01)
    x = (150, x[min(range(len(y)), key=lambda i: abs(y[i]-86))])
    a = (86,86)
    ax2.plot(x,a, 'r--', color='black',  linewidth=1)

    x = (50,50)
    a = (0,86)
    ax2.plot(x,a, 'r--', color='black', linewidth=1)

    #*************** FIRST LEGEND ********************

    labels = ['SS - Nick', 'SS - Auto + Clean', 'SS - Auto Only', 'SS - Catalin',  'Other - James (HHMI)', 'KK - Dan (Allen)']

    leg = plt.legend([nick_patch, clean_patch, dirty_patch, cat_patch, james_patch, dan_patch, dash_patch], labels, 
    loc=1, prop={'size':13}) #


    #*************************************** SECOND PLOT **********************************************

    ax = plt.subplot(1, 2, 2)
    title('Undersplit Errors + Noise (Units)' , fontsize=20)

    ax.set_xlim(0,150)
    ax.set_ylim((0,100))
    ax.set_xlabel('Average Peak-to-Peak (PTP) of Spikes in Unit (uV)',fontsize=15)

    #ax.set_ylabel('% Undersplit - Unit spikes from multiple cells',multialignment='center', fontsize=17)
    #ax.yaxis.tick_left()
    #ax.yaxis.set_label_position("left")

    all_cell='james9/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1], color='green', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    james_n_cells = float(len(set(bestmatch)))/float(n_cells)*100.0
    #print james_n_cells

    all_cell='nick/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1], color='blue', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    nick_n_cells = float(len(set(bestmatch)))/n_cells*100.0
    #print nick_n_cells

    all_cell='dan/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1], color='orange', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    dan_n_cells = float(len(set(bestmatch)))/n_cells*100.0

    all_cell='cat/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1],  color='red', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    cat_n_cells = float(len(set(bestmatch)))/n_cells*100.0

    all_cell='dirty_auto/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1],  color='cyan', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    dirty_auto_n_cells = float(len(set(bestmatch)))/n_cells*100.0

    all_cell='clean_auto/'
    undersplit = np.genfromtxt(sim_dir+all_cell+'undersplit.csv', dtype='float32', delimiter=",")
    dash_patch, = ax.plot(undersplit[0], undersplit[1],  color='magenta', linewidth=3)
    bestmatch = np.genfromtxt(sim_dir+all_cell+'bestmatch.csv', dtype='float32', delimiter=",")
    clean_auto_n_cells = float(len(set(bestmatch)))/n_cells*100.0


    #****************** INSERTED PLOT ***********
    ax2 = ax.twinx()

    ax2.set_xlim(0,150)
    ax2.set_ylim(0, 250)

    plt.yticks([0,43,50,86,100])

    x = np.arange(7.,150.,.01)
    y = (3500/(x-7)+5).tolist()

    ax2.plot(x,y, '--k', color='black',  linewidth=3)

    x = (150, x[min(range(len(y)), key=lambda i: abs(y[i]-43))])
    a = (43,43)
    ax2.plot(x,a, 'r--', color='black',linewidth=1)

    x = (100,100)
    a = (0,43)
    ax2.plot(x,a, 'r--', color='black', linewidth=1)

    x = np.arange(7.,150.,.01)
    x = (150, x[min(range(len(y)), key=lambda i: abs(y[i]-86))])
    a = (86,86)
    ax2.plot(x,a, 'r--', color='black',  linewidth=1)

    x = (50,50)
    a = (0,86)
    ax2.plot(x,a, 'r--', color='black', linewidth=1)

    subplots_adjust(hspace = 0.6)
    subplots_adjust(wspace = 0.22)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

    #**************************************************************************************************************************************
    #******************************************* 4 PLOT COMPARISONS *****************************************************************
    #**************************************************************************************************************************************

    #ax1 = plt.subplot2grid((2,2),(0,0))

    ax = plt.subplot(2, 2, 1)

    plt.suptitle('Comparison Metrics (In-Vitro Dataset)', fontsize = 18)
    ax.get_xaxis().set_visible(False)

    ax.set_xlim(0, 200)
    ax.set_ylim(0, 100)

    x = (0,200)
    a = (80,80)
    ax.plot(x,a, 'r--', color='black', linewidth=1)

    x = (0,200)
    a = (90,90)
    ax.plot(x,a, 'r--', color='black', linewidth=1)


    print "Total ground truth spikes: ", tsf.n_cell_spikes

    title('% Average (Per Unit): Purity, Undersplit, Noise' , fontsize=15)

    all_cell='nick/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent=0
    plt.bar(4+indent, ave_purity, 6, color='blue')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='blue', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    all_cell='clean_auto/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent+=8
    plt.bar(4+indent, ave_purity, 6, color='magenta')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='magenta', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    all_cell='dirty_auto/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent+=8
    plt.bar(4+indent, ave_purity, 6, color='cyan')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='cyan', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    all_cell='cat/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent+=8
    plt.bar(4+indent, ave_purity, 6, color='red')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='red', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    all_cell='james9/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent+=8
    plt.bar(4+indent, ave_purity, 6, color='green')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='green', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    all_cell='dan/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity) #WHY NEED THIS?
    ave_purity=sum(purity)/len(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.array(sortedspikes_fromcells) 
    other_cells = (sortedspikes_fromcells - n_sorted_spikes*purity/100.)/(n_sorted_spikes)*100.
    ave_undersplit = sum(other_cells)/len(other_cells)
    noise = (n_sorted_spikes-sortedspikes_fromcells)/(n_sorted_spikes)*100
    ave_noise = sum(noise)/len(noise)

    indent+=8
    plt.bar(4+indent, ave_purity, 6, color='orange')
    indent+=6
    plt.bar(4+indent, ave_undersplit, 6, color='orange', hatch='//')
    indent+=6
    plt.bar(4+indent, ave_noise, 6, color='black')

    #*** LEGEND ***

    #labels = ['Solid Color', 'SS - Auto + Clean', 'SS - Auto Only', 'SS - Catalin',  'Other - HHMI', 'KK - Dan']

    labels = ['Correct', 'Undersplit', 'Noise']
    solid = mpatches.Patch(color = 'yellow', label = 'Correct')
    hashed = mpatches.Patch(color = 'yellow', hatch='//', label = 'Undersplit')
    black = mpatches.Patch(color = 'black', label = 'Noise')

    ax.legend([solid, hashed, black], labels, loc=1, prop={'size':13}) 

    #*************** LEGEND ********************

    labels = ['SS - Nick', 'SS - Auto + Clean', 'SS - Auto Only', 'SS - Catalin',  'Sort_9', 'Ephys_Sort']

    leg = ax.legend([nick_patch, clean_patch, dirty_patch, cat_patch, james_patch, dan_patch, dash_patch], labels, 
    loc=4, prop={'size':13}) #



    #****************************************************** SECOND PLOT ***************************************

    ax = plt.subplot(2, 2, 3)
    ax.get_xaxis().set_visible(False)

    ax.set_xlim(0, 200)
    ax.set_ylim(0, 100)
    x = (0,200)
    a = (80,80)
    ax.plot(x,a, 'r--', color='black', linewidth=1)

    x = (0,200)
    a = (90,90)
    ax.plot(x,a, 'r--', color='black', linewidth=1)

    title('% Total (Sum Units): Purity, Undersplit, Noise' , fontsize=15)

    all_cell='nick/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    nick_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    #indent=0
    #plt.bar(4+indent, ave_purity, 6, color='blue')
    #indent+=6
    #plt.bar(4+indent, ave_undersplit, 6, color='blue', hatch='//')
    #indent+=6
    #plt.bar(4+indent, ave_noise, 6, color='black')

    indent=0
    plt.bar(4+indent, correct*scale, 6, color='blue')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='blue', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')
    #plt.text(4+indent+4, 102*scale, 'P: %d%%'%int(nick_purity), ha='center', va='bottom')

    all_cell='clean_auto/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    cleanauto_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    indent+=8
    plt.bar(4+indent, correct*scale, 6, color='magenta')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='magenta', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')

    all_cell='dirty_auto/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    dirtyauto_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    indent+=8
    plt.bar(4+indent, correct*scale, 6, color='cyan')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='cyan', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')

    all_cell='cat/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    cat_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    indent+=8
    plt.bar(4+indent, correct*scale, 6, color='red')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='red', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')

    all_cell='james9/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    james_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    indent+=8
    plt.bar(4+indent, correct*scale, 6, color='green')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='green', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')

    all_cell='dan/'
    purity = np.genfromtxt(sim_dir+all_cell+'purity.csv', dtype='float32', delimiter=",")
    purity=np.array(purity)
    n_sorted_spikes = np.genfromtxt(sim_dir+all_cell+'size.csv', dtype='float32', delimiter=",")
    scale = 1 #sum(n_sorted_spikes)/tsf.n_cell_spikes
    correct = sum(n_sorted_spikes*purity/100)/sum(n_sorted_spikes)*100
    sortedspikes_fromcells = np.genfromtxt(sim_dir+all_cell+'sortedspikes_fromcells.csv', dtype='float32', delimiter=",")
    noise = sum(n_sorted_spikes-sortedspikes_fromcells)/sum(n_sorted_spikes)*100
    other_cells = sum(sortedspikes_fromcells - n_sorted_spikes*purity/100.)/sum(n_sorted_spikes)*100.
    dan_purity = float(sum(n_sorted_spikes*purity))/float(sum(n_sorted_spikes))

    indent+=8
    plt.bar(4+indent, correct*scale, 6, color='orange')
    indent+=6
    plt.bar(4+indent, other_cells*scale, 6, color='orange', hatch='//')
    indent+=6
    plt.bar(4+indent, noise*scale, 6, color='black')

    #*************** LEGEND ********************

    #labels = ['Solid Color', 'SS - Auto + Clean', 'SS - Auto Only', 'SS - Catalin',  'Other - HHMI', 'KK - Dan']

    labels = ['Correct', 'Undersplit', 'Noise']
    solid = mpatches.Patch(color = 'pink', label = 'Correct')
    hashed = mpatches.Patch(color = 'pink', hatch='//', label = 'Undersplit')
    black = mpatches.Patch(color = 'black', label = 'Noise')

    ax.legend([solid, hashed, black], labels, loc=1, prop={'size':13}) 
    #leg.get_frame().set_alpha(1.0)

    ##**************************** THIRD PLOT ******************************

    ax = plt.subplot(2, 2, 4)
    subplots_adjust(hspace = 0.35)
    #subplots_adjust(wspace = 0.22)
    title('% Correctly Assigned Spikes (Normalized to Total Cell Spikes)' , fontsize=14)

    ax.get_xaxis().set_visible(False)

    # this is another inset axes over the main axes
    #b = axes([0.36, 0.5, .11, .12])
    #setp(b, xlim=(0,135),ylim=(0,100), xticks=[], yticks=[0,25,50,75,100])

    nick_time = 30
    nick_cores = 4
    cat_time = 60
    cat_cores = 1
    autoclean_time = 6
    autoclean_cores = 1
    autodirty_time = 1
    autodirty_cores = 1
    james_time = 200
    james_cores = 16
    dan_time = 85
    dan_cores = 8
    
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 120)

    indent=0
    ax.bar(4, nick_time, 8, color='blue')
    ax.text(8, 1.05*nick_time, '%dm'%int(nick_time), ha='center', va='bottom')

    indent+=16
    ax.bar(24, autoclean_time, 8, color='magenta')
    ax.text(28, 1.05*autoclean_time, '%dm'%int(autoclean_time), ha='center', va='bottom')

    ax.bar(44, autodirty_time, 8, color='cyan')
    ax.text(48, 1.05*autodirty_time, '%dm'%int(autodirty_time), ha='center', va='bottom')

    ax.bar(64, cat_time, 8, color='red')
    ax.text(68, 1.05*cat_time, '%dm'%int(cat_time), ha='center', va='bottom')

    ax.bar(84, james_time, 8, color='green')
    #ax.text(88, 1.05, '???', ha='center', va='bottom')

    ax.bar(104, dan_time, 8, color='orange')
    ax.text(108, 1.05*dan_time, '%dm'%int(dan_time), ha='center', va='bottom')

    title('Sorting Time    &     No. Cores Used', fontsize=15)
    plt.ylabel('Sorting Time (mins)',fontsize=15)


    ax2 = ax.twinx()
    #color=[black]

    #plt.ylabel('# Cores used in sort',color='black', fontsize=15)

    ax2.set_ylim(0, 12)
    ax2.set_xlim(0, 200)

    ax2.bar(12, nick_cores, 8, color='blue', hatch='//')
    ax2.text(16, 1.05*nick_cores, '%dC'%int(nick_cores), ha='center', va='bottom')

    ax2.bar(32, autoclean_cores, 8, color='magenta',hatch='//')
    ax2.text(36, 1.05*autoclean_cores, '%dC'%int(autoclean_cores), ha='center', va='bottom')
    ax2.bar(52, autodirty_cores, 8, color='cyan', hatch='//')
    ax2.text(56, 1.05*autodirty_cores, '%dC'%int(autodirty_cores), ha='center', va='bottom')
    ax2.bar(72, cat_cores, 8, color='red', hatch='//')
    ax2.text(76, 1.05*cat_cores, '%dC'%int(cat_cores), ha='center', va='bottom')

    ax2.bar(92, james_cores, 8, color='green')
    #ax2.text(96, .05, '???', ha='center', va='bottom')
    ax2.bar(112, dan_cores, 8, color='orange', hatch='//')
    ax2.text(116, 1.05*dan_cores, '8', ha='center', va='bottom')

    #ax2.plot(t, s2, 'r.')
    ax2.set_ylabel('# Cores Used in Sorting', color='black')


    #*************** LEGEND ********************

    labels = ['SS - Nick', 'SS - Auto + Clean', 'SS - Auto Only', 'SS - Catalin',  'Sort_9', 'Ephys_Sort']

    leg = plt.legend([nick_patch, clean_patch, dirty_patch, cat_patch, james_patch, dan_patch, dash_patch], labels, 
    loc=4, prop={'size':13}) #

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

    #********************PLOT AVERAGE UNIT RELIABILITY
    ax = plt.subplot(2, 2, 2)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 100)

    title('Average Unit Reliability (= Unit Purity * Cell Completeness)' , fontsize=14)

    indent=0

    all_cell='nick/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='blue')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='blue',  hatch='//')
    indent+=20

    all_cell='clean_auto/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='magenta')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='magenta',  hatch='//')
    indent+=20

    all_cell='dirty_auto/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='cyan')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='cyan',  hatch='//')
    indent+=20

    all_cell='cat/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='red')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='red',  hatch='//')
    indent+=20

    all_cell='james9/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='green')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='green',  hatch='//')
    indent+=20

    all_cell='dan/'
    ptp = np.genfromtxt(sim_dir+all_cell+'ptps.csv', dtype='float32', delimiter=",")
    purity = np.genfromtxt(sim_dir+all_cell+'/purity.csv', dtype='float32', delimiter=",")
    purity = np.array(purity)
    completeness = np.genfromtxt(sim_dir+all_cell+'/completeness.csv', dtype='float32', delimiter=",")
    completenss = np.array(completeness)
    a = sum(completeness/100*purity/100.0)/len(purity)*100
    ax.bar(8+indent, a, 8, color='orange')
    a = 0
    counter = 1
    for i in range(len(purity)):
        if ptp[i]>50.0:
            a+=(completeness[i]/100*purity[i]/100.0)
            counter+=1
    a = a/counter*100
    ax.bar(16+indent, a, 8, color='orange',  hatch='//')

    labels = ['All Units', 'Units > 50uv PTP']
    solid = mpatches.Patch(color = 'pink', label = 'All Units')
    hashed = mpatches.Patch(color = 'pink', hatch='//', label = 'Units > 50uv PTP')

    ax.legend([solid, hashed], labels, loc=1, prop={'size':13}) 

    plt.show()

def Plot_errors(tsf, Sort1):

    print "Plotting sorting errors"

    full_sort = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/full_sort.csv', dtype='float32', delimiter=",")
    bestmatch = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/bestmatch.csv', dtype='int32', delimiter=",")
    ptps = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/ptps.csv', dtype='float32', delimiter=",")

    purity = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/purity.csv', dtype='float32', delimiter=",")
    completeness = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/completeness.csv', dtype='float32', delimiter=",")
    n_sorted_spikes = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/size.csv', dtype='float32', delimiter=",")
    sortedspikes_fromcells = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/sortedspikes_fromcells.csv', dtype='float32', delimiter=",")

    spike_membership=[]
    for i in range(len(bestmatch)):
        spike_membership.append([])
    
    parser = csv.reader(open(tsf.sim_dir+Sort1.filename+'/spike_membership.csv'), delimiter=',')
    counter=0
    for l in parser: 
        for i in range(len(l)):
            spike_membership[counter].append(int(l[i]))
        counter+=1

    x = np.zeros((tsf.n_electrodes,30),dtype=np.float32)
    for i in range(tsf.n_electrodes):
        x[i]= tsf.Siteloc[i*2] + np.array(arange(0,30,1))

    #Plot 6 plots: ground truth, main unit, noise, + additionally mis-identified units
    for k in range(Sort1.n_units):

        #First locate the 2nd and 3rd best matching units - if available
        counts = []
        counts.append(list(count_unique(spike_membership[k])[0]))
        counts.append(list(count_unique(spike_membership[k])[1]))

        print "PTP of Unit: ", ptps[k]
        print counts
        #Remove bestmatch and noise from the 2D list
        for i in range(len(counts[0])):
            if (counts[0][i] == bestmatch[k]):
                del counts[0][i]
                del counts[1][i]
                break
        for i in range(len(counts[0])):
            if (counts[0][i] == 999):
                del counts[0][i]
                del counts[1][i]
                break

        second_match_index=999
        if (len(counts[1])>0):
            max_index, max_value = max(enumerate(counts[1]), key=operator.itemgetter(1))
            if (max_index<>None):
                print "Second unit index: ",  counts[0][max_index]
                second_match_index = counts[0][max_index]
                del counts[0][max_index]
                del counts[1][max_index]

        third_match_index=999
        if (len(counts[1])>0):
            max_index, max_value = max(enumerate(counts[1]), key=operator.itemgetter(1))
            if (max_index<>None):
                print "Third unit index: ",  counts[0][max_index]
                third_match_index =  counts[0][max_index]

        plt.suptitle('Analysis for: "' + Sort1.filename + '" (In-Vitro Data)    Unit: '+str(k+1) + ' / '+str(Sort1.n_units) 
        + ' , No. of spikes: '+str(len(Sort1.units[k])) + ", PTP: " + str(ptps[k]) + 'uV', fontsize = 18)


        # PLOT TEMPLATES
        #**************************************
        ax1 = plt.subplot2grid((2,6),(0,0))
        y1 = np.zeros((tsf.n_electrodes,30), dtype=np.float32)
        counter=0
        for j in range(len(Sort1.units[k])):
            if(spike_membership[k][j]==bestmatch[k]):
                print "Spike: ", j, " of ", len(Sort1.units[k]), " belongs to main unit: ", bestmatch[k]
                counter+=1
                for i in range(tsf.n_electrodes):
                    plt.plot(x[i], tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2-60*tsf.Siteloc[i*2+1], color='blue', alpha=0.2)
                    y1[i]+=tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2
                    #print y1[i]
                    #time.sleep(10)

        y1=y1/float(counter)
        for i in range(tsf.n_electrodes):
            plt.plot(x[i],y1[i]-60*tsf.Siteloc[i*2+1], color='black',linewidth=3,alpha=1.0)

        title('Best Match Cell: ' +str(bestmatch[k]) + '\n No. Spikes: '+str(counter)+" / "+str(len(Sort1.units[k])) , fontsize=13)
        main_spikes = counter
        #Cat: Plot legends; the scale has to be multiplied by V_HP_Scale
        #ax1.text(-25*1E+4/tsf.SampleFrequency,-12500, '100uV', rotation=90)
        #ax1.text(-19*1E+4/tsf.SampleFrequency,-14650, '1.0ms', rotation=0)
        #x_local=[-18.0,-18.0]
        #y_local=[-12000.0,-12000.0-100*2*10]
        #ax1.plot(x_local,y_local, linewidth=5, color='black') 

        #x_local=[-18.0,-8.0]
        #y_local=[-14000.0,-14000]
        #ax1.plot(x_local,y_local, linewidth=5, color='black')

        ax1.set_ylim(-15000,1000)
        ax1.xaxis.set_major_formatter(plt.NullFormatter())
        ax1.yaxis.set_major_formatter(plt.NullFormatter())


        #**************************************
        ax4 = plt.subplot2grid((2,6),(0,1))
        y3 = np.zeros((tsf.n_electrodes,30), dtype=np.float32)
        counter=0
        for j in range(len(Sort1.units[k])):
            if(spike_membership[k][j]==second_match_index and second_match_index<>999):
                counter+=1
                print "Spike: ", j, " of ", len(Sort1.units[k]), " belongs to second unit"
                for i in range(tsf.n_electrodes):
                    plt.plot(x[i], tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2-60*tsf.Siteloc[i*2+1], color='magenta', alpha=0.2)
                    y3[i]+=tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2

        if(second_match_index==999):
            title('No additional units', fontsize=15)
        else:
            title('2nd Best Match cell: '+  str(second_match_index)+ '\nNo. of Spikes: '+str(counter)+ ' / '+str(len(Sort1.units[k])), fontsize=13)
            y3=y3/float(counter)
            for i in range(tsf.n_electrodes):
                plt.plot(x[i],y3[i]-60*tsf.Siteloc[i*2+1], color='black',linewidth=3,alpha=1.0)

        second_spikes = counter
        plt.ylim(-15000,1000)
        ax4.xaxis.set_major_formatter(plt.NullFormatter())
        ax4.yaxis.set_major_formatter(plt.NullFormatter())

        #**************************************
        ax5 = plt.subplot2grid((2,6),(0,2))
        y4 = np.zeros((tsf.n_electrodes,30), dtype=np.float32)
        counter=0
        for j in range(len(Sort1.units[k])):
            if(spike_membership[k][j]==third_match_index and third_match_index<>999):
                counter+=1
                print "Spike: ", j, " of ", len(Sort1.units[k]), " belongs to third unit"
                for i in range(tsf.n_electrodes):
                    plt.plot(x[i], tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2-60*tsf.Siteloc[i*2+1], color='green', alpha=0.2)
                    y4[i]+=tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2

        if(third_match_index==999):
            title('No additional units', fontsize=15)
        else:
            title('3rd Best Match cell: '+  str(third_match_index)+ '\nNo. of Spikes: '+str(counter)+ ' / '+str(len(Sort1.units[k])), fontsize=13)

            y4=y4/float(counter)
            for i in range(tsf.n_electrodes):
                plt.plot(x[i],y4[i]-60*tsf.Siteloc[i*2+1], color='black',linewidth=3,alpha=1.0)

        third_spikes = counter

        plt.ylim(-15000,1000)
        ax5.xaxis.set_major_formatter(plt.NullFormatter())
        ax5.yaxis.set_major_formatter(plt.NullFormatter())

        #**************************************
        ax2 = plt.subplot2grid((2,6),(1,0))
        title('Matching Cell: ' +str(bestmatch[k])+ "\n Total Spikes: "+str(len(tsf.cell_spikes[bestmatch[k]])), fontsize=13)
        y0 = np.zeros((tsf.n_electrodes,30), dtype=np.float32)
        counter=0
        for j in range(len(tsf.cell_spikes[bestmatch[k]])):
            print "Cell: ", bestmatch[k], " Spike: ", j, " of ", len(tsf.cell_spikes[bestmatch[k]])
            counter+=1
            if(counter<1000):
                for i in range(tsf.n_electrodes):
                    plt.plot(x[i], tsf.ec_traces[i][tsf.cell_spikes[bestmatch[k]][j]-15:tsf.cell_spikes[bestmatch[k]][j]+15]*2-60*tsf.Siteloc[i*2+1], 
                    color='black', alpha=0.2)
                    y0[i]+=tsf.ec_traces[i][tsf.cell_spikes[bestmatch[k]][j]-15:tsf.cell_spikes[bestmatch[k]][j]+15]*2

        y0=y0/float(min(1000,counter))
        for i in range(tsf.n_electrodes):
            plt.plot(x[i],y0[i]-60*tsf.Siteloc[i*2+1], linewidth=3,color='white')

        plt.ylim(-15000,1000)
        ax2.xaxis.set_major_formatter(plt.NullFormatter())
        ax2.yaxis.set_major_formatter(plt.NullFormatter())

        #**************************************
        ax3 = plt.subplot2grid((2,6),(1,1))
        y2 = np.zeros((tsf.n_electrodes,30), dtype=np.float32)
        counter=0
        for j in range(len(Sort1.units[k])):
            if(spike_membership[k][j]==999):
                print "Spike: ", j, " of ", len(Sort1.units[k]), " is noise."
                counter+=1
                if(counter<1000): #Plot only first 1000 noise spikes;
                    for i in range(tsf.n_electrodes):
                        if Sort1.units[k][j]>15: #Skip spikes that are in the first 15 timesteps in the record;
                            plt.plot(x[i], tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2-60*tsf.Siteloc[i*2+1], color='red', alpha=0.2)
                            y2[i]+=tsf.ec_traces[i][Sort1.units[k][j]-15:Sort1.units[k][j]+15]*2

        y2=y2/float(min(1000,counter))
        for i in range(tsf.n_electrodes):
            plt.plot(x[i],y2[i]-60*tsf.Siteloc[i*2+1], linewidth=3,color='black')
        noise_spikes=counter

        plt.ylim(-15000,1000)
        title('Noise (or In-Vivo Unit) \n No of Spikes: ' + str(counter)+" / "+str(len(Sort1.units[k])), fontsize=13)
        ax3.xaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_formatter(plt.NullFormatter())



        #Plot purity and completnes pie_charts
        #**************************************
        ax6 = plt.subplot2grid((2,6),(1,2))

        ##PURITY
        #a = purity[k]/100.0 
        #b = (sortedspikes_fromcells[k] - n_sorted_spikes[k]*purity[k]/100.0)/n_sorted_spikes[k] #undersplit
        #c = 1-a-b
        #colors = ['green','red', 'black']
        #draw_pie(ax6,[a, b, c], 500,3000 , size=2000, colors=colors)

        ##COMPLETENESS
        #colors = ['blue', 'black']
        #a = completeness[k]/100.0
        #b = 1-a
        #draw_pie(ax6,[a,b], 500,1000, size=2000, colors=colors)

        #title('Purity & Completeness', fontsize=13)
        title('Unit Composition (%)', fontsize=13)
        #colors = ['blue', 'magenta', 'green', 'red']
        a = main_spikes/ float(len(Sort1.units[k]))
        b = second_spikes/float(len(Sort1.units[k]))
        c = third_spikes/float(len(Sort1.units[k]))
        e = noise_spikes/float(len(Sort1.units[k]))
        d = 1-a-b-c-e
        #draw_pie(ax6,[a,b,c,d], 500,1000, size=2000, colors=colors)

        ax6.bar(25, 100, 50, color='red')
        ax6.bar(25, (a+b+c+d)*100, 50, color='pink')
        ax6.bar(25, (a+b+c)*100, 50, color='green')
        ax6.bar(25, (a+b)*100, 50, color='magenta')
        ax6.bar(25, a*100, 50, color='blue')

        plt.ylim(0,100)
        plt.xlim(0,100)
        ax6.xaxis.set_major_formatter(plt.NullFormatter())
        #ax6.yaxis.set_major_formatter(plt.NullFormatter())

        #**************************************
        ax7 = plt.subplot2grid((2,6),(0,3), rowspan=2)
        for i in range(tsf.n_electrodes):
            plt.plot(x[i],y1[i]-60*tsf.Siteloc[i*2+1], color='blue',linewidth=3,alpha=0.6)
            plt.plot(x[i],y2[i]-60*tsf.Siteloc[i*2+1], color='red',linewidth=3,alpha=1)

        title('Correct vs. Noise', fontsize=15)
        plt.ylim(-15000,1000)
        ax7.xaxis.set_major_formatter(plt.NullFormatter())
        ax7.yaxis.set_major_formatter(plt.NullFormatter())

        #ax7.text(-25,-12500, '100uV', rotation=90)
        #ax7.plot([-18.0,-18.0],[-12000.0,-12000.0-100*2*10], linewidth=5, color='black') #Cat: the scale has to be multiplied by V_HP_Scale

        #ax7.text(-19,-14650, '1.0ms', rotation=0)
        #ax7.plot([-18.0,-8.0],[-14000.0,-14000.0], linewidth=5, color='black')

        #**************************************
        ax8 = plt.subplot2grid((2,6),(0,4), rowspan=2)
        if(second_match_index<>999):
            for i in range(tsf.n_electrodes):
                plt.plot(x[i],y1[i]-60*tsf.Siteloc[i*2+1], color='blue',linewidth=3,alpha=.6)
                plt.plot(x[i],y3[i]-60*tsf.Siteloc[i*2+1], color='magenta',linewidth=3,alpha=1)
            title('Correct vs. 2nd Largest', fontsize=15)
        else:
            title('No additional matches')
        plt.ylim(-15000,1000)
        ax8.xaxis.set_major_formatter(plt.NullFormatter())
        ax8.yaxis.set_major_formatter(plt.NullFormatter())

        #**************************************
        ax9 = plt.subplot2grid((2,6),(0,5), rowspan=2)
        if(third_match_index<>999):
            for i in range(tsf.n_electrodes):
                plt.plot(x[i],y1[i]-60*tsf.Siteloc[i*2+1], color='blue',linewidth=3,alpha=.6)
                plt.plot(x[i],y4[i]-60*tsf.Siteloc[i*2+1], color='green',linewidth=3, alpha=1)
            title('Correct vs. 3rd Largest', fontsize=15)
        else:
            title('No additional matches')

        plt.ylim(-15000,1000)
        ax9.xaxis.set_major_formatter(plt.NullFormatter())
        ax9.yaxis.set_major_formatter(plt.NullFormatter())

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        #pylab.savefig(tsf.sim_dir+Sort1.filename+'/Unit_'+str(k)+'.png')
        #plt.get_current_fig_manager().window.showMaximized()
        #plt.savefig(tsf.sim_dir+Sort1.filename+'/Unit_'+str(k)+'.png', bbox_inches='tight', dpi=600)
        plt.show()

def Make_tsf(tsf, Sort1, file_out, unit_no, index_1, index_2):
    "Making TSF: use only iformat = 1002"

    iformat=1002
    f = open(file_out, 'wb')
    f.write(tsf.header)
    f.write(struct.pack('i', iformat))
    f.write(struct.pack('i', tsf.SampleFrequency))
    f.write(struct.pack('i', tsf.n_electrodes))
    f.write(struct.pack('i', tsf.n_vd_samples))
    f.write(struct.pack('f', tsf.vscale_HP))

    for i in range (tsf.n_electrodes):
        f.write(struct.pack('h', tsf.Siteloc[i*2]))
        f.write(struct.pack('h', tsf.Siteloc[i*2+1]))
        f.write(struct.pack('i', tsf.Readloc[i]))

    print "Writing data"
    for i in range (tsf.n_electrodes):
        print "Writing ch: ", i
        tsf.ec_traces[i].tofile(f)

    #Determine spike sets to save as ground truth
    maxchan = np.genfromtxt(tsf.sim_dir+Sort1.filename+'/maxchan.csv', dtype='float32', delimiter=",")

    spike_membership=[]
    for i in range(len(maxchan)):
        spike_membership.append([])
    
    parser = csv.reader(open(tsf.sim_dir+Sort1.filename+'/spike_membership.csv'), delimiter=',')
    counter=0
    for l in parser: 
        for i in range(len(l)):
            spike_membership[counter].append(int(l[i]))
        counter+=1

    print spike_membership[unit_no]

    fake_spike_times=[]
    fake_spike_assignment=[]
    fake_spike_channels=[]
    n_spikes=0
    for i in range(len(spike_membership[unit_no])):
        if(spike_membership[unit_no][i]==index_1):
            fake_spike_times.append(Sort1.units[unit_no][i])
            fake_spike_assignment.append(1)
            fake_spike_channels.append(maxchan[unit_no+1])
            n_spikes+=1
        if(spike_membership[unit_no][i]==index_2):
            fake_spike_times.append(Sort1.units[unit_no][i])
            fake_spike_assignment.append(2)
            fake_spike_channels.append(maxchan[unit_no+1])
            n_spikes+=1
    fake_spike_times=np.array(fake_spike_times, dtype=np.int32)
    fake_spike_assignment=np.array(fake_spike_assignment, dtype=np.int32)
    fake_spike_channels=np.array(fake_spike_channels, dtype=np.int32)

    print fake_spike_times[0:10]
    print fake_spike_assignment[0:10]
    print fake_spike_channels[0:10]


    f.write(struct.pack('i', n_spikes)) #Write # of fake spikes

    print "No. cell spikes: ", n_spikes
    if (n_spikes>0):
        if (iformat==1001):
            f.write(struct.pack('i', 30)) # = struct.unpack('i',fin.read(4))[0] 
            f.write(struct.pack('i', n_spikes)) #Write # of fake spikes
        fake_spike_times.tofile(f)
        fake_spike_assignment.tofile(f) 
        fake_spike_channels.tofile(f) 
    f.close()


def find_nearest(array,value,dist):
    #print "Looking for: ", value
    #print "Diff to closest value: ", array[np.abs(array-value).argmin()]-value
    if abs(array[np.abs(array-value).argmin()]-value) < dist:
        idx = (np.abs(array-value)).argmin()
        return idx
    else:
        return None

def fit_log(x,a,b,c,d):
    return a + b*np.log(100-1/x)

def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)

def exponential(x,a,b,c,d):
    return a*np.exp(b*x+c)+d

def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def sigmoid2(x, x0, k):
    y = 100 / (1 + np.exp(-k*(x-0.5)))
    return y 

#def sigmoid3(x, a, b,c,d):
    #y = 1-np.exp((-(x-a)/b)**c) #  + b/x
    #return y

def asymptoticx(x, a, b,c,d):
    y =a+b*x/(x+d) #+ /x
    return y

def asymptoticlogx(x, a, b,c,d,e,f,g):
    y =a+ b*np.log(x+d) #/(1+log(x+f)) #+ /x
    return y


def sigmoid_function(xdata, x0, k):
    y = np.exp(-k*(xdata-x0)) / (1 + np.exp(-k*(xdata-x0)))
    return y

def residuals(p,x,y):
    return y - sigmoid(p,x)

def resize(arr,lower=0.0,upper=1.0):
    arr=arr.copy()
    if lower>upper: lower,upper=upper,lower
    arr -= arr.min()
    arr *= (upper-lower)/arr.max()
    arr += lower
    return arr

def erf(x, a1, a2, a3, a4, a5, p):

    # constants
    #a1 =  0.254829592
    #a2 = -0.284496736
    #a3 =  1.421413741
    #a4 = -1.453152027
    #a5 =  1.061405429
    #p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

def draw_pie(ax,ratios=[0.4,0.3,0.3], X=0, Y=0, size = 1000, colors=colors):
    N = len(ratios)
 
    xy = []
 
    start = 0.
    for ratio in ratios:
        x = [0] + np.cos(np.linspace(2*math.pi*start,2*math.pi*(start+ratio), 30)).tolist()
        y = [0] + np.sin(np.linspace(2*math.pi*start,2*math.pi*(start+ratio), 30)).tolist()
        xy1 = zip(x,y)
        xy.append(xy1)
        start += ratio
 
    for i, xyi in enumerate(xy):
        ax.scatter([X],[Y] , marker=(xyi,0), s=size, facecolor=colors[i] )

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def gauss(x, *p):
    A, mu, S = p
    return A*np.exp(-(x-mu)**2/(2.*S**2))

