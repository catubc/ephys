'''Sergey's CSD code 
NB: Possible new version available; 
TODO: units are not clear but will be required for publication;
'''

import numpy as np
import matplotlib.pyplot as plt
from math import *
import h5py
import imp
import struct


#LOAD ALLEN INST DATA from hdf5 files
#datadir = '/media/cat/4TB/in_vivo/dan/M181420/181420_2015-04-07_18-10-46_flash_vis/'
#fname = 'M181420_181420_2015-04-07_18-10-46_flash_vis'
#f5 = h5py.File(datadir+fname+'.hdf5','r')


#LOAD Nick lab Data
#file_name = '/media/cat/4TB/in_vivo/nick/ptc15/87 - track 7c spontaneous craziness/87 - track 7c spontaneous craziness'
#file_name = '/media/cat/4TB/in_vivo/nick/ptc18/05-tr1-blankscreen/05-tr1-blankscreen'
#file_name = '/media/cat/4TB/in_vivo/nick/ptc18/05-tr1-blankscreen/05-tr1-blankscreen'


#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/63-tr5c-blankscreen/63-tr5c-blankscreen'
#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/61-tr5c-blankscreen/61-tr5c-blankscreen'

#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/02-tr1-blankscreen_while_decr_propofol/02-tr1-blankscreen_while_decr_propofol'
#file_name = '/media/cat/4TB/in_vivo/nick/ptc20/13-tr1-blankscreen/13-tr1-blankscreen'
#file_name = '/media/cat/4TB/in_vivo/nick/ptc22/27-tr2-blankscreen/27-tr2-blankscreen'



#LOAD TIM lab Data
#file_name = '/media/cat/4TB/in_vivo/tim/2015-7-22-13-HL'
#file_name = '/media/cat/4TB/in_vivo/tim/2015-7-23-11-FL'

#data_dir = '/media/cat/4TB/in_vivo/tim/2015-7-22/'
#file_name = [
#'2015-7-22-9-W1',
#'2015-7-22-10-w2',
#'2015-7-22-11-w3',
#'2015-7-22-12-v1',
#'2015-7-22-13-V2',
#'2015-7-22-14-FL',
#'2015-7-22-16-HL',
#'2015-7-22-17-HL',
#'2015-7-22-18-A1',
#'2015-7-22-19-A2'
#]

data_dir = '/media/cat/4TB/in_vivo/tim/2015-7-23/'
file_name = [
'2015-7-23-5-w1',
'2015-7-23-6-w2',
'2015-7-23-7-v1',
'2015-7-23-8-v2',
'2015-7-23-9-a1',
'2015-7-23-10-A2',
'2015-7-23-11-FL',
'2015-7-23-12-FL',
'2015-7-23-13-HL',
'2015-7-23-14-HL',
]


for file_no in range(len(file_name)):
    tsf_name=data_dir + file_name[file_no]+'/'+file_name[file_no]+'_lp.tsf'
    print "Load extracellular file: ", tsf_name
        
    fin = open(tsf_name, "rb")
    header = fin.read(16)
    iformat = struct.unpack('i',fin.read(4))[0] 
    SampleFrequency = struct.unpack('i',fin.read(4))[0] 
    n_electrodes = struct.unpack('i',fin.read(4))[0] 
    n_vd_samples = struct.unpack('i',fin.read(4))[0] 
    vscale_HP = struct.unpack('f',fin.read(4))[0] 

    if iformat==1002:
        Siteloc = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ecp =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)
    ecp.shape = n_electrodes, n_vd_samples

    #**************** COMPUTE STIM_TRIGGERED LFP AVERAGES ***************
    lfp_total=[]
    no_lfps = 1

    for s in range(no_lfps):  #Loop over population spikes
        stim_file = data_dir + file_name[file_no]+'/'+file_name[file_no]+'_stim_times.csv'
        stim_array= np.loadtxt(stim_file,dtype=np.float32) #Load data in seconds
        print stim_array

        #Compute stim array triggered averages;
        window = 0.500 #+/- secs to average 
        print "Computing cumulative LFP: ", s
        SampleFrequency = 20000 #Samples per sec 
        print "SampleFrequency: ", SampleFrequency
        print n_electrodes

        lfp_cumulative = np.zeros((n_electrodes,2*window*SampleFrequency), dtype=np.float32)
        print lfp_cumulative

        for i in range(n_electrodes):
            counter=0
            for j in stim_array:

                start=int(j*SampleFrequency-window*SampleFrequency)
                end=min(start+2*window*SampleFrequency, n_vd_samples) #only search to end of recording

                if (end-start)==2*window*SampleFrequency: 
                    lfp_cumulative[i] += ecp[i][start:end]
                    counter+=1

        print "counter:",counter
        lfp_cumulative /= counter

        lfp_total.append(lfp_cumulative)

    #ecp_trial_avgCSD = f5[fname+'/stim_triggered_ecp_bright'][...] #f5['/DeepCSD4/stim_triggered_ecp'][...]
    ecp_trial_avgCSD=[]
    for i in range(no_lfps):
        ecp_trial_avgCSD.append(lfp_total[i])

    #probe_layout = f5[fname+'/probe_layout'][...]  #HDF5 files only
    probe_layout = Siteloc[1::2]
    print "probe_layout: ", probe_layout

    #vscale = f5[fname+'/vscale'][...]  #HDF5 files only
    vscale = vscale_HP


    #ecp_trial_avgCSD=ecp_trial_avgCSD[::2]*vscale #Grabs every other electrode - required for IMEC
    #ecp_trial_avgCSD=np.flipud(ecp_trial_avgCSD)


    z = probe_layout*1E-3 #Convert to mm size
    csdmod = imp.load_source('csd_est_funds','csd_est_funcs.py')
    

    sigma = 0.3  # extracellular conductivity (mS/mm)
    b = 1.0   # assumed radius of the column (mm)
    SNR = 6.0    # average ratio of signal to noise on each channel

    [A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
    [W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
    [W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check

    #[fig100,fig101] = csdmod.show_operators(A,A0,W,W0,R,R0)    # display forward and inverse operators 
    #plt.show()


    for i in range(no_lfps):
        csd3=np.dot(W,ecp_trial_avgCSD[i])   # units: mS/mm^3*mV = uA/mm^3
        
        lfp3_0mean = np.mean(ecp_trial_avgCSD[i], axis=0)
        lfp3 = ecp_trial_avgCSD[i]-lfp3_0mean
        
        tlims =[-0.1,0.1]
        
        fig1=plt.figure(1)
        
        CSD=csd3
        plt.subplot(1,1,1);

        dt=1./ecp_trial_avgCSD[i].shape[1]
        t=np.arange(-ecp_trial_avgCSD[i].shape[1]/2.,ecp_trial_avgCSD[i].shape[1]/2.)*dt

        #plt.subplot(1,1,1)
        VspanLFP = np.array([-1, 1])*0.8*abs(CSD).max()
        plt.imshow(CSD, vmin = -1000, vmax = 1000,extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto',
                   interpolation='lanczos'); 
        plt.xlim(tlims)
        plt.ylabel("Distance from top of electrode (mm)")
        plt.xlabel("Time from stimulus onset (seconds)")
        
        plt.ylim(max(probe_layout)*1E-3,min(probe_layout)*1E-3)
        plt.colorbar();
        plt.plot([0,0],[0,1500], 'r--', linewidth=2, color='black')
        plt.title('CSD for recording: '+ file_name[file_no])

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.savefig(data_dir + file_name[file_no]+'/'+file_name[file_no]+'.png', bbox_inches='tight', dpi=100)
    plt.show()
