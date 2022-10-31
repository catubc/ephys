from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

from scipy import stats
import numpy as np
import time, math
import sys
import os.path
from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab
import re
import glob
import filter

global colors

colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']


#************************** LOAD PTCS FILES ****************************

np.random.seed(12345678)  #fix random seed to get the same result
np.random.seed(seed=None)

#***************************************************************************************
#************************************NICK DATA******************************************
#***************************************************************************************

ptcs_flag = 1

#PTC15
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc15/'
#sorted_file = '87 - track 7c spontaneous craziness'
#sorted_file = '86 - track 7c slope-b'

#PTC17 - Mostly desynchronized states with some mid-recording transitions
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'
#sorted_file = '02-tr1-blankscreen'     #~40mins desynch; transition to SYNCH in last ~5mins
#sorted_file = '03-tr1-mseq32_40ms'     #43mins desynch; mid-transition to SYNCH
#sorted_file = '04-tr1-MVI_1400_5s'      #40mins desynch 20-30 band; no synch
#sorted_file = '06-tr1-MVI_1406-1410'    #40mins desynch 20-30 fluctuating band; no synch
#sorted_file = '07-tr1-driftbar_longbar' #8min desynch 20-30Hz band; no synch
#sorted_file = '08-tr1-blankscreen'      #21mins desynch; ~9mins SYNCH period in middle
#sorted_file = '32-tr1-blankscreen'      #40mins desynch; short 8 min ~SYNCH period
#sorted_file = '41-tr1-blankscreen'     #8mins desynch;
#sorted_file = '44-tr2b-blankscreen'     #9mins desynch w. very STRONG multiple DISTINCT BANDS; 
#sorted_file = '69-tr2b-blankscreen'     #40mins desynch; Gamma x 2 bands! very interesting
#sorted_file = '70-tr2b-blankscreen_66Hz' #15mins desynch; single ~25hz band;

#PTC18 - Mostly desynch states 
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc18/'
#sorted_file = '05-tr1-blankscreen' #desynchronized ~35Hz strong band
#sorted_file = '08-tr1-blankscreen_66Hz' #desynchronized !35hz
#sorted_file = '09-tr1-blankscreen_200Hz' #desynchronized 35-40Hz band
#sorted_file = '36-tr1-blankscreen_66Hz' #desynchronized 
#sorted_file = '37-tr1-blankscreen_200Hz' #desynchronized 
#sorted_file = '70-tr2c-blankscreen_200Hz' #Mixed state sync/desynchronized  58 UNITS - interesting locking 4 units     SORTED - NG
#sorted_file = '71-tr2c-blankscreen_66Hz' #looks up-state synchronized; 21 LFP units; 1, 2, 7, 9, 11, 13, 14, 17 show some locking; 
                                            #LFP 19 lots of locking by most cells; also 20 a bit muckier
#sorted_file = '79-tr2c-blankscreen' #5mins recording; mid section shows some upphases? not clear       SORTED - NG
                                        #no locking; 

#PTC20
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc20/'
#sorted_file = '04-tr1-blankscreen' #
#sorted_file = '13-tr1-blankscreen' #~15 units, majority show narrow locking;                                      SORTED - GOOD
#sorted_file = '17-tr1-blankscreen'  #Cat: no big pops; specgram: bimodal 2 & 8Hz; INHIBITION IN a couple of cells; 
#sorted_file = '19-tr1-blankscreen' #Cat: few pop spikes; 2 cells Oscillate;
#sorted_file = '20-tr1-blankscreen_66Hz' #Loads of oscillating spikes
#sorted_file = '30-tr1-blankscreen'  #Cat: many oscillations
#sorted_file = '34-tr2-blankscreen'  #not many large pop spikes; < 50% phase locked units;
#sorted_file = '39-tr2-blankscreen' #Mainly oscillations;
#sorted_file = '50-tr2-blankscreen'  #5 of 8 units locked at some level
#sorted_file = '52-tr2-blankscreen'  #4 of 7 -8 units locked
#sorted_file = '54-tr2-blankscreen_66Hz'  # 6 of 7units oscillating
#sorted_file = '58-tr3-spontaneous' #2 of 8 show locking; despite several large pop -spikes
#sorted_file = '64-tr3-blankscreen' #Cat doesn't load properly; 10 of 24 units show locking
#sorted_file = '66-tr3-blankscreen_66Hz' #doesn't load properly; majority lock; lots of power in lfp                 SORTED - Decent, BUT OSCILLATORY NEAR LFP EVENT 
#sorted_file = '71-tr3-blankscreen' #13 of 26 lock; interesting inhibition in some units                             SORTED - GOOD only short time scale...

#PTC21
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#sorted_file = '02-tr1-blankscreen_while_decr_propofol' #30+ units almost all units lock across all pop spike ranges
#sorted_file = '06-tr1-blankscreen' #5 units only; some locking; many pop spikes/events - 
#sorted_file = '10-tr1b-blankscreen' #6 of 7 units some lockking, esp for larger pop-spikes
#sorted_file = '13-tr2-spontaneous' #12units, very small though; 4 show some locking; 
#sorted_file = '19-tr2-blankscreen' #20units, all show gradual increase in locking; cat's sort; ~13 have broad peaks
#sorted_file = '21-tr2-blankscreen' #14units, ~8 some locking; cat's sort
#sorted_file = '27-tr2-blankscreen_66Hz' #20+units, some locking (<50%)
#sorted_file = '37-tr2-blankscreen' #SOMETHING WRONG WITH nchans values;
#sorted_file = '42-tr3-blankscreen' #Some locking < 50%
#sorted_file = '45-tr4-spontaneous' #some locking < 50%; lots of high amplitude pop spikes
#sorted_file = '50-tr5-spontaneous' #15-16 small units, no locking
#sorted_file = '61-tr5c-blankscreen' #30 of 35 units - VERY STRONG LOCKING across all amplitudes; Large pop spikes      SORTED - SUPER SYNCHRONIZED RECORDING
#sorted_file = '63-tr5c-blankscreen' #30 units - very strong locking;                                                  SORTED - SUPER SYNC
#sorted_file = '65-tr5c-blankscreen' #30+units, very little locking; low pop spike amplitudes                        SORTED - MIXED - Units 28, 35 show some locking
#sorted_file = '67-tr5c-blankscreen' #40+units; Verty strong locking across most; High amplitude POPs                SORTED - SUPER SYNC; many interesting units
#sorted_file = '68-tr5c-blankscreen_66Hz' #35 of 43 units Very strong locking; high POPS distribution + tail
#sorted_file = '77-tr6-spontaneous' #55 units almost no locking at all;
#sorted_file = '85-t6b-blankscreen' #6units, 5 of them strong locking
#sorted_file = '87-t6b-blankscreen' #7units, all show locking
#sorted_file = '89-t6b-blankscreen' #8units, all show locking;

#PTC22
#sim_dir = '/media/cat/4TB/in_vivo/nick/ptc22/'
#sorted_file = '07-tr1-blankscreen' #desynchronized
#sorted_file = '09-tr1-blankscreen' #desynchronized - but to small freq band ~appears synchornized on specgram 
#sorted_file = '11-tr1-blankscreen' #desyncrhonized
#sorted_file = '21-tr1-blankscreen' #desynchronized
#sorted_file = '27-tr2-blankscreen'#; ptcs_flag = 0 #synchronized; Nick's Sort is great
#sorted_file = '32-tr2-blankscreen'#; #synchronized - Nick's sort shows some structure
#sorted_file = '36-tr2-blankscreen'#; #synchronized - Nick's sort shows some structure

#***************************************************************************************
#************************************TIM DATA*******************************************
#***************************************************************************************

#2015-5-6
#sim_dir = '/media/cat/4TB/in_vivo/tim/dongsheng/2015-5-6/'
#sorted_file = '2015-5-6-2'
#sorted_file = '2015-5-6-5'

#**************************** NICK'S CAT DATA ***************************
#track_name = 'ptc17-tr2b'     
#track_name = 'ptc18-tr2c'
#track_name = 'ptc20_tr1'
#track_name = 'ptc21-tr5c'
#track_name = 'ptc21-tr5c-57'
#track_name = 'ptc22-tr1'
#track_name = 'ptc22-tr2'

#********* DONGSHENG'S MICE ************
#track_name = '2015-11-27-2'    #Anesthetized cortex
#track_name = '2015-11-27-4'    #Awake cortex
#track_name = '2015-11-27-14'   #Anesthetized subcortical
#track_name = '2015-11-27-16'   #Awake subcortical
#track_name = '2015-12-1-1'   #cortex: anesth + awake; NB: HAS DIFFERENT 10 x faster TIME COMPRESSION THAN other mice;


#********* CAT'S MICE ***********
track_name = '2016_05_03/tr2'    #Anesthetized cortex




lfp_resplit = True  #Resplit the lfp only when loading entire tracks

if track_name == 'ptc17-tr2b':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc17/'
    #lfp_sort_file = 'track2b_all_notch_1_240Hz_5subsample'  #Ptc17 tr2b
    lfp_sort_file = 'track2b_all_notch_5subsample'  #Ptc17 tr2b
    #recordings = ['44','45','46','47','48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73']
    #recordings = ['51']
    #recordings = ['48','51','58','59','60','67','68']; lfp_resplit=False #Recs containing sync periods
    recordings = ['54'] #Recs containing sync periods
    recs = recordings; lfp_resplit=False #['50']#'47', '48', '49', '50', '51','52','53']

    #recs = ['51'];  lfp_resplit=False

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit()

if track_name == 'ptc18-tr2c':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc18/'
    lfp_sort_file = 'track2c_all_notch_5subsample'  #Ptc17 tr2b
    recordings = ['47','48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94']
    #recordings = ['70', '71', '79']; lfp_resplit = False #Spont recs
    recs = recordings; lfp_resplit=False
    #recs = ['59']

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit() #concatenate track wide LFP recs and zero out si < threshold for LFP sorting

if track_name == 'ptc20_tr1':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc20/'
    lfp_sort_file = 'track1_all_notch_5subsample'  #Ptc17 tr2b
    recordings = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19']
    #recordings = ['70', '71', '79']; lfp_resplit = False #Spont recs
    recs = recordings;  lfp_resplit=False

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit()


if track_name == 'ptc21-tr5c':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc21/'
    lfp_sort_file = 'track5c_59-75_notch_1_240Hz_5subsample'    #Ptc21 tr5c
    #recordings = ['59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75']
    #recordings = ['61', '63', '65', '68']; lfp_resplit = False #Spontaneous recs
    #recordings = ['59','60','62','64','66']; lfp_resplit = False #nat scenes
    recordings = ['60']
    recs =  recordings;  lfp_resplit=False

if track_name == 'ptc21-tr5c-57':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc21/'
    lfp_sort_file = '57-tr5c-mseq32_40ms_all_notch_5subsample'
    recordings = ['57']
    recs =  recordings;  lfp_resplit=True

if track_name == 'ptc22-tr1':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc22/'
    lfp_sort_file = 'track1_all_notch_5subsample'    #Ptc21 tr5c
    #recordings = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22']
    recordings = ['09']
    #recordings = ['61', '63', '65', '68']; lfp_resplit = False #Spontaneous recs
    #recordings = ['05','06','08','10','19','20']; lfp_resplit = False #nat scenes
    recs =  recordings #['60']#, '60', '61']

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit()

if track_name == 'ptc22-tr2':
    sim_dir = '/media/cat/12TB/in_vivo/nick/ptc22/'
    lfp_sort_file = 'track2_all_notch_5subsample'    #Ptc21 tr5c
    #recordings = ['23','24','25','26','27','28','29','30','31','32','33','34','35','36']
    recordings = ['26']
    #recordings = ['27', '32', '36']; lfp_resplit = False #Spontaneous recs
    #recordings = ['28','33']#,'08','10','19','20']; lfp_resplit = False #nat scenes
    recs =  recordings; lfp_resplit=False #['60']#, '60', '61']

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit()

if track_name == 'ptc22-tr2-26':
    pass


#********** DONGSHENG'S MICE *************

if track_name == '2015-11-27-2':     #Anesthetized cortex
    sim_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/'
    lfp_sort_file = '2015-11-27-2-10electrodein-iso1'   
    recordings = [lfp_sort_file] 
    recs = recordings 

if track_name == '2015-11-27-4':
    sim_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/'
    lfp_sort_file = '2015-11-27-4-10electrodein-iso0'   
    recordings = [lfp_sort_file] 
    recs = recordings 

if track_name == '2015-11-27-16':
    sim_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/'
    lfp_sort_file = '2015-11-27-16-deep-iso0'
    recordings = [lfp_sort_file] 
    recs = recordings 

if track_name == '2015-12-1-1':
    sim_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-1/'
    lfp_sort_file = '2015-12-1-1-10electrodeiniso1'   
    recordings = [lfp_sort_file]
    recs = recordings

    if False: Concatenate_LFP_sort(sim_dir, recordings, track_name); quit()


#********CAT'S DATA **********

if track_name == '2016_05_03/tr2':
    sim_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_03/tr2/'
    lfp_sort_file = '2016_05_03_tr2_spont_afterinsertion_160503_162252'   
    recordings = [lfp_sort_file]
    recs = recordings


#*************************************************************************
#**************************** LOAD .lfp.zip FILE  ************************
#*************************************************************************
#NB: Load original .lfp.zip files from Martin; Do not load compressed/uncompressed concatenated LFPs
if True:
    lfp = Object_temp()
    lfp.data = []           #LFP data
    lfp.t_start = []        #Start of lfp data relative high-pass record
    lfp.t_end = []          #end of lfp data
    lfp.rec = []            #recording index (string) 
    lfp.time_index = []     #running time index of data - need in order to index into LFP event data file
    lfp.rec_name = []       #holds names or recordings
    lfp.sampfreq = 1000
    counter=0
    for recording in recordings:
        print "Loading lfp zip: ", recording
        print sim_dir+recording+"*"
        temp_name = glob.glob(sim_dir+recording+"*")[0]
        temp_name = temp_name.replace('_lfp',''); temp_name = temp_name.replace('_mua','')  #Dongsheng directories are duplicated
        file_name = temp_name.replace(sim_dir, '')
        
        data_in  = np.load(sim_dir+file_name+'/'+file_name+'.lfp.zip')
        lfp.rec_name.append(file_name)
        lfp.rec.append(recording)
        lfp.t_start.append(data_in['t0']*1E-6)      #Convert to seconds
        lfp.t_end.append(data_in['t1']*1E-6)        #convert to seconds
        #print lfp.t_start, lfp.t_end
        start_offset = np.zeros((len(data_in['data']), int(lfp.t_start[counter]*1E3)), dtype=np.int16)
        
        #Remove 60Hz frequency
        temp_data = data_in['data']
        for k in range(len(temp_data)):
            temp_data[k] = np.array(filter.notch(temp_data[k])[0], dtype=np.int16) # remove 60 Hz mains noise
        data_out = temp_data

        lfp.data.append(np.concatenate((start_offset, data_out), axis=1))

        lfp.time_index.append(lfp.t_end[counter])
        counter+=1


#***** LOAD SUA SORT ****
print "LOADING SUA SPIKES"
Sorts_sua = Concatenate_ptcs(sim_dir, track_name, recordings)


#******* LOAD LFP SORT*****
#Load whole track sorted LFP events
if 'tr' in track_name: #Save lfp sort file as time corrected and in each recording directory
    print "Loading whole track sorted lfp events"

    work_dir = sim_dir + lfp_sort_file + "/"
    Sort_lfp = Loadptcs(lfp_sort_file, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
    Sort_lfp.name=file_name
    Sort_lfp.filename=lfp_sort_file
    Sort_lfp.directory=sim_dir

    if lfp_resplit: Split_lfp_ptcs(Sort_lfp, lfp, recordings) #Convert track_wide lfp to individual rec lfp sorts; RUN ONCE!
    
    Sorts_lfp = Concatenate_lfp_events(sim_dir, recordings)

else:
    print "Loading single lfp ptcs file" #Tim's data

    work_dir = sim_dir + lfp_sort_file + "/"
    lfp_ptcs = lfp_sort_file+'_lp_compressed'
    Sort_lfp = Loadptcs(lfp_ptcs, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
    Sort_lfp.name=file_name
    Sort_lfp.filename=lfp_sort_file
    Sort_lfp.directory=sim_dir
    
    Save_lfp_ptcs(Sort_lfp, lfp, recordings, track_name)
    
    Sorts_lfp = Concatenate_lfp_events(sim_dir, recordings)

##Load single recording LFP events
#else: 
    #print "Loading single track sorted lfp events"
    #file_name = file_name+"_lp_compressed"
    #work_dir = sim_dir + file_name + "/"
    #print "Loading LFP events: ", file_name
    #Sort2 = Loadptcs(file_name, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
    #Sort2.name=file_name
    #Sort2.filename=file_name
    #Sort2.directory=work_dir
    

#***************************************************************************************
#********************************** ANALYSIS METHODS ************************************
#***************************************************************************************


##Overlayed multi-taper specgram / SUA / LFP events / Synchrony index; Fig 3 top;
Plot_MUA_vs_LFP_events(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs[0])    #Plot specgram + LFP events + SUA events + MUA 

#Statistics on msl piece-wise for particular lfp cluster and recording; Fig 4 top; Fig 6;
#Plot_msl(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Control on msl; Fig 4 middle;
#Plot_msl_controls(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Plt msl comparison across lfp clusters; Fig 4 bottom
#Plot_msl_compare_across_lfp_clusters(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Plt msl comparison across time clusters; Fig 5;
#Plot_msl_compare_across_time(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Statistics on msl only for synchronized period of recording: Fig 7
#Plot_percentlock_synchronized_state(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Statistics of msl for syn periods - colapse across all recs    #Fig 7; also Fig 8 output (non_lock and lock arrays)
#Plot_percentlock_synchronized_state_allrec(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name)  

#Plot lfp cluster templates w. standard deviations; Fig 3 bottom;
#Plot_templates(sim_dir, Sorts_lfp, lfp, recs, track_name)  


quit()
#*************** OLD FUNCTIONS **************

#Plot_rasters(Sorts, Sort2)

#Plot_rasters_SUA_vs_LFP(Sorts, Sort2)

#Plot_firingrate(Sort1)
#Compute_kuebler(Sort1, Sort2)
#Compute_upstate(Sort1)

#Sort2=False
#Plot_LFP(lfp) #,Sort2)

#Plot_CSD(lfp)

#Plot_upstate(Sort1,Sort2,lfp)

#Plot_specgram(lfp)

#Plot_PSD(lfp,'PSD')

#Plot_triggered_activity(Sort1,Sort2,lfp)

#Plot_LFP_triggered(Sort1, Sort2, lfp)

#Plot_SD_vs_specgram(Sort1, tsf, lfp)

#Plot_trackwide_LFP_triggered_activity(Sorts, Sort2, tsf)


#Plot_multitaper_specgram(time_series, sampfreq)





#***************************************************************************************
#************************************** OLD CODE ***********************************
#***************************************************************************************

##********Compare Nick's export LFP w. Martin's

#print "Loading LFP"
#fname = sim_dir+file_name+'/'+file_name+'.lfp'
#fin = open(fname, "rb")

#n_electrodes = 10
#header = fin.read(24)
#surf_filename = fin.read(256)

#n_lfp_samples = struct.unpack('i',fin.read(4))[0] 
#n_lfp_records = struct.unpack('i',fin.read(4))[0] 
#n_vd_samples = n_lfp_samples

#lfp_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_lfp_samples)
#lfp_traces.shape = n_electrodes, n_lfp_samples                

##Careful: Nick's LFP may not align with original data; Needs separate loading/consideration
#lfp_time_stamps =  np.fromfile(fin, dtype=np.int64, count=n_lfp_records)
#lfp_time_stamps.shape= n_lfp_records/n_electrodes, n_electrodes
#lfp_time_stamps=lfp_time_stamps.T/1E+6/60. #Convert from usec to minutes???

#lfp_sites = np.fromfile(fin, dtype=np.int32, count=10)
#lfp_sites_positions = np.fromfile(fin, dtype=np.int32, count=10)

#xx = np.arange(0, len(lfp_traces[0]),1)
#plt.plot(xx, lfp_traces[0]*0.488)
#data_plot = np.concatenate((np.zeros(lfp.t_start[0]*1E3, dtype=np.int16), data_in['data'][0]))
#print len(xx), len(data_plot)
#martin_lfp = data_plot*data_in['uVperAD']
#plt.plot(xx, martin_lfp)

#print martin_lfp/lfp_traces[0]
#plt.show()
#quit()



#***************************************************************************************
#******************************** LOAD LFP & GENERATE SPECGRAM *************************
#***************************************************************************************

if True:
    
    #Use.tsf lfp file to generate specgram directly
    if False: 
        if (os.path.exists(tsf.fname[:-4]+'specgram.npy')==False):
            spec_ch=9 #Choose specgram channel, lower is usually prefered
            show_plot = True
            Plot_specgram_tsf(tsf, show_plot=show_plot, spec_ch=spec_ch)
            if show_plot: quit()
        else:
            tsf.specgram = np.load(tsf.fname[:-4]+'specgram.npy')

    if False: #TEMPORARY FOR plotting synchrony index
        
        sim_dir = '/media/cat/4TB/in_vivo/tim/2015-11-18/'
        sorted_file = '2015-11-18-2-9electrodein'
                    
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_lp.tsf', "rb")

        print "Loading ", sim_dir + sorted_file+'/' + sorted_file + '.tsf'

        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
        
        temp_data = np.load(sim_dir+file_name+'/'+file_name+'.lfp.zip')

        compressed_lfp.vscale = temp_data['uVperAD']
        for i in range(10):
            data = temp_data['data'][i]*compressed_lfp.vscale

        
        spec_ch=9 #Choose specgram channel, lower is usually prefered
        show_plot = True
        Plot_specgram_tsf(tsf, show_plot=show_plot, spec_ch=spec_ch)
        if show_plot: quit()
        

    #Load single LFP from .lfp files made from .srf files
    if False:
        sorted_file = '2015-8-6-1'
        fname = sim_dir + sorted_file+'/' + sorted_file + '.lfp'
        lfp = Load_lfp(fname)

        lfp.header = sorted_file
        lfp.fname = sim_dir + sorted_file+'/' + sorted_file + '_lp.tsf'    

        #Plot_PSD(lfp,'PSD')
        spec_ch=7 #Choose specgram channel, lower is usually prefered
        show_plot = True
        Plot_specgram(lfp, show_plot=show_plot, spec_ch=spec_ch)
        if show_plot: quit()

    if False:
        #Load multiple LFP
        for i in range(len(sorted_files)):
            fname = sim_dir + sorted_files[i]+'/' + sorted_files[i] + '_lp.tsf'
            lfp = Load_lfp(fname)

            #lfp.header = sorted_files[i]
            lfp.fname = sim_dir + sorted_files[i]+'/' + sorted_files[i] + '_lp.tsf'

            #Plot_PSD(lfp,'PSD')
            spec_ch=7 #Choose specgram channel, lower is usually prefered
            show_plot = True
            Plot_specgram(lfp, show_plot=show_plot, spec_ch=spec_ch)

        if show_plot: quit()




#**** COMPRESS LFP FOR LFP SORTING *****************************

if False:
    Compress_LFP(lfp)

#*** CONCATENATE TRACK WIDE LFP TRACES ************************

if False: 
    #recordings = ['59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75']
    recordings = ['all']
    Concatenate_LFPtraces(sim_dir, track_name, recordings)
    quit()

#*** PLOT TRACK WIDE SPECGRAMS *************************

if False: 
    track_name = 'tr3'
    Plot_trackwide_specgram(sim_dir, track_name)


