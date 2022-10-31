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


#************************** CODE START ****************************
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

np.random.seed(12345678)  #fix random seed to get the same result

#***************************** PICK SIM_DIR **************************
#***** IN VITRO ******
#sim_dir = '/media/cat/Data1/in_vitro/all_cell_random_depth/'

#***** IN VIVO *****
#NICK
#sim_dir = '/media/cat/Data1/in_vivo/nick/87 - track 7c spontaneous craziness/'
#sim_dir = '/media/cat/Data1/in_vivo/nick/ptc17/08-tr1-blankscreen/'
#sim_dir = '/media/cat/Data1/in_vivo/nick/ptc21/19-tr2-blankscreen/'
#sim_dir = '/media/cat/Data1/in_vivo/nick/ptc21/21-tr2-blankscreen/'
#sim_dir = '/media/cat/Data1/in_vivo/nick/saline/2012-01-31_test1/'
#sim_dir = '/media/cat/Data1/in_vivo/nick/saline/911a_newlab_test/'

#DAN
#sim_dir = '/media/cat/Data1/in_vivo/dan/openephysLGN_128chs/'

#PETER
#sim_dir = '/media/cat/Data1/in_vivo/peter/'


#****** IN SILICO ******
#sim_dir = '/media/cat/Data1/in_silico/scinet/June29/ECP/'
#sim_dir = '/media/cat/4TB/in_silico/ucsd_July5_rat_5k_20khz/ECP_0_noise/'
#sim_dir = '/media/cat/4TB/in_silico/ucsd_Sep6_rat_3k_20Khz/ECP_0/'
#sim_dir = '/media/cat/4TB/in_silico/ucsd_July5_rat_5k_20khz/ECP_0_noise/'
#sim_dir = '/media/cat/4TB/in_vivo/peter/baseline.h5/'
#sim_dir = '/media/cat/4TB/in_vivo/sev/M71/DeepCSD5/'
#sim_dir = '/media/cat/4TB/in_vivo/dan/M150566/2015-04-01_16-41-14_flashes_vis/'





#***************************** PICK TSF_NAME **************************
#*****IN VITRO*****
#tsf_name = '20Khz_10cells_truth'
#tsf_name = 'Traces/James_in-vitro_extracellular_traces'


#*****IN VIVO***** 
#NICK
#tsf_name = '87 - track 7c spontaneous craziness' 
#tsf_name = '2012-01-31_test1'
#tsf_name = '911a_newlab_test'
#tsf_name = '08-tr1-blankscreen'
#tsf_name = '19-tr2-blankscreen'
#tsf_name = '21-tr2-blankscreen'

#DAN
#tsf_name = 'openephysLGN_128chs_filtered' 

#PETER
#tsf_name = 'baseline.h5'

#TIM
#sim_dir = '/media/cat/4TB/in_vivo/tim/'


#************ SFN 2015 poster results ****************
sim_dir = '/media/cat/4TB/SFN_2015_poster/'
data_dir = 'data/'

anonymous = True

sorted_dirs = [ 'cat_sort/', 'nicksw_sort_pp1/', 'dan_sort/'] #,'peter_sort/', 'nickst_sort/', 'pierre_sort/', 'martin_sort/']
#sorted_dirs = [ 'cat_sort2/']

sorted_file = 'vivo_1'


#********************* READ TSF FILES *******************
tsf_name = data_dir + '/'+ sorted_file + '_truth_filtered'

f = open(sim_dir + tsf_name + '.tsf', "rb")
print "Loading ", sim_dir + tsf_name + '.tsf'
tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
tsf.sim_dir = sim_dir
tsf.tsf_name = tsf_name
f.close()

if False:
    x = np.arange(0, tsf.n_vd_samples, dtype=np.float32) #Time in samplesize timesteps
    x = x/tsf.SampleFrequency
    font_size=40
    for i in range(8):
        plt.plot(x[102000:105500],-np.float32(tsf.ec_traces[i][102000:105500])/4.+i*30, color='blue', linewidth=4)
    plt.ylim(240,-20)
    plt.xlim(5.1,5.15)
    plt.xlabel("Time (seconds)", fontsize=font_size)
    plt.ylabel("Depth ($\mu$m)", fontsize=font_size)

    plt.tick_params(axis='both', which='major', labelsize=font_size)
    plt.tick_params(axis='both', which='minor', labelsize=font_size)
    plt.show()
    quit()

##*************** READ GROUND TRUTH FILE : SORT1 ************************

ptcs_flag = 1

if True:
    Sort1 = Loadptcs(sorted_file+'_truth_filtered', sim_dir+data_dir, ptcs_flag, save_timestamps=False)
    Sort1.name=sorted_file+'_truth_filtered'
    Sort1.filename=sorted_file+'_truth_filtered'
    Sort1.rawfile= tsf_name
    Sort1.directory=sim_dir+data_dir
    Sort1.sim_dir = sim_dir
    Sort1.data_dir = data_dir

    ##Dark neuron plots; load sort2 for everyone
    #Sort1 = Loadptcs('ECP_0_allcells', sim_dir+sorted_dirs[0], ptcs_flag, save_timestamps=False)
    #Sort1.name='ECP_0_allcells'
    #Sort1.filename='ECP_0_allcells'
    #Sort1.rawfile= tsf_name
    #Sort1.directory=sim_dir+sorted_dirs[0]
    #Sort1.sim_dir = sim_dir
    #Sort1.data_dir = sorted_dirs[0]

    if False: Plot_histogram(Sort1)
    
##************** READ 2ND SORT DATA **********

for sorted_dir in sorted_dirs:

    print "******* LOADING: ", sorted_dir, " ************* "

    #*** READ PTCS FILES *****
    if ('cat' in sorted_dir) or ('nicksw' in sorted_dir):
        Sort2 = Loadptcs(sorted_file, sim_dir+sorted_dir, ptcs_flag, save_timestamps=False)
        Sort2.name=sorted_file
        Sort2.filename=sorted_file
        Sort2.rawfile= tsf_name
        Sort2.directory=sim_dir+sorted_dir
        Sort2.sim_dir = sim_dir

    #*** READ CSV FILES *****
    if ('dan' in sorted_dir) or ('nickst' in sorted_dir) or ('pierre' in sorted_dir) or ('peter' in sorted_dir) or ('martin' in sorted_dir):
        fname = sim_dir + sorted_dir + sorted_file + '.csv'

        #flag = 1 for James' .csv; flag = 2 for Dan .csv; flag = 4 and 5 for CK's row-wise data;
        Sort2 = Loadcsv(fname, tsf, sim_dir + sorted_file + '/',2)  
        Sort2.directory=sim_dir + sorted_dir
        Sort2.name = sorted_file
        Sort2.filename=sorted_file
        Sort2.sim_dir = sim_dir
        Sort2.chanpos=[999]
        f.close()

    #Plot firing rate distribution of ground truth cells:
    if False:
        firing_rate = []
        for i in range(len(Sort2.units)):
            if len(Sort2.units[i])/240>1:
                #for j in range(30):
                    firing_rate.append(len(Sort1.units[i])/325)
            
        font_size=80
        #Plot histogram from cumulative trains
        #plt.title("Firing rate distributions", fontsize=font_size)
        plt.ylabel("No. of cells", fontsize=font_size)
        plt.xlabel("Frequency", fontsize=font_size)
        plt.xlim(0,20)
        plt.ylim(0,10)
        bin_width = 1
        x = np.arange(0,20,bin_width)
        print x
        y = np.histogram(firing_rate, bins = x)
        #plt.bar(y[0], color='blue')#, linewidth=4)
        plt.bar(y[1][:-1], y[0], bin_width-.1, color='blue', alpha=0.65)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.tick_params(axis='both', which='major', labelsize=font_size)
        plt.tick_params(axis='both', which='minor', labelsize=font_size)
        plt.show()    

        quit()

    ##Plot cross-correlograms between pop spikes and each unit:
    ##Synchronized state from: 300sec -> 750sec
    #temp2 = np.array(Sort2.units[0],dtype=np.float32)/Sort1.samplerate

    #for i in range(Sort1.n_units):
        #print "plot: ", i
        #temp1 = np.array(Sort1.units[i], dtype=np.float32)/Sort1.samplerate
        ##temp2 = np.array(Sort1.units[i],dtype=float32)/float(Sort1.samplerate)
        ##x1 = temp2[np.where(np.logical_and(temp2>=300., temp2<=750.))]
        ##if len(x1)==0: continue
        ##y = scipy.signal.correlate(temp1,temp2, mode='full')
        #y = np.correlate(temp1,temp2, mode='full')
        ##y/=min(y)
        #x = np.arange(len(y))
        #plt.bar(x,y)
        ##plt.xlim(-500.,+500.)
        ##figManager = plt.get_current_fig_manager()
        ##figManager.window.showMaximized()
        #plt.show()

    Sort2.n_spikes=0
    for i in range(len(Sort2.units)):
        Sort2.n_spikes+= len(Sort2.units[i])
    #print "Total Sorted spikes: ", Sort2.n_spikes
    #print ""

    #Run the compare_sorts algorithm on the data if not already done
    Compare_sorts(Sort1, Sort2, tsf)
    #Compare_sorts(Sort2, Sort1, tsf)

    #Compute the FP_histograms needed for some of the plots - NOT ALWYAS NEEDED
    #Compute_FP_histograms(Sort1, Sort2, tsf)

    #Plot_Composite_metrics(sim_dir, sorted_dirs, sorted_file, anonymous)

    #Compute_darkneuron(Sort2)
    #quit()


#Plot_fpfn(sim_dir, sorted_dirs, sorted_file, anonymous)

Plot_Composite_metrics(sim_dir, sorted_dirs, sorted_file, anonymous)


##************ SAVE .BIN FILE
#file_name = sim_dir+'Traces/james_spikes_timesteps.txt'

#with open(file_name, "w") as file:
    #writer = csv.writer(file, delimiter=',')
    #for i in range(Sort2.n_units):
        #writer.writerow(np.array(Sort2.units[i], dtype=np.float32))
#quit()

#****************** READ CSV SORTED DATA ********************

#file_name = '/media/cat/Data1/in_vivo/87 - track 7c spontaneous craziness.ptcs'
#f = open(file_name, "r")
#csv_file=Loadcsv()

##********** CHECK DISTRIBUTIONS ***********

#fname = sim_dir + sorted_file + '.ptcs'
#f = open(fname, "rb")

#work_dir = sim_dir
#Sort1 = Loadptcs(f, work_dir)
#Sort1.name=sorted_file
#Sort1.filename=sorted_file
#Sort1.rawfile= tsf_name
#Sort1.directory=work_dir

##Load the locations of each channel - required for distribution plots
#depth=[]
#for k in range(Sort1.nneurons):
    #depth.append(Sort1.chanpos[Sort1.maxchan[k]][1])

#distributions(tsf, depth, Sort1)
#quit()

#******************** CHECK CELLS AGAINST UNITS *******************
####Check to see if data has already been checked already

#print "Checking sort for: ", sorted_file 
#check_dir = sim_dir + sorted_file + '/'
#NOT USED ANYMORE
#Check_sort(tsf, Sort2, check_dir) #Do the comparisons OR load from disk of already done; 
#ptcs_name = sorted_file+'/'

#Check_sort_plot(tsf, Sort1, sim_dir , ptcs_name) #Simplify arguments to just Sort1.directory

#******************* MAKE SINGLE UNIT TSF FILES ***********************

#file_out = work_dir + 'Unit1.tsf'

#unit_no = 1
#index_1= 1
#index_2=2
#Make_tsf(tsf, Sort1, file_out, unit_no, index_1, index_2)

#****************** CHECK TWO SORTS OF SAME DATA AGAINST EACH OTHER **********************

#Compare_sorts(Sort1, Sort2, tsf) 
#Compare_sorts(Sort2, Sort1, tsf) 
#Compare_sorts_plot(Sort1, Sort2, tsf)

#Plot_confusion(Sort1, Sort2, tsf)

#Plot_fpfn(sim_dir, sorted_dirs, sorted_file, anonymous)


Composite_metrics(sim_dir, sorted_dirs, sorted_file, anonymous)

#****************** MATCH UNITS OF DIFFERENT SORTED TIMES AGAINST EACH OTHER *************

#Compare_sorts(Sort1, Sort2, tsf) 
#Compare_sorts(Sort2, Sort1, tsf) 

#Plot_rasters(Sort1)

#Match_units(Sort1, Sort2, tsf, tsf2)


#********************* CONFUSION MATRIX PLOTS **************************

#Plot_specificity(Sort1, Sort2, tsf)
#Plot_sensitivity(Sort1, Sort2, tsf)
#Plot_table(Sort1, Sort2, tsf)

#Plot_confusion(Sort1, Sort2, tsf)


#Plot_errors(tsf, Sort1)
