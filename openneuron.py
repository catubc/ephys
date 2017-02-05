#Main file for doing analysis on simultaneous ephys and ophys data

from animal import *
from animal_analysis import *

import pprint
import numpy as np


#Initialize pyqtgraph - needs to be done at beginning - or even better in a main conditional below...
#from pyqtgraph.Qt import QtCore, QtGui

#app = QtGui.QApplication([])    #Starts up the QtGui; makes gl object that needs to be passed to graph routines...


#******************************************
#************* DATA FILES *****************
#******************************************
home_dir = '/media/cat/500GB/in_vivo/tim/cat/'
#home_dir = '/media/cat/12TB/in_vivo/tim/cat/'
#animal_names = ['2016_07_11_vsd']
#animal_names = ['2016_07_12_vsd']
#animal_names = ['2016_07_15_vsd']
animal_names = ['2016_07_20_vsd']

#******************************************
#********** LOAD ANIMALS  ********************
#******************************************

animals = []
for animal_name in animal_names:
    if False:        
        animal = Animal(animal_name, home_dir)

        animal.rhd_to_tsf()
        animal.tsf_to_lfp()
        animal.lfp_compress()
        animal.bin_to_npy()
        animal.rhd_digital_save()

    else:
        animal = Animal(animal_name, home_dir)

    #List processed files
    #pprint.pprint(animal.filenames)

    #Compute STA maps
    if False: 
        ptcs_files = glob.glob(animal.home_dir+animal.name+'/tsf_files/*.ptcs')
        for ptcs_file in ptcs_files:
            if ('track2' in ptcs_file) and ('2nd' in ptcs_file) and ('compressed' in ptcs_file):
            #if ('track2' in ptcs_file) and ('1st' not in ptcs_file) and ('compressed' in ptcs_file):
            #if ('track2' in ptcs_file) and ('compressed' in ptcs_file):
                compute_sta(animal, ptcs_file)

    #Plot STA maps
    if False: 
        #prefix = ''; 
        #prefix = '_hp'
        prefix = '_lp_compressed'
        for k in range(0,len(animal.filenames),1):
            if ('track2' in animal.filenames[k]) and ('2nd' not in animal.filenames[k]):# and ('compressed' in animal.filenames[k]):
                sta_maps(animal, k, prefix)

    #Make STA movies
    if False: 
        pprint.pprint(animal.filenames)
        #prefix = ''
        #prefix = '_hp'
        prefix = '_lp_compressed'
        for k in range(0,len(animal.filenames),1):
            if ('track1' in animal.filenames[k]):# and ('3rd' in animal.filenames[k]):
            #if ('track2' in animal.filenames[k])  and ('1st' not in animal.filenames[k]):
                sta_movies(animal, k, prefix)

    #Generate specgrams
    if False: 
        for k in range(0,len(animal.filenames),1):
            animal.load_channel(animal.filenames[k].replace('rhd_files','tsf_files')+'_lp.tsf', channel=55) #Loads single channel as animal.tsf
            Specgram_syncindex(animal.tsf)


    #Generate MSL plots
    if True:
        for rec_index in range(0,len(animal.filenames),1):
            if ('track2' in animal.filenames[rec_index]):# and ('3rd' in animal.filenames[rec_index]):
                msl_plots(animal, rec_index)
    

    quit()




##Dimensionality reduction
#dim_red_data = multi_dim_scaling(all_traces, 4, home_dir, mouse.name)      #Methods: 0: MDS; 1:tSNE; 2: PCA; 3: Barnes-Hut tSNE


##Clustering: Kmeans vs. Meanshift
#if True:
#    n_clusters = 25
#    cluster_labels = KMEANS(dim_red_data, n_clusters)
#else:
#    cluster_labels = Meanshift(dim_red_data)
#    n_clusters = np.max(cluster_labels)+1

##Plot pyqtgraph 3d data
#if False: plot_pyqt(app, dim_red_data, cluster_labels)


##Plot trace averages and SDs
#selected_cluster, n_traces = plot_traces_pyqt(app, mouse, all_traces, cluster_labels, n_clusters)

##Plot DFF for individual trials
#plot_selected_cluster_DFF(mouse, cluster_labels, selected_cluster, n_traces)


print "Clean exit..."

quit()





