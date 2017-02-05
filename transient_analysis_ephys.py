from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(formatter={'float': '{: 0.6f}'.format})

#************************* NAT SCENE 5 SEC REPEATS ***************
#work_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'
#file_name = '04-tr1-MVI_1400_5s'

#work_dir = '/media/cat/12TB/in_vivo/nick/ptc17/'
#file_name = '36-tr1-MVI_1419_5s_66Hz'

#work_dir = '/media/cat/12TB/in_vivo/nick/ptc17/'
#file_name = '37-tr1-MVI_1419_5s_200Hz'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'
#file_name = '51-tr2b-MVI_1406-1410'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'
#file_name = '57-tr2b-MVI_1400_5s'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'
#file_name = '58-tr2b-MVI_1400_5s_66Hz'
#state_change = 0 #Seconds into recording that state change occurs

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '60-tr5c-MVI_1419_5s'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '62-tr5c-MVI_1403_5s'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '64-tr5c-MVI_1419_5s'

work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
file_name = '66-tr5c-MVI_1403_5s'

#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '88-t6b-MVI_1419_5s'

#*********************** NAT SCENE LONG MOVIES *********************
#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '58-tr5c-MVI_1403-1405'

#*********************** DRIFT BAR *********************
#work_dir = '/media/cat/12TB/in_vivo/nick/ptc21/'
#file_name = '56-tr5c-driftbar_shortbar'

##*********************** SPONTANEOUS ACTIVITY *******************
#work_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#file_name = '61-tr5c-blankscreen'

#*******************************************************************************************************************************
#*************************************************** BEGIN PROCESSING **********************************************************
#*******************************************************************************************************************************
#Exp type:
nat_scene_5sec = False
nat_scene_long = False
spontaneous = False
drift_bar = False

if '5s' in file_name:
    nat_scene_5sec = True
elif 'MVI' in file_name:
    nat_scene_5sec = True
    nat_scene_long = True
if 'blankscreen' in file_name:
    spontaneous = True
if 'drift' in file_name:
    drift_bar = True

rasters = True
lfp = False

#*********************** LOAD .DIN DATA STIMULUS *********************
#Load .din from nat scenes files or just assign random trigger time in spont activity
#NB: Here can set spontaenous start/end times as well;
if nat_scene_5sec:
    stim_times = Load_din_5sec(work_dir, file_name)  #TODO!!!! Find zero-frame manually!!!
    #movie_frames = [100,101,102,103,105]
    movie_frames = np.arange(0,400,1)
    spont_activity_length=1.0
elif nat_scene_long:
    stim_times = Load_din_5sec(work_dir, file_name)  #TODO!!!! Find zero-frame manually!!!
    print "stim_times[0]: ", stim_times*1E6
    movie_frames = [100,101]
    spont_activity_length=0.0
elif spontaneous:
    stim_times = [[100.0, 150.0]]  #Pick a period of time to look at in secs
    movie_frames= [0] #spont recs only have 1 trial
    spont_activity_length=0.0
elif drift_bar:
    experiment = 7
    stim_times = Load_din_driftbar(work_dir, file_name, experiment)  #TODO!!!! Find zero-frame manually!!!
    print "stim_times[0]: ", stim_times*1E6
    movie_frames = [0]
    spont_activity_length=0.0

#*********************** LOAD TSF DATA *********************
#tsf_name = work_dir + file_name+'/'+file_name+'.tsf'
#f = open(tsf_name, "rb")
#tsf = Tsf_file(f, work_dir+file_name+'/')  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
#f.close()

#******************** PARAMETERS *********************
millisecs_load = stim_times[0][1]*1E3-stim_times[0][0]*1E3 + spont_activity_length*1E3 #How much time to consider from stim_start

#******************** LOAD SINGLE UNIT DATA ******************
ptcs_dir = work_dir + file_name + '/'
print "Loading: ", file_name
Sort = Loadptcs(file_name, ptcs_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
Sort.name=file_name
Sort.filename=file_name
Sort.directory=work_dir
print "No. units: ", len(Sort.units)

#*******************************************************************************
#********************************* LFP DATA LOAD *******************************
#*******************************************************************************
if lfp:
    bin_width = 1 #bin width for computing vectors from lfp

    #Load Martin's .npy ziped lfp data
    print "Loading lfp ..."
    temp_data = np.load(work_dir+file_name+'/'+file_name+'.lfp.zip')
    #print temp_data.keys()
    
    vscale = temp_data['uVperAD']
    lfp_data=[]
    for i in range(len(temp_data['chans'])): #Loop over # of lfp chans
        lfp_data.append(filter.notch(temp_data['data'][i]*vscale)[0])
        #plt.plot(lfp_data[i]+1000*i)#, linewidth=2)
    #plt.show()

    def round_down(num, divisor):
        return num - (num%divisor)
    
    #Bin LFP in a particular band
    lowcut=30
    highcut=100
    vector_array = []
    lfp_matrix = []
    for i in movie_frames:  #NB: FIX THIS; 1) append multiple vectors; 2) remove length_seg variable - use only start_stim
        print ""
        print "LFP for movie frame: ", i
        temp_vector = Compute_binned_lfp(lfp_data, stim_times[i], spont_activity_length, bin_width, lowcut, highcut)
        temp_vector = temp_vector[:round_down(len(temp_vector),10)]
        vector_array.extend(temp_vector)
        lfp_matrix.append(temp_vector)
    
    vector_data = vector_array
    print "All movie array shape lfp: ", np.array(vector_data).shape
    print "LFP matrix ", np.array(lfp_matrix).shape
        
    #Plot LFP 
    #LFP_PSTH_plot(lfp_matrix, bin_width, work_dir, file_name, movie_frames, lowcut, highcut, spont_activity_length)
    

#*************************************************************************************
#********************************* RASTER DATA LOAD  *********************************
#************************************************************************************* 
if rasters: 
    bin_width = 25 #bin width for computing firing rates

    #Plot SUA rasters for repeat movies
    Plot_stimulus_triggered_rasters(Sort, stim_times)

    #Computes spike times stim_start..stim_end for all repeats in movies; or for entire length of spont rec.
    save_response_array_individual = Compute_responses(Sort, stim_times, spont_activity_length) 

    if True:
        np.save(work_dir+file_name[:2]+'_binned_raster',save_response_array_individual)
    
    if True:
        file_load = [64, 88]
        temp_data = []
        for p in file_load:
            load_matrix = np.load(work_dir+str(p)+'_binned_raster.npy')
            print load_matrix.shape
            temp_data.append(load_matrix)
        
        out_data = temp_data[0]
        for i in range(len(temp_data)-1):
            out_data = np.vstack((out_data,temp_data[i+1]))
        temp_data = out_data

    print temp_data.shape

    binned_array = np.array(save_response_array_individual)
    print binned_array.shape
    
    joined_array = np.concatenate((binned_array, temp_data), axis=0)
    print joined_array.shape
    #quit()

    #Plot SUA PSTH
    PSTH_plot(millisecs_load, bin_width, joined_array, work_dir, file_name, movie_frames)

    #Vectorize Rasters
    vector_data = Compute_firing_rates(joined_array, millisecs_load, bin_width)
    
    #print vector_matrix.shape
    #quit()

    vector_array = []
    vector_matrix = []
    for i in movie_frames:
        vector_array.extend(vector_data[i])
        vector_matrix.append(vector_data[i])
    vector_data = vector_array

    print "All movie array shape: ", np.array(vector_data).shape
    #quit()

    #Display rasters
    if False:
        vector_img=[]
        for i in [42,  48,  54,  58,  79, 115]:#range(10,30,1):#len(vector_array)):
            vector_img.append(vector_array[i])
            print vector_array[i]
            
        plt.imshow(np.array(vector_img).T)
        plt.show()
        quit()

#*************************************************************************************
#********************************* CLASSIFIERS  **************************************
#************************************************************************************* 

if True:
    print "Running classifiers..."

    #Training data:
    X = np.array(vector_array)
    
    print "Vector matrix size: ", np.array(vector_matrix).shape

    #Training classes:
    Y = []
    for i in range(len(vector_matrix)):
        Y.extend(np.arange(0,len(vector_matrix[i]),1))
    Y=np.array(Y)
    print X.shape, Y.shape

    print "Training classifier..."
    if False:
        #AdaBoosted decision tree
        from sklearn.ensemble import AdaBoostClassifier
        from sklearn.tree import DecisionTreeClassifier
    
        print "Adaboost classifier"
        clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=3),
                                 algorithm="SAMME",
                                 n_estimators=500)
        print "fitting data"
        clf.fit(X, Y)

    if False:
        #AdaBoosted decision tree
        from sklearn.ensemble import AdaBoostRegressor
        from sklearn.tree import DecisionTreeRegressor
    
        print "Adaboost regressor"
        rng = np.random.RandomState(1)
        clf = AdaBoostRegressor(DecisionTreeRegressor(max_depth=30),
                                n_estimators=300, random_state=rng)
        print "fitting data"
        clf.fit(X, Y)
        
    if True:
        from sklearn import linear_model

        print "BayesianRidge"
        clf = linear_model.BayesianRidge()
        clf.fit(X, Y)

    if False:
        from sklearn import linear_model

        print "Perceptron"
        clf = linear_model.Perceptron()
        clf.fit(X, Y)

    #Plot results as bar graph
    ax = plt.subplot(111)
    for frame_number in range(len(vector_matrix[0])):
        frame_index = np.arange(frame_number,len(Y),len(movie_frames))
        frame_result = clf.predict(X[frame_index])
        counter = 0
        for i in range(len(frame_index)):
            if frame_index[i]==frame_result[i]: counter+=1
        
        temp_float = counter/float(len(frame_index))
        plt.bar(frame_number*bin_width,temp_float, width=1, color='green')
    
    plt.xlim(0,len(vector_matrix[i]))
    plt.ylim(0,0.1)
    plt.show()

    quit()



#*************** SELECT AND CONCATENATE MOVIE FRAME RATE VECTORS (+ FIND DUPLICATES) ************
#concatenate firing rates into long stack; Not sure required, MDS might return same absolute coords for independently split data

#index_duplicates,duplicates = Find_duplicates(cumulative_array)

#for k in range(len(index_duplicates)):
#    print index_duplicates[k]
#    print duplicates[k]

#*****************************************************************************
#*********************** DIMENSIONALITY REDUCTION  ***************************
#*****************************************************************************
if True:
    #Method: 0 = MDS SMACOF;  1 = t-SNE; 2 = PCA; 3 = Sammon
    method = 0

    #Save output to file (or load from file if existing)
    data_type = 'rasters'
    if lfp: data_type = 'lfp'
    fname = work_dir+file_name+'/'+data_type+'_dimred_method_'+str(method)+"_startframe_"+str(movie_frames[0])+'_nframes_'+str(len(movie_frames))+'_binwidth_'+str(bin_width)+"ms"

    if os.path.isfile(fname+'.npy'):
        multi_trial = np.load(fname+'.npy') 
        single_trial= multi_trial[0]

        print "Loaded saved dim reduction file (shape: ", multi_trial.shape, ")"
    else:
        multi_trial = []
        single_trial = []
        #for i in range(len(vector_data)):
        
        #Do not loop over movie repeats; concatenate all movies and do single dim reduction.
        data_out = multi_dim_scaling(vector_data, method)
        
        #Chunk data back out to represent 
        multi_trial = np.split(data_out, len(movie_frames))
        single_trial= multi_trial[0]

        #Save data for future processing
        np.save(fname, multi_trial)


#*****************************************************************************
#******************** CLUSTER TIME POINTS (MPL) AND VISUALIZE ****************
#*****************************************************************************
if False:
    loaded_list = cluster_matplotlib(multi_trial[0], fname)
    
    ##if loading from file
    #import cPickle
    #file_name = fname+"_clusters.txt"
    #loaded_list = cPickle.load(open(file_name, 'rb'))
   
    clusters = []
    for i in range(len(loaded_list)):
        clusters.append(loaded_list[i])

    #Display clusters as cluster in image 
    if False:
        for k in range(len(clusters)):
            vector_img=[]
            for i in clusters[k]:#len(vector_array)):
                vector_img.append(vector_array[i])
           
            plt.title("Cluster: "+str(k)+ "  #pts: " + str(len(clusters[k])))
            plt.imshow(np.array(vector_img).T)
            plt.show()

    #Display clusters as shaded bin in SUA rasters
    if True: 
        for i in range(1,len(clusters),1): #Skip 0 cluster that contains Zero firing rate vector
            ax = plt.subplot(1, 1, 1)
            #print len(clusters[i])
            #print clusters[i][0]
            #print clusters[i][1]
            #quit()

            #Plot rasters during stim_times (e.g. [[100.0, 150.0]])
            for j in range(len(Sort.units)):
                temp_array = np.array(Sort.units[j])
                t0= stim_times[0][0]*Sort.samplerate
                t1= stim_times[0][1]*Sort.samplerate
                temp_index = np.where(np.logical_and(temp_array>=t0, temp_array<=t1))[0]  #Convert to sample units

                plt.scatter((temp_array[temp_index])/Sort.samplerate*1E3,[j]*len(temp_index), s = 10, color='blue') #units in ms

            for k in range(len(clusters[i])):
                p = ax.axvspan(clusters[i][k]*bin_width +stim_times[0][0]*1E3, clusters[i][k]*bin_width+bin_width+stim_times[0][0]*1E3, facecolor='red', alpha=0.15)
                
            #plt.xlim(0, stim_times[0][1]-stim_times[0][0])
            #plt.ylim(n_repeats,-0.5)
            plt.show()


    quit()
    

#*****************************************************************************
#************************** OPENGL ROUTINES **********************************
#*****************************************************************************
if True:
    from opengl_pyqt_classes import *
    if __name__ == '__main__':
        
        # prevents "The event loop is already running" messages when calling ipshell():
        QtCore.pyqtRemoveInputHook()
        app = QtGui.QApplication(sys.argv) #Cat: app start

        mainwindow = MainWindow()
        mainwindow.show()
        print "number of time segs: ", len(movie_frames)

        mainwindow.coords = multi_trial[0]
        mainwindow.single_trial = single_trial
        mainwindow.multi_trial = multi_trial
        mainwindow.repeats = len(movie_frames)
        mainwindow.file_name = file_name#+" len: " + str(length_segment*1E-3)+"sec  bin: " + str(bin_width)+"ms "
        
        sys.exit(app.exec_()) #loop here using .exec_...


#*********************************************************************


#Matplotlib 3d 
#fig = plt.figure()
#ax = fig.gca(projection='3d')

#multi_trials=[]

##Compute MDS reduction individually for each movie repeat
#if False:
    #for i in range(3):
    ##for i in range(len(cell_rate_trans)):

        #print "Movie repeat: ", i
        #dists = sklearn.metrics.pairwise.pairwise_distances(cell_rate_trans[i])

        #adist = np.array(dists)
        #amax = np.amax(adist)
        #adist /= amax

        #mds = manifold.MDS(n_components=3, dissimilarity="euclidean", random_state=6)
        #results = mds.fit(adist)

        #coords = results.embedding_

        #if plotting:
            #for i in range(len(coords)-1):
                #loc_1 = [coords[i,0],coords[i,1],coords[i,2]]
                #loc_2 = [coords[i+1,0],coords[i+1,1],coords[i+1,2]]

                #ax.plot([loc_1[0],loc_2[0]],[loc_1[1],loc_2[1]],[loc_1[2],loc_2[2]],color=(float(i)/len(coords),0,float(i)/len(coords)))
                #ax.scatter(coords[i, 0],coords[i, 1],coords[i, 2], s=20, color=(float(i)/len(coords),0,float(i)/len(coords)))
        #multi_trials.append(coords*1E3)

    #if plotting:
        #plt.show()

