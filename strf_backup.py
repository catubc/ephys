import numpy as np

from tsf_ptcs_classes import *
from strf_utils import *
from movie import *

import multiprocessing as mp

global din


#***************
def Parallel_average(data):
    output = np.average(data,axis=0)
    return output
    
def find_previous_frame_index(din, target):
    index = np.argmin(np.abs(din[:,0] - target))
    if target< din[:,0][index]:
        return index-1
    else:    
        return index

def Parallel_find_previous_frame_index((targets)):
    global din
    indexes = []
    for i in range(len(targets)):
        index = np.argmin(np.abs(din[:,0] - targets[i]))
        
        if (targets[i] < din[:,0][index]):
            indexes.append(din[:,1][index-1])
        else:    
            indexes.append(din[:,1][index])
    return indexes

#*******************

n_procs = 10

main_dir = '/media/cat/12TB/in_vivo/nick/ptc21/'

#file_name = '60-tr5c-MVI_1419_5s'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-25/MVI_1419'

#file_name = '58-tr5c-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'
#Examples of good reverse correlation right out of nat scenes
#Unit :  7  No. spikes:  3651
#Unit :  10  No. spikes:  1195

#file_name = '25-tr2-MVI_1406-1410'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1406-1410'

file_name = '75-tr5c-mseq32_40ms'
#Unit 12 has 29k spikes and good single phase strf
#file_name = '93-t6b-mseq32_40ms'
movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

din_file = main_dir + file_name + '/' + file_name+'.din'

#**** Load spiketimes from .ptcs file ****
work_dir = main_dir + file_name + "/"
print "Loading: ", file_name
Sort = Loadptcs(file_name, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
Sort.name=file_name
Sort.filename=file_name
Sort.directory=work_dir
print "No. units: ", len(Sort.units)

#**** Load movie frames ****
movie = Movie()
movie.load(movie_file)
print "Length of movie (in frames): ", len(movie.frames)
#for i in range(0,len(movie.frames),10):
#    plt.imshow(movie.frames[i], cmap=cm.gray)
#    plt.show()
#quit()

#**** Subselect movie area ****
if 'MSEQ' in movie_file:
    temp_movie = np.zeros((len(movie.frames),10,10),dtype=np.float32)
    for i in range(len(movie.frames)):
        print "Subselecting image frame: ", i
        temp_movie[i] = movie.frames[i][16:26,14:24]
    movie.frames = temp_movie

else:
    temp_movie = np.zeros((len(movie.frames),40,40),dtype=np.float32)
    for i in range(len(movie.frames)):
        print "Subselecting frame: ", i
        temp_movie[i] = movie.frames[i][100:140,140:180]

    movie.frames = temp_movie

#**** Load .din stimulus timestamps
f = open(din_file, 'rb')
din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols
f.close()
print "Len of din file: ", len(din)

#plt.plot(sort(din[:,1],axis=0))
#plt.show()
#quit()

revcorr = False
strf_ridge = True

#Compute average image 
if revcorr:
    print "Computing average frame - parallel"
    chunks = int(len(movie.frames)/n_procs) #Break up the array into n_procs that are "chunk" long each; throw out residue!
    temp4=[]
    average_movie = []
    for i in range(n_procs):
        temp4.append(movie.frames[i*chunks:(i+1)*chunks])

    pool = mp.Pool(n_procs)
    average_movie.extend(pool.map(Parallel_average, temp4))
    pool.close()

    average_movie = np.average(average_movie, axis=0)


#Loop over all units and compute STRFs - spatial correlation only
#for unit in range(0,len(Sort.units),1):
for unit in [12]:#range(0,len(Sort.units),1):
    print "Unit : ", unit, " No. spikes: ", len(Sort.units[unit])
    
    frame_indexes = []
    print "Finding nearest frame indexes"
    spike_times = np.array(Sort.units[unit])*40 #Convert spike times to usec to match din times;
    correct_indexes = np.where(np.array(spike_times)<1E14)  #REmove weird 1E+19 values from the spiketimes; SS bug!
    spike_times = spike_times[correct_indexes]

    if len(spike_times)<n_procs: continue  #Skip if too few spikes - cant' be parallelized; should eventually process these units also
    
    if revcorr: 
        ##parallel computation frame_indexes
        chunks = int(len(spike_times)/n_procs) #Break up the array into n_procs that are "chunk" long each
        temp4=[]
        for i in range(n_procs):
            temp4.append(spike_times[i*chunks:(i+1)*chunks])
        
        pool = mp.Pool(n_procs)
        frame_indexes.extend(pool.map(Parallel_find_previous_frame_index, temp4))
        pool.close()
        
        #Remove blank frames not present in movie
        blank_indexes = np.where(np.array(frame_indexes)==65535)
        frame_indexes = np.delete(frame_indexes, blank_indexes)

        images_total=[]
        images_triggered = []
        images_ave = []
        lags = 20
        for j in range(lags):
            images_triggered.append([])
            images_ave.append([])
        
        print "Computing triggered images"
        for j in range(0,lags,1):
            print "Processing lag: ", j, " of total: ", lags
            for i in range(len(frame_indexes)):
                images_triggered[j].append(movie.frames[max(0,min(frame_indexes[i]+j-lags/2,len(movie.frames)-1))])
            
            #print "Averaging triggered frames"
            chunks = int(len(images_triggered[j])/n_procs) #Break up the array into n_procs that are "chunk" long each; discard residue!!
            temp4=[]
            for i in range(n_procs):
                temp4.append(images_triggered[j][i*chunks:(i+1)*chunks])

            pool_average = []
            pool = mp.Pool(n_procs)
            pool_average.extend(pool.map(Parallel_average, temp4))
            pool.close()
            
            images_ave[j] = np.average(pool_average,axis=0)
            images_ave[j] = images_ave[j] - average_movie
            images_total.append(images_ave[j])

        images_plot=images_total[0]
        blank_image = np.zeros((len(images_total[0]),3),dtype=np.uint8)+np.max(images_total)
        for i in range(1,len(images_total),1):
            images_plot=np.hstack((images_plot,blank_image))
            images_plot=np.hstack((images_plot,images_total[i]))
            
            
        plt.imshow(images_plot, cmap=cm.gray,interpolation='none')
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.title(file_name + " unit : " + str(unit) + " no. spikes: " + str(len(spike_times)))
        plt.tight_layout()
        plt.show()

    if strf_ridge:
        '''
        :param input: A numpy array of shape (num_time_points, num_spatial_channels).
        :param output: A numpy array of shape (num_time_points).
        :param lags: An array of integers that specify which time points should be included. For example,
                     to fit a STRF that uses only the stimulus information at time t to predict the output
                     at time t, specify lags=[0]. To fit a STRF that uses the stimulus information at
                     time t and the previous 10 time points, specify lags=range(0, 11).
        :param alpha: The regularization parameter for ridge regression. Try a bunch!
        :param verbose:
        :return: strf,bias: strf is a numpy array of shape (num_spatial_channels,len(lags)), which is the
                receptive field. bias is a scalar.'''

        #Downsample movies:
        if False:
            from skimage.measure import block_reduce
            reduced_movies = []
            for i in range(len(movie.frames)):
                reduced_movies.append(block_reduce(movie.frames[i], block_size=(2,2), func=np.mean))
            movie.frames = np.array(reduced_movies)

        #Make input_array from movie.frames
        print "Flattening movies"
        input_array = np.zeros((len(movie.frames), len(movie.frames[0].ravel())))
        for i in range(len(movie.frames)):
            input_array[i] = movie.frames[i].ravel()

        #Make output_array from spike_times; first compute frame-indexes for each spike
        chunks = int(len(spike_times)/n_procs) #Break up the array into n_procs that are "chunk" long each
        temp4=[]
        for i in range(n_procs):
            temp4.append(spike_times[i*chunks:(i+1)*chunks])

        frame_indexes = []
        pool = mp.Pool(n_procs)
        frame_indexes.extend(pool.map(Parallel_find_previous_frame_index, temp4))
        pool.close()

        #Remove blank frames not present in movie
        blank_indexes = np.where(np.array(frame_indexes)==65535)
        print "No. of blank_indexes: ", len(blank_indexes)
        frame_indexes = np.delete(frame_indexes, blank_indexes)
        
        bin_count = np.zeros(len(movie.frames),dtype=int32)
        for i in range(len(movie.frames)):
            bin_count[i]=len(np.where(frame_indexes==i)[0])
       
        output_array = np.array(bin_count)        

        #Reduce matrix sizes
        if False:
            input_array = input_array[0:4000]
            output_array = output_array[0:4000]
        
        print input_array.shape
        
        lags = [-1,0,1]
        alpha = [1.0, 6.5]
        
        #Serial ridge regression over alphas
        #strf, bias = fit_strf_ridge(input_array, output_array, lags=lags, alpha=alpha[0], verbose=True)
        #print shape(strf)
        #print np.array(strf).shape

        #Parallel ridge regression over alphas
        strf_array = []
        input_array = [input_array]*len(alpha)
        output_array = [output_array]*len(alpha)
        lags = [lags]*len(alpha)
        pool = mp.Pool(n_procs)
        strf_array.extend(pool.map(Parallel_fit_strf_ridge, zip(input_array, output_array, lags, alpha)))
        pool.close()

        #Loop over computed alphas
        for a in range(len(alpha)):
            print shape(strf_array)
            strf = np.array(strf_array)[a].T
            
            print strf.shape

            #Get strf product with original data for each lag
            strf_temp = []
            for i in range(len(lags)):
                print input_array[a][i].shape
                print strf[i].shape
                strf_temp.append(np.multiply(input_array[a][i],strf[i])
            
            strf_temp = np.array(strf_temp)
            print strf_temp.shape
            
            #Average image for each lag ?!
            image1d_average = []
            for i in range(len(lags)):
                image1d_average.append(np.average(strf_temp[i], axis=0))
                #image1d_average.append(strf_temp[i])

            print shape(image1d_average)

            #Convert 1d image back to 2d for display
            image2d = []
            for i in range(len(lags)):
                image2d.append(np.reshape(image1d_average[i], movie.frames[0].shape))
            image2d = np.array(image2d)
            
            #Add white separators between lags
            images_plot=image2d[0]
            blank_image = np.zeros((len(image2d[0]),3),dtype=np.uint8)+np.max(image2d)
            for i in range(1,len(image2d),1):
                images_plot=np.hstack((images_plot,blank_image))
                images_plot=np.hstack((images_plot,image2d[i]))
            
            plt.imshow(images_plot, cmap=cm.gray, interpolation='none')
            plt.title(file_name + " unit : " + str(unit) + " no. spikes: " + str(len(spike_times)) + " lags: " + str(lags[0])+ " alpha: " + str(alpha))

            plt.show()
        
        

        



















































