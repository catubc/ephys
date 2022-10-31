import numpy as np
import multiprocessing as mp

from tsf_ptcs_classes import *
from strf_utils import *
from skimage import feature

global din, lum_contrast_movie, movie_frames, frame_indexes


n_procs = 6

def Parallel_find_previous_frame_index((spikes)):
    global din
    indexes = []
    for i in range(len(spikes)):
        index = np.argmin(np.abs(din[:,0] - spikes[i]))
        
        if (spikes[i] < din[:,0][index]):
            indexes.append(din[:,1][index-1])
        else:    
            indexes.append(din[:,1][index])
    return indexes

def Find_previous_frame_index(spikes):
    global din
    indexes = []
    for i in range(len(spikes)):
        index = np.argmin(np.abs(din[:,0] - spikes[i]))
        
        if (spikes[i] < din[:,0][index]):
            indexes.append(din[:,1][index-1])
        else:    
            indexes.append(din[:,1][index])
            
    return indexes

def Parallel_average(data):
    output = np.average(data,axis=0)
    return output

def Luminance_contrast(frame_in, frame_out):
    frame_ave = np.average(frame_in) 
    #for i in range(0,len(frame_in),1):
        #for j in range(0,len(frame_in[0]),1):
            #frame_out[i,j] = frame_in[i][j]/frame_ave   #Definition of luminance-contrast from Ringach 2002
                                                        ##This should be pythonized!!!
    frame_out=frame_in/frame_ave
    
    return frame_out

def Edge_detection(frame_in, frame_out):
    frame_out = feature.canny(frame_in)

    return frame_out

def Parallel_recleastsq((tau)):
    global lum_contrast_movie, movie_frames, frame_indexes
    
    n_frames = len(lum_contrast_movie)
    delta = 1E6
    P = delta * (np.zeros(len(movie_frames[0]), dtype=np.float32)+1.0)
    w = np.zeros(len(movie_frames[0]), dtype=np.float32)
    #N = len(movie_frames[0])*len(movie_frames[0])
    mu = 0
    beta = 0.99
    
    frame_spikes = []       #Keep track of the number of spikes fired in each frame (+ tau)
    #frame_out = np.zeros((len(movie_frames[0]), len(movie_frames[0][0])),dtype=np.float32) #Initialize once and zero every frame
    for i in range(0,n_frames,1): #Loop over movie frames
        #print "Procesing frame : ", i
        I_ijn = lum_contrast_movie[i]                                                       #Line 5 
        #I_ijn = movie_frames[i]
        
        u_n = I_ijn                                                                         #Line 6: no need to subsample unless using entire frame to start
        #Compute number of spikes that fall in the frame
        #r_n = len(np.where(np.logical_and(frame_indexes>=i+tau, frame_indexes <= i+tau+1))[0])      #Line 7
        r_n = len(np.where(frame_indexes==i+tau)[0])      #Line 7
        #print "No. spikes: ", r_n
        frame_spikes.append(r_n)
        mu = beta*mu + (1-beta)*r_n                                                         #Line 8
        r_n = r_n - mu                                                                      #Line 9
        
        I_n = u_n * P                                                                       #Line 10
        
        kappa_n = 1 + I_n*u_n                                                               #Line 11
        
        k_n = (P*u_n)/kappa_n                                                               #Line 12
        
        #alpha_n = r_n -w * u_n                                                              #Line 13
        #w = w+ k_n * alpha_n                                                                #Line 15
        
        T = k_n * I_n                                                                       #Line 16
        
        P = P - T
    
    #Decorrelate the stimulus
    decorrelated_frames = []
    for i in range(0,n_frames,1): #Loop over movie frames
        temp = P * movie_frames[i]
        decorrelated_frames.append(temp)
    
    strf = np.zeros((len(movie_frames[0]), len(movie_frames[0][0])),dtype=np.float32)
    for i in range(0,n_frames,1):
        strf+= decorrelated_frames[i] * frame_spikes[i]
    
    return strf


class nat_scene_movie():
    
    def load(self, file_name, asarray=False, flip=False):
        """Load movie frames"""
        # figure out the local path to the same movie file:
        #pathparts = core.pathdecomp(self.static.fname) # as it existed on the stim computer
        #movi = pathparts.index('mov')
        #tail = os.path.join(pathparts[movi+1:]) # everything after 'mov' folder
        #MOVIEPATH = get_ipython().user_ns['MOVIEPATH']
        #fullfname = os.path.join(MOVIEPATH, *tail) # full fname with local MOVIEPATH
        self.f = file(file_name, 'rb') # open the movie file for reading in binary format
        headerstring = self.f.read(5)
        if headerstring == 'movie': # a header has been added to the start of the file
            self.ncellswide, = struct.unpack('H', self.f.read(2)) # 'H'== unsigned short int
            self.ncellshigh, = struct.unpack('H', self.f.read(2))
            self.nframes, = struct.unpack('H', self.f.read(2))
            if self.nframes == 0:
                # this was used in ptc15 mseq movies to indicate 2**16 frames, shouldn't
                # really worry about this, cuz we're using slightly modified mseq movies now
                # that don't have the extra frame at the end that the ptc15 movies had (see
                # comment in Experiment module), and therefore never have a need to indicate
                # 2**16 frames
                self.nframes = 2**16
            self.offset = self.f.tell() # header is 11 bytes long
        else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
            self.f.seek(0)
            self.ncellswide = mseq.ncellshigh = 64
            self.nframes = 6000
            self.offset = self.f.tell() # header is 0 bytes long
        self.framesize = self.ncellshigh*self.ncellswide

        # read in all of the frames. Maybe check first to see if file is > 1GB. If so,
        # _loadaslist() to prevent trying to allocate one huge piece of contiguous memory and
        # raising a MemoryError, or worse, segfaulting
        if asarray:
            self._loadasarray(flip)
        else:
            self._loadaslist(flip)
        leftover = self.f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            pprint(leftover)
            print(self.ncellswide, self.ncellshigh, self.nframes)
            raise RuntimeError('There are unread bytes in movie file %r. Width, height, or nframes is incorrect in the movie file header.' % self.static.fname)
        self.f.close() # close the movie file

    def _loadasarray(self, flip):
        self.frames = np.fromfile(self.f, dtype=np.uint8, count=self.nframes*self.framesize)
        self.frames.shape = (self.nframes, self.ncellshigh, self.ncellswide)
        if flip:
            self.frames = self.frames[::, ::-1, ::] # flip all frames vertically for OpenGL's bottom left origin

    def _loadaslist(self, flip):
        self.frames = []
        for framei in xrange(self.nframes): # one frame at a time...
            frame = np.fromfile(self.f, dtype=np.uint8, count=self.framesize) # load the next frame
            frame.shape = (self.ncellshigh, self.ncellswide)
            if flip:
                frame = frame[::-1, ::] # flip all frames vertically for OpenGL's bottom left origin
            self.frames.append(frame)
            

#***********************************************************************************************
#*************************************** PTC17 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/12TB/in_vivo/nick/ptc17/'

#file_name = '03-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '06-tr1-MVI_1406-1410'; movie_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/2007-11-24/MVI_1406-1410'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'

#file_name = '12-tr1-MVI_1413-1416_128x128x1024'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/f-1/MVI_1413-1416_128x128x1024'

#file_name = '25-tr1-MVI_1398-1399_200Hz'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1398-1399'

#file_name = '31-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#***********************************************************************************************
#*************************************** PTC18 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/12TB/in_vivo/nick/ptc18/'

#file_name = '33-tr1-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'

#***********************************************************************************************
#*************************************** PTC21 *************************************************
#***********************************************************************************************

#file_name = '60-tr5c-MVI_1419_5s'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-25/MVI_1419'

#file_name = '58-tr5c-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'
#Examples of good reverse correlation right out of nat scenes
#Unit :  7  No. spikes:  3651
#Unit :  10  No. spikes:  1195

#file_name = '25-tr2-MVI_1406-1410'; movie_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/2007-11-24/MVI_1406-1410';  main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

#file_name = '75-tr5c-mseq32_40ms'
#Unit 12 has 29k spikes and good single phase strf

#file_name = '93-t6b-mseq32_40ms';  movie_file = '/media/cat/4TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '57-tr5c-mseq32_40ms'; movie_file = '/media/cat/4TB/mseq/MSEQ32'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

#file_name = '60-tr5c-MVI_1419_5s'; movie_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/2007-11-25/MVI_1419';  main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

file_name = '59-tr5c-MVI_1403-1405_66Hz'; movie_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/2007-11-24/MVI_1403-1405';  main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'


#***********************************************************************************************
#*************************************** PTC22 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/12TB/in_vivo/nick/ptc22/'

#file_name = '04-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '05-tr1-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'

#file_name = '17-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '19-tr1-MVI_1406-1410'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1406-1410'


#************MAIN CODE***********

#**** Load spiketimes from .ptcs file ****
work_dir = main_dir + file_name + "/"
print "Loading: ", file_name
Sort = Loadptcs(file_name, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
Sort.name=file_name
Sort.filename=file_name
Sort.directory=work_dir
print "No. units: ", len(Sort.units)

#**** Load .din stimulus timestamps
din_file = main_dir + file_name + '/' + file_name+'.din'
f = open(din_file, 'rb')
din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols containing [time_stamp in usec : frame number]
f.close()
print "No of frames in din file: ", len(din)
print "Length of recording: ", float(din[-1][0])*1E-6/60., " mins." 

#**** Load movie frames ****
movie = nat_scene_movie()
movie.load(movie_file)
print "Length of movie (in frames): ", len(movie.frames)
#for i in range(0,len(movie.frames),10):
#    plt.imshow(movie.frames[i], cmap=cm.gray)
#    plt.show()
#quit()

#**** Subselect movie area ****
if True:
    width = 50
    print "Selecting smaller visual area : ", width, " x ", width 
    if 'MSEQ' in movie_file:
        pass
    #    temp_movie = np.zeros((len(movie.frames),10,10),dtype=np.float32)
    #    for i in range(len(movie.frames)):
    #        #print "Subselecting image frame: ", i
    #        temp_movie[i] = movie.frames[i][16:26,14:24]
    #    movie.frames = temp_movie

    else:
        start_index_w = (len(movie.frames[0])-width)/2
        end_index_w = (len(movie.frames[0])-width)/2+width
        start_index_h = (len(movie.frames[0][0])-width)/2
        end_index_h = (len(movie.frames[0][0])-width)/2+width
        temp_movie = np.zeros((len(movie.frames),width,width),dtype=np.float32)
        for i in range(len(movie.frames)):
            #print "Subselecting frame: ", i
            temp_movie[i] = movie.frames[i][start_index_w:end_index_w,start_index_h:end_index_h]
        movie.frames = temp_movie


#*******************************************************************************************************
#*************************************** COMPUTE STRF **************************************************
#*******************************************************************************************************

if 'MSEQ' in movie_file:
    revcorr = True
    recursive_leastsq = False
else:
    revcorr = False
    recursive_leastsq = True
#strf_ridge = False

#unit = 12  #Must pick correct cell from sorted data to look at


#Make luminance-contrast movies for natural scenes:
if 'MSEQ' in movie_file:
    pass
else:
    #Make Luminance Contrast Movies
    frame_out = np.zeros((len(movie.frames[0]), len(movie.frames[0][0])),dtype=np.float32) #Initialize once and zero every frame
    if (os.path.exists(movie_file+"_luminance_contrast_"+str(width)+"pixels.npy")==False):
        lum_contrast_movie = np.zeros((len(movie.frames),len(movie.frames[0]),len(movie.frames[0][0])),dtype=np.float32)
        for i in range(len(movie.frames)):
            print "Luminance contrast conversion frame: ", i, " / ", len(movie.frames), " video size: ", width, " x ", width
            lum_contrast_movie[i] = Luminance_contrast(movie.frames[i], frame_out)
            
        np.save(movie_file+"_luminance_contrast_"+str(width)+"pixels.npy", lum_contrast_movie)
        
    #Load luminance-contrast movie from file
    else:
        lum_contrast_movie = np.load(movie_file+"_luminance_contrast_"+str(width)+"pixels.npy")
        print "Length of lum_constrast saved movie: ", len(lum_contrast_movie)

    #Make Edge Map Movies
    frame_out = np.zeros((len(movie.frames[0]), len(movie.frames[0][0])),dtype=np.float32) #Initialize once and zero every frame
    if (os.path.exists(movie_file+"_edge_"+str(width)+"pixels.npy")==False):
        edge_map_movie = np.zeros((len(movie.frames),len(movie.frames[0]),len(movie.frames[0][0])),dtype=np.float32)
        for i in range(len(movie.frames)):
            print "Edge conversion frame: ", i, " / ", len(movie.frames), " video size: ", width, " x ", width
            edge_map_movie[i] = Edge_detection(movie.frames[i], frame_out)
            
        np.save(movie_file+"_edge_"+str(width)+"pixels.npy", edge_map_movie)
        
    #Load luminance-contrast movie from file
    else:
        edge_map_movie = np.load(movie_file+"_edge_"+str(width)+"pixels.npy")
        print "Length of edge movie: ", len(edge_map_movie)

#lum_contrast_movie = edge_map_movie

#ax = plt.subplot(121)
#plt.imshow(edge_map_movie[0])
#ax = plt.subplot(122)
#plt.imshow(movie.frames[0])
#plt.show()


if False:
    non_lock_spikes_array = np.load(main_dir+file_name+'/'+file_name+'_non_locked_spikes.npy')
    lock_spikes_array = np.load(main_dir+file_name+'/'+file_name+'_locked_spikes.npy')
    lfp_events_array = np.load(main_dir+file_name+'/'+file_name+'.lfp.events.npy')

#****************************************************************************************************************************
#****************************************************************************************************************************
#****************************************************************************************************************************
min_spikes = n_procs
tau = np.arange(-8,20,1)  #This sets the lags for reclsq and also routines below

#***************************** RECURSIVE LEAST SQUARES
if recursive_leastsq:         
    
    movie_frames = movie.frames

    #unit_list = [11, 13]
    unit_list = np.arange(10,20,1) #len(Sort.units),1)
    total_img = []
    unit_name = []
    unit_nspikes = []
    for unit in unit_list: #range(10, 15 ,1): #len(Sort.units),1):

        #Find nearest frame index
        spike_times = np.array(Sort.units[unit])*40 #Convert spike times to usec to match din times;
        spike_times = spike_times[spike_times<1E12]   #Remove weird 1E+19 values from the spiketimes; SS bug!
        #non_lock_spikes = non_lock_spikes_array[unit]*1E6
        #non_lock_spikes = non_lock_spikes[non_lock_spikes<1E12]
        #lock_spikes = lock_spikes_array[unit]*1E6
        #lock_spikes = lock_spikes[lock_spikes<1E12]
        
        if len(spike_times)<10: continue
        
        print "Processing unit: ", unit, " unique unit id: ", Sort.uid[unit], " no. spikes: ", len(spike_times) # \
        #, " non_lock_spikes: ", len(non_lock_spikes), " lock_spikes: ", len(lock_spikes)
        
        #test_spikes = [all_spikes, non_lock_spikes, lock_spikes]
        #for spike_times in test_spikes:
        frame_indexes = []
        if len(spike_times)>min_spikes:
            ##parallel computation frame_indexes
            chunks = int(len(spike_times)/n_procs) #Break up the array into n_procs that are "chunk" long each
            temp4=[]
            for i in range(n_procs):
                temp4.append(spike_times[i*chunks:(i+1)*chunks])
            
            pool = mp.Pool(n_procs)
            frame_indexes.extend(pool.map(Parallel_find_previous_frame_index, temp4))
            pool.close()
        else:
            frame_indexes = Find_previous_frame_index(spike_times)
            
        frame_indexes = np.sort(np.array(frame_indexes).ravel(),axis=0)

        
        #n_frames = len(lum_contrast_movie)
        multi_frame_img = []
        pool = mp.Pool(n_procs)
        multi_frame_img.extend(pool.map(Parallel_recleastsq,tau))
        pool.close()
        
        #for k in range(len(multi_frame_img)):
            #temp1_img = np.array(multi_frame_img[k])-127
            #v_abs=np.max(np.abs(temp1_img))
            #multi_frame_img[k] = temp1_img/v_abs

        img_plot = multi_frame_img[0]
        for i in range(1,len(multi_frame_img),1):
            img_plot = np.hstack((img_plot,multi_frame_img[i]))

        #img_plot = (img_plot - np.min(img_plot))/(np.max(img_plot)-np.min(img_plot))

        total_img.append(img_plot)
        unit_name.append(Sort.uid[unit])
        unit_nspikes.append(len(spike_times))

    #Vstack all the cells into large image; 
    img_out = total_img[0]
    for j in range(1, len(total_img),1):
        img_out = np.vstack((img_out, total_img[j]))
    
    #Average every 4 frames:
    img_out2 = []
    for k in range(7):
        temp_array = []
        for p in range(4):
            temp_array.append(img_out[:, k*200+p*50:k*200+(p+1)*50])
        img_out2.append(np.average(np.array(temp_array), axis=0))
    
    #Horizontal restack of the averaged chunks
    img_out = img_out2[0]
    for j in range(1, len(img_out2),1):
        img_out = np.hstack((img_out, img_out2[j]))
    
    #Invert the total_img
    img_out = np.fliplr(img_out)
    
    #Add separators to data:
    for k in range(8):
        img_out[:,k*50-1:k*50+1]=0.
    for j in range(img_out.shape[0]/50):
        img_out[j*50-1:j*50+1,:]=0.

    ax = plt.subplot(111)
    ax.get_yaxis().set_ticks([])
    x_original = np.arange(0,len(img_out[0])+50,len(img_out)/len(unit_nspikes))
    x_label = np.int32(np.arange(-100,41,20))
    plt.xticks(x_original, x_label, fontsize=20)
    plt.xlabel("Time course of RF (ms)", fontsize=20)
    
    y_original = np.arange(25,len(img_out),50)
    y_label = []
    for k in range(len(unit_nspikes)):
        y_label.append("Unit: "+str(unit_name[k])+"\n#: "+str(unit_nspikes[k]))
    plt.yticks(y_original, y_label, rotation=90, fontsize=20)
    
    plt.title(file_name+ "  Unit: "+str(unit), fontsize=25)
    
    plt.imshow(img_out, interpolation='none')
    plt.show()

    
    quit()
    plt.suptitle(file_name)
    plt.draw()
    fig1.savefig("/home/cat/"+file_name+".png", dpi=2400)#, strf)        
    plt.close()
    quit()

#***************************** REVERSE CORRELATION **************
if revcorr:
    
    total_img = []
    unit_name = []
    unit_nspikes = []
    
    print "Computing average frame - parallel"
    chunks = int(len(movie.frames)/n_procs) #Break up the array into n_procs that are "chunk" long each; throw out residue!
    temp4=[]
    for i in range(n_procs):
        temp4.append(movie.frames[i*chunks:(i+1)*chunks])

    average_movie = []
    pool = mp.Pool(n_procs)
    average_movie.extend(pool.map(Parallel_average, temp4))
    pool.close()

    average_movie = np.average(average_movie, axis=0)

    for unit in range(0,len(Sort.units),1):
        print "Unit : ", unit, " No. spikes: ", len(Sort.units[unit])
        
        frame_indexes = []
        print "Finding nearest frame indexes"
        spike_times = np.array(Sort.units[unit])*40 #Convert spike times to usec to match din times;
        correct_indexes = np.where(np.array(spike_times)<1E14)  #Remove weird 1E+19 values from the spiketimes; SS bug!
        spike_times = spike_times[correct_indexes]

        #if len(spike_times)<n_procs: continue  #Skip if too few spikes - cant' be parallelized; should eventually process these units also
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
        for j in range(len(tau)):
            images_triggered.append([])
            images_ave.append([])
        
        print "Computing triggered images"
        for j in tau:
            print "Processing lag: ", j, " of total: ", len(tau)
            for i in range(len(frame_indexes)):
                images_triggered[j].append(movie.frames[max(0,min(frame_indexes[i]+j,len(movie.frames)-1))])
            
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

        #Load images in horizontal stack
        images_plot=images_total[0]
        blank_image = np.zeros((len(images_total[0]),3),dtype=np.uint8)+np.max(images_total)
        for i in range(1,len(images_total),1):
            images_plot=np.hstack((images_plot,blank_image))
            images_plot=np.hstack((images_plot,images_total[i]))
            

        total_img.append(images_plot)
        unit_name.append(Sort.uid[unit])
        unit_nspikes.append(len(spike_times))
        
    fig1 = plt.gcf()
    n_units = len(total_img)
    for i in range(n_units):
        ax=plt.subplot(max(Sort.uid)+1,1,Sort.uid[i]+1)
        plt.imshow(total_img[i], cmap=cm.jet, interpolation='none')
        if unit_nspikes<min_spikes:
            plt.ylabel(str(unit_name[i]), fontsize = 1, color='red')
        else:
            plt.ylabel(str(unit_name[i]), fontsize = 1)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

    plt.suptitle(file_name)
    plt.draw()
    #fig1.savefig("/home/cat/"+file_name+".png", dpi=2400)#, strf)        
    fig1.savefig("/home/cat/"+file_name+".png", dpi=2400)#, strf)        
    plt.close()
    quit()

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
    for unit in range(0,len(Sort.units),1):
        print "Unit : ", unit, " No. spikes: ", len(Sort.units[unit])
        
        frame_indexes = []
        print "Finding nearest frame indexes"
        spike_times = np.array(Sort.units[unit])*40 #Convert spike times to usec to match din times;
        correct_indexes = np.where(np.array(spike_times)<1E14)  #Remove weird 1E+19 values from the spiketimes; SS bug!
        spike_times = spike_times[correct_indexes]


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
        
        print "Original input_array shape: ", input_array.shape
        
        lags = [-2,-1,0,1,2]
        alpha = [[-1E3], [-1E4], [-650], [-50], [-10.0]]
        
        #Serial ridge regression over alphas
        #strf, bias = fit_strf_ridge(input_array, output_array, lags=lags, alpha=alpha[0], verbose=True)

        #Parallel ridge regression over alphas
        strf_array = []
        input_array = [input_array]*len(alpha)
        output_array = [output_array]*len(alpha)
        lags = [lags]*len(alpha)
        pool = mp.Pool(n_procs)
        strf_array.extend(pool.map(Parallel_fit_strf_ridge, zip(input_array, output_array, lags, alpha)))
        pool.close()

        print "Lags: ", lags
        print "Alpha: ", alpha

        #Loop over computed alphas
        images_alphas = []
        for a in range(len(alpha)):
            strf = np.array(strf_array)[a].T
            
            #Get strf product with original data for each lag
            strf_temp = []
            for i in range(len(lags[0])):
                strf_temp.append(np.multiply(input_array[a],strf[i]))
            
            strf_temp = np.array(strf_temp)
            
            #Average image for each lag ?!
            image1d_average = []
            for i in range(len(lags[0])):
                image1d_average.append(np.average(strf_temp[i], axis=0))

            #Convert 1d image back to 2d for display
            image2d = []
            for i in range(len(lags[0])):
                image2d.append(np.reshape(image1d_average[i], movie.frames[0].shape))
            image2d = np.array(image2d)
            
            #Add white separators between lags
            images_plot=image2d[0]
            blank_image = np.zeros((len(image2d[0]),3),dtype=np.uint8)+np.max(image2d)
            for i in range(1,len(lags[0]),1):
                images_plot=np.hstack((images_plot,blank_image))
                images_plot=np.hstack((images_plot,image2d[i]))
        
            images_alphas.append(images_plot)
        
        images_plot=images_alphas[0]
        for i in range(1,len(images_alphas),1):
            images_plot=np.vstack((images_plot,images_alphas[i]))
        
        plt.imshow(images_plot, cmap=cm.gray, interpolation='none')
        plt.title(file_name + " unit : " + str(unit) + " no. spikes: " + str(len(spike_times)) + " lags: " + str(lags[0])+ " alpha: " + str(alpha))

        plt.show()
        
        

        



















































