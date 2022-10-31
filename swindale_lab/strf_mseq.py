import numpy as np
import multiprocessing as mp

from tsf_ptcs_classes import *
#from strf_utils import *
#from movie import *

import operator

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

#********** Load simple mseq movie file ****************

class mseq_movie():
    
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
            

def Parallel_average(data):
    output = np.average(data,axis=0)
    return output
    

def Parallel_find_previous_frame_index((targets)):
    global din
    indexes = []
    for i in range(len(targets)):
        #index = np.argmin(np.abs(din[:,0] - targets[i]))
        print targets[i]
        quit()
        index = np.argmin(np.abs(din[:,0] - targets[i]))
        
        if (targets[i] < din[:,0][index]):
            indexes.append(din[:,1][index-1])
        else:    
            indexes.append(din[:,1][index])
    return indexes

def Find_previous_frame_index((targets)):
    global din
    indexes = []
    for i in range(len(targets)):
        #index = np.argmin(np.abs(din[:,0] - targets[i]))
        
        index = np.argmin(np.abs(din[:,0] - targets[i]))
        #print index
        #print targets[i]
        #if (targets[i] < din[:,0][index]):
        #    print din[:,0][index-1]
        indexes.append(din[:,1][index-1])
        #else:    
        #    print din[:,0][index]
        #    indexes.append(din[:,1][index])

        #print ""

    #quit()

    return np.array(indexes)


n_procs = 6

#***********************************************************************************************
#*************************************** PTC15 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/12TB/in_vivo/nick/ptc15/'

#movie_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/mseq/MSEQ32'



#***********************************************************************************************
#*************************************** PTC17 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'

#file_name = '03-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '06-tr1-MVI_1406-1410'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1406-1410'

#file_name = '12-tr1-MVI_1413-1416_128x128x1024'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/f-1/MVI_1413-1416_128x128x1024'

#file_name = '25-tr1-MVI_1398-1399_200Hz'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1398-1399'

#file_name = '31-tr1-mseq32_40ms'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'

#file_name = '48-tr2b-mseq32_40ms'; mseq_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/mseq/MSEQ32'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'

#file_name = '60-tr2b-mseq16_40ms'; mseq_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/mseq/MSEQ16'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc17/'


#***********************************************************************************************
#*************************************** PTC18 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/12TB/in_vivo/nick/ptc18/'

#file_name = '33-tr1-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'

#***********************************************************************************************
#*************************************** PTC21 *************************************************
#***********************************************************************************************

#main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

#file_name = '60-tr5c-MVI_1419_5s'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-25/MVI_1419'

#file_name = '58-tr5c-MVI_1403-1405'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1403-1405'
#Examples of good reverse correlation right out of nat scenes
#Unit :  7  No. spikes:  3651
#Unit :  10  No. spikes:  1195

#file_name = '25-tr2-MVI_1406-1410'
#movie_file = '/media/cat/12TB/in_vivo/nick/stimuli/Polytrode_Data/movies/2007-11-24/MVI_1406-1410'


#********* MSEQ DATA ********
#file_name = '57-tr5c-mseq32_40ms'; mseq_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/mseq/MSEQ32'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'

file_name = '75-tr5c-mseq32_40ms'; mseq_file = '/media/cat/4TB/in_vivo/nick/Stimulus/Polytrode_Data/movies/mseq/MSEQ32'; main_dir = '/media/cat/4TB/in_vivo/nick/ptc21/'
#Unit 12 has 29k spikes and good single phase strf

#file_name = '93-t6b-mseq32_40ms'
#movie_file = '/media/cat/4TB/in_vivo/nick/stimuli/Polytrode_Data/movies/mseq/MSEQ32'



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

#Load non_locking_spikes from file:
if True:
    non_lock_spikes_array = np.load(main_dir+file_name+'/'+file_name+'_non_locked_spikes.npy')
    lock_spikes_array = np.load(main_dir+file_name+'/'+file_name+'_locked_spikes.npy')
    lfp_events_array = np.load(main_dir+file_name+'/'+file_name+'.lfp.events.npy')


#**** Load .din stimulus timestamps
din_file = main_dir + file_name + '/' + file_name+'.din'
f = open(din_file, 'rb')
din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols containing [time_stamp in usec : frame number]
f.close()

#Remove blank (black) screen frames from stimulus
din = din[din[:,1]!=65535]

print "No of frames in din file: ", len(din)
rec_length = float(din[-1][0])*1E-6     #rec length in seconds
print "Length of recording: ", rec_length/60., " mins." 

#**** Load movie frames ****
#movie = Movie()
movie = mseq_movie()
movie.load(mseq_file)
print "Length of movie (in frames): ", len(movie.frames)


#**************************************************************
#*********************** COMPUTE STRF *************************

#Single unit RFs
if True:
    unit_list = np.arange(0,len(Sort.units),1)
    for unit in unit_list:
        #Find # of spikes following each frame

        all_spikes = np.array(Sort.units[unit])*40.         #From sample time to usec to match din times.
        non_lock_spikes = non_lock_spikes_array[unit]*1E6   #From sample time to usec
        lock_spikes = lock_spikes_array[unit]*1E6

        #if len(all_spikes)!=1461: continue
        print unit
        print "unique unit: ", Sort.uid[unit]
        
        all_spikes=all_spikes[all_spikes<5E10]
        non_lock_spikes=non_lock_spikes[non_lock_spikes<5E10]
        lock_spikes=lock_spikes[lock_spikes<5E10]

        #Look only at higher firing rate cells;
        fire_rate = len(all_spikes)/rec_length 
        if fire_rate<0.250: continue

        print "Unit: ", unit, " f_rate: ", fire_rate, " # all_spikes: ", len(all_spikes), "  # non_lock: ", len(non_lock_spikes), "  # lock: ", len(lock_spikes)

        #Loop over spikes:
        n_frames = 28
        all_spikes_indexes = []
        for i in range(len(all_spikes)):
            index = np.argmin(np.abs(din[:,0] - all_spikes[i]))     #Find index of nearest frame for each spike time
            if din[:,0][index]>all_spikes[i]: index-=1              #Ensure index is to previuos frame; so if time > spike time, use previous frames
            if index>=n_frames:                                     #Ensure not going into negative frames for very early spikes;
                all_spikes_indexes.append(din[:,1][index-n_frames:index])
        all_spikes_indexes = np.array(all_spikes_indexes).T #Transpose to convert each row to a time index (e.g. row 0 contains all indexes at t=0)

        non_lock_spikes_indexes = []
        for i in range(len(non_lock_spikes)):
            index = np.argmin(np.abs(din[:,0] - non_lock_spikes[i]))
            if din[:,0][index]>non_lock_spikes[i]: index-=1  
            if index>=n_frames:
                non_lock_spikes_indexes.append(din[:,1][index-n_frames:index])
        non_lock_spikes_indexes = np.array(non_lock_spikes_indexes).T #Transpose to convert each row to a time index (e.g. row 1 contains all indexes at t=0)

        lock_spikes_indexes = []
        for i in range(len(lock_spikes)):
            index = np.argmin(np.abs(din[:,0] - lock_spikes[i]))
            if din[:,0][index]>lock_spikes[i]: index-=1  
            if index>=n_frames:
                lock_spikes_indexes.append(din[:,1][index-n_frames:index])
        lock_spikes_indexes = np.array(lock_spikes_indexes).T #Transpose to convert each row to a time index (e.g. row 1 contains all indexes at t=0)


        time_averages = np.zeros((int(n_frames/4.),32,32), dtype=np.float32)
        counter=1
        for k in range(0,n_frames,4):
            ax=plt.subplot(3,int(n_frames/4.),counter)
            temp_img = []
            for f in range(4):
                temp_img.append(np.average(np.array(movie.frames)[all_spikes_indexes[k+f]], axis=0))
            temp_img = np.average(temp_img, axis=0)
            v_abs=np.max(np.abs(temp_img))-127
            plt.imshow(temp_img, vmin=127-v_abs, vmax=127+v_abs, interpolation='none')
            if k ==0:
                plt.ylabel("all spikes\n#: "+str(len(all_spikes)))
                ax.xaxis.set_visible(False)
            ax.yaxis.set_ticklabels([])    
            ax.xaxis.set_ticklabels([])    
            
            if len(non_lock_spikes)>0:
                ax=plt.subplot(3,int(n_frames/4.),counter+int(n_frames/4.))
                temp_img = []
                for f in range(4):
                    temp_img.append(np.average(np.array(movie.frames)[non_lock_spikes_indexes[k+f]], axis=0))
                temp_img = np.average(temp_img, axis=0)
                v_abs=np.max(np.abs(temp_img))-127
                plt.imshow(temp_img, vmin=127-v_abs, vmax=127+v_abs, interpolation='none')
                if k ==0:
                    plt.ylabel("non_lock spikes\n#: "+str(len(non_lock_spikes)))
                ax.yaxis.set_ticklabels([])    
                ax.xaxis.set_ticklabels([])    

            if len(lock_spikes)>0:
                ax=plt.subplot(3,int(n_frames/4.),counter+2*int(n_frames/4.))
                temp_img = []
                for f in range(4):
                    temp_img.append(np.average(np.array(movie.frames)[lock_spikes_indexes[k+f]], axis=0))
                temp_img = np.average(temp_img, axis=0)
                v_abs=np.max(np.abs(temp_img))-127
                plt.imshow(temp_img, vmin=127-v_abs, vmax=127+v_abs, interpolation='none')
                if k ==0:
                    plt.ylabel("lock_spikes\n#: "+str(len(lock_spikes)))
                    ax.xaxis.set_visible(False)
                ax.yaxis.set_ticklabels([])    
                ax.xaxis.set_ticklabels([])    

            counter+=1

        plt.suptitle(file_name+ "  Unit: " + str(unit), fontsize=30)
        plt.show()

#********LFP Event based RFs
if True:
    unit_list = np.arange(0,len(lfp_events_array),1)
    for unit in unit_list:
        #Find # of spikes following each frame

        all_spikes = np.array(lfp_events_array[unit])*1E3         #From sample time to usec to match din times.
        all_spikes=all_spikes[all_spikes<5E10]

        print "LFP Unit: ", unit, " # all_spikes: ", len(all_spikes)

        #Loop over spikes:
        n_frames = 28
        all_spikes_indexes = []
        for i in range(len(all_spikes)):
            index = np.argmin(np.abs(din[:,0] - all_spikes[i]))     #Find index of nearest frame for each spike time
            if din[:,0][index]>all_spikes[i]: index-=1              #Ensure index is to previuos frame; so if time > spike time, use previous frames
            if index>=n_frames:                                     #Ensure not going into negative frames for very early spikes;
                all_spikes_indexes.append(din[:,1][index-n_frames:index])
        all_spikes_indexes = np.array(all_spikes_indexes).T #Transpose to convert each row to a time index (e.g. row 0 contains all indexes at t=0)

        time_averages = np.zeros((int(n_frames/4.),32,32), dtype=np.float32)
        counter=1
        for k in range(0,n_frames,4):
            ax=plt.subplot(3,int(n_frames/4.),counter)
            temp_img = []
            for f in range(4):
                temp_img.append(np.average(np.array(movie.frames)[all_spikes_indexes[k+f]], axis=0))
            temp_img = np.average(temp_img, axis=0)
            v_abs=np.max(np.abs(temp_img))-127
            plt.imshow(temp_img, vmin=127-v_abs, vmax=127+v_abs, interpolation='none')
            if k ==0:
                plt.ylabel("all spikes\n#: "+str(len(all_spikes)))
                ax.xaxis.set_visible(False)
            ax.yaxis.set_ticklabels([])    
            ax.xaxis.set_ticklabels([])    
            
            counter+=1

        plt.suptitle(file_name+ "  Unit: " + str(unit), fontsize=30)
        plt.show()
    




quit()
    





#Loop over frames:
if False:
    img_bins = []
    for k in range(4):
        img_bins.append([])
        
    for i in range(65535):
        print i
        frame_times = din[din[:,1]==i][:,0]
        
        for k in range(0,80000,20000):
            img_bins[k/20000].extend([i]*len(np.where(np.logical_and(spike_times>=frame_times[0]+k, spike_times<=frame_times[0]+k+20000))[0]))

    #Make average frame for each time segment
    images_total=np.zeros((4,32,32), dtype=np.float32)
    for k in range(4):
        temp_array = np.array(movie.frames)[img_bins[k]]
        images_total[k] = np.average(temp_array,axis=0)














