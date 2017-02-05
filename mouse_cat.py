#Mouse class and subclasses for processing lever pull data projec Greg Silasi and Tim Murphy lab
#Code: Cat Mitelut

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
import cPickle as pickle
import dill
import gc
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import shutil
from imreg import *

class Mouse(object):      
    
    def __init__(self, mouse_name, home_dir, n_sec, window, n_pixels):
        
        self.name = mouse_name
        self.home_dir = home_dir
        self.window = window
        self.n_sec = n_sec
        self.n_pixels = n_pixels
        
    def process(self):

        self.load_filenames()   #Loads tif files, event_files, lever_files
        
        self.load_sessions()    #Loads session data

        self.stroke = Stroke(self.sessions, self.home_dir, self.name)      #Load stroke information; need to have processed session info first to place location of stroke

        self.save_mouse()       #Save mouse file names and trace data to pickle object (NB: no DFF data saved)

        

    def load_filenames(self):
        '''Load event files, tiffs and lever position files'''
        
        #LOAD EVENT FILE NAMES
        event_files = os.listdir(self.home_dir + self.name + '/event_files')

        #ORDER EVENT FILES BY DATE
        file_order = []
        for event_file in event_files:
            if '2015' in event_file: 
                temp = event_file.replace(self.name+'_','')
                temp = temp[:temp.find('2015')-1]
                month = temp[:temp.find('-')]
                day = temp[temp.find('-')+1:]
                #print month,day
                file_order.append('2015'+month.zfill(2)+day.zfill(2))
            else: 
                temp = event_file.replace(self.name+'_','')
                temp = temp[:temp.find('2016')-1]
                month = temp[:temp.find('-')]
                day = temp[temp.find('-')+1:]
                #print month,day
                file_order.append('2016'+month.zfill(2)+day.zfill(2))

        indexes= np.argsort(np.int32(file_order)) #Sorted indexes for input files.
        event_files = list(np.array(event_files)[indexes])
            
        #LOAD LEVER FILES NAMES
        lever_files = []
        del_file_counter = []
        counter=0
        for event_file in event_files:
            if (os.path.exists(self.home_dir+self.name+'/leverPos/'+event_file[:-6]+'leverPos'+event_file[-6:])==True):
                lever_files.append(self.home_dir+self.name+'/leverPos/'+event_file[:-6]+'leverPos'+event_file[-6:])
            else:
                print "Missing lever pos file: ", self.home_dir+self.name+'/leverPos/'+event_file[:-6]+'leverPos'+event_file[-6:]
                del_file_counter.append(counter)
            counter+=1
            
        #del_file_counter = np.array(del_file_counter)
        for i in reversed(del_file_counter):
            del event_files[i]
            
        #LOAD TIF FILE NAMES; need to convert to letter dates;
        tif_files = []
        months = {'1':'Jan', '2':'Feb', '3':'Mar','4':'Apr','5':'May','6':'Jun',
                  '7':'Jul','8':'Aug','9':'Sep','10':'Oct','11':'Nov','12':'Dec'}
        counter=0
        while True:
            event_file=event_files[counter].replace(self.name+'_','')
            month = event_file[:2].replace('-','')
            day = event_file[len(month):len(month)+3].replace('-','')
            #print event_file, "  month: ", month,  "day: ", day
            
            if 'AM' in event_file: am_pm = 'am'
            else: am_pm = 'pm'
            
            dir_name = glob.glob(self.home_dir+self.name+'/tif_files/'+self.name+am_pm+'_'+months[month]+day+"_*")   #use .tif extension otherwise will load .npy files
            
            if len(dir_name)==1: #Found tif file
                make_dir = dir_name[0]
                fname = glob.glob(make_dir+"*")[0]
                new_name = make_dir +'/'+ fname.replace(self.home_dir+self.name+'/tif_files/','')+'.tif'

                        
                #fnames = glob.glob(make_dir+"*")   #use .tif extension otherwise will load .npy files
                #if not os.path.exists(make_dir):        #Make directory
                    #os.makedirs(make_dir)

                    #for fname in fnames:
                        #new_name = make_dir +'/'+ fname.replace(self.home_dir+self.name+'/tif_files/','')
                        #print new_name
                        #shutil.move(fname, new_name)    

                tif_files.append(new_name)
                counter+=1
                
            elif len(dir_name)>1:
                print dir_name
                print "---------TOO MANY TIF FILES-------"
                quit()
                
            elif len(dir_name)==0:
                print "Missing tif file: ", self.home_dir+self.name+'/tif_files/'+self.name+am_pm+'_'+months[month]+day+"*"
                del lever_files[counter]
                del event_files[counter]
                
            if counter==len(event_files): break

        for k in range(len(event_files)):
            event_files[k] = self.home_dir+self.name+'/event_files/'+event_files[k]
    
        self.event_files = event_files
        self.lever_files = lever_files
        self.tif_files = tif_files

    def load_sessions(self):

        self.sessions = []

        #if os.path.exists(home_dir+mouse+"/"+mouse+"_imaging.npy")==False:
        counter=-1
        for tif_file, event_file, lever_file in zip(self.tif_files, self.event_files, self.lever_files):
            print ""
            print "******************************************************"
            print "Session: ", len(self.sessions), " / ", len(self.tif_files)
            print ".......: ", tif_file
            print ".......: ", event_file
            print ".......: ", lever_file
    
            #counter+=1; print counter
            #if counter!=6: continue
            
            self.sessions.append(Session(tif_file, event_file, lever_file, self.window, len(self.sessions),self.home_dir, self.name))
            #if len(self.sessions)>15: return
            
            
    def save_mouse(self):
        print "Saving mouse to disk..."
        
        f = open(self.home_dir+self.name+'/'+self.name+'.pkl', 'wb')
        pickle.dump(self, f)


    def load_mouse(self, mouse_name, home_dir):

        print "Loading mouse from disk..."
        
        f = open(home_dir+mouse_name+'/'+mouse_name+'.pkl', 'rb')
        gc.disable()                    # disable garbage collector
        data = pickle.load(f)
        gc.enable()             # enable garbage collector again

        #data.stroke = Stroke(data.home_dir, data.name, data.sessions)      #Need to reload this as not saved in all mice; CAN REMOVE LATER
        
        return data


    def move_tifs(self):
        import shutil
        
        for session in self.sessions:
            
            fnames = glob.glob(session.tif_file[:-4]+"*")   #use .tif extension otherwise will load .npy files

            dir_name = session.tif_file.replace('.tif','/')
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)

                for fname in fnames:
                    new_name = dir_name + fname.replace(self.home_dir+self.name+'/tif_files/','')
                    print fname
                    print new_name
                    shutil.move(fname, new_name)    
            
        #quit()


    def save_DFF(self):
        ''' NOT USED CURRENTLY '''
        
        print "Loading all DFF - saving to single file (unit8)"
        
        DFF_array = []
        if (os.path.exists(self.home_dir+self.name+'/'+self.name+'_DFF.npy')==False):
            for session in self.sessions:
                print session.tif_file
                DFF_array.extend(np.load(session.tif_file[:-4]+'_3sec.npy'))
        
            #Convert list to np array
            DFF_array = np.array(DFF_array)
            print DFF_array.shape

            print "Writing DFF_array to disk..."
            np.save(self.home_dir+self.name+'/'+self.name+'_DFF.npy', DFF_array)

            #CODE FOR CONVERTING TO 64pixel by averaging... not great, but may be ok for debugging.
            #DFF_array_64 = block_reduce(DFF_array, block_size=(1,1,2,2), func=np.mean)
            #np.save(self.home_dir+self.name+'/'+self.name+'_DFF_64.npy', DFF_array_64)
     
        else:
            print "DFF already saved..."
            
            DFF_array = np.load(self.home_dir+self.name+'/'+self.name+'_DFF.npy')

        print DFF_array.shape
        quit()
        
        

class Session(object):
    
    def __init__(self, tif_file, event_file, lever_file, window, index, home_dir, name):
    
        self.home_dir = home_dir
        self.name = name
        self.tif_file = tif_file
        self.event_file = event_file
        self.lever_file = lever_file
        self.window = window
        self.index = index                  #Save session index for access inside session object
        self.aligned_images = []            #Make empty list to ensure
        
        self.load_lever()                   #Load lever positions

        if self.reclength == 0: return
    
        self.convert_tif()                  #Convert .tif -> .npy format; save to disk

        self.align_images()                 #Align raw_images to first session frame 1000

        self.compute_DFF()                  #Compute DFF on aligned images

        self.make_trials()                  #Make and populate individual reward trials for each session


    def load_lever(self):

        plotting = False

        print " ... loading lever positions..."
        #Make lever object
        self.trials = []       #List for each trial 
        self.positions = []    #All lever positions 
        self.threshold = []    #Holds '44' values which indicates threshold crossing; '77' otherwise
        self.times = []        #Corresponding times for each lever position
        self.code = []         #Code for trial
        
        #****************LOAD EVENT FILE****************
        text_file = open(self.event_file, "r")
        lines = text_file.read().splitlines()
        event_text = []
        for line in lines:
            event_text.append(re.split(r'\t+',line))
        
        #Delete false starts from event file
        for k in range(len(event_text)-1,-1,-1):        #Search backwards for the 1st occurence of "date" indicating last imaging start
                                                        #NB: There can be multiple false starts WHICH DON"T LINE UP - NEED TO IGNORE SUCH SESSIONS
            if event_text[k][0]=='date': 
                event_text = event_text[k+2:]         #Remove first 2 lines
                break
        
        if len(event_text)==0:
            print "empty session"
            self.reclength = 0
            return 
        else:
            self.reclength = float(event_text[-1][3])
        #print self.reclength
        #quit()

        #**************** LOAD LEVER POSTIION FILE **********
        text_file = open(self.lever_file, "r")
        lines = text_file.read().splitlines()
        lever_text = []
        for line in lines:
            lever_text.append(re.split(r'\t+',line))
        lever_text = np.array(lever_text)

        #Delete false starts from lever position file
        for k in range(0,len(lever_text),1):     #Search for 1st occurence of first correct value in event file
            if lever_text[k][0]==event_text[0][0]: 
                lever_text = lever_text[k:]        
                break
        
        #Convert lever_text array into data for the lever object
        trial_counter=-1
        lever_counter=0
        event_counter=-1
        for event in event_text[:-1]:       #SHITY CODE... NEED TO REDO THESE LOOPS
            event_counter+=1
            if lever_counter==len(lever_text): break
            if (event[0]==lever_text[lever_counter][0]) and (event[1]==lever_text[lever_counter][2]):
                trial_counter+=1
                #Save code from file
                self.code.append(lever_text[lever_counter][2])
                
                #Save relative time
                self.trials.append(float(event_text[event_counter][3]))

                #Save threshold value, positions and times
                self.threshold.append([])
                self.positions.append([])
                self.times.append([])

               # print "*", event[0], event[1]
               # print "*", lever_text[lever_counter][0], lever_text[lever_counter][2]
               # print "*", self.trials[trial_counter], '\n\n'


                while True: 
                    lever_counter+=1
                    if lever_counter==len(lever_text): break
                    #if (lever_text[k][0][:5]=='2016-') or (lever_text[k][0][:5]=='2015-'):   #CAREFUL THIS MAY BE FOOLED BY OTHER  DATES
                    if (lever_text[lever_counter][1]=='None'):   break 
                    
                    #Save lever positions
                    self.times[trial_counter].append(float(lever_text[lever_counter][0]))       #Time values from leverpos file; e.g. 2.91, 2.99, 3.07...5.02
                    self.positions[trial_counter].append(float(lever_text[lever_counter][1]))   #Position vals from levpos file; e.g. 0, 0, 1, 11, 34, 44, ...
                    self.threshold[trial_counter].append(float(lever_text[lever_counter][2]))   #Threshold vals from levps file; e.g. 77, 77, 77, 44, 77, 77...
          
        #quit()
        
        #Generate continuous time stream array
        self.abstimes = []

        #Absolute positions
        self.abspositions = []
        
        #Trial code for each recorded position
        self.abscodes = []
        
        #Successful codes: '44' value;
        self.absthresholds = []

        colors = { '04': 'blue', '02': 'red', '07': 'black'}
        #Load absolute positions 
        
        #plotting=True
        '''  NB: This code saves data in chunks and not in real time; the indexes of lever position data (abspositions) do match real time data (abstimes);
                These time values (abstimes) can be reliably used to match imaging data;
        '''
        for k in range(len(self.threshold)):
            if 44 in self.threshold[k]:  #Search trial for 44 value to match to lever.trials times
                if self.code[k]=='04':
                    #print k
                    thresh_index = self.threshold[k].index(44)
                    #print thresh_index
                    ##print self.times[k]
                    #print self.times[k][thresh_index]
                    #print self.trials[k]
                    
                    #Save absolute times to nearest milisecond
                    times = np.int32((np.array(self.times[k]) - self.times[k][thresh_index] +self.trials[k])*1000)
                    #lever.abstimes.extend([times[0]/1000.-0.001]) #Add extra vals before and after 
                    self.abstimes.extend(times/1000.)
                    #lever.abstimes.extend([times[-1]/1000.+0.001])
                    
                    #Save absolute lever positions
                    #lever.abspositions.extend([0])
                    self.abspositions.extend(self.positions[k])
                    #lever.abspositions.extend([0])

                    #Save trial code for each position
                    self.abscodes.extend([self.code[k]]*(len(self.positions[k])))

                    #Save threshold values as stream
                    self.absthresholds.extend(self.threshold[k])

                    if plotting:
                        temp_time= []
                        temp_time.extend(times/1000.)
                        temp_lever = []
                        temp_lever.extend(self.positions[k])
                        plt.plot(temp_time, temp_lever, linewidth=1.5, color=colors[self.code[k]])

        if plotting:
            #plt.plot(lever.abstimes, lever.abspositions, linewidth=1.5)
            plt.ylabel("Lever Position", fontsize=25)
            plt.xlabel("Time (sec)", fontsize=25)
            plt.plot([0,self.abstimes[-1]],[11,11], 'r--', linewidth = 2, color='black', alpha =0.5)
            plt.plot([0,self.abstimes[-1]],[59,59], 'r--', linewidth = 2, color='black', alpha =0.5)
            plt.show()
        
        print "Event Codes:  # 02: ", len(np.where(np.array(self.code)=='02')[0]), ",   # 04: ", len(np.where(np.array(self.code)=='04')[0]), \
        ",   # 07: ", len(np.where(np.array(self.code)=='07')[0]), ",   # 08: ", len(np.where(np.array(self.code)=='08')[0]), \
        ",   # 00: ", len(np.where(np.array(self.code)=='00')[0])
        
        
        #SET THRESHOLDS AND REWARD CODES
        self.trigger = 44
        self.trigger_indexes = np.where(np.array(self.absthresholds)==self.trigger)[0]
        #print self.trigger_indexes 

        #SELECT CODE 
        self.reward_code = '04'  #Choose codes to process data;
        indexes_04 = np.where(np.array(self.abscodes)[self.trigger_indexes]==self.reward_code)[0]
        self.trigger_indexes = self.trigger_indexes[indexes_04]
        
        #print self.trigger_indexes
        #print np.array(self.abstimes)[self.trigger_indexes]
        #print "# trials for code: ", self.reward_code, " ('44' index): ", len(self.trigger_indexes)

        #quit()

    
    def convert_tif(self):
    
        if (os.path.exists(self.tif_file[:-4] +'.npy')==False):
            print "...read: ", self.tif_file
            images_raw = tiff.imread(self.tif_file)

            print "... saving .npy"
            np.save(self.tif_file[:-4], images_raw)


    def compute_DFF(self):

        if os.path.exists(self.tif_file[:-4]+"_"+str(self.window/30)+"sec_traces.npy"): 
            self.traces=np.load(self.tif_file[:-4]+"_"+str(self.window/30)+"sec_traces.npy")
            return

        #Load rotated imaging data
        if len(self.aligned_images)==0: 
            print "...loading aligned img data..."
            self.aligned_images = np.load(self.tif_file[:-4]+"_aligned.npy")
        
        self.n_images=len(self.aligned_images)
        
        #COMPUTE absolute times of triggers from lever pull data
        trigger_times = np.array(self.abstimes)[self.trigger_indexes]

        #FIND NEAREST PREVIOUS IMAGING FRAME TIME FOR EACH LEVER POS THRESHOLD
        frame_times = np.linspace(0, self.reclength, self.n_images)
        img_frame_triggers = []
        for i in range(len(trigger_times)):
            img_frame_triggers.append(self.find_previous(frame_times, trigger_times[i]))
        print "Number of triggers: ", len(img_frame_triggers)

        #CHECK img_rate 
        img_rate = self.n_images/self.reclength
        print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", img_rate
        
        print "...computing DF/F..."
        
        #data_array = np.zeros((window*2,128,128), dtype=np.float32)
        data = []
        traces = []
        counter=-1
        plotting=True
        for trigger in img_frame_triggers:
            counter+=1
            print "...trial trigger frame: ", trigger
            #NB: Skip first 100 frames or so as the camera is not on and DF/F averaging will not work
            if trigger < (2*self.window+300) or trigger>(self.n_images-self.window): 
                print "..skip: too close to start/end ..."
                continue  #Skip if too close to start/end

            #***PROCESS IMAGING; load data before and after pull
            data_chunk = self.aligned_images[trigger-self.window:trigger+self.window]
            
            #DF/F computation; uses baseline -2*window .. -window; e.g. for 3sec windows, baseline: average(-6sec..-3sec)
            baseline = np.average(self.aligned_images[trigger-2*self.window:trigger-self.window], axis=0)
            data_chunk = (data_chunk-baseline)/baseline 
            data.append(data_chunk)
            
            #***PROCESS TRACES - WORKING IN DIFFERENT TIME SCALE
            #lever_window = 50
            lever_window = 120*3    #NB: Lever window is computing in real time steps @120Hz; and the trigger_indexes are discontinuous;
            t = np.linspace(-lever_window*0.0082,lever_window*0.0082, lever_window*2)
            lever_trace = self.abspositions[self.trigger_indexes[counter]-lever_window:self.trigger_indexes[counter]+lever_window]
            
            
            if len(lever_trace)!=len(t): 
                print "...missing lever trace data ... extrapolating..."
                lever_trace = np.zeros(lever_window*2,dtype=np.float32)
                for k in range(-lever_window,lever_window,1):
                    lever_trace[k+lever_window] = self.abspositions[k+lever_window]     #Double check this...
                if plotting: 
                    plt.plot(t, lever_trace, color='red')
                    plt.show()
            else:
                if plotting: 
                    plt.plot(t, lever_trace, color='black')
                    plt.show()

            traces.append(lever_trace)

        np.save(self.tif_file[:-4]+'_'+str(self.window/30)+"sec_traces", traces)
        self.traces = traces

        #Save individual trial time dynamics
        data = np.float16(data)
        for k in range(len(data)):
            np.save(self.tif_file[:-4]+'_'+str(k).zfill(4), data[k])
            
        #Set aligned_images to empty 
        self.aligned_images = []

    def align_images(self):

        if os.path.exists(self.tif_file[:-4]+"_aligned.npy")==False:
            
            raw_images=np.load(self.tif_file[:-4]+'.npy')
            np.save(self.tif_file[:-4]+'_frame1000', raw_images[1000])
            
            #Save 1st session frame 1000 to be used by all other code
            print self.index
            if self.index==0: np.save(self.home_dir+self.name+'/'+self.name+'_align_frame', raw_images[1000])
            
            #Load 1st session frame1000:
            first_img = np.load(self.home_dir+self.name+'/'+self.name+"_align_frame.npy")

            #Load current session frame1000
            current_img = raw_images[1000]
            
            #Find translation between frame1000 current session and 1st session
            #im2, scale, angle, t = similarity(im0, plot_img)
            t = translation(first_img, current_img)
            
            #Load raw .npy imaging data
            aligned =np.load(self.tif_file[:-4]+".npy")
            for k in range(len(aligned)):
                aligned[k] = np.roll(aligned[k], t[0], axis=0)
                aligned[k] = np.roll(aligned[k], t[1], axis=1)
            
            #Save rotated imaging data
            np.save(self.tif_file[:-4]+"_aligned", aligned)
            self.aligned_images = aligned

    def make_trials(self):
        ''' Saves rewarded events from each session into individual trial objects for later analysis'''
        
        self.trials = []
        for k in range(len(self.traces)):
            trial = Trial(self.traces[k])
            self.trials.append(trial)

    def find_previous(self, array, value):
        temp = (np.abs(array-value)).argmin()
        if array[temp]>value: return temp-1
        else: return temp
        
    #def load_DFF(self):
        #''' NOT USED ''' 
        #print " session: ", self.event_file

        #self.DFF = np.load(self.DFF_fname+'.npy')                
        
        #for k in range(len(self.DFF)):
            #self.trials[k].DFF = self.DFF[k]
      


class Trial(object):
    ''' Make individual trial swithin each session representing individual DF/F stacks and lever position traces'''
    def __init__(self, traces):

        #self.DFF = DFF
        self.traces = traces


class Stroke(object):
    ''' Make individual trial swithin each session representing individual DF/F stacks and lever position traces'''
    
    def __init__(self, sessions, home_dir, name):

        stroke_data = np.loadtxt(home_dir+name+'/stroke.txt', dtype=str)
        
        if len(stroke_data)==3: 
            self.ampm = stroke_data[1]
            self.kind = stroke_data[2]
            self.day = stroke_data[0][3:5]
            self.month = stroke_data[0][0:2]
            self.year = '20'+stroke_data[0][6:8]

        else:
            print "-----Stroke file problem"
            quit()
        

        #Find relative location of stroke in rewarded trials for later referencing/indexing
        trial_counter = 0
        for session in sessions:
            date = session.event_file.replace(home_dir,'').replace(name+'/event_files/','').replace('.txt','').replace(name+'_','')
            month = date[0:date.index('-')];    date = date[date.index('-')+1:]
            day = date[0:date.index('-')];    date = date[date.index('-')+1:]
            year = date[0:date.index('_')]; 
            
            if year == self.year:
                if (int(month) >= int(self.month)) and (int(day)>=int(self.day)):
                    #print self.month, self.day, self.year
                    #print session.event_file
                    self.trial = trial_counter
                    #print "Stroke trial: ", self.trial
                    break

            trial_counter+=len(session.trials)
    
       




