#Created on Mon Jan 20 14:48:45 2014
#Function to get all the analog entities from the *.mcd files
#created by the MCRack software from multichannelsystems.
#It takes the complete file path as input and returns a
#dictionary containing the name of each channel as a key and
#its contents as values (numpy arrays)
#for this to work, the following libraries for python must be installed:
#neuroshare bindings from http://pythonhosted.org/neuroshare/
#numpy
#This code is distributed under:
#creative commons attribution-shareAlike 4.0 international (CC BY-SA 4.0) license.
#@author: andre maia chagas -
#find this and more open source tools @
# www.openeuroscience.wordpress.com
#function to get the recording of the digital line from
 
import numpy as np
import matplotlib.pyplot as plt
import struct
from tsf_ptcs_classes import *

def MCD_read(MCDFilePath):
 
    #import necessary libraries
 
    import neuroshare as ns
    import numpy as np
 
    #open file using the neuroshare bindings
 
    fd = ns.File(MCDFilePath)
 
    #create index
 
    indx = 0
 
    #create empty dictionary
 
    data = dict()
 
    #loop through data and find all analog entities
 
    for entity in fd.list_entities():
        print "looping over entities: ", entity
        analog = 2
 
        #if entity is analog
 
        if entity.entity_type == analog:
 
            #store entity in temporary variable
 
            dummie = fd.entities[indx]
 
            #get data from temporary variable
 
            data1,time,count=dummie.get_data()
 
            #create channel names
 
            channelName = entity.label[0:4]+entity.label[23:]
 
            #store data with name in the dictionary
 
            data[channelName] = np.array(data1)
 
        #increase index
 
        indx = indx + 1
 
    #return dictionary
 
    return data


data_dir = '/media/cat/4TB/in_vivo/tim/dongsheng/2015-11-27/'
file_name = ['2015-11-27-4-10electrodein-iso0.mcd']


#data_dir = '/media/cat/4TB/in_vivo/tim/dongsheng/2015-12-1/'
#file_name = ['2015-12-1-1-10electrodeiniso1.mcd']



for i in range(len(file_name)):
    file_name[i]=data_dir + file_name[i][:-4]+'/'+file_name[i]

for file_ in file_name:
    print "LOADING FILE: ", file_
    data = MCD_read(file_)

    data_keys = data.keys()
    print data_keys
    
    order = np.array([15, 10, 16, 11, 14, 12, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13])
    #order = np.array([15, 10, 16, 11, 14, 666, 12, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13]) #Some recs have extra dummy value

    ecp = np.zeros((16,len(data[data_keys[0]])), dtype=np.int16)
    counter=0
    for key in data_keys:
        #print " data: ", key
        #print "saving order: ", order[counter], " data: ", key
        if order[counter] != 666:
            print "saving electrode: ", order[counter]-1, " from data: ", key
            ecp[order[counter]-1]=data[key]*1E+7 #Use 1E6 conversion to uV and then multiply by 10 to gain more resolution using vscale_HP parameter below
            counter+=1
        else:
            counter+=1

    ##*************************** SAVE COMPRESSED LOW PASS .TSF FILE *************************************

    if True:
        header = 'Test spike file '
        iformat = 1002
        n_electrodes = 16
        Compress_factor = 20
        SampleFrequency = 1000 * Compress_factor
        vscale_HP = .1
        print "original data length: ", len(ecp[0])
        n_vd_samples = int(len(ecp[0])/200.)
        print "n_vd_samples: ", n_vd_samples

        probe = [16, 1, 15, 2, 12, 5, 13, 4, 10, 7, 9, 8, 11, 6, 14, 3]

        SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (n_electrodes):
            SiteLoc[i][0]=0
            SiteLoc[i][1]=i*100
            
        file_name1 = file_[:-4] + '_lp_compressed.tsf'
        f1 = open(file_name1, 'wb')

        f1.write(header)
        f1.write(struct.pack('i', iformat))
        f1.write(struct.pack('i', SampleFrequency))
        f1.write(struct.pack('i', n_electrodes))
        f1.write(struct.pack('i', n_vd_samples))
        f1.write(struct.pack('f', vscale_HP))
        for i in range (n_electrodes):
            f1.write(struct.pack('h', SiteLoc[i][0]))
            f1.write(struct.pack('h', SiteLoc[i][1]))
            f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

        ecp_lp=ecp.copy()
        ecp_temp = ecp.copy()

        ecp_lp = ecp_temp - wavelet(ecp_lp, wname="db4", maxlevel=6)
        
        ecp_temp = []
        for i in range(len(ecp_lp)):
            #ecp_temp.append(ecp_lp[i][:n_vd_samples])
            ecp_temp.append(ecp_lp[i])
        ecp_temp = np.int16(ecp_temp)
        
        print "Writing data"
        for i in range(len(ecp_lp)):
            #ecp[i] = filter.notch(ecp[i])[0]
            #ecp.tofile(f1)
            ecp_temp[probe[i]-1][::200].tofile(f1)
            #ecp[np.where(order==i+1)[0][0]].tofile(f2)

        f1.write(struct.pack('i', 0)) #Write # of fake spikes
        f1.close()
        print "DONE"
        #quit()

#************************************* SAVE Low pass filter .tsf *******************
    if True:
        header = 'Test spike file '
        iformat = 1002
        n_electrodes = 16
        SampleFrequency = 1000
        vscale_HP = .1
        print "original data length: ", len(ecp[0])
        n_vd_samples = int(len(ecp[0])/20.)
        print "n_vd_samples: ", n_vd_samples

        probe = [16, 1, 15, 2, 12, 5, 13, 4, 10, 7, 9, 8, 11, 6, 14, 3]

        SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (n_electrodes):
            SiteLoc[i][0]=0
            SiteLoc[i][1]=i*100
            
        file_name1 = file_[:-4] + '_lp.tsf'
        f1 = open(file_name1, 'wb')

        f1.write(header)
        f1.write(struct.pack('i', iformat))
        f1.write(struct.pack('i', SampleFrequency))
        f1.write(struct.pack('i', n_electrodes))
        f1.write(struct.pack('i', n_vd_samples))
        f1.write(struct.pack('f', vscale_HP))
        for i in range (n_electrodes):
            f1.write(struct.pack('h', SiteLoc[i][0]))
            f1.write(struct.pack('h', SiteLoc[i][1]))
            f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

        ecp_lp=ecp.copy()
        ecp_temp = ecp.copy()

        ecp_lp = ecp_temp - wavelet(ecp_lp, wname="db4", maxlevel=6)
        
        ecp_temp = []
        for i in range(len(ecp_lp)):
            #ecp_temp.append(ecp_lp[i][:n_vd_samples])
            ecp_temp.append(ecp_lp[i])
        ecp_temp = np.int16(ecp_temp)
        
        print "Writing data"
        for i in range(len(ecp_lp)):
            #ecp[i] = filter.notch(ecp[i])[0]
            #ecp.tofile(f1)
            ecp_temp[probe[i]-1][::20].tofile(f1)
            #ecp[np.where(order==i+1)[0][0]].tofile(f2)

        f1.write(struct.pack('i', 0)) #Write # of fake spikes
        f1.close()
        print "DONE"
        
    ##***************************************** SAVE RAW .TSF FILE *************************************
    if True:
        header = 'Test spike file '
        iformat = 1002
        n_electrodes = 16
        SampleFrequency = 20000
        vscale_HP = .1
        n_vd_samples = len(ecp[0])
        print "n_vd_samples: ", n_vd_samples

        probe = [16, 1, 15, 2, 12, 5, 13, 4, 10, 7, 9, 8, 11, 6, 14, 3]

        SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (n_electrodes):
            SiteLoc[i][0]=0
            SiteLoc[i][1]=i*100
            
        file_name1 = file_[:-4] + '_raw.tsf'
        f1 = open(file_name1, 'wb')

        f1.write(header)
        f1.write(struct.pack('i', iformat))
        f1.write(struct.pack('i', SampleFrequency))
        f1.write(struct.pack('i', n_electrodes))
        f1.write(struct.pack('i', n_vd_samples))
        f1.write(struct.pack('f', vscale_HP))
        for i in range (n_electrodes):
                f1.write(struct.pack('h', SiteLoc[i][0]))
                f1.write(struct.pack('h', SiteLoc[i][1]))
                f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays
        
        print "Writing data"
        for i in range(len(ecp)):
            ecp[probe[i]-1].tofile(f1)

        f1.write(struct.pack('i', 0)) #Write # of fake spikes
        f1.close()
        print "DONE"




    ##***************************************** SAVE HIGH PASS .TSF FILE *************************************
    if True:
        header = 'Test spike file '
        iformat = 1002
        n_electrodes = 16
        SampleFrequency = 20000
        vscale_HP = .1
        n_vd_samples = len(ecp[0])
        print "n_vd_samples: ", n_vd_samples

        probe = [16, 1, 15, 2, 12, 5, 13, 4, 10, 7, 9, 8, 11, 6, 14, 3]

        SiteLoc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (n_electrodes):
            SiteLoc[i][0]=0
            SiteLoc[i][1]=i*100
            
        file_name1 = file_[:-4] + '_hp.tsf'
        f1 = open(file_name1, 'wb')

        f1.write(header)
        f1.write(struct.pack('i', iformat))
        f1.write(struct.pack('i', SampleFrequency))
        f1.write(struct.pack('i', n_electrodes))
        f1.write(struct.pack('i', n_vd_samples))
        f1.write(struct.pack('f', vscale_HP))
        for i in range (n_electrodes):
                f1.write(struct.pack('h', SiteLoc[i][0]))
                f1.write(struct.pack('h', SiteLoc[i][1]))
                f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

        print "Wavelet filtering data"
        ecp_hp = wavelet(ecp, wname="db4", maxlevel=6)
        
        print "Writing data"
        for i in range(len(ecp)):
            ecp_hp[probe[i]-1].tofile(f1)


        f1.write(struct.pack('i', 0)) #Write # of fake spikes
        f1.close()
        print "DONE"

