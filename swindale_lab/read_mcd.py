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


file_name = '/media/cat/4TB/in_vivo/tim/2015-7-22/2015-7-22-9-W1/2015-7-22-9-W1.mcd'

data = MCD_read(file_name)

print data
print dir(data)
print data[0]

for entity in data:
    print entity.label, entity.entity_type

