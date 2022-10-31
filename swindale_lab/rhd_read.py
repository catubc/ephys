'''Read .rhd files, convert to correct electrode mapping and save to .tsf file
NB: There are 2 possible mapping depending on the insertion of the AD converter 
TODO: implement a wavelet high pass filter directly to avoid SpikeSorter Butterworth filter artifacts
'''

from load_intan_rhd_format import *
import matplotlib.pyplot as plt

connector_map = []

#A64 Omnetics adaptor
adaptor_map = []
adaptor_map.append([34,35,62,33,60,54,57,55,10,8,11,5,32,3,30,31])
adaptor_map.append([64,58,63,56,61,59,52,50,15,13,6,4,9,2,7,1])
adaptor_map.append([53,51,49,47,45,36,37,38,27,28,29,20,18,16,14,12])
adaptor_map.append([48,46,44,42,40,39,43,41,24,22,26,25,23,21,19,17])

adaptor_layout1=[]      #Concatenated rows
for maps in adaptor_map:
    adaptor_layout1.extend(maps)

#Intan adapter - if inserted right-side up
intan_map = []
intan_map.append(list(reversed([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])))     #NB: need to reverse these arrays:  list(reversed(...))
intan_map.append(list(reversed([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])))
intan_map.append(list(reversed([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])))
intan_map.append(list(reversed([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])))

intan_layout1=[]
for maps in intan_map:
    intan_layout1.extend(maps)


#Intan adapter - if inserted upside-down; no need to reverse
intan_map = []
intan_map.append([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])
intan_map.append([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])
intan_map.append([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])
intan_map.append([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])

intan_layout3=[]
for maps in intan_map:
    intan_layout3.extend(maps)


#A1x64 probe layout
a = [27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,28,29,30,31,32]
b = [37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,36,35,34,33]
probe_map = a+b
probe_map[::2] = a
probe_map[1::2] = b


print "channel    x_loc    y_loc"
for i in range(64):
    print intan_layout1[adaptor_layout1.index(probe_map[i])], '    ', 30*(i%2), '    ', i*23
    
quit()

#print intan_layout3

ch_no = 0
#print probe_map[ch_no]
#print adaptor_layout1.index(probe_map[ch_no])
#print intan_layout3[adaptor_layout1.index(probe_map[ch_no])]

#quit()


#*****2016_04_26
#Track 1
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_1/2016_04_26_cortex_initial_insert_probe9093_160426_174828.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_1/2016_04_26_cortex_initial_insert_probe9093_160426_174231.rhd'

#Track 2
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_2/2016_04_26_cortex_initial_insert_2ndprobe_160426_180145.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_2/2016_04_26_cortex_initial_insert_2ndprobe_160426_180245.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_2/2016_04_26_cortex_initial_insert_2ndprobe_continuous_recording_160426_180425.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_2/2016_04_26_cortex_initial_insert_2ndprobe_continuous_recording_deep_160426_181346.rhd'

#Track 3
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_26/track_3/2016_04_26_cortex_initial_insert_2ndprobe_continuous_recording_3rd_track_deep_160426_185358.rhd'

#*****#2016-04-28 Saline tests
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_28/2015_04_28_test_3chprobe_160428_205917.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2015_04_28/2015_04_28_test_2ndprobe_160428_210514.rhd'

#*****#2016-05-03
#Track 1
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/2016_05_03_tr1_spont_afterinsertion_noisyamp_160503_145951.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr1/2016_05_03_tr1_spont_deep_noisyamp_160503_152854.rhd'

#Track 2
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_spont_afterinsertion_160503_162252.tsf'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_laser225_160503_163505.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_laser535_160503_165846.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_deep_spont_160503_175515.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_deep_laser535_160503_171740.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_deep_laser535_2_160503_172757.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_deep_laser535_3_160503_173756.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/2016_05_03_tr2_deep_laser535_4_160503_174635.rhd'

#Track 3
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/2016_05_03_tr3_spont_160503_182605.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/2016_05_03_tr3_deep_spont_160503_184142.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/2016_05_03_tr3_deep_laser535_160503_185155.rhd'
#file_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/2016_05_03_tr3_deep_laser535_2_160503_192518.rhd'

#*****2016-05-24
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spont_tr1_rec2_160524_215834.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spont_tr1_rec3_160524_222211.rhd'

#*****2016-05-25
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr1_deep_laser535_10ms_5sec_5rep_160525_174613.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr1_spont_afterinsertion_160525_170616.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr2_deep_laser535_10ms_5sec_5rep_160525_193958.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr2_spont_afterinsertion_160525_185339.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr2_transition_to_deep_160525_192816.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr32_deep_laser535_7ms_1sec_15rep_160525_224253.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr31_deep_laser535_7ms_5sec_5rep_160525_214419.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr3_spont_cortex_160525_205132.rhd'
#file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_05_25_chr2/rhd_files/2016_05_25_tr3_transition_to_deep_160525_212701.rhd'

#*****2016-05-26
file_names = [
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr1_spontaneous_cortex_01_160527_150948.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr1_spontaneous_cortex_02_160527_152937.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr1_spontaneous_cortex_03_160527_154808.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr1_spontaneous_cortex_04_160527_163101.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr2_spontaneous_cortex_05_160527_165921.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr2_spontaneous_cortex_06_160527_171746.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr2_spontaneous_cortex_07_160527_173500.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr2_spontaneous_cortex_08_160527_175843.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr3_spontaneous_cortex_09_160527_183014.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr3_spontaneous_cortex_10_160527_184615.rhd',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr3_spontaneous_cortex_11_160527_190528.rhd',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/rhd_files/2016_05_27_tr3_spontaneous_cortex_12_160527_192127.rhd',
]

#**************************** PROCESS DATA **************************
load_rhd = True
save_tsf = True
load_meta = True

for file_name in file_names:
    print ""; print ""
    print "Processing: \n", file_name
    #Load track name
    start = file_name.find('tr')
    end = file_name[start:].find('_')
    track = file_name[start:start+end]

    if load_rhd:
        if False:       #DON"T ERASE - NEED IT CONCATENATE 1MIN RECS OUT OF INTANT....
            file_names = [
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213221.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213321.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213421.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213521.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213621.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213721.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213821.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_213921.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214021.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214121.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214221.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214321.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214421.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214521.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214621.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214721.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214821.rhd',
            '/media/cat/12TB/in_vivo/tim/cat/2016-05-24/spontaneous_tr1_160524_214921.rhd',
            ]
            
            ec_traces = []
            for k in range(64):
                ec_traces.append([])
            for f_name in file_names:
                print "Loading: ", f_name
                data = read_data(f_name)
                for k in range(64):
                    ec_traces[k].extend(data['amplifier_data'][k])
            file_name = f_name
            
        else: 
            data = read_data(file_name)
            ec_traces = data['amplifier_data']

    #*************Save data to .tsf
    if save_tsf:

        SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
        print "SampleFrequency: ", SampleFrequency

        header = 'Test spike file '
        iformat = 1002
        n_vd_samples = len(ec_traces[0])
        n_electrodes = 64
        vscale_HP = 1
        n_cell_spikes = 0

        Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (n_electrodes):
            Siteloc[i][0]=30*(i%2)
            Siteloc[i][1]=i*23
        
        #ec_traces = np.float32(ec_traces)/vscale_HP
        print "Converting data to int16..."
        ec_traces = np.array(ec_traces, dtype=np.int16)
            
        #Set file name to correct directory
        file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:][:-4]+'.tsf'
        fout = open(file_out, 'wb')
        fout.write(header)
        fout.write(struct.pack('i', 1002))
        fout.write(struct.pack('i', SampleFrequency))
        fout.write(struct.pack('i', n_electrodes))
        fout.write(struct.pack('i', n_vd_samples))
        fout.write(struct.pack('f', vscale_HP))
        
        for i in range (n_electrodes):
            #print i+1, Siteloc[i]
            fout.write(struct.pack('h', Siteloc[i][0]))
            fout.write(struct.pack('h', Siteloc[i][1]))
            fout.write(struct.pack('i', i+1))
       
        for i in range(n_electrodes):
            #ec_traces[intan_layout3[adaptor_layout1.index(probe_map[i])]].tofile(fout)  #Backside adaptor insertion
            ec_traces[intan_layout1[adaptor_layout1.index(probe_map[i])]].tofile(fout)  #Frontside

        fout.write(struct.pack('i', n_cell_spikes))
    
    #Read and save digital input data
    laser_filename = file_name[0:file_name.find('rhd_files')]+'laser_files/'+file_name[file_name.find('rhd_files')+10:][:-4]+'_laser_times'
    meta_filename = file_name[0:file_name.find('rhd_files')]+'laser_files/'+file_name[file_name.find('rhd_files')+10:][:-4]+'_meta_data'
    camera_pulses_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:][:-4]+'_camera_pulses'
    camera_onoff_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:][:-4]+'_camera_onoff'

    SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
    print "SampleFrequency: ", SampleFrequency

    counter=0
    #laser_pulses = data['board_dig_in_data'][counter]; np.save(laser_filename, laser_pulses); counter+=1
    #meta_data = data['board_dig_in_data'][counter]; np.save(meta_filename, meta_data); counter+=1
    camera_frames = data['board_dig_in_data'][counter]; np.save(camera_frames_filename, camera_frames); counter+=1
    camera_onoff = data['board_dig_in_data'][counter]; np.save(camera_onoff_filename, camera_onoff); counter+=1
    #plt.plot(camera_onoff)
    #plt.show()



quit()
#**************OLD CONVERSION OF META DATA FROM SERIAL PORT ***** DON"T ERASE!!!!!!
if False:
    ##Quick plot meta/laser data
    #ax1 = plt.subplot(211)
    #plt.plot(board_dig_in_data[0])  #Laser pulse info
    #ax2 = plt.subplot(212)          #Meta data
    #plt.plot(board_dig_in_data[1])
    #plt.show()
    
    print len(laser_pulses)
    
    #Load laser on times from file
    laser_times = np.where(laser_pulses==1)[0]/float(SampleFrequency)

    #Parse laser times to find beginning of each chunk 
    print laser_times[0:2580]
    laser_start_times = []
    laser_start_times.append(laser_times[0])
    for k in range(len(laser_times)):
        if (laser_times[k]-laser_start_times[-1])>1.0:      #Search for next laser time starts excluding 1.0 seconds of time
            laser_start_times.append(laser_times[k])
            
    laser_start_times = np.array(laser_start_times)
    print laser_start_times
    print len(laser_start_times)
    
    #area_repeats = []
    #for k in range(24):
        #area_repeats.append(laser_start_times[k::24])
    #np.save(file_name[0:-4]+"_area_repeats", area_repeats)


    #Look at meta data
    n = meta_data
    print meta_data
    
    #plt.plot(meta_data[0:3000000])
    #plt.plot(meta_data[170160:172000])
    #plt.show()
    
    #parsing = True
    #for i in range(len(meta_data)):
    #    if (meta_data[i]!=1) and parsing:
    #        parsing = False

    print np.argmax(n==0)

    index_scale = SampleFrequency/9600.  #Divide SampleFrequency by baud rate... can't recall why;
    
    for i in range (170160,172000,1):
        print "i: ", i
        if meta_data[i]!=1:
            byte = ''
            for j in range(8):
                print 1.5*index_scale+i
                print n[1.5*index_scale+i]
                byte=byte + str(n[(1.5+j)*index_scale])
            
            print byte
            
            i+=8
            quit()




    for i in range(255):
        n=n[np.argmax(n==0):]
        #print n[0:35]
        byte = ''
        for j in range(8):
            #print 1.5*index_scale+i
            #print n[1.5*index_scale+i]
            byte=byte + str(n[(1.5+j)*index_scale])

        byte = byte[::-1]
        #print byte
        print int(byte, 2)
       # print chr (128+int(byte, 2))
        n=n[8*index_scale+3:]

    quit()
    



plt.plot(board_dig_in_data[0])
plt.show()
quit()
print binascii.unhexlify('%x' %n)

#for i in range(len(board_dig_in_data[0])):
    







amp_data = data['amplifier_data']

print amp_data

for i in range(len(amp_data)):
    plt.plot(amp_data[i][10000:20000]/5000. + i+1)

plt.ylabel("Channels")
plt.xlabel("timestep")
plt.show()


