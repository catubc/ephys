#Gcamp6 and ephys concatenation code

from PIL import Image
import struct
import numpy as np
import matplotlib.pyplot as plt
from tsf_ptcs_classes import *

def load_tsf(ephys_file_name):
    with open(ephys_file_name,'rb') as fin:
        header = fin.read(16)                                   #Header info; not used
        iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
        SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
        n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
        n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples 
        vscale_HP1 = struct.unpack('f',fin.read(4))[0]          #Scaling of int2 values below to save space, currently 0.1
        
        Siteloc1 = np.zeros((2*n_electrodes), dtype=np.int16)
        Readloc = np.zeros((n_electrodes), dtype=np.int32)
        for i in range(n_electrodes):
            Siteloc1[i*2] = struct.unpack('h', fin.read(2))[0]
            Siteloc1[i*2+1] = struct.unpack('h', fin.read(2))[0]
            Readloc[i] = struct.unpack('i', fin.read(4))[0]

        ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples) 
        
        #Reshape traces
        ec_traces.shape = n_electrodes, n_vd_samples
        fin.close()
    
    return ec_traces

def save_tsf(file_name, track, all_traces):
    SampleFrequency = 25000
    print "SampleFrequency: ", SampleFrequency

    header = 'Test spike file '
    iformat = 1002; n_electrodes = 64; vscale_HP = 1.0; n_cell_spikes = 0
    n_vd_samples = len(all_traces[1])
    print n_vd_samples

    Siteloc = np.zeros((n_electrodes,2), dtype=np.int16) #Read as 1D array
    for i in range (n_electrodes): 
        Siteloc[i][0]=30*(i%2)
        Siteloc[i][1]=i*23
        
    #Set file name to correct directory
    file_out = file_name[:file_name.find('tsf_files')+10]+track+'.tsf'

    fout = open(file_out, 'wb')
    fout.write(header)
    fout.write(struct.pack('i', 1002))
    fout.write(struct.pack('i', SampleFrequency))
    fout.write(struct.pack('i', n_electrodes))
    fout.write(struct.pack('i', n_vd_samples))
    fout.write(struct.pack('f', vscale_HP))
    
    for i in range (n_electrodes): 
        fout.write(struct.pack('h', Siteloc[i][0]))
        fout.write(struct.pack('h', Siteloc[i][1]))
        fout.write(struct.pack('i', i+1))
   
    for i in range(n_electrodes):
        all_traces[i].tofile(fout)  #Frontside

    fout.write(struct.pack('i', n_cell_spikes))

def load_tif(img_file):
    if (os.path.exists(img_file[:-4]+'.npy')==False):
        img = Image.open(img_file)
        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1
                print counter
        img.seek(0)

        n_pixels = img.size[0]
        n_frames = counter
        images_raw = np.zeros((n_frames, 128, 128), dtype = np.float16)
        for i in range(0, n_frames,1): 
            try:
                img.seek(i)
                print "Loading frame: ", i
                #images_raw [i] = np.flipud(np.fliplr(np.float16(img))) #FLIP IMAGES FOR Experiments Nov and Dec 2015
                images_raw[i] = np.float16(img) #2016-1-11 2016-1-14 experiment no flipping needed
            except EOFError:
                break

        np.save(img_file[:-4], images_raw)
        return images_raw

    else:
        return np.load(img_file[:-4]+'.npy')
        


ephys_file_names = [
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr1_spontaneous_cortex_01_160527_150948.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr1_spontaneous_cortex_02_160527_152937.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr1_spontaneous_cortex_03_160527_154808.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr1_spontaneous_cortex_04_160527_163101.tsf'
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr2_spontaneous_cortex_05_160527_165921.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr2_spontaneous_cortex_06_160527_171746.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr2_spontaneous_cortex_07_160527_173500.tsf',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr2_spontaneous_cortex_08_160527_175843.tsf'
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr3_spontaneous_cortex_09_160527_183014.tsf',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr3_spontaneous_cortex_10_160527_184615.tsf',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr3_spontaneous_cortex_11_160527_190528.tsf',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/2016_05_27_tr3_spontaneous_cortex_12_160527_192127.tsf'
]

camera_onoff_names = [
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr1_spontaneous_cortex_01_160527_150948_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr1_spontaneous_cortex_02_160527_152937_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr1_spontaneous_cortex_03_160527_154808_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr1_spontaneous_cortex_04_160527_163101_camera_onoff.npy'
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr2_spontaneous_cortex_05_160527_165921_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr2_spontaneous_cortex_06_160527_171746_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr2_spontaneous_cortex_07_160527_173500_camera_onoff.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr2_spontaneous_cortex_08_160527_175843_camera_onoff.npy'
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr3_spontaneous_cortex_09_160527_183014_camera_onoff.npy',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr3_spontaneous_cortex_10_160527_184615_camera_onoff.npy',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr3_spontaneous_cortex_11_160527_190528_camera_onoff.npy',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/camera_files/2016_05_27_tr3_spontaneous_cortex_12_160527_192127_camera_onoff.npy'
]

img_files = [
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/01_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/02_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/03_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/04_sponGCAMP6f_30Hz_9000fr_iso1.tif'
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/05_sponGCAMP6f_30Hz_9000fr_iso1.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/06_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/07_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
#'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/08_sponGCAMP6f_30Hz_9000fr_iso1.25.tif'
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/09_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/10_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/11_sponGCAMP6f_30Hz_9000fr_iso1.25.tif',
'/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/tif_files/12_sponGCAMP6f_30Hz_9000fr_iso1.25.tif'
]

#Load track name
start = ephys_file_names[0].find('tr')
end = ephys_file_names[0][start:].find('_')
track = ephys_file_names[0][start:start+end]
print "Track: ", track

#****** LOAD DATA ? CONVERT TO SINGLE TRACK: .TSF AND .NPY IMAGES
if False:
    all_traces = []
    all_img = []
    for ephys_file_name, camera_onoff_name, img_file in zip(ephys_file_names, camera_onoff_names, img_files):
        print ""
        print ephys_file_name
        #print camera_onoff_name

        #Load ephys
        ec_traces = load_tsf(ephys_file_name)

        #Load camera onoff data
        camera_onoff = np.load(camera_onoff_name)
        camera_on = np.where(camera_onoff==1)[0][0]
        camera_off = np.where(camera_onoff==1)[0][-1]
        
        ec_traces = ec_traces[:,camera_on:camera_off]
        all_traces.append(ec_traces)
        
        #Load imaging data
        img = load_tif(img_file)
        all_img.append(img)

    all_traces = np.hstack((all_traces))
    print all_traces.shape
    save_tsf(ephys_file_names[0], track, all_traces)

    all_img = np.concatenate(np.array(all_img), axis=0)
    print all_img.shape
    file_out = img_files[0][:img_files[0].find('tif_files')+10]+track
    np.save(file_out, all_img)

#****** SAVE CAMERA TIMES
#Need to save image times for concatenated recording as the sample rate slightly varies between recs
#THIS IS A BIT TRICKY... SEEMS STILL LOOSING ~50ms for a 1hr recording...
times=[]
offset=0
for camera_onoff_name, img_file in zip(camera_onoff_names, img_files):
    print img_file
    camera_onoff = np.load(camera_onoff_name)
    camera_on = np.where(camera_onoff==1)[0][0]
    camera_off = np.where(camera_onoff==1)[0][-1]
    
    #Load imaging data
    img = load_tif(img_file)
    t = np.linspace(0, camera_off-camera_on, len(img)+1)+offset
    times.extend(t[:-1])

    offset = t[-1]
    print t[0], t[-2]

times = np.array(times)
file_out = camera_onoff_names[0][:camera_onoff_names[0].find('camera_files')+13]+track
np.save(file_out, times)





