import matplotlib.pyplot as plt
import numpy as np
from scipy import misc
import cPickle as pickle
from libtiff import TIFF


main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/'

data_file = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit16_ch09_all_3sec_window_02036_spikes.npy'

#data = np.nan_to_num(np.load(data_file))

data = np.load(data_file)

n_pixels = len(data[0])
generic_mask_file = main_dir +'genericmask.txt'
coords_generic = np.loadtxt(generic_mask_file)
v_min = np.nanmin(data)
#Ablate generic map
for i in range(len(data)):
    for j in range(len(coords_generic)):
        data[i][min(n_pixels-1,int(coords_generic[j][0]))][min(n_pixels-1,int(coords_generic[j][1]))]=v_min
                

#print np.min(data), np.max(data)
#for k in range(len(data)):
#    plt.imshow(data[90+k], vmin=np.min(data), vmax=np.max(data), interpolation='none')
#    plt.show()

print data.shape
print data.dtype
data = (data-np.min(data))/(np.max(data)-np.min(data))*255

data = np.uint8(data)


print data.dtype

tiff = TIFF.open(data_file+'.tif', mode='w')
tiff.write_image(data)
tiff.close()

print data[0]
