'''Compute Firing rate distrbutions from ptcs files using last spike time across all units; necessary without loading .tsf files
'''

from tsf_ptcs_classes import *
import matplotlib.pyplot as plt

#Track 1
#dir_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr1/'
#file_name = '2016_05_03_tr1_spont_afterinsertion_noisyamp_sorted_160503_145951'
#file_name = '2016_05_03_tr1_spont_deep_noisyamp_sorted_160503_152854'

#Track 2
dir_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr2/'
file_name = '2016_05_03_tr2_deep_laser535_sorted_160503_171740'
#file_name = '2016_05_03_tr2_deep_laser535_2_sorted_160503_172757'
#file_name = '2016_05_03_tr2_deep_laser535_3_sorted_160503_173756'
#file_name = '2016_05_03_tr2_deep_laser535_4_sorted_160503_174635'

#Track 3
#dir_name = '/media/cat/12TB/in_vivo/cat/2016_05_03/tr3/'
#file_name = '2016_05_03_tr3_deep_laser535_sorted_160503_185155'
#file_name = '2016_05_03_tr3_deep_laser535_2_160503_192518'


work_dir = dir_name
sort_file = file_name
Sort = Loadptcs(sort_file, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
Sort.name=file_name
Sort.filename=sort_file
Sort.directory=work_dir

print "# units: ", len(Sort.units)

#Compute approximate length of recording; THIS IS NOT ACCURATE
rec_length = []
for k in range(len(Sort.units)):
    rec_length.append(np.max(Sort.units[k]))
rec_length = np.max(rec_length)/Sort.samplerate
print rec_length/60.

#Firing rate distributions
fire_rate = []
for k in range(len(Sort.units)):
    fire_rate.append(len(Sort.units[k])/rec_length)

print fire_rate

plt.title("Fire Rate Distributions")
plt.tick_params(axis='both', which='major', labelsize=30)
x = np.linspace(0,50,100)
y = np.histogram(fire_rate, bins = x)
plt.plot(x[:-1], y[0], color='blue', linewidth=5)
plt.suptitle(file_name+"\n# units: "+str(len(Sort.units)), fontsize=30)
plt.show()
