#Sergey Gratiy's method for computing smooth CSD

import numpy as np
from csd_est_funcs import *
import pylab as plt

# units:

# LFP should be in (mV)

# the CSD will be in (uA/mm^3)



nch=32  # number of channels should correspond to the dimensions of the data
channels = np.arange(0,nch)
dz=0.02 # interchannel distance in (mm)
zs = channels*dz    # channel locations
sigma =0.3  # extracellular conductivity (mS/mm)
b=0.5   # radius of the column (mm)
SNR=20.0    # average ration of signla to noise on each channel


[A,A0] = forward_operator(zs,b,sigma) # A -dimensional (mm^3/mS) and A0 -dimensionless operators 
[W,R] = inverse_tikhonov(A,SNR) # mS/mm^3
[W0,R0] = inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check

[fig100,fig101] = show_operators(A,A0,W,W0,R,R0)    # display forward and ivnerse operators 



fdir='/home/sergeyg/ephys/fsf_lfp_oblique_sev/'
fname='M72.1CDS1.lfp.npy'
lfp72 = np.load(fdir+fname)

csd72=np.dot(W,lfp72)   # units: mS/mm^3*mV = uA/mm^3

LFP=lfp72
CSD=csd72


fig1 = plt.figure(1)
plt.subplot(2,3,1);
VspanLFP = np.array([-1, 1])*abs(LFP).max()
plt.imshow(LFP, vmin = VspanLFP[0], vmax = VspanLFP[1],cmap='jet',aspect='auto'); 
plt.colorbar();
plt.title('LFP')


plt.subplot(2,3,4);
VspanCSD = np.array([-1, 1])*abs(CSD).max()
plt.imshow(CSD, vmin = VspanCSD[0], vmax = VspanCSD[1],cmap='jet',aspect='auto'); 
plt.colorbar();
plt.title('CSD')




#extent=[time[0],time[-1],z[-1],z[0]],cmap='jet',origin='upper',aspect='auto'    
plt.show()
