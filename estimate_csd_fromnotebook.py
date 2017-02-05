
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from math import *
import h5py

# reload external librariesthat have changed during execution
get_ipython().magic(u'load_ext autoreload')

# check on every execution call
get_ipython().magic(u'autoreload 2')


# In[4]:

datadir = '/data/mat/Cat/Sergey_data/'
fname = 'M71.hdf5'
f5 = h5py.File(datadir+fname,'r')
print f5.keys()
print f5['/CSD3'].keys()


# In[5]:

ecp_trial_avgCSD3 = f5['/DeepCSD4/stim_triggered_ecp'][...]
ecp_trial_avgCSD4 = f5['/DeepCSD5/stim_triggered_ecp'][...]
probe_layout = f5['/CSD3/probe_layout'][...]
vscale = f5['/CSD3/vscale'][...]


# In[6]:

z = probe_layout[:,1]*1E-3
dt=1./ecp_trial_avgCSD3.shape[1]
t=np.arange(0,ecp_trial_avgCSD3.shape[1])*dt
print t[1:3], t[-1]
print vscale


# In[7]:

import imp
csdmod = imp.load_source('csd_est_funds','/home/sergeyg/mypy/csd_estimation/csd_est_funcs.py')


sigma =0.3  # extracellular conductivity (mS/mm)
b=1.0   # assumed radius of the column (mm)
SNR=20.0    # average ratio of signal to noise on each channel


[A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
[W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
[W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check

#[fig100,fig101] = csdmod.show_operators(A,A0,W,W0,R,R0)    # display forward and inverse operators 
#plt.show()


# In[ ]:

csd3=np.dot(W,ecp_trial_avgCSD3)   # units: mS/mm^3*mV = uA/mm^3


lfp3_0mean = np.mean(ecp_trial_avgCSD3, axis=0)
lfp3 = ecp_trial_avgCSD3-lfp3_0mean



tlims =[0,0.5]

fig1=plt.figure(1)
# csd3
LFP = ecp_trial_avgCSD3
plt.subplot(2,3,1);
VspanLFP = np.array([-1, 1])*abs(LFP).max()
plt.imshow(LFP, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.xlim(tlims)

plt.colorbar();
plt.title('LFP')

LFP = lfp3
plt.subplot(2,3,2);
VspanLFP = np.array([-1, 1])*abs(LFP).max()
plt.imshow(LFP, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.xlim(tlims)
plt.colorbar();
plt.title('LFP zero mean')

CSD=csd3
plt.subplot(2,3,3);
VspanLFP = np.array([-1, 1])*0.8*abs(CSD).max()
plt.imshow(CSD, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.xlim(tlims)
plt.colorbar();
plt.title('CSD')



#csd4

csd4=np.dot(W,ecp_trial_avgCSD4)   # units: mS/mm^3*mV = uA/mm^3

lfp4_0mean = np.mean(ecp_trial_avgCSD3, axis=0)
lfp4 = ecp_trial_avgCSD4-lfp4_0mean

LFP = ecp_trial_avgCSD4
plt.subplot(2,3,4);
VspanLFP = np.array([-1, 1])*abs(LFP).max()
plt.imshow(LFP, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.colorbar();
plt.title('LFP')
plt.xlim(tlims)

LFP = lfp4
plt.subplot(2,3,5);
VspanLFP = np.array([-1, 1])*abs(LFP).max()
plt.imshow(LFP, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.colorbar();
plt.title('LFP zero mean')
plt.xlabel('time (s)')
plt.xlim(tlims)

CSD=csd4
plt.subplot(2,3,6);
VspanLFP = np.array([-1, 1])*0.8*abs(CSD).max()
plt.imshow(CSD, vmin = VspanLFP[0], vmax = VspanLFP[1],extent=[t[0],t[-1],z[-1],z[0]],cmap='jet',aspect='auto'); 
plt.colorbar();
plt.title('CSD')
plt.xlim(tlims)





fig2 = plt.figure(2)

plt.subplot(2,1,1); plt.plot(t,lfp3_0mean)
plt.xlim(tlims)

plt.subplot(2,1,2); plt.plot(t,lfp3_0mean)
plt.xlim(tlims)

plt.xlabel('time (s)')



plt.show()


# In[ ]:



