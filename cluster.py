from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *
from random import randrange

from scipy import stats
import numpy as np
import time, math
import sys
import os.path
from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab

import pyximport
pyximport.install(build_in_temp=False, inplace=True)


#colors = np.array(['blue', 'red', 'green', 'orange', 'magenta', 'cyan'])
colors = np.array(['orange', 'green', 'magenta', 'cyan'])

#file_name = '/media/cat/4TB/in_vivo/nick/ptc21/67-tr5c-blankscreen/LFP_clusters.dat'

file_name = '/media/cat/4TB/in_vivo/nick/ptc21/track5c_59-75_notch_1_240Hz_5subsample/3_vs_6.dat'

pca_data = np.loadtxt(file_name)
pca_data = pca_data.astype(np.float)

pca_data = pca_data[:,0:2]
#pca_data = np.array([pca_data[:,0], pca_data[:,1]])
print pca_data.shape

#Cluster the log-log ISI plots
#K Means

if True:
    n_clusters = 2
    #from sklearn import cluster, datasets
    #clusters = cluster.KMeans(n_clusters)
    #clusters.fit(pca_data)
    
    labels_ = KMEANS(pca_data, n_clusters)
    
    

#Mean shift
if False:
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from sklearn.datasets.samples_generator import make_blobs

    ###############################################################################
    # Generate sample data
    #centers = [[1, 1], [-1, -1], [1, -1]]
    #X, _ = make_blobs(n_samples=10000, centers=centers, cluster_std=0.6)

    X=pca_data
        
    ###############################################################################
    # Compute clustering with MeanShift

    # The following bandwidth can be automatically detected using
    bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=len(X))

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)
    #print labels_unique
    #print n_clusters_
    #print labels

    print("number of estimated clusters : %d" % n_clusters)
    #quit()

ax = plt.subplot(1,1,1)

#plt.scatter(np.array(pca_data[0]), np.array(pca_data[1]), s=5, color=colors[clusters.labels_])
#plt.scatter(pca_data.T[0], pca_data.T[1], s=40, color=colors[clusters.labels_])
plt.scatter(pca_data.T[0], pca_data.T[1], s=40, color=colors[labels_])

#mng = plt.get_current_fig_manager()
#mng.resize(*mng.window.maxsize())
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.show()
