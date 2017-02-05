import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse, io
import cPickle as pickle


data = '/media/cat/12TB/in_vivo/seamans/barak/Front_Beh_2014-04-26_12-27-45/CSC9.npy'

data = np.load(data)
print len(data)
plt.plot(data[:1000000])
plt.show()


with open(data, 'rb') as infile:
    data = pickle.load(infile)

print data
data = sparse.csr_matrix(data)
data = data.toarray()

print data.shape
#data = sparse.coo_matrix((y['data'],(y['row'],y['col'])),shape=y['shape'])

plt.imshow(data, interpolation='none')
plt.show()
