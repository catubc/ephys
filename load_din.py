import numpy as np
import matplotlib.pyplot as plt
import struct

file_name = '/home/cat/neuron/in_vivo/nick/ptc17/04-tr1-MVI_1400_5s.din'

din_data = np.fromfile(file_name,dtype=np.int64)
print len(din_data)
times = din_data[::2]
print times[0:1000]
print len(times)
frames = din_data[1::2]
print frames[1220:2220]

print len(frames)
