import matplotlib.pyplot as plt

import numpy as np
x = [1,2,3,1,2]; 
y = np.random.rand(10)
out1 = np.correlate(x,y,mode='full')
out2 = np.correlate(y,x,mode='full')
plt.plot(out1)
plt.plot(out2)
plt.show()
