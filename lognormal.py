import numpy as np
import matplotlib.pyplot as plt

ax1 = plt.subplot(1,1,1)

mu, sigma = 1, .5 # mean and standard deviation
s = np.random.lognormal(mu, sigma, 1000)
print float(len(s[s>2.0]))/len(s)
quit()
count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')

ax1.set_xscale('symlog', linthreshy=0.0001, basex=10)
plt.xlim(0.001,10)

plt.axis('tight')

plt.show()

