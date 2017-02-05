from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import time
from pylab import *
from scipy.interpolate import interp1d
import matplotlib.mlab as mlab
import math

np.random.seed(12345678)  #fix random seed to get the same result

sim_dir = '/media/cat/Data1/in_vitro/'

cell_name=[None]*12
cell_name[1]='040912_3_4_2' 
cell_name[2]='040912_3_2_2' 
cell_name[3]='040912_3_4_1' 
cell_name[4]='040912_3_2_1' 
cell_name[5]='060912_3_2_1' 
cell_name[6]='080612_2_2_1' 
cell_name[7]='080612_3_1_1' 
cell_name[8]='290512_1_2_1' 
cell_name[9]='300412_2_2_1' 
cell_name[10]='300412_2_2_2' 
cell_name[11]='310812_2_4_1' 

numBins = 50
binwidth = 1

KS = [None]*11
SAT = [None]*11
for i in range(11):      

    file_name = sim_dir + cell_name[i+1] + '/' + 'PCA_values.dat'
    #electrode_array = np.genfromtxt(file_name,delimiter="\n")
    f = open(file_name,'r')
    lines = f.readlines()
    f.close()
    x=[]
    y=[]
    for line in lines:
        p = line.split()
        x.append(float(p[0]))
        y.append(float(p[1]))

    x = np.asarray(x)
    y = np.asarray(y)

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    x = x-sum(x)/len(x)
    y = y-sum(y)/len(y)
    
    i+=1    
    fig6 = plt.figure(6)
    ax = plt.subplot(6, 2, i)
    
    bb = axes([0.12+((i+1)%2)*0.43, 0.92-((i+1)/2)*.137, .11, .11], frameon=False)
    #t = arange(0.0, 10.0, 0.001)
    #r = exp(-t[:1000]/0.05)      
    plt.scatter(x,y,s=.01, color='blue')
    setp(bb, xticks=[], yticks=[])
    #title('PCA Plot'+str(i))
    if(i==1): 
        ax.text(0.10,0.85, 'PCA',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes,
            color='black', fontsize=11)

    x1 = ax.hist(x, bins=range(int(min(x)), int(max(x)) + binwidth, binwidth),color='blue',alpha=0.4)

    n_len = len(x)
    qk = np.random.normal(0, np.std(x), n_len)
    KS[i-1] = stats.ks_2samp(x, qk)[1]
    a = np.histogram(qk, bins=range(int(min(qk)), int(max(qk)) + binwidth, binwidth))
    
    mean = 0
    variance = max(a[1])-min(a[1])
    sigma = math.sqrt(variance)
    x = np.arange(min(a[1]),max(a[1]),.1)
    ax.plot(x,(mlab.normpdf(x,mean,sigma)/max(mlab.normpdf(x,mean,sigma)))*max(a[0]), linewidth=2, color='red')
    

    ax.text(0.34,0.83, "Cell: " + cell_name[i], 
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=11,fontweight='bold')

    ax.text(0.80,0.85, 'KS: ' + str(round(KS[i-1]*100,2))+'%',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=11)

    ax.set_yticklabels([])
    ax.set_xlim((-100,40))

    #setp(a, xlim=(0,.2), xticks=[], yticks=[])


    #plt.plot(x1[0]) #, markersize=2, linewidth=2)
    #plt.plot(y1) #, markersize=2, linewidth=2)
    #plt.title(group.name);    plt.xlim(visparams.tlim)

    #if ig>6:            plt.xlabel('time (s)'); 

    file_name = sim_dir + cell_name[i] + '/' + 'SAT_values.dat'
    #electrode_array = np.genfromtxt(file_name,delimiter="\n")
    f = open(file_name,'r')
    lines = f.readlines()
    f.close()
    x=[]
    y=[]
    for line in lines:
        p = line.split()
        y.append(float(p[1]))

    y = np.asarray(y)
    SAT[i-1] = np.average(y)

plt.show()



#plt.scatter(KS,SAT,s=300, color='blue')
#plt.show()


x = KS
y = SAT

coefficients = polyfit(x, y, 2)
polynomial = poly1d(coefficients)
xs = arange(0.,max(x), 0.001)
ys = polynomial(xs)

f = figure()
ax = f.add_subplot(111)

ax.scatter(x, y, s=200, color='blue')
ax.plot(xs, ys, color='blue')
ax.set_xlim((-0.02,max(x)))
ax.set_ylim((0,100))


ax.text(0.35,0.25, "Red: Spike Detection Threshold (5.0 x Noise)", 
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='red', fontsize=13,fontweight='bold')

ax.text(0.8,0.43, "Blue: poly fit (2nd order)", 
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='blue', fontsize=13,fontweight='bold')

ax.text(0.20,0.83, "Dashed Black Line: 10% KS null hypothesis", 
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='black', fontsize=13,fontweight='bold')

ax.text(0.54,0.96, "More Gaussian", 
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='black', fontsize=13,fontweight='bold')

ax.text(0.02,0.96, "Less Gaussian", 
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='black', fontsize=13,fontweight='bold')

x = [0,1]
y = [30,30]
ax.plot(x,y,color='red', linestyle='--', lw=3)


x = [0.1,0.1]
y = [0,100]
ax.plot(x,y,color='black', linestyle='--', lw=3)


title('KS value vs. PTP (Noise ~6uV) ')

ylabel('PTP Amplitude (uV)')
xlabel('Kolmogorov-Smirnov P-value')
plt.show()



#print len(x1), x1
#print len(y1),y1
#n1=len(x1)
#n2=len(y1)

#rvs1 = stats.norm.rvs(size=n1, loc=0., scale=0)
#rvs2 = stats.norm.rvs(size=n2, loc=0., scale=0)



