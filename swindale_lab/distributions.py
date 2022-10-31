#Code to plot the distributions of amplitude values across all channels of an electrode
#Nick and Cat's work on skewed gaussian distributions of data
import numpy as np
import scipy.optimize as optimization
from scipy.integrate import simps, trapz
import matplotlib.pyplot as plt

from tsf_ptcs_classes import *

class distributions(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, tsf, depth, Sort1):

        #Run distribution plots; can further modularize this  as needed;
        self.plot_distributions(tsf, depth, Sort1)

    def plot_distributions(self,tsf, depth, Sort1):
        #tsf.ec_distribution.shape = tsf.n_electrodes * tsf.n_vd_samples, 

        print depth

        #hist, bin_edges = np.histogram(tsf.ec_distribution, density=True)
        #bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
        #p0 = [1., 0., 1.]
        graph11=[]
        graph12=[]
        graph13=[]
        graph21=[]
        graph22=[]
        graph23=[]

        n_electrodes = tsf.n_electrodes
        for i in range(n_electrodes):
            print "Ch: ", i 
            y = np.histogram(tsf.ec_distribution[i], bins=np.arange(-250, 250, 1))
            #yy = hist(tsf.ec_distribution[i], bins=np.arange(-250, 250, 1.0))[0]

            #plt.plot(y[1][:246],y[0][:246], color='black', linewidth=2, alpha=1)
            #plt.plot(y[1][255:-1],y[0][255:], color='black', linewidth=2, alpha=1)
            #plt.xlim((-240,240))

            #plt.axvspan(-5.0,5.0, color='black', alpha=0.3)

            #xx=[-5.0, -5.0]
            #yy=[0,y[0][245]]
            #print xx, yy
            #plt.plot(xx,yy, 'r--', color='black', linewidth=1, alpha=1)

            #xx=[5.0, 5.0]
            #yy=[0,y[0][255]]
            #print xx, yy
            #plt.plot(xx,yy, 'r--', color='black', linewidth=1, alpha=1)

            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.show()


            #PLOT BLUE PLOTS **********************************


            #plt.axvspan(-5.0,5.0, color='black', alpha=0.3)
            data=np.zeros((245), dtype=np.int32)
            data=y[0][0:245]
            x=np.zeros((245), dtype=np.int32)
            x[0:245] = np.arange(5, 250, 1.0)

            p0 = [-1000., 0, 1000]
            sigma = np.zeros((len(x)), dtype=np.int32)+1
            #data+=1 #Must add 1 or get singularities in fitting curve
            A, mu, S = optimization.curve_fit(gauss, x, data, p0, sigma)[0]
            xx=np.arange(-250,500,1.0)
            Gaussian2 = A*np.exp(-(xx-mu)**2/(2.*S**2))
            a1= simps(Gaussian2, dx=5)


            #plt.plot(y[1][:246],y[0][:246], color='black', linewidth=2, alpha=1)
            #plt.plot(y[1][255:-1],y[0][255:], color='black', linewidth=2, alpha=1)
            #plt.xlim((-240,240))
            #plt.plot(xx-255, Gaussian2, color='blue', linewidth=4, alpha='0.5')
            #plt.plot(x-255, data, 'r--', color='blue', linewidth=2, alpha='0.5')
            #plt.ylim((0,10000))

            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.show()
            #*********************MAGENTA CURVE******************

            data=np.zeros((245), dtype=np.int32)
            print len(y[251:])
            data[0:245]=y[0][251:251+245]
            x=np.zeros((245), dtype=np.int32)
            x[0:245] = np.arange(-250, -5, 1.0)

            p0 = [-1000., 0, 1000]
            sigma = np.zeros((len(x)), dtype=np.int32)+1
            A, mu, S = optimization.curve_fit(gauss, x, data, p0, sigma)[0]
            xx=np.arange(-500,250,1.0)
            Gaussian2 = A*np.exp(-(xx-mu)**2/(2.*S**2))
            a2 = simps(Gaussian2, dx=5)

            #plt.plot(y[1][:246],y[0][:246], color='black', linewidth=2, alpha=1)
            #plt.plot(y[1][255:-1],y[0][255:], color='black', linewidth=2, alpha=1)
            #plt.xlim((-240,240))
            #plt.axvspan(-5.0,5.0, color='black', alpha=0.3)
            #plt.plot(xx+251, Gaussian2, color='red', linewidth=4, alpha='0.5')
            #plt.plot(x+251,data, 'r--', color='red', linewidth=2, alpha='0.5')
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.show()

            #graph21.append(sum(data))
            #graph22.append(mu)
            #graph23.append(S)
            graph11.append(sum(data))
            graph12.append(mu)
            graph13.append(S)
            print A, mu, S

            #print A, mu, S

        x=[]

        for i in range(n_electrodes):
            x.append(tsf.Siteloc[i*2+1])

        print "Siteloc: ", x

        #graph1.sort()
        #graph2.sort()
        #graph3.sort()
        #print graph1

        #plt.scatter (x,graph1)
        #plt.plot (x,graph1, color='red', linewidth=3)
        #plt.ylim((0,1.5))
        #print max(x)
        #print len(x)

        #xs = linspace(0, max(x), len(x))
        #plt.scatter(x,graph1, color='red', alpha=1.0)
        #plt.plot(x,graph1, color='black', alpha=1.0)
        #plt.scatter(x,graph1, color='black', alpha=0.7) 
        #plt.show()

        inds = np.array(x).argsort()
        #y_sorted = np.array(graph11)[inds]
        #plt.scatter(np.sort(x), y_sorted, color='red', alpha=1.0, s=30)
        #plt.plot(np.sort(x), y_sorted, color='red', alpha=1.0)

        #y_sorted = np.array(graph21)[inds]
        #plt.scatter(np.sort(x), y_sorted, color='blue', alpha=1.0, s=30)
        #plt.plot(np.sort(x), y_sorted, color='blue', alpha=1.0)
        #plt.show()

        #y_sorted = np.array(graph12)[inds]
        #plt.scatter(np.sort(x), y_sorted, color='red', alpha=1.0, s=30)
        #plt.plot(np.sort(x), y_sorted, color='red', alpha=1.0)

        ax = plt.subplot(1, 1, 1)
        ax.set_ylabel('Asymetry index value', fontsize=16)

        y_sorted = np.array(graph12)[inds]/min(graph12)
        ax.set_ylim(0,max(y_sorted)) #Max number of units per channel
        print "Sorted: ", y_sorted
        plt.scatter(np.sort(x), y_sorted, color='blue', alpha=1.0, s=30)
        plt.plot(np.sort(x), y_sorted, color='blue', alpha=1.0)

#        plt.bar(24, autoclean_time, 8, color='magenta')

        ax2 = ax.twinx()

        ax2.set_ylabel('Number of Units on max channel', fontsize=16)
        #ax2.set_xlim(0,150)
        #ax2.set_xlabel('PTP (uV)',fontsize=20)

        #depth_hist = np.histogram(depth, bins=np.arange(0,1500,22))[0]
        #depth_hist = np.histogram(depth, bins=range(0,max(tsf.Siteloc),33))
        ##print depth_hist[0]
        ##print max(depth_hist[0])

        #depth=depth_hist[0]
        #print depth
        ##depth=(depth)/float(max(depth_hist[0])) 
        ##print depth
        ##depth_hist/max(depth_hist)
        #print "Units depth_histogram", depth_hist[1]
        #for i in range (len(depth)):
            #print i, "depth: ", depth_hist[1][i], " # cells: ", depth[i]
            #plt.bar(depth_hist[1][i]*Sort, depth[i],10, color='red', alpha=0.5)

        ##for i in range(len(depth)):
        ##    plt.bar(depth[i], 10, 10, color='green')
        #ax2.set_ylim(0, max(depth_hist[0])) #Max number of units per channel

        #VERSION 2: Scale by # of spikes in unit
        depth=[]
        metric=[]
        for k in range(Sort1.nneurons):
            metric.append([Sort1.chanpos[Sort1.maxchan[k]][1], len(Sort1.units[k])*Sort1.ptp[k]])

        #for k in range(len(metric)):
            #for j in range(len(metric)):
                #print metric[k][0], metric[j][0]
                #if metric[k][0]==metric[j][0]:
                    #temp = metric[k][0]
                    #metric[k][0] = 0 
                    #metric[j][1]+= temp #Cat: pool all the similar depth indexes together;

        #depth_hist = np.histogram(depth, bins=range(0,max(tsf.Siteloc),22)) #**********THIS STEP SIZE NEEDS TO CHANGE!!!!
        #depth=depth_hist[0]

        print metric
        #print "Units depth_histogram", depth_hist[1]
        for i in range (len(metric)):
            #print i, "depth: ", depth_hist[1][i], " # cells: ", depth[i]
            plt.bar(metric[i][0], metric[i][1] ,10, color='red', alpha=0.5)

        #for i in range(len(depth)):
        #    plt.bar(depth[i], 10, 10, color='green')
        ax2.set_ylim(bottom=0) #Max number of units per channel

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        quit()

        #***************************************** SECOND PLOTS ??? *********************

        print "inds ", inds
        print y_sorted

        y_sorted = np.array(graph13)[inds]
        plt.scatter(np.sort(x), y_sorted, color='red', alpha=1.0, s=30)
        plt.plot(np.sort(x), y_sorted, color='red', alpha=1.0)

        print "inds ", inds
        print y_sorted

        #y_sorted = np.array(graph23)[inds]
        #plt.scatter(np.sort(x), y_sorted, color='blue', alpha=1.0, s=30)
        #plt.plot(np.sort(x), y_sorted, color='blue', alpha=1.0)
        plt.show()

        #plt.yscale('log')

    #tsf.make_DCoffset_tsf()

