################################################################################
# program: psd_gui.py
# author: Tom Irvine
# version: 1.9
# date: September 19, 2013
# description:  
#    
#  Determine the power spectral density of a signal.
#  The file must have two columns: time(sec) & amplitude.
#              
################################################################################

from __future__ import print_function

import sys

if sys.version_info[0] == 2:
    print ("Python 2.x")
    import Tkinter as tk
    from tkFileDialog import asksaveasfilename
    import tkMessageBox
           
if sys.version_info[0] == 3:
    print ("Python 3.x")    
    import tkinter as tk 
    from tkinter.filedialog import asksaveasfilename 
    import tkinter.messagebox as tkMessageBox    
    

from tompy import read_two_columns_from_dialog,sample_rate_check



from math import sqrt,log
from scipy.fftpack import fft

import numpy as np

import matplotlib.pyplot as plt

################################################################################

class tk_PSD:

    def __init__(self,parent): 
        self.master=parent        # store the parent
        top = tk.Frame(parent)    # frame for all class widgets
        top.pack(side='top')      # pack frame in parent's window
        
        self.sr=0
        self.dt=0
        self.b=[]
        self.a=[]
        self.num=0
        
        self.ss=np.zeros(100,'int')
        self.seg=np.zeros(100,'int')
        self.i_seg=np.zeros(100,'int')
        self.ddf=np.zeros(100,'float') 
        
        self.r_ss=[]
        self.r_seg=[]
        self.r_i_seg=[]
        self.r_ddf=[]         
        
        self.r_ss=[]
        self.r_seg=[]
        self.r_i_seg=[]
        self.r_ddf=[]
        self.rms=0
        self.freq=[]
        self.full=[]
        self.mH=0
        self.NW=0
        self.mmm=0

################################################################################

        crow=1

        self.hwtext1=tk.Label(top,text='Power Spectral Density Function')
        self.hwtext1.grid(row=crow, column=0, columnspan=6, pady=8,sticky=tk.W)

################################################################################

        crow=crow+1

        self.hwtext2=tk.Label(top,\
            text='The input file must have two columns:  time(sec) & amplitude')
        self.hwtext2.grid(row=crow, column=0, columnspan=6, pady=8,sticky=tk.W)

################################################################################

        crow=crow+1

        self.button_read = tk.Button(top, text="Read Input File",command=self.read_data)
        self.button_read.config( height = 2, width = 15 )
        self.button_read.grid(row=crow, column=0,columnspan=1, pady=9,sticky=tk.W)  

################################################################################

        crow=crow+1

        self.hwtextadv=tk.Label(top,text='Select Analysis Option')
        self.hwtextadv.grid(row=crow, column=0, pady=10)
        
        
        self.hwtextat=tk.Label(top,text='Select Amplitude Type')
        self.hwtextat.grid(row=crow, column=1, padx=5, pady=10)
        
        self.hwtextat=tk.Label(top,text='Enter Time History Amplitude Unit')
        self.hwtextat.grid(row=crow, column=2, padx=5, pady=10)

################################################################################

        crow=crow+1
        
        myframe=tk.Frame(top)
        myframe.grid(row=crow, column=0,padx=3)
        scrollbar = tk.Scrollbar(myframe) 
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.Lba = tk.Listbox(myframe, width=35, yscrollcommand=scrollbar.set) 
        self.Lba.pack()
        scrollbar.config(command=self.Lba.yview)
        
        self.Lbat = tk.Listbox(top,exportselection=0)
        self.Lbat.insert(1, "Acceleration")
        self.Lbat.insert(2, "Velocity")
        self.Lbat.insert(3, "Displacement")
        self.Lbat.insert(4, "Force")        
        self.Lbat.insert(5, "Pressure")
        self.Lbat.insert(6, "Stress")          
        self.Lbat.grid(row=crow, column=1, padx=5)
        self.Lbat.select_set(0) 
             
        self.Lbat.bind('<<ListboxSelect>>', self.unit_change)
        
        self.iur=tk.StringVar()       
        self.iur.set('G')  
        self.iu_entry=tk.Entry(top, width = 8,textvariable=self.iur)
        self.iu_entry.grid(row=crow, column=2,padx=5, pady=1,sticky=tk.N)        

################################################################################

        crow=crow+1

        self.hwtext4=tk.Label(top,text='Mean Removal')
        self.hwtext4.grid(row=crow, column=0, columnspan=1, pady=10)

        self.hwtext5=tk.Label(top,text='Window')
        self.hwtext5.grid(row=crow, column=1, columnspan=1, pady=10)

################################################################################

        crow=crow+1

        self.Lb_mr = tk.Listbox(top,height=2,exportselection=0)
        self.Lb_mr.insert(1, "Yes")
        self.Lb_mr.insert(2, "No")
        self.Lb_mr.grid(row=crow, column=0, pady=4)
        self.Lb_mr.select_set(0) 

        self.Lb_win = tk.Listbox(top,height=2,exportselection=0)
        self.Lb_win.insert(1, "Rectangular")
        self.Lb_win.insert(2, "Hanning")
        self.Lb_win.grid(row=crow, column=1, pady=4)
        self.Lb_win.select_set(0) 
        
        crow=crow+1        
        
        self.hwtextf1=tk.Label(top,text='Min Freq (Hz)')
        self.hwtextf1.grid(row=crow, column=0,padx=5, pady=8)
        self.hwtextf1.config(state='disabled')

        self.hwtextf2=tk.Label(top,text='Max Freq (Hz)')
        self.hwtextf2.grid(row=crow, column=1,padx=5, pady=8)
        self.hwtextf2.config(state='disabled')

################################################################################

        crow=crow+1    
        
        self.f1r=tk.StringVar()  
        self.f1r.set('')  
        self.f1_entry=tk.Entry(top, width = 9,textvariable=self.f1r)
        self.f1_entry.grid(row=crow, column=0,padx=5, pady=1)
        self.f1_entry.config(state='disabled')

        self.f2r=tk.StringVar()  
        self.f2r.set('')  
        self.f2_entry=tk.Entry(top, width = 9,textvariable=self.f2r)
        self.f2_entry.grid(row=crow, column=1,padx=5, pady=1)     
        self.f2_entry.config(state='disabled')
        
        crow=crow+1  

        self.button_calculate = tk.Button(top, text="Calculate PSD", command=self.psd_calculation)
        self.button_calculate.config( height = 2, width = 15,state = 'disabled')
        self.button_calculate.grid(row=crow, column=0, pady=20) 
        
        self.button_save = tk.Button(top, text="Save PSD", command=self.export_psd)
        self.button_save.config( height = 2, width = 15,state = 'disabled')
        self.button_save.grid(row=crow, column=1, pady=20) 
        
        self.button_quit=tk.Button(top, text="Quit", command=lambda root=root:quit(root))
        self.button_quit.config( height = 2, width = 15 )
        self.button_quit.grid(row=crow, column=2, padx=10,pady=20)
        
        
################################################################################

    def export_psd(self):
        output_file_path = asksaveasfilename(parent=root,title="Enter the PSD filename")           
        output_file = output_file_path.rstrip('\n')    
        WriteData2(len(self.freq),self.freq,self.full,output_file)
    
    def unit_change(self,val):
        sender=val.widget
        nat= int(sender.curselection()[0])
        
#        print 'nat=%d' %nat

        if(nat==0):
           self.iur.set('G') 
        else:
           self.iur.set('') 
           
    @classmethod     
    def Hanning_initial(cls,mmm):
        H=np.zeros(mmm,'f')
        tpi=2*np.pi    
        alpha=np.linspace(0,tpi,mmm)
        ae=sqrt(8./3.)
        H=ae*0.5*(1.-np.cos(alpha))                
        return H               

    def psd_calculation(self):
        n= int(self.Lba.curselection()[0])
        
        mr_choice=int(self.Lb_mr.curselection()[0])
        h_choice =int(self.Lb_win.curselection()[0])        
        

        self.NW=self.r_i_seg[n]
        self.mmm = 2**int(log(float(self.num)/float(self.NW))/log(2))
        
        print (n)
        print (self.r_ddf[n])
        
        self.df=1./(self.mmm*self.dt)
        self.mH=((self.mmm/2)-1)
        
        if h_choice==1:
            H=tk_PSD.Hanning_initial(self.mmm)

        print (" ")
        print ("     number of segments   NW= %d " %self.NW)
        print ("       samples/segments  mmm= %d " %self.mmm)
        print (" half samples/segment-1   mH=%d  " %self.mH)
        print (" ")
        print ("        df=%6.3f Hz" %self.df)

        self.full=np.zeros(self.mH,'f')
        mag_seg=np.zeros(self.mH,'f') 

        amp_seg=np.zeros(self.mmm,'f')

        nov=0

        for ijk in range (1,2*self.NW):
            
        
            amp_seg[0:self.mmm]=self.b[(0+nov):(self.mmm+nov)]    

            nov=nov+int(self.mmm/2)
    
            if mr_choice==0 or h_choice==1:
                amp_seg-=np.mean(amp_seg)
        
            if h_choice==1:
                amp_seg*=H

            Y = fft(amp_seg)
 
            mag_seg = tk_PSD.magnitude_resolve(self.mmm,self.mH,Y)    
   
            self.full+=mag_seg
    
      
        den=self.df*(2*self.NW-1)

        self.full/=den

        ms=sum(self.full)

        self.freq=np.zeros(self.mH,'f')

        maxf=(self.mH-1)*self.df

        self.freq=np.linspace(0,maxf,self.mH)
    
        
    
        tempf=self.freq[0:self.mH-1]    
        tempa=self.full[0:self.mH-1]
        self.freq=tempf
        self.full=tempa      
    
        self.rms=sqrt(ms*self.df)
        

      
#        tk_PSD.psd_plots
  
#    @classmethod        
#    def psd_plots(cls,self):

        print ('plotssss')

        f1=float(self.f1r.get())
        f2=float(self.f2r.get()) 

        ys=1.0e+20
        
        for i in range(0,len(self.full)):
            if(self.freq[i]>=f1 and self.freq[i]<=f2):
                if(ys>self.full[i]):
                    ys=self.full[i]
            
        
        ymax=10**np.ceil(np.log10(max(self.full)))
        ymin=10**np.floor(np.log10(ys))

        if(ymin<(ymax/10000)):
            ymin=ymax/10000


        nat= int(self.Lbat.curselection()[0])
        
        str=self.iur.get()
        
        if(nat==0):
            out1='Accel (' +str+ '^2/Hz)'
        if(nat==1):
            out1='Velocity (' +str+ '^2/Hz)'
        if(nat==2):
            out1='Disp (' +str+ '^2/Hz)'
        if(nat==3):
            out1='Force (' +str+ '^2/Hz)'    
        if(nat==4):
            out1='Pressure (' +str+ '^2/Hz)'
        if(nat==5):
            out1='Stress (' +str+ '^2/Hz)'
          

        print (" ")
        print (" view plots ")   
        
        plt.gca().set_autoscale_on(False)
        
        out2="%6.3g" %self.rms

        plt.close(2)
        plt.figure(2)     
        plt.plot(self.freq,self.full)
        title_string='Power Spectral Density   '+out2+' GRMS Overall '
        plt.title(title_string)
        plt.xlim([f1,f2])
        plt.ylim([ymin,ymax])
        plt.ylabel(out1)     
        plt.xlabel(' Frequency (Hz) ')
        plt.grid(True)
        plt.savefig('power_spectral_density')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()       
        
        self.button_save.config(state = 'normal')
       
    @classmethod           
    def magnitude_resolve(cls,mmm,mH,Y):
#
       mHm1=mH-1
       z=np.zeros(mH,'f')
       mag_seg=np.zeros(mH,'f')
#      
       z=abs(Y)/float(mmm)
#
       mag_seg[0]=z[0]**2
#
       mag_seg[1:mHm1]=((2*z[1:mHm1])**2)/2  
#
       return mag_seg          


    def read_data(self):            
            
        self.a,self.b,self.num=read_two_columns_from_dialog('Select Input File')
        
        dur=self.a[self.num-1]-self.a[0]
        self.dt=dur/float(self.num)
        
        self.sr=1./self.dt
        
        self.sr,self.dt=sample_rate_check(self.a,self.b,self.num,self.sr,self.dt)
        
        plt.ion()
        plt.clf()
        plt.figure(1)

        plt.plot(self.a, self.b, linewidth=1.0,color='b')        # disregard error
       
        plt.grid(True)
        plt.xlabel('Time (sec)')
##        plt.ylabel(self.y_string.get())  
        plt.title('Input Time History')
    
        nat= int(self.Lbat.curselection()[0])
        
        str=self.iur.get()
        
        if(nat==0):
            out1='Accel (' +str+ ')'
        if(nat==1):
            out1='Velocity (' +str+ ')'
        if(nat==2):
            out1='Disp (' +str+ ')'
        if(nat==3):
            out1='Force (' +str+ ')'    
        if(nat==4):
            out1='Pressure (' +str+ ')'
        if(nat==5):
            out1='Stress (' +str+ ')'
        
        plt.ylabel(out1)       
         
        
        plt.draw()

        print ("\n samples = %d " % self.num)
        
        self.button_calculate.config(state = 'normal')     
        
        self.advise(self)
        
        self.hwtextf1.config(state='normal')
        self.hwtextf2.config(state='normal')
        self.f1_entry.config(state='normal')
        self.f2_entry.config(state='normal')
        
        a=10**np.ceil(np.log10(1./dur))        
        
        out1='%8.4g' %(a)
        out2='%8.4g' %(self.sr/2.)        
        
        self.f1r.set(out1)  
        self.f2r.set(out2)               
    
        
    @classmethod      
    def advise(cls,self):
        
        n=self.num
        
        NC=0
        
        for i in range(1,1000):
    
            nmp = 2**(i-1)
   
            if(nmp <= n ):
                self.ss[i] = 2**(i-1)
                self.seg[i] = np.round(n/self.ss[i])
                self.i_seg[i] = np.floor(self.seg[i])
                tseg=self.dt*self.ss[i]
                self.ddf[i]=1./tseg                
                NC=NC+1
            else:
                break


        t_ss= self.ss[NC+1:0:-1]
        t_seg= self.seg[NC+1:0:-1]
        t_i_seg= self.i_seg[NC+1:0:-1]        
        t_ddf= self.ddf[NC+1:0:-1]          
        
        k=0
        
        for i in range(0,NC):
            if( t_seg[i]>0 ):
                out1='sps=%d,  df=%6.3g Hz,  sdof=%d' \
                                       %(t_ss[i],t_ddf[i],2*t_i_seg[i])
                self.Lba.insert(i, out1)

                self.r_ss.append(t_ss[i])
                self.r_seg.append(t_ddf[i])
                self.r_i_seg.append(t_i_seg[i])        
                self.r_ddf.append(t_ddf[i]) 

                k=k+1
                    
            if(i==12):
                break

        self.Lba.select_set(0) 
        
################################################################################

def WriteData2(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    for i in range (0, nn):
        outfile.write(' %10.6e \t %8.4e \n' %  (aa[i],bb[i]))
    outfile.close()
root = tk.Tk()

################################################################################

root.minsize(500,500)
root.geometry("650x620")

root.title("psd_gui.py ver 1.9  by Tom Irvine") 

tk_PSD(root)
root.mainloop()  