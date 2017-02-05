import numpy as np

class Electrode(object):
    
    def __init__(self):
        
        pass
        
    def load(self,  electrode_no):

        electrode_array=np.load('/media/cat/4TB/in_silico/ucsd_Sep20_rat_3k_20Khz/modelsdir/input/lam_electrode_rat.npy')

        e=[]
        e.append(electrode_array.T[0:30])
        e.append(electrode_array.T[30:38])
        e.append(electrode_array.T[38:48])
        e.append(electrode_array.T[48:80])
        e.append(electrode_array.T[80:170])

        #Add 5th layout - full column 1000um depth probe
        temp=[]
        for i in range(50):
            temp.append([0,-i*20,0])
        e.append(temp) 

        #Save probe layouts for later use as part of electrode objects
        self.probes=np.array(e)

        print "Electrode #: ", electrode_no
       
        alpha=1.
        electrode_colours = [0, .5, 0, alpha]
        
        #LOAD QUAD ELECTRODE SURFACES
        if True:
            n_electrodes = len(e[electrode_no])

            site_width = 12.  #WIdth of site in um
            quad_faces = np.zeros((4,3),dtype=np.float32)

            quad_faces[0] = (0, -site_width/2., -site_width/2.)
            quad_faces[1] = (0, site_width/2., -site_width/2.)
            quad_faces[2] = (0, site_width/2., site_width/2.)
            quad_faces[3] = (0, -site_width/2., site_width/2.)

            #Switch y<->z and invert z 
            temp = e[electrode_no][:,1].copy()
            e[electrode_no][:,1] = - e[electrode_no][:,2]
            e[electrode_no][:,2] = temp

            no_back_faces = 1
            if electrode_no == 4: no_back_faces = 6

            self.electrode = np.zeros(((n_electrodes+no_back_faces)*4, 3), dtype=np.float32)
            self.electrode_colors = np.zeros(((n_electrodes+no_back_faces)*4, 4), dtype=np.float32)+electrode_colours

            for i in range(n_electrodes): 
                #FULL PROBE
                self.electrode[i*4] = quad_faces[0]+e[electrode_no][i]
                self.electrode[i*4+1] = quad_faces[1]+e[electrode_no][i]
                self.electrode[i*4+2] = quad_faces[2]+e[electrode_no][i]
                self.electrode[i*4+3] = quad_faces[3]+e[electrode_no][i]

                ##ONLY HALF PROBE
                #self.electrode[i*4] = quad_faces[0]+e[electrode_no][i*2+1] 
                #self.electrode[i*4+1] = quad_faces[1]+e[electrode_no][i*2+1]
                #self.electrode[i*4+2] = quad_faces[2]+e[electrode_no][i*2+1]
                #self.electrode[i*4+3] = quad_faces[3]+e[electrode_no][i*2+1]

            #ADD BACKING SURFACES to electrodes
            if electrode_no==0:
                self.electrode[n_electrodes*4]=[-1,-1035,-11]
                self.electrode[n_electrodes*4+1]=[-1,-1035,+37]
                self.electrode[n_electrodes*4+3]=[-1,-1375,-11]
                self.electrode[n_electrodes*4+2]=[-1,-1375,+37]
                
                for i in range(4):
                    self.electrode_colors[n_electrodes*4+i]=[.1, .1, .1, alpha]

            if electrode_no==1:
                self.electrode[n_electrodes*4]=[-1,-1080,-37]
                self.electrode[n_electrodes*4+1]=[-1,-1080,+37]
                self.electrode[n_electrodes*4+3]=[-1,-1335,-10]
                self.electrode[n_electrodes*4+2]=[-1,-1335,10]

                for i in range(4):
                    self.electrode_colors[n_electrodes*4+i]=[.1, .1, .1, alpha]

            if electrode_no==2:
                self.electrode[n_electrodes*4]=[-1,-1020,-37]
                self.electrode[n_electrodes*4+1]=[-1,-1020,+97]
                self.electrode[n_electrodes*4+3]=[-1,-1385,-37]
                self.electrode[n_electrodes*4+2]=[-1,-1385,97]
                
                for i in range(4):
                    self.electrode_colors[n_electrodes*4+i]=[.1, .1, .1, alpha]      

            if electrode_no==4:

                #print e[electrode_no]

                self.electrode[n_electrodes*4]=[-55,-1025,-115]
                self.electrode[n_electrodes*4+1]=[-55,-1380,-115]
                self.electrode[n_electrodes*4+2]=[-55,-1380,-85]
                self.electrode[n_electrodes*4+3]=[-55,-1025,-85]

                self.electrode[n_electrodes*4+4]=[-55,-1025,-15]
                self.electrode[n_electrodes*4+5]=[-55,-1380,-15]
                self.electrode[n_electrodes*4+6]=[-55,-1380,+15]
                self.electrode[n_electrodes*4+7]=[-55,-1025,+15]
                
                self.electrode[n_electrodes*4+8]=[-55,-1025,115]
                self.electrode[n_electrodes*4+9]=[-55,-1380,115]
                self.electrode[n_electrodes*4+10]=[-55,-1380,85]
                self.electrode[n_electrodes*4+11]=[-55,-1025,85]
                

                self.electrode[n_electrodes*4+12]=[45,-1025,-115]
                self.electrode[n_electrodes*4+13]=[45,-1380,-115]
                self.electrode[n_electrodes*4+14]=[45,-1380,-85]
                self.electrode[n_electrodes*4+15]=[45,-1025,-85]

                self.electrode[n_electrodes*4+16]=[45,-1025,-15]
                self.electrode[n_electrodes*4+17]=[45,-1380,-15]
                self.electrode[n_electrodes*4+18]=[45,-1380,+15]
                self.electrode[n_electrodes*4+19]=[45,-1025,+15]
                
                self.electrode[n_electrodes*4+20]=[45,-1025,115]
                self.electrode[n_electrodes*4+21]=[45,-1380,115]
                self.electrode[n_electrodes*4+22]=[45,-1380,85]
                self.electrode[n_electrodes*4+23]=[45,-1025,85]                
                
                for i in range(4*no_back_faces):
                    self.electrode_colors[n_electrodes*4+i]=[.1, .1, .1, alpha]      

            if electrode_no==5:

                #print e[electrode_no]

                self.electrode[n_electrodes*4]=[-55,-1025,-115]
                self.electrode[n_electrodes*4+1]=[-55,-1380,-115]
                self.electrode[n_electrodes*4+2]=[-55,-1380,-85]
                self.electrode[n_electrodes*4+3]=[-55,-1025,-85]

                self.electrode[n_electrodes*4+4]=[-55,-1025,-15]
                self.electrode[n_electrodes*4+5]=[-55,-1380,-15]
                self.electrode[n_electrodes*4+6]=[-55,-1380,+15]
                self.electrode[n_electrodes*4+7]=[-55,-1025,+15]
                
                self.electrode[n_electrodes*4+8]=[-55,-1025,115]
                self.electrode[n_electrodes*4+9]=[-55,-1380,115]
                self.electrode[n_electrodes*4+10]=[-55,-1380,85]
                self.electrode[n_electrodes*4+11]=[-55,-1025,85]
                

                self.electrode[n_electrodes*4+12]=[45,-1025,-115]
                self.electrode[n_electrodes*4+13]=[45,-1380,-115]
                self.electrode[n_electrodes*4+14]=[45,-1380,-85]
                self.electrode[n_electrodes*4+15]=[45,-1025,-85]

                self.electrode[n_electrodes*4+16]=[45,-1025,-15]
                self.electrode[n_electrodes*4+17]=[45,-1380,-15]
                self.electrode[n_electrodes*4+18]=[45,-1380,+15]
                self.electrode[n_electrodes*4+19]=[45,-1025,+15]
                
                self.electrode[n_electrodes*4+20]=[45,-1025,115]
                self.electrode[n_electrodes*4+21]=[45,-1380,115]
                self.electrode[n_electrodes*4+22]=[45,-1380,85]
                self.electrode[n_electrodes*4+23]=[45,-1025,85]                
                
                for i in range(4*no_back_faces):
                    self.electrode_colors[n_electrodes*4+i]=[.1, .1, .1, alpha]      



        probe_name = str(electrode_no)
        return probe_name
        
