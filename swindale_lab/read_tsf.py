file_name = 'ECP_1.tsf'

with open(file_name,'rb') as fin:
    header = fin.read(16)                                   #Header info; not used
    iformat = struct.unpack('i',fin.read(4))[0]             #Default: '1002'
    SampleFrequency = struct.unpack('i',fin.read(4))[0]     #Sample frequency
    n_electrodes = struct.unpack('i',fin.read(4))[0]        #No. of electrodes
    n_vd_samples = struct.unpack('i',fin.read(4))[0]        #No. of samples
    vscale_HP = struct.unpack('f',fin.read(4))[0]           #Scaling of int2 values below to save space, currently 0.1

    Siteloc1 = np.zeros((2*n_electrodes), dtype=np.int16)
    Readloc = np.zeros((n_electrodes), dtype=np.int32)
    for i in range(n_electrodes):
        Siteloc1[i*2] = struct.unpack('h', fin.read(2))[0]
        Siteloc1[i*2+1] = struct.unpack('h', fin.read(2))[0]
        Readloc[i] = struct.unpack('i', fin.read(4))[0]

    ec_traces =  np.fromfile(fin, dtype=np.int16, count=n_electrodes*n_vd_samples)
    ec_traces.shape = n_electrodes, n_vd_samples

    n_cell_spikes = struct.unpack('i',fin.read(4))[0]
    print "No. cell spikes: ", n_cell_spikes
    if (n_cell_spikes>0):
        fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
        fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=n_cell_spikes)
