from neuron import h
import numpy as np
from scipy.signal import butter, lfilter, freqz

h.load_file('stdgui.hoc')
h.load_file('interpxyz.hoc')

# units are mv and mA

dedx = 1 # dedx  gradient in v/m

def init_model():
    h.load_file('MRGaxon.hoc')          # import the MRG model. In this file the paramaters of the neuron can be altered.

    for sec in h.allsec():
        sec.insert('xtra')              # Insert 'xtra' in de sections. This couples the location of the extracellular stimulation to the potential at the nerve section

    h.define_shape()        # crestes cordinates of the nodes in 3D
    h.grindaway()           # in interpxyz.hoc, determines interpolated locations of nodes

    for sec in h.allsec():
        for seg in sec:
            #seg.xtra._ref_ex = seg._ref_e_extracellular
            h.setpointer(sec(seg.x)._ref_e_extracellular, 'ex', sec(seg.x).xtra)        # Pointer so extracellular can call correctly on xtra

    v = {}
    for sec in h.allsec():
        v[str(sec)] = h.Vector().record(sec(0.5)._ref_v)                                # Setup recording vector
        # es[sec] = h.Vector().record(sec(0.5)._ref_e_extracellular)
    tvec = h.Vector().record(h._ref_t)                                                  # Setup time vector

    print("Init_model complete")

    return v, tvec

def calcesI(x1, y1, amp1, x2, y2, amp2, sigma_e = 2.76e-07):                                # x,y are electrode cordinates. amp = stimulation amplitude as found in SD curve
    for sec in h.allsec():
        for seg in sec:
            r1 = np.sqrt((x1 - seg.x_xtra)**2 + y1**2)                                     # calculate the distance from electrode to the segment
            seg.es1_xtra = (1e-3*amp1)/(4*np.pi*sigma_e*r1)                                      #calculate the extracellular portential for unit current
            r2 = np.sqrt((x2 - seg.x_xtra)**2 + y2**2)                                     # calculate the distance from electrode to the segment
            seg.es2_xtra = (1e-3*amp2)/(4*np.pi*sigma_e*r2)
    
    print("calcesI complete")

def Stimulation(v, Dur):
   #Dur = 0.5 # duration
   dt = Dur/1e3
   Del = 1 # delay until start stimulation

   dat = {}
   t0 = 0

   tstop = max([5, Del+Dur+4])                  # make sure stumulation is at least 5 seconds
   h.dt = dt
   
   #if Dur < 1:
   #     h.dt = dt
   #else:
   #     h.dt = dt*2
    
    # Create stimulation wave
   i = int((tstop-t0)/dt)
   j = int(Del/dt)
   unit_wave = np.zeros(i)
   while (j)<((((Dur/dt)/2)+Del/dt)+1):                   # create a unit squire wave
       unit_wave[j] = 1
       j = j+1
   while (j)<((Dur/dt+Del/dt)+1):
       unit_wave[j] = 2
       j = j+1

   n = int((tstop-t0)/dt)
   t = np.linspace(t0, tstop, n) # Vector with time stamp
   #wave = np.multiply(unit_wave, amp)   #Multiply unit wave with aplitude so wave has correct aplitude
   #wave = butter_lowpass_filter(wave, 10, dt)

   setStim(unit_wave, t, dt, tstop)
   h.run()
       
   return dat
   
  

def butter_lowpass(cutoff, dt, order=1):
    nyq = 0.5 / dt
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, dt, order=1):
    b, a = butter_lowpass(cutoff, dt, order=order)
    y = lfilter(b, a, data)
    return y

def setStim(data, t, dt, tstop): #Dur, Del, amp, w, dt, t0, tstop, tau, sigma, biphasic, ipd, ratio):
    if 'stim_amp' not in globals():
        global stim_amp 
        stim_amp = h.Vector() 
    if 'stim_time' not in globals():
        global stim_time
        stim_time = h.Vector()
    # data, t = waveform(Dur, Del, amp, w, dt=dt, t0=t0, tstop=tstop, tau=tau, sigma=sigma, biphasic=biphasic, ipd=ipd, ratio=ratio)
    h.dt = dt
    h.tstop = tstop
    stim_amp.resize(np.size(data))
    stim_amp.fill(0)
    stim_time.resize(np.size(t))
    stim_time.fill(0)
    for i in range(np.size(data)):
        stim_amp[i] = data[i]
        stim_time[i] = t[i]
    attach_stim()
    
    # return data, t

def attach_stim():
    # # now drive h.is_xtra
    stim_amp.play(h._ref_stim_xtra, stim_time, 1)


def findExLoc(rec): 
    t = []
    loc = []
    v = {}
    for vec in rec:
        v[vec] = rec[vec].to_python()
        if np.where(np.diff(np.sign(v[vec])))[0].size != 0:
            t.append(np.where(np.diff(np.sign(v[vec])))[0][0])
    t_ex = min(t)
    for vec in rec:
        if np.where(np.diff(np.sign(v[vec])))[0].size != 0:
            if np.where(np.diff(np.sign(v[vec])))[0][0] == t_ex:
                loc.append(vec)
    
    return v, t, loc, t_ex
  

