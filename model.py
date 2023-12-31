from neuron import h
import numpy as np
from scipy.signal import butter, lfilter, freqz

h.load_file('stdgui.hoc')
h.load_file('interpxyz.hoc')

# units are mv and uA

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

def calcesI(x, y, sigma_e = 2.76e-07, sigma_t = (1/1200), sigma_l = (1/300)):           # sigma[ohm/cm]
    for sec in h.allsec():
        for seg in sec:
            #r = np.sqrt((x - seg.x_xtra)**2 + y**2)                                     # calculate the distance from electrode to the segment
            #seg.es_xtra = 1e-3/(4*np.pi*sigma_e*r)                                      #calculate the extracellular portential for unit current
            # The unit of current is mA so mV is in mV. However uA is used for the amplitude becuase otherwise the input values are to low therefore it is multiplied with 1e-3 so Ve is in mV
            seg.es_xtra = 1e-3/(4*np.pi*(np.sqrt((sigma_l*sigma_t)*(y/10000)**2 + sigma_t**2*((x - seg.x_xtra)/10000)**2)))       # Calculation of Ve according to Grill 1999 es[mV] I[mA]
    
    print("calcesI complete")

def Stimulation(v, Dur, amp = -1, depth = -1,loc_find = 1): # Dur[ms] amp[uA]
   #Dur = 0.5 # duration
   dt = Dur/1e3
   Del = 1 # delay until start stimulation

   dat = {}
   t0 = 0

   tstop = max([5, Del+Dur+4]) # make sure stumulation is at least 5 seconds
   h.dt = dt
   
   #if Dur < 1:
   #     h.dt = dt
   #else:
   #     h.dt = dt*2
    
    # Create stimulation wave
   i = int((tstop-t0)/dt)
   j = int(Del/dt)
   unit_wave = np.zeros(i)
   while (j)<((Dur/dt+Del/dt)+1):   # create a unit squire wave
       unit_wave[j] = 1
       j = j+1

   n = int((tstop-t0)/dt)
   t = np.linspace(t0, tstop, n) # Vector with time stamp
   #wave = np.multiply(unit_wave, amp)   #Multiply unit wave with aplitude so wave has correct aplitude
   #wave = butter_lowpass_filter(wave, 10, dt)

###########################################################################################################################
   excite = 0
   ex = 0 #last amplitude causing excitation
   nex = 0 #last amplitude failed to cause excitation

    # for i in np.logspace(1, depth, (2-depth)):
    # print(i)
   while abs(amp)<500000:                     # While loop that searches for lowest stimulation amplitude for given pulse length
       print('amp = ' + str(amp))
       wave = np.multiply(unit_wave, amp)   #Multiply unit wave with aplitude so wave has correct aplitude
       setStim(wave, t, dt, tstop)          # define stimulation
       h.run()                              
       recdat = np.array(v['node[10]'])     #put results of node 10 in array.
       if np.max(recdat) > 0:               # Check if a action potential has been generated at node 10
           ex = amp
           amp = (nex + ex)/2
       elif np.max(recdat) <= 0:
           if ex == 0:
               nex = amp
               amp = amp*2
           else:
               nex = amp
               amp = (nex + ex)/2

       if np.abs(ex - nex) < (10**depth) and ex != 0:
            # print(amp)
            # print(np.max(recdat))
           print('excitation =' + str(ex))
           excite = 1
           break  
        

   pwr = np.multiply(wave, wave)
 
   if excite == 0:
       print('NO th found for wave' + ' at duration ' + str(Dur))
   else:
      print('th found for wave' + ' at duration ' + str(Dur))
    # plt.plot(t_vec, v)

   dat['wave'] = unit_wave*ex
   dat['pwr'] = pwr
   dat['t'] = t
   dat['dur'] = Dur
   dat['th'] = ex
   dat['ex'] = excite
   dat['dt'] = h.dt
 
   
   if loc_find == 1 and excite == 1:                # Create data set for all nodes for the given stimualtion amplitude that gives activation
       _ = setStim(unit_wave*ex, t, dt, tstop)
       h.run()
       v, t_ex_lst, loc, t_ex = findExLoc(v)
       dat['ex_loc'] = loc
       dat['t_onset'] = t_ex - Del
       dat['vm_ex'] = {}
       for l in loc:
           dat['vm_ex'][l] = v[l]
       dat['t_ex'] = t_ex_lst
    
   return dat

def SD_Curve(v):
    data = {}
    SD_Data = {}
    SD_Threshold = []
    SD_Duration = []
    i = 0.1
    m = 0
    Amp = -1
    while i < 2:
        data = Stimulation(v, Dur=i, amp = Amp)
        SD_Threshold.append(data['th'])
        SD_Duration.append(i)
        i = i+0.1
        m = m+1
        if data['th'] != 0:
            Amp = round(data['th'])
        
    SD_Data['SD_Threshold'] = SD_Threshold
    SD_Data['SD_Duration'] = SD_Duration
    return SD_Data


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
    stim_amp.play(h._ref_stim_xtra, stim_time, 1)  #run stim_apm into "stim" in xtra.mod with time vector
  

