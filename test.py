import numpy as np

def testje():

   Dur = 0.5 # duration
   dt = Dur/1e3
   Del = 1 # delay until start stimulation

   dat = {}
   t0 = 0

   tstop = max([5, Del+Dur+4])                  # make sure stumulation is at least 5 seconds
   
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
   return unit_wave
