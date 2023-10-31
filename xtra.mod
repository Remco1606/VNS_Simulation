: $Id: xtra.mod,v 1.4 2014/08/18 23:15:25 ted Exp ted $
: 2018/05/20 Modified by Aman Aberra 

NEURON {
	SUFFIX xtra
	RANGE es1, es2 : (es = max amplitude of the potential)		
	RANGE x, y, z, type, order
	GLOBAL stim : (stim = active electrode, amp = stimulation amplitude)
	POINTER ex 
}

PARAMETER {	
	es1 = 0 (mV)
	es2 = 0 (mV)
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)		
	type = 0 (1) : numbering system for morphological category of section - unassigned is 0
	order = 0 (1) : order of branch/collateral. 
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	stim (unitless) 		
	area (micron2)
}

INITIAL {
	if (stim == 1){
		ex = 1*es2
	}
	else if (stim == 2){
		ex = 1*es2
	} else{
		ex = 0
	}	
}


BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  if (stim == 1){
		ex = 1*es2
	}
	else if (stim == 2){
		ex = 1*es2
	} else{
		ex = 0
	}
}

