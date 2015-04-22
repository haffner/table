#!/usr/bin/env python

import numpy as np
import h5py
import sys

d=h5py.File(sys.argv[2],"w")

with h5py.File(sys.argv[1],"r") as s:
	# copy
	s.copy("ControlFileContents", d)
	s.copy("RingAdded", d)
	s.copy("SensorZenith", d)
	s.copy("SolarZenith", d)          
	s.copy("TotalO3ProfileName", d)      
	s.copy("c0", d) 
	s.copy("c1", d)
	s.copy("c2", d)                       
	s.copy("pressure", d) 
	s.copy("wlen", d) 
 
	# slice
	d["f0flux"]=s["f0flux"].value[1,:]  
	d["knr2"]  =s["knr2"].value[1,:]      
	d["lgi0"]  =s["lgi0"].value[1,:]      
	d["sb"]    =s["sb"].value[1,:]          
	d["tr"]    =s["tr"].value[1,:]          
	d["z1"]    =s["z1"].value[1,:]          
	d["z2"]    =s["z2"].value[1,:]          
