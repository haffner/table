#!/usr/bin/env python

import numpy as np
import h5py
import sys

fout=h5py.File(sys.argv[2],"w")

with h5py.File(sys.argv[1],"r") as f:
	fout["ControlFileContents"]=f["ControlFileContents"].value
	fout["SensorZenith"]=f["SensorZenith"].value          
	fout["SolarZenith"]=f["SolarZenith"].value         
	fout["TotalO3ProfileName"]=f["TotalO3ProfileName"].value       
	fout["c0"]=f["c0"].value
	fout["c1"]=f["c1"].value                         
	fout["c2"]=f["c2"].value
	fout["pressure"]=f["pressure"].value
	fout["wlen"]=f["wlen"].value         
	# slice these
	fout["f0flux"]=f["f0flux"].value[1,:]  
	fout["knr2"]=f["knr2"].value[1,:]      
	fout["lgi0"]=f["lgi0"].value[1,:]      
	fout["sb"]=f["sb"].value[1,:]          
	fout["tr"]=f["tr"].value[1,:]          
	fout["z1"]=f["z1"].value[1,:]          
	fout["z2"]=f["z2"].value[1,:]          
