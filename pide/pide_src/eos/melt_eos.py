#!/usr/bin/env python3

import numpy as np

def Howell_Green_Powell_2018_ds633_MeltEOS(T,P,sio2,al2o3,mgo,feo,cao,na2o,k2o,tio2,mno,p2o5,cr2o3,h2o):

	from burnman.minerals.HGP_2018_ds633 import make_melt_class, q4L, sl1L, wo1L, fo2L, fa2L, jdL, hmL, ekL, tiL, kjL, ctL, h2o1L
	
	#defining the thermodynamic environment
	endmembers = [q4L, sl1L, wo1L, fo2L, fa2L, jdL, hmL, ekL, tiL, kjL, ctL, h2o1L]
	melt_class = make_melt_class(endmembers)
	
	vp_melt = np.zeros(len(T))
	density_melt = np.zeros(len(T))
	bulk_modulus_melt = np.zeros(len(T))
	
	for i in range(T):
	
		melt = MyMeltClass([sio2[i],al2o3[i],mgo[i],feo[i],cao[i],
		na2o[i],k2o[i],tio2[i],mno[i],p2o5[i],cr2o3[i],h2o[i]])
		
		melt.set_state(P[i]*1e9,T[i]) #Pa and K
		
		vp_melt[i] = melt.vp
		density_melt[i] = melt.density
		bulk_modulus_melt[i] = melt.density * (melt.vp**2.0)
		
	return density_melt,vp_melt,bulk_modulus_melt

