#!/usr/bin/env python3

import numpy as np
R_const = 8.3144621

def Kohlstedt1996_RwdSol(T,P,depth,h2o_fug,o2_fug,fe_rwd_wds,method):

	#ppm value is calculated using the Table 2 in Bolfan-Casanova et al. (2000, EPSL)

	return 28777.0 * np.ones(len(T))
	
def Kohlstedt1996_WdsSol(T,P,depth,h2o_fug,o2_fug,fe_rwd_wds,method):

	#ppm value is calculated using the Table 2 in Bolfan-Casanova et al. (2000, EPSL)

	return 25580.0 * np.ones(len(T))
	
def Inoue1994_WdsSol(T,P,depth,h2o_fug,o2_fug,fe_rwd_wds,method):

	#ppm value is calculated using the Table 2 in Bolfan-Casanova et al. (2000, EPSL)
	
	return 31975.0 * np.ones(len(T))
	

