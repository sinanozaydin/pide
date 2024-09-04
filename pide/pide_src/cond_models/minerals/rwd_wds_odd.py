#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Yoshino2012_DryRingwoodite_Xfe(T, P, water, xFe ,param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A = 1885.0
	E = 193000.0 #J/mol
	alpha = 146000.0 #J/mol
	dv = -0.33e-6 #m3/mol
	beta = 0.72e-6 #m3/mol
	P = P * 1e9 #converting to Pa - J/m3
	
	cond = A * xFe * np.exp(-((E - (alpha * xFe**(1.0/3.0))) + (P*(dv - beta*xFe))) / (R_const*T))
	
	if mechanism == 'proton':
	
		raise ValueError('Proton conduction is not included in electrical conductivity model: Yoshino2012_DryRingwoodite_Xfe')
		
	else:
			
		return cond
	
def Yoshino2009_DryRingwoodite_Xfe(T, P, water, xFe ,param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A = 10042.0
	E = 208000.0
	alpha = 155000.0
	
	cond = A * xFe * np.exp(-((E - (alpha * xFe**(1.0/3.0)))) / (R_const*T))
	
	if mechanism == 'proton':
	
		raise ValueError('Proton and ionic conduction is not included in electrical conductivity model: Yoshino2009_DryRingwoodite_Xfe')
		
	else:
			
		return cond
	
def Yoshino2012b_DryWad(T, P, water, xFe ,param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A = 5.51 #S/m
	E = 94000.0 #j/mol
	dv = -0.62e-6 #converting go m3/mol
	P = P * 1e9 #converting to Pa / J/m3
	
	cond = A * np.exp(-(E + (P*dv)) / (R_const*T))
	
	if mechanism == 'proton':
	
		raise ValueError('Proton and ionic conduction is not included in electrical conductivity model: Yoshino2012b_DryWad')
		
	else:
			
		return cond
	