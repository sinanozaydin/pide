#!/usr/bin/env python3

import numpy as np

def Ozaydin2020_Opx1(al_opx, P, cpx_opx, p_change, d_opx_ol, method):

	part = (1.393 * al_opx) + 1.947 
	
	return part
	
def Ozaydin2020_Opx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	if method == 'array':
	
		al_opx = al_opx[0]
		p = p[0]
		part = np.zeros(len(al_opx))
		
		for i in range(0, len(al_opx)):
		
			if p[i] < p_change:
			
				part[i] = (1.393 * al_opx[i]) + 1.947
			
			else:
				part[i] = (0.088 * al_opx[i]) + 0.881
				
	else:
	
		if p < p_change:
			
			part = (1.393 * al_opx) + 1.947
			
		else:
			
			part = (0.088 * al_opx) + 0.881
				
	return part
	
def Ferot2012_Opx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (0.0127 * (p**2)) - (0.2942 * p) + 2.9146
	
	return part
	
def Kovacs2012b_Opx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (-1.343 * p) + 10.048
	
	return part
	
def Sakurai2014_Opx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (3.438 * al_opx) + 0.743
	
	return part
	
def Withers2011_Opx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (-0.151 * p) + 1.542
	
	return part
	
def Ozaydin2020_Cpx1(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = d_opx_ol * ((0.077 * p) + 1.137)
	
	return part
	
def Ozaydin2020_Cpx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = d_opx_ol * 1.88
	
	return part
	
def Aubaud2004_Cpx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 1.4
	
def Demouchy2016_Cpx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 2.39
	
def Demouchy2017_Cpx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 1.9
	
def Kovacs2012_Cpx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 3.35
	
def Novella2014_Cpx2(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 1.65
	
def Tenner2009_Cpx(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	return d_opx_ol * 1.35
	
def Hirschmann2009_OpxMelt(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (0.0015 * al_opx) + 0.0082
	
	return part
	
def Hirschmann2009_CpxMelt(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (0.0015 * al_opx) + 0.0082
	
	return part
	
def Novella2014_OpxMelt(al_opx, p, cpx_opx, p_change, d_opx_ol, method):

	part = (0.0028 * al_opx) + 0.0033
	
	return part
