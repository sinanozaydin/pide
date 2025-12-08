#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Hu2024_Apatite(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E = 205000.0
	dv = 9.31e3
	sigma = 6.32
	
	cond = 10**sigma * np.exp(-(E + (dv*P)) / (R_const*T))
	
	return cond

def Dobson2000_Magnesiowustite_xfe_014(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	cond =10**(3.13-(0.4484-0.002694*P)*0.4343/(T*8.617e-5))

	return cond