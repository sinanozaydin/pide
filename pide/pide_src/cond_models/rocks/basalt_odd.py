#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2015_HydrousBasalt(T,P,water, param1, fo2 = None, fo2_ref = None, method = None):

	A1 = 1.49
	A2 = 1.81
	E = 114e3
	r = 1.65
	alpha = 0.26
	beta = 0.2

	cond = (A1 + (A2 * water**r)) * np.exp(-(E + (alpha*(water**beta)))/ (R_const*T))

	return cond

def Hui2015_QuartzAndesiteGrainBoundary(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	sigma = 10.0**2.36
	E = 90e3
	dv = 0.56e3

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond


def Hui2015_QuartzAndesiteGrainInterior(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	sigma = 10.0**1.63
	E = 76e3
	dv = 4.96e3

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

def Hui2017_PyroxeneAndesite(T, P, water, param1, fo2 = None, fo2_ref = None, method = None):

	sigma = 10.0
	E = 48000.0
	dv = -6.75e3

	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond

