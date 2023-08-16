#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Shen2021(T, P, water, method):

	E = 52210.0
	dv = 0.33e3
	sigma = 10.0**1.93
	
	cond = sigma * np.exp(-(E + (P * dv)) / (R_const * T))

	return cond
