#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2014_DryandWetOlivine_param1_fo2_param2_fo2ref(T, P, water, param1, param2, fo2, fo2ref, method):

	dv_dai2014 = -0.86 * 1e3 # cm^3 /mol Taken from Dai2014b-PEPI, Error is insignificant, since the effect itself is insignificant...
	p_ref = 4.0
	cw_ref = 460.0 / 1e4
	e1_dai = 74000.0
	e1_di_err = 3000
	e2_dai = 115000
	r_dai = 0.8
	q_dai = -0.066 #Taken from Dai2014c-P"EPI, Error is insignificant.
	q_dry = 0.16666 #Taken from Constable (2006)
	h2o = water / (1e4)
	
	#Changing fo2 to complex to avoid Runtime errors
	fo2 = np.array(fo2, dtype = np.complex)
	fo2_ref = np.array(fo2ref, dtype = np.complex)

	dai_ref = ((10.0**0.48) * np.exp(-((e1_dai) + (p_ref*dv_dai2014)) / (R_const * T))) + (10.0**2.84 * np.exp(-((e2_dai) + (p_ref*dv_dai2014)) / (R_const * T)))
	cond_wet = (dai_ref * (h2o / (cw_ref))**(r_dai) * (fo2/fo2_ref)**(q_dai)) *  np.exp(- ((P - p_ref) * dv_dai2014) / (R_const*T))
	dai_dry =  10**2.4 * ((fo2/fo2_ref)**(q_dry)) * np.exp(-(154000.0 + (P * dv_dai2014)) / (R_const * T))
	
	cond = np.real(dai_dry + cond_wet)

	return cond
	
def Dai2020_WetOlivine_200ppmTi(T,P, water, param1, param2, fo2, fo2ref, method):

	dv_dai2014 = -0.86 * 1e3 # m^3 /mol Taken from Dai2014b-PEPI, Error is insignificant, since the effect itself is insignificant...
	q_dai = -0.066 #Taken from Dai2014c-P"EPI, Error is insignificant.
	q_dry = 0.16666 #Taken from Constable (2006)

	h2o = water / (1e4)

	A_dai_200 = 3.01

	r_dai_200 = 0.51

	E_dai_200 = 87000.0

	cond_wet = (10.0 ** A_dai_200) * ((h2o)**r_dai_200) * np.exp(-(E_dai_200) / (R_const * T))

	dai_dry =  10**2.4 * ((fo2/fo2ref)**(q_dry)) * np.exp(-(154000.0 + (P * dv_dai2014)) / (R_const * T))

	cond = dai_dry + cond_wet

	return cond
	
def Poe2010_DryandWetOlivine(T, P, water, param1, param2, fo2, fo2ref, method):

	E1 = [146e3,126e3]
	E2 = [112e3,150e3]
	E3 = [129e3,81.2e3]
	sigma_1 = [2.52,2.59]
	sigma_2 = [1.139,3.46]
	sigma_3 = [1.99,1.02]
	alpha_1 = 1180
	alpha_2 = 1430
	alpha_3 = 700
	
	water = water / 1e4
	
	cond_1 = (sigma_1[0] * water * np.exp(-(E1[0] + (alpha_1 * water))) / (R_const*T)) + (sigma_1·[1] * water * np.exp(-(E1[1] + (alpha_1 * water)) / (R_const*T)))
	cond_2 = (sigma_2[0] * water * np.exp(-(E2[0] + (alpha_2 * water))) / (R_const*T)) + (sigma_1·[1] * water * np.exp(-(E1[1] + (alpha_1 * water)) / (R_const*T)))
	cond_3 = (sigma_3[0] * water * np.exp(-(E3[0] + (alpha_3 * water))) / (R_const*T)) + (sigma_1·[1] * water * np.exp(-(E1[1] + (alpha_1 * water)) / (R_const*T)))
		
	cond = (cond_1*cond_2*cond_3)**(1.0/3.0)
	
	return cond
	
def Constable2006_dryOlivine_fo2(T, P, water, param1, param2, fo2, fo2ref, method):

	e = 1.602e-19 #charge of the electron in coulombs
	k = 8.617e-5
	kT = T * k
	bfe = (5.06e24) * np.exp(-0.357 / kT)
	bmg = (4.58e26) * np.exp(-0.752 / kT)
	ufe = (12.2e-6) * np.exp(-1.05 / kT)
	umg = (2.72e-6) * np.exp(-1.09 / kT)
	condfe = bfe + (3.33e24 * np.exp(-0.02/kT) * (fo2*1e5)**(0.16666))
	condmg = bmg + (6.21e30 * np.exp(-1.83/kT) * (fo2*1e5)**(0.16666))

	cond = (condfe * ufe * e) + (2.0*condmg* umg * e)

	return cond

def Dai2014_DryOlivine_param1_xfe(T, P, water, param1, param2, fo2, fo2ref, method):

	sigma_dai2014c = 10**(2.77 + (-1.19 * param1))
	e_dai2014c = 162000.0 + (-63000 * param1)

	cond = sigma_dai2014c * np.exp(-(e_dai2014c) / (R_const * T))

	return cond

def Fullea2011_DryOlivine_param1_xfe(T, P, water, param1, param2, fo2, fo2ref, method):

	sigma_i_fullea = 10**4.73
	sigma_pol_fullea = 10**2.7 #average value of used models Fullea et al. (2011)

	e_i_fullea = 231000.0

	a_ful = [1.642,0.246,-4.85,3.259] #polynomial coefficients that calculates the fe-dependency of act. enthalpies from Omura et al (1989)
	dv_ful = 0.68 * 1e-6
	fe_pol_fullea = (a_ful[0] + (a_ful[1] * param1) + (a_ful[2] * (param1**2.0)) + (a_ful[3]* (param1**3.0)) + (P * dv_ful)) * 1e5

	cond = (sigma_i_fullea * np.exp(-e_i_fullea / (R_const * T))) + (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))

	return cond
	
def Pommier2018_ShearedDryOlivine(T, P, water, param1, param2, fo2, fo2ref, method):

	A = [284.0,780.0,261.0]

	E = [126400.0,122700.0,114500.0]

	sigma = np.zeros((3,len(T)))

	for i in range(0,3):
		sigma[i] = A[i] * np.exp(-E[i] / (R_const*T))

	cond = (sigma[0] * sigma[1] * sigma[2])**(1.0/3.0)

	return cond
	
def Yoshino2012_DryOlivine_param1_xfe(T, P, water, param1, param2, fo2, fo2ref, method):

	#Conductivity model from Yoshino et al. (2012, JGR:SE)
	#Function does not contain the effects of pressure since it has almost no effect.

	R_const = 8.3144621

	A0_yosh = 10.0**2.72
	H_yosh = 196000.0
	alpha_yosh = 149000.0

	cond = (A0_yosh) * param1 * np.exp(- ((H_yosh) - ((alpha_yosh) * param1**(0.33))) / (R_const * T))

	return cond
	
def Fei2020_WetOlivineIonic_Isotropic(T, P, water, param1, param2, fo2, fo2ref, method):

	sigma1 = 10**11.1
	E1 = 372e3
	dv1 = 3.8e3
	sigma2 = 10**2.3
	E2 = 139e3
	dv2 = 0.3e3
	r = 1.3
	
	water = water / 1e4 #converting to wt%
	
	cond = ((sigma1 / T) * (water ** r) * np.exp(-(E1 + (P*dv1)) / (R_const*T))) +\
		((sigma2) * np.exp(-(E2 + (P*dv2)) / (R_const*T)))
		
	return cond
		