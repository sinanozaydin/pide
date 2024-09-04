#!/usr/bin/env python3

import numpy as np

R_const = 8.3144621

def Dai2014_DryandWetOlivine_fo2(T, P, water, xFe ,param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	dv_dai2014 = -0.86 * 1e3 # cm^3 /mol Taken from Dai2014b-PEPI, Error is insignificant, since the effect itself is insignificant...
	p_ref = 4.0
	cw_ref = 460.0 / 1e4
	e1_dai = 74000.0
	e1_di_err = 3000
	e2_dai = 115000
	r_dai = 0.8
	q_dai = -0.066 #Taken from Dai2014c-P"EPI, Error is insignificant.
	q_dry = 0.16666 #Taken from Constable (2006)
		
	dai_ref = ((10.0**0.48) * np.exp(-((e1_dai) + (p_ref*dv_dai2014)) / (R_const * T))) + (10.0**2.84 * np.exp(-((e2_dai) + (p_ref*dv_dai2014)) / (R_const * T)))
	cond_wet = (dai_ref * (water / (cw_ref))**(r_dai) * (fo2/fo2_ref)**(q_dai)) *  np.exp(- ((P - p_ref) * dv_dai2014) / (R_const*T))
	dai_dry =  10**2.4 * ((fo2/fo2_ref)**(q_dry)) * np.exp(-(154000.0 + (P * dv_dai2014)) / (R_const * T))
	
	if mechanism == None:
		cond = np.real(dai_dry + cond_wet)
	elif mechanism == 'proton':
		cond = np.real(cond_wet)
	elif mechanism == 'dry':
		cond = np.real(dai_dry)
	else:
		cond = np.real(dai_dry)
		print('There is no polaron & ionic distinction for Dai2014_DryandWetOlivine_fo2.  Using Dry conductivity output...')

	return cond
	
def Dai2020_WetOlivine_200ppmTi_fo2(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	dv_dai2014 = -0.86 * 1e3 # m^3 /mol Taken from Dai2014b-PEPI, Error is insignificant, since the effect itself is insignificant...
	q_dai = -0.066 #Taken from Dai2014c-P"EPI, Error is insignificant.
	q_dry = 0.16666 #Taken from Constable (2006)

	A_dai_200 = 3.01

	r_dai_200 = 0.51

	E_dai_200 = 87000.0

	cond_wet = (10.0 ** A_dai_200) * ((water)**r_dai_200) * np.exp(-(E_dai_200) / (R_const * T))

	dai_dry =  10**2.4 * ((fo2/fo2_ref)**(q_dry)) * np.exp(-(154000.0 + (P * dv_dai2014)) / (R_const * T))
	
	if mechanism == None:
		cond = dai_dry + cond_wet
	elif mechanism == 'proton':
		cond = cond_wet
	elif mechanism == 'dry':
		cond = dai_dry
	else:
		cond = dai_dry
		print('Warning: There is no polaron & ionic distinction for Dai2020_WetOlivine_200ppmTi_fo2. Using Dry conductivity output...')

	return cond
	
def Dai2020_WetOlivine_683ppmTi_fo2(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	dv_dai2014 = -0.86 * 1e3 # m^3 /mol Taken from Dai2014b-PEPI, Error is insignificant, since the effect itself is insignificant...
	q_dai = -0.066 #Taken from Dai2014c-P"EPI, Error is insignificant.
	q_dry = 0.16666 #Taken from Constable (2006)

	A_dai_200 = 3.01

	r_dai_200 = 0.51

	E_dai_200 = 87000.0

	cond_wet = (10.0 ** A_dai_200) * ((water)**r_dai_200) * np.exp(-(E_dai_200) / (R_const * T))

	dai_dry =  10**2.4 * ((fo2/fo2_ref)**(q_dry)) * np.exp(-(154000.0 + (P * dv_dai2014)) / (R_const * T))

	if mechanism == None:
		cond = dai_dry + cond_wet
	elif mechanism == 'proton':
		cond = cond_wet
	elif mechanism == 'dry':
		cond = dai_dry
	else:
		cond = dai_dry
		print('Warning: There is no polaron & ionic distinction for Dai2020_WetOlivine_683ppmTi_fo2. Using Dry conductivity output...')

	return cond
	
def Poe2010_DryandWetOlivine(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	E1 = [146e3,126e3]
	E2 = [112e3,150e3]
	E3 = [129e3,81.2e3]
	sigma_1 = [2.52,2.59]
	sigma_2 = [1.139,3.46]
	sigma_3 = [1.99,1.02]
	alpha_1 = 1180
	alpha_2 = 1430
	alpha_3 = 700
	
	if (mechanism == None):
		cond_1 = (10**sigma_1[0] * np.exp(-(E1[0])) / (R_const*T)) + (10**sigma_1[1] * water * np.exp(-(E1[1] + (alpha_1 * water)) / (R_const*T)))
		cond_2 = (10**sigma_2[0] * np.exp(-(E2[0])) / (R_const*T)) + (10**sigma_2[1] * water * np.exp(-(E2[1] + (alpha_2 * water)) / (R_const*T)))
		cond_3 = (10**sigma_3[0] * np.exp(-(E3[0])) / (R_const*T)) + (10**sigma_3[1] * water * np.exp(-(E3[1] + (alpha_3 * water)) / (R_const*T)))
			
	elif mechanism == 'proton':
		cond_1 = (10**sigma_1[1] * water * np.exp(-(E1[1] + (alpha_1 * water)) / (R_const*T)))
		cond_2 = (10**sigma_2[1] * water * np.exp(-(E2[1] + (alpha_2 * water)) / (R_const*T)))
		cond_3 = (10**sigma_3[1] * water * np.exp(-(E3[1] + (alpha_3 * water)) / (R_const*T)))
		
	elif (mechanism == 'polaron') or (mechanism == 'dry'):
		cond_1 = (10**sigma_1[0] * np.exp(-(E1[0])) / (R_const*T))
		cond_2 = (10**sigma_2[0] * np.exp(-(E2[0])) / (R_const*T))
		cond_3 = (10**sigma_3[0] * np.exp(-(E3[0])) / (R_const*T))
		
	cond = (cond_1*cond_2*cond_3)**(1.0/3.0)
	return cond
	
def Constable2006_dryOlivine_fo2(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	e = 1.602e-19 #charge of the electron in coulombs
	k = 8.617e-5
	kT = T * k
	bfe = (5.06e24) * np.exp(-0.357 / kT)
	bmg = (4.58e26) * np.exp(-0.752 / kT)
	ufe = (12.2e-6) * np.exp(-1.05 / kT)
	umg = (2.72e-6) * np.exp(-1.09 / kT)
	condfe = bfe + (3.33e24 * np.exp(-0.02/kT) * (fo2*1e5)**(0.16666))
	condmg = bmg + (6.21e30 * np.exp(-1.83/kT) * (fo2*1e5)**(0.16666))

	if (mechanism == 'proton') or (mechanism == 'ionic'):
		raise ValueError('Proton conduction is not included in electrical conductivity model: Constable2006_dryOlivine_fo2')
	elif (mechanism == 'polaron') or (mechanism == 'dry') or (mechanism == None):
		cond = (condfe * ufe * e) + (2.0*condmg* umg * e)
	else:
		raise ValueError('Unknown mechanism: ' + mechanism)
	return cond

def Dai2014_DryOlivine_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_dai2014c = 10**(2.77 + (-1.19 * xFe))
	e_dai2014c = 162000.0 + (-63000 * xFe)

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Dai2014_DryOlivine_xFe')
	else:
		cond = sigma_dai2014c * np.exp(-(e_dai2014c) / (R_const * T))

	return cond

def Fullea2011_DryOlivine_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma_i_fullea = 10**4.73
	sigma_pol_fullea = 10**2.7 #average value of used models Fullea et al. (2011)

	e_i_fullea = 231000.0

	a_ful = [1.642,0.246,-4.85,3.259] #polynomial coefficients that calculates the fe-dependency of act. enthalpies from Omura et al (1989)
	dv_ful = 0.68 * 1e-6
	fe_pol_fullea = (a_ful[0] + (a_ful[1] * xFe) + (a_ful[2] * (xFe**2.0)) + (a_ful[3]* (xFe**3.0)) + (P * dv_ful)) * 1e5
	
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Fullea2011_DryOlivine_xFe')
	else:
		cond = (sigma_i_fullea * np.exp(-e_i_fullea / (R_const * T))) + (sigma_pol_fullea * np.exp(-fe_pol_fullea / (R_const * T)))

	return cond
	
def Pommier2018_ShearedDryOlivine(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A = [284.0,780.0,261.0]

	E = [126400.0,122700.0,114500.0]

	sigma0 = A[0] * np.exp(-E[0] / (R_const*T))
	sigma1 = A[1] * np.exp(-E[1] / (R_const*T))
	sigma2 = A[2] * np.exp(-E[2] / (R_const*T))

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Fullea2011_DryOlivine_xFe')
	else:
		cond = (sigma0 * sigma1 * sigma2)**(1.0/3.0)

	return cond
	
def Yoshino2012_DryOlivine_xFe(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	#Conductivity model from Yoshino et al. (2012, JGR:SE)
	#Function does not contain the effects of pressure since it has almost no effect.

	R_const = 8.3144621

	A0_yosh = 10.0**2.72
	H_yosh = 196000.0
	alpha_yosh = 149000.0

	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Yoshino2012_DryOlivine_xFe')
	else:
		cond = (A0_yosh) * xFe * np.exp(- ((H_yosh) - ((alpha_yosh) * xFe**(0.33))) / (R_const * T))

	return cond
	
def Fei2020_WetOlivineIonic_Isotropic(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	sigma1 = 10**11.1
	E1 = 372e3
	dv1 = 3.8e3
	sigma2 = 10**2.3
	E2 = 139e3
	dv2 = 0.3e3
	r = 1.3
		
	if mechanism == 'proton':
		raise ValueError('Proton conduction is not included in electrical conductivity model: Yoshino2012_DryOlivine_xFe')
	else:
		cond = ((sigma1 / T) * (water ** r) * np.exp(-(E1 + (P*dv1)) / (R_const*T))) +\
			((sigma2) * np.exp(-(E2 + (P*dv2)) / (R_const*T)))
		
	return cond

def Novella2017_HDiffusion(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None, mechanism = None):

	A_H = [-0.7,-5,-3.5]

	EH = [229000.0,172000.0,188000.0]

	if method == 'array':
		T = T[0]
		DH_0 = np.zeros((3,len(T)))

	elif method == 'index':
		T = T
		DH_0 = np.zeros((3,1))

	DH_0 = np.zeros((3,len(T)))

	for i in range(0,3):

		DH_0[i] = 10**(A_H[i]) * np.exp(-EH[i] / (R_const * T))


	DH = (DH_0[0] * DH_0[1] * DH_0[2])**(1.0/3.0)

	return DH

def Sun2019_HDiffusion(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None,mechanism = None):

	A_H = -7.4
	EH = 130000.0

	r_sun = 0.41

	DH = (10.0**A_H) * (water**r_sun) * np.exp(-EH / (R_const * T))

	return DH

def DuFrane2012_HDiffusion(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None,mechanism = None):

	A_H = [-4.9,-5.4,-8.4]
	EH = [140000.0,170000.0,100000.0]

	if method == 'array':
		T = T[0]
		DH_0 = np.zeros((3,len(T)))

	elif method == 'index':
		T = T
		DH_0 = np.zeros((3,1))

	DH_0 = np.zeros((3,len(T)))

	for i in range(0,3):

		DH_0[i] = 10**(A_H[i]) * np.exp(-EH[i] / (R_const * T))

	DH = (DH_0[0] * DH_0[1] * DH_0[2])**(1.0/3.0)

	return DH

def Kohlstedt1998_HDiffusion(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None,mechanism = None):

	A_H = [-3.85,-3.82,-6.83]

	EH = [145000.0,180000.0,110000.0]

	if method == 'array':
		T = T[0]
		DH_0 = np.zeros((3,len(T)))
	elif method == 'index':
		T = T
		DH_0 = np.zeros((3,1))

	DH_0 = np.zeros((3,len(T)))

	for i in range(0,3):
		DH_0[i] = 10**(A_H[i]) * np.exp(-EH[i] / (R_const * T))

	DH = (DH_0[0] * DH_0[1] * DH_0[2])**(1.0/3.0)

	return DH

def Demouchy2006_HDiffusion(T, P, water, xFe, param1, fo2 = None, fo2_ref = None, method = None,mechanism = None):

	A_H = [-4.5,-4.5,-1.4]

	EH = [204000.0,204000.0,258000.0]

	if method == 'array':
		T = T[0]
		DH_0 = np.zeros((3,len(T)))

	elif method == 'index':
		T = T
		DH_0 = np.zeros((3,1))

	for i in range(0,3):
		DH_0[i] = 10**(A_H[i]) * np.exp(-EH[i] / (R_const * T))

	DH = (DH_0[0] * DH_0[1] * DH_0[2])**(1.0/3.0)

	return DH
		