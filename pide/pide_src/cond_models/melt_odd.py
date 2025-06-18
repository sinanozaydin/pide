#!/usr/bin/env python3

import numpy as np
from pide.utils.utils import text_color

R_const = 8.3144621

def Sifre2014_Wet_Carbonated(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	Melt_CO2 = Melt_CO2 * 1e-4 #converting ppm to wt percent

	E_CO2 = 789166.0 * np.exp(-0.1808 * Melt_CO2) + 32820.0
	E_H2O = 88744.0 * np.exp(-0.388 * Melt_H2O) + 73029.0

	sigma_co2 = np.exp((5.5 * (1e-5) * E_CO2) + 5.7956)
	sigma_h2o = np.exp((4.54 * (1e-5) * E_H2O) + 5.5607)

	cond = (sigma_h2o * np.exp(-E_H2O / (R_const * T))) +\
	 	(sigma_co2 * np.exp(-E_CO2 / (R_const * T)))

	return cond

def Pommier2011_SIGMELTS_WetSilicate_Melt(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	"""
	converted to python from webscraped javascript code hosted at SIGMELTS website: https://calcul-isto.cnrs-orleans.fr/apps/sigmelts/#
	"""

	if method == 'array':
		T = T[0]
		P = P[0]*1e3 #converting into MPa
		Melt_H2O = Melt_H2O[0]
		Melt_Na2O = Melt_Na2O[0]
		Melt_SiO2 = Melt_SiO2[0]
		cond = np.zeros(len(T))
	else:
		P = P * 1e3

	def _get_params(sio2):
		# Select model coefficients based on SiO2 content
		if sio2 < 50.0:
			c_intercept   =  28.1072140868
			c_sio2        = -0.30953383
			c_temperature =  0.0004035372
			c_pressure    =  0.0016782775
			c_ln_water    =  0.1189143906
			c_na2o        = -1.0946261697
			c_water       = -0.2161950362
			c_q_rt        = -0.916865944
		elif sio2 < 70.0:
			c_intercept   =  17.85482826
			c_sio2        = -0.052183931
			c_temperature = -0.001171623
			c_pressure    =  0.001009267
			c_ln_water    =  0.08745967
			c_na2o        = -0.707477243
			c_water       = -0.433092118
			c_q_rt        = -0.919252205
		else:
			c_intercept   = -19.08557442
			c_sio2        =  0.489734869
			c_temperature = -0.003553936
			c_pressure    =  0.001306084
			c_ln_water    =  0.290884415
			c_na2o        =  0.0
			c_water       = -0.512471336
			c_q_rt        = -1.143904969

		return c_intercept,c_sio2,c_temperature,c_pressure,c_pressure,c_ln_water,c_na2o,c_water,c_q_rt

	# Activation energy for Pommier model
	ea_pom = (
		7.091674171904e-15 / (5.73165e-21 * Melt_Na2O + 4.042e-20)
		- 1.08426212889391294371e3 * Melt_Na2O
		+ 6.67436156229428138586e3
		+ P * 20.70882392
		- 1000.0 * Melt_H2O ** 2)
	
	# Activation energy for Gaillard model
	
	ea_gai = -2925.0 * np.log(Melt_H2O) + 64132.0 + 20.0 * P
	
	# Natural log conductivity using Gaillard model
	lns_gai = np.log(-78.9 * np.log(Melt_H2O) + 4000.0) - ea_gai / (R_const * T)

	if method == 'index':
		c_intercept,c_sio2,c_temperature,c_pressure,c_pressure,c_ln_water,c_na2o,c_water,c_q_rt = _get_params(sio2=Melt_SiO2)
		lns_pom = (
		c_intercept
		+ c_sio2 * Melt_SiO2
		+ c_temperature * T
		+ c_pressure * P
		+ c_ln_water * np.log(Melt_H2O)
		+ c_na2o * Melt_Na2O
		+ c_water * Melt_H2O
		+ c_q_rt * ea_pom / (8.314472 * T))
		# Return sigma based on comparison
		if Melt_SiO2 > 6 and lns_pom > lns_gai:
			cond =  np.exp(lns_gai)
		else:
			cond =  np.exp(lns_pom)
	else:
		
		lns_pom = np.zeros(len(T))
		for i in range(0,len(T)):
			c_intercept,c_sio2,c_temperature,c_pressure,c_pressure,c_ln_water,c_na2o,c_water,c_q_rt = _get_params(sio2=Melt_SiO2[i])

			lns_pom[i] = (
				c_intercept
				+ c_sio2 * Melt_SiO2[i]
				+ c_temperature * T[i]
				+ c_pressure * P[i]
				+ c_ln_water * np.log(Melt_H2O[i])
				+ c_na2o * Melt_Na2O[i]
				+ c_water * Melt_H2O[i]
				+ c_q_rt * ea_pom[i] / (8.314472 * T[i]))	

			if Melt_SiO2[i] > 6 and lns_pom[i] > lns_gai[i]:
				cond[i] =  np.exp(lns_gai[i])
			else:
				cond[i] =  np.exp(lns_pom[i])

	return cond

def Ni2011_WetBasalt(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):
	
	if method == 'array':
		P = P[0]
		T = T[0]
		Melt_H2O = Melt_H2O[0]

		cond = np.zeros(len(T))

		for i in range(len(T)):
			if T[i] > 1146.8:
				cond[i] = 10**(2.172 - ((860.82 - (204.46*np.sqrt(Melt_H2O[i]))) / (T[i] - 1146.8)))

			else:
				cond[i] = np.nan
	else:

		if T > 1146.8:
			cond = 10**(2.172 - ((860.82 - (204.46*np.sqrt(Melt_H2O))) / (T - 1146.8)))
		else:
			cond = np.nan
	# import ipdb
	# ipdb.set_trace()
	return cond


def TyburczyWaff1983_DryTholeiite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	sigma_0_low = 1.12e5
	E_low = 112000.0
	dv_low = 4.6

	sigma_0_high = 2.15e5
	E_high = 153000.0
	dv_high = -0.06
	
	if method == 'array':
		P = P[0]
		T = T[0]

		cond = np.zeros(len(T))

		for i in range(0,len(T)):
	
			if P[i] < 0.9:
				cond[i] = sigma_0_low * np.exp(-(E_low + (dv_low * P[i] * 1e3)) / (R_const * T[i]))
			else:
				cond[i] = sigma_0_high * np.exp(-(E_high + (dv_high * P[i] * 1e3)) / (R_const * T[i]))
				
	elif method == 'index':
	
		if P < 0.9:
			cond = sigma_0_low * np.exp(-(E_low + (dv_low * P * 1e3)) / (R_const * T))
		else:
			cond = sigma_0_high * np.exp(-(E_high + (dv_high * P * 1e3)) / (R_const * T))

	return cond

def TyburczyWaff1983_DryAndesite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	sigma_0_low = 1.01e3
	E_low = 78000.0
	dv_low = 17.9

	sigma_0_high = 6.61e3
	E_high = 117000.0
	dv_high = 3.25
	
	if method == 'array':

		T = T[0]
		P = P[0]
		
		cond = np.zeros(len(T))
		for i in range(0,len(T)):
	
			if P[i] < 0.9:
				cond[i] = sigma_0_low * np.exp(-(E_low + (dv_low * P[i] * 1e3)) / (R_const * T[i]))
			else:
				cond[i] = sigma_0_high * np.exp(-(E_high + (dv_high * P[i] * 1e3)) / (R_const * T[i]))
				
	elif method == 'index':
	
		if P[i] < 0.9:
			cond = sigma_0_low * np.exp(-(E_low + (dv_low * P * 1e3)) / (R_const * T))
		else:
			cond = sigma_0_high * np.exp(-(E_high + (dv_high * P * 1e3)) / (R_const * T))

	return cond

def Guo2017_WetAndesite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	cond = 10**(5.23 - (0.56 * (Melt_H2O**0.6)) - ((8130.4 - (1462.7*(Melt_H2O**0.6)) + ((581.3 - (12.7*Melt_H2O**2)) * P)) / T))

	return cond

def Laumonier2017_WetAndesite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	P = P * 1e5 #converting gpa to bars

	a = -0.34
	b = 8.96
	c = -8.07e-6
	d = 1.67e-4
	e = -9627
	f = 1.25e5
	g = -1.46e-1
	h = 2.462

	sigma_0 = np.exp(((a*Melt_H2O) + b) + (P * ((c*Melt_H2O) + d)))
	Ea = (e * Melt_H2O) + f
	dV = (g * Melt_H2O) + h

	cond = sigma_0 * np.exp(-(Ea + (P * dV)) / (R_const * T))

	return cond

def Laumonier2015_WetDacite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	P = P * 1e5 #converting gpa to bars

	a = -0.064
	b = 5.96
	c = 1.06e-5
	d = 2.49e-5
	e = -6146.0
	f = 88440.0
	g = 0.176
	h = 0.388

	sigma_0 = np.exp(((a*Melt_H2O) + b) + (P * ((c*Melt_H2O) + d)))
	Ea = (e * Melt_H2O) + f
	dV = (g * Melt_H2O) + h

	cond = sigma_0 * np.exp(-(Ea + (P * dV)) / (R_const * T))

	return cond

def Gaillard2004_WetRhyolite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):
	
	if method == 'array':
		if np.any(Melt_H2O) == 0:
			print(text_color.RED + 'Gaillard2004_WetRhyolite cannot accomodate any dry material.' + text_color.END)
	else:
		if Melt_H2O == 0.0:
			print(text_color.RED + 'Gaillard2004_WetRhyolite cannot accomodate any dry material. + text_color.END')
	
	P = P * 1e3 #to MPa

	Ea = (-2925.0 * np.log(Melt_H2O)) + 64132.0
	sigma_0 = (-78.9 * np.log(Melt_H2O)) + 754.0

	cond = sigma_0 * np.exp((-Ea - (20*P)) / (R_const * T))

	return cond

def Guo2016_WetRhyolite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	cond = 10**(2.983 - (0.0732*Melt_H2O) -\
		((3528 - (233.8*Melt_H2O) + ((763 - 7.5*Melt_H2O**2)*P)) / T))

	return cond

def Chen2018_WetGranite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	sigma_0 = 6.673 - (0.491*P)
	Ea = 58987.0 - (2200 * np.log(Melt_H2O + 0.046))
	dV = 10927.0 - (1981 * Melt_H2O)

	cond = sigma_0 * np.exp(-(Ea + (P * dV)) / (R_const * T))

	return cond

def Guo2018_WetGranite(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):


	cond = 10.0**(3.205 - (0.102 * Melt_H2O) -\
		(-(4228.5 - (354.7 * Melt_H2O) + (693.6 * P)) / T))
		
	return cond

def Poe2008_Phonotephrite_Average(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):

	Ea = 116000.0
	sigma_cond = 34700.0
	
	cond = sigma_cond * np.exp((-Ea) / (R_const*T))
		
	return cond
	
def Poe2012_WetPanteleriteGlass(T, P, Melt_H2O, Melt_CO2, Melt_Na2O, Melt_K2O, Melt_SiO2, method):
	
	if np.mean(Melt_H2O) == 0.0:
		raise ValueError('The electrical conductivity model Poe2012_WetPanteleriteGlass does not support dry melt assemblages.')
		
	sigma_cond = 1260.0
	E = 96485.3

	cond = sigma_cond*(Melt_H2O**(1.25))*np.exp((-0.812*E)/(R_const*T))
	
	return cond
	
