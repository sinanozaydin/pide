import numpy as np
from scipy.special import erf

def Q1_Sobolev1996(P,T):
	#Attenuation model from Sobolev et al. (1996)
	A = 0.148
	f = 0.02
	omega = 2.*np.pi*f
	a = 0.15
	H = 500000.
	V = 0.000020
	#universal gas constant
	R = 8.31446
	E = H + P*1e9*V # Gpa to Pa
	Q = A*(omega**a)*np.exp((a*E)/R/T)
	v_anelasticity_b = 1. - 2./Q/np.tan(a*np.pi/2.)/9
	v_anelasticity_p = 1. - 2./Q/np.tan(a*np.pi/2.)/9
	v_anelasticity_s = 1. - 1./Q/np.tan(a*np.pi/2.)/2
	
	return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s
	
def Q2_Berckhemer1982(P,T):

	#Attenuation model from Berckhemer et al. (1982)
	A = 0.0002
	f = 0.02
	omega = 2.*np.pi*f
	a = 0.25
	H = 584000.
	V = 0.000021
	#universal gas constant
	R = 8.31446
	E = H + P*1e9*V # Gpa to Pa
	Q = A*(omega**a)*np.exp((a*E)/R/T)
	v_anelasticity_b = 1. - 2./Q/np.tan(a*np.pi/2.)/9
	v_anelasticity_p = 1. - 2./Q/np.tan(a*np.pi/2.)/9
	v_anelasticity_s = 1. - 1./Q/np.tan(a*np.pi/2.)/2
	
	return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s
	
def JacksonFaul2010(P,T):

	#Attenuation model from Jackson & Faul, 2010
	A = 816
	t = 75 # time period (seconds)
	grain_size = 10 # grain size (mm)
	a = 0.36
	E = 293e3 # Activation energy
	V = 1.2e-5 # Activation volumn
	#universal gas constant
	R = 8.314472
	H = E + P*1e9*V # Gpa to Pa
	Qs_inv = A*(t/grain_size/1e3*np.exp(-H/R/T))**a
	Qp_inv = Qs_inv * 4/9
	v_anelasticity_b = 1. - Qp_inv/np.tan(a*np.pi/2.)/2.
	v_anelasticity_p = 1. - Qp_inv/np.tan(a*np.pi/2.)/2.
	v_anelasticity_s = 1. - Qs_inv/np.tan(a*np.pi/2.)/2.
	
	return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s
	
def YamauchiTakei2016(P,T):

	#Attenuation model from Yamauchi & Takei, 2016
	R=8.3145
	Ab=0.664
	alpha=0.38
	tauP=6.e-5
	Teta=0.94
	beta=0.
	delphi=0.
	gamma=5.
	lambdaphi=0
	TKr=1473.
	Pr=1.5e9
	
	mu0 = 72.45
	dmudT = -0.01094
	dmudP = 1.987
	eta0 = 6.22*10**21
	E = 462.5e3
	Va = 7.913e-6
	dTdz = 2.25
	sol50 = 1326 # That might change due to composition
	dep = P * 33
	Pg=(dep/33.)
	P=Pg*1.e9 # Gpa to P
	TK = T # T in K
	Tsol = sol50 + (dTdz * (dep - 50.0))
	Tn = TK / (Tsol + 273.0)
	# ---- Aeta ----
	Aeta = np.where(
		Tn < Teta,
		1.0,
		np.where(
			Tn < 1.0,
			np.exp((-1.0 * ((Tn - Teta) / (Tn - (Tn * Teta)))) * np.log(gamma)),
			(1.0 / gamma) * np.exp(-delphi),
		),
	)
	# ---- Viscosity ----
	eta = (
		eta0
		* np.exp((E / R) * (1.0 / TK - 1.0 / TKr))
		* np.exp((Va / R) * (P / TK - Pr / TKr))
		* Aeta
	)
	# ---- Unrelaxed compliance ----
	Ju = 1.0 / (1.0e9 * (mu0 + (dmudP * Pg) + (dmudT * T)))
	# ---- Maxwell time and scaled period ----
	tauM = eta * Ju
	tau = (3.0 * dep) / 4.2
	tauS = tau / (2.0 * np.pi * tauM)
	# ---- Anelastic amplitude (Ap) ----
	Ap = np.where(
		Tn < 0.91,
		0.01,
		np.where(
			Tn < 0.96,
			0.01 + (0.4 * (Tn - 0.91)),
			np.where(Tn < 1.0, 0.03, 0.03 + beta),
		),
	)
	# ---- Stress parameter (sigmap) ----
	sigmap = np.where(
		Tn < 0.92,
		4.0,
		np.where(Tn < 1.0, 4.0 + (37.5 * (Tn - 0.92)), 7.0),
	)
	factor = (1.+((Ab*(tauS**alpha))/alpha)+((np.sqrt(2.*np.pi)/2.)*Ap*sigmap*(1.-erf((np.log(tauP/tauS))/(np.sqrt(2.)*sigmap)))))
	v_anelasticity_b = 5/9+4/np.sqrt(factor)/9
	v_anelasticity_p = 5/9+4/np.sqrt(factor)/9
	v_anelasticity_s = 1/np.sqrt(factor)
	
	return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s