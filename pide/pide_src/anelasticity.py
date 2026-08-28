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
	
def YamauchiTakei2016_old(P,T):

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

def YamauchiTakei2016(
	P,
	T,
	t_sol_C = None,
	period_s=None,
	depth_km=None,
	grain_size_m=None,
	ref_grain_size_m=None,
	grain_size_exponent=3.0,
	melt_fraction=0.0,
	lambda_phi=0.0,
	beta_phi=0.0,
	delta_poro=0.0,
	vp_vs_ratio=np.sqrt(3.0),
	t_sol_50km_C=1326.0,
	dtsol_dz=1.0,
	km_per_GPa=33.0,
	exact=True,
	full_output=False,
	return_extras = False):
	"""
	Near-solidus anelasticity model of Yamauchi & Takei (2016).

	Reference
	---------
	Yamauchi, H. & Takei, Y. (2016), Polycrystal anelasticity at near-solidus
	temperatures, J. Geophys. Res. Solid Earth, 121, 7790-7820,
	doi:10.1002/2016JB013316.

	What the model does
	-------------------
	Converts a pressure/temperature point into a relaxed (anelastic) shear
	modulus, and from there into velocity reduction factors and attenuation.
	The chain is:

		P, T
		  -> solidus Tm(z), homologous temperature Tn = T / Tm
		  -> viscosity eta(T, P, Tn)                         [eq. 17, 18]
		  -> unrelaxed (anharmonic) modulus mu_U, Ju = 1/mu_U
		  -> Maxwell time tauM = eta * Ju
		  -> normalised seismic period tauS_n = tauS / (2*pi*tauM)
		  -> relaxation spectrum X = XB + XP                 [eq. 11-14]
		  -> complex compliance J* = J1 + i*J2               [eq. 15]
		  -> Vs, Qs                                          [eq. 16]

	The point of YT2016, versus the older McCarthy et al. (2011) master
	curve, is that the amplitude A_P and width sigma_P of the high-frequency
	dissipation peak both grow as Tn approaches 1. That breaks the Maxwell
	frequency scaling and makes Vs drop steeply just below the solidus, with
	no melt required.

	Parameters
	----------
	P : float or ndarray
		Pressure in GPa.
	T : float or ndarray
		Temperature in KELVIN.
	period_s : float or ndarray, optional
		Seismic wave period tau_S in seconds. If None, falls back to the
		heuristic tau_S = 3*z/4.2, i.e. a wave of wavelength 3*z travelling
		at 4.2 km/s. That heuristic is NOT from YT2016. Pass the period of
		the data you are actually modelling if you know it.
	depth_km : float or ndarray, optional
		Depth. If None, computed as P * km_per_GPa.
	grain_size_m, ref_grain_size_m, grain_size_exponent
		Optional grain size term (d/dr)**m in eq. 17. YT2016 fitted the
		reference viscosity with dr = d, so leaving these as None reproduces
		the paper and implicitly assumes the mantle grain size is the one
		their inversion absorbed into eta_r.
	melt_fraction, lambda_phi, beta_phi, delta_poro
		Direct melt effects: viscosity reduction exp(-lambda*phi), peak
		amplitude offset beta(phi), and poroelastic softening of Ju. All
		zero reproduces YT2016 Table 4, where melt is neglected because
		phi << 1% in the oceanic upper mantle.
	vp_vs_ratio : float
		Unrelaxed Vp/Vs, used to split Vp into its bulk and shear parts.
		sqrt(3) is a Poisson solid.
	t_sol_50km_C, dtsol_dz
		Solidus at 50 km depth (degC) and its gradient (degC/km). YT2016
		fixed Tm(50 km) = 1326 degC and inverted for Tm(75 km) = 1351 degC,
		so their gradient is exactly 1.0 degC/km. The rheology constants
		hard-coded below came out of that same joint inversion, so changing
		the gradient breaks their internal consistency. It matters: at
		100 km and 1650 K, a gradient of 1.0 gives a 3.1% Vs reduction and
		a gradient of 2.25 gives 1.1%.
	exact : bool
		If True, use the full eq. 16 with the (J2/J1) correction. If False,
		use the Vs ~ 1/sqrt(rho*J1) approximation. The difference is O(Q^-2)
		and stays below 0.1% for anything above Q ~ 30.
	full_output : bool
		If True, also return the intermediate quantities as a dict.

	Returns
	-------
	(v_b, v_p, v_s) : tuple
		Multiplicative velocity factors V_anelastic / V_unrelaxed for bulk
		sound speed, P wave and S wave. These are RATIOS, not velocities.
		Multiply them by the anharmonic velocities computed elsewhere, and
		make sure that elsewhere uses the same mu_U convention as below,
		otherwise the Maxwell time here is inconsistent with the velocities
		you apply it to.
	extras : dict
		Only if full_output=True. Holds Tn, eta, mu_U_GPa, tauM, tauS_n,
		J1_over_Ju, J2_over_Ju, Qs_inv and Qp_inv.
	"""

	# ==================================================================
	# Constants
	# ==================================================================
	R_GAS = 8.3145  # gas constant, J / (mol K)

	# YT2016 Table 4: shape of the relaxation spectrum, fitted to the
	# borneol experiments and treated as universal for polycrystalline
	# aggregates. Do not retune these casually.
	A_B = 0.664       # amplitude of the high-temperature background XB  [eq. 8]
	ALPHA = 0.38      # slope, XB = A_B * tau_n**alpha                   [eq. 8]
	TAU_P_N = 6.0e-5  # normalised time scale of the high-frequency peak [eq. 9]
	T_ETA_N = 0.94    # Tn above which the activation energy increases   [eq. 18]
	GAMMA = 5.0       # extra viscosity reduction factor at Tn = 1       [eq. 18]

	# YT2016 Figure 20: rheology and elasticity from inverting the model
	# against the Pacific Vs(T) profiles of Priestley & McKenzie (2013) at
	# 50 and 75 km depth. These are the output of ONE joint inversion, so
	# they are consistent with each other and with the solidus assumed
	# there. Changing one in isolation breaks that.
	MU0_GPA = 72.45    # unrelaxed shear modulus at 0 degC and 0 Pa   [GPa]
	DMU_DT = -0.01094  # d(mu_U)/dT                                   [GPa/K]
	DMU_DP = 1.987     # d(mu_U)/dP                                   [-]
	ETA_R = 6.22e21    # reference viscosity at Tr, Pr, dr            [Pa s]
	E_ACT = 462.5e3    # activation energy H                          [J/mol]
	V_ACT = 7.913e-6   # activation volume V                          [m3/mol]
	T_REF_K = 1473.0   # reference temperature Tr = 1200 degC         [K]
	P_REF_PA = 1.5e9   # reference pressure Pr = 1.5 GPa              [Pa]

	# MU0_GPA is defined at 0 degC (YT2016 Fig. 20 caption), so DMU_DT has
	# to multiply temperature in CELSIUS, not kelvin. Getting this wrong
	# puts mu_U about 3 GPa (5%) too low.
	T_MODULUS_REF_K = 273.15

	P = np.asarray(P, dtype=float)
	T_K = np.asarray(T, dtype=float)

	if depth_km is None:
		depth_km = P * km_per_GPa
	P_Pa = P * 1.0e9

	# ==================================================================
	# Homologous temperature Tn = T / Tm, both in kelvin
	# ==================================================================
	# Tn is the master variable of the whole model. Everything that makes
	# YT2016 different from a plain master-curve model is a function of Tn:
	# the viscosity drop A_eta, the peak amplitude A_P, the peak width.
	if t_sol_C is None:
		T_sol_C = t_sol_50km_C + dtsol_dz * (depth_km - 50.0)
	else:
		T_sol_C = np.asarray(t_sol_C, dtype= float)
	
	Tn = T_K / (T_sol_C + 273.15)

	# ==================================================================
	# Extra viscosity reduction A_eta(Tn), eq. 18
	# ==================================================================
	#	Tn <  0.94       A_eta = 1                  (normal Arrhenius creep)
	#	0.94 <= Tn < 1   drops smoothly 1 -> 1/gamma
	#	Tn >= 1          A_eta = (1/gamma)*exp(-lambda*phi)
	#
	# The middle branch is equivalent to switching to a higher activation
	# energy H + dH, with dH = R*Tm*T_eta_n/(1 - T_eta_n)*ln(gamma). It
	# encodes the observation that grain boundaries soften BEFORE melting
	# starts ("premelting"), not at the solidus itself.
	#
	# np.where evaluates every branch, so clip the argument first: the
	# middle branch divides by Tn*(1 - T_eta_n) and would blow up at Tn = 0
	# even where the result is discarded.
	Tn_safe = np.clip(Tn, 1e-6, None)
	A_eta_middle = np.exp(
		-((Tn_safe - T_ETA_N) / (Tn_safe * (1.0 - T_ETA_N))) * np.log(GAMMA)
	)
	A_eta = np.where(
		Tn < T_ETA_N,
		1.0,
		np.where(
			Tn < 1.0,
			A_eta_middle,
			(1.0 / GAMMA) * np.exp(-lambda_phi * melt_fraction),
		),
	)

	# ==================================================================
	# Viscosity, eq. 17
	# ==================================================================
	# eta = eta_r * (d/dr)^m * exp[H/R*(1/T - 1/Tr)]
	#             * exp[V/R*(P/T - Pr/Tr)] * A_eta(Tn)
	if grain_size_m is not None and ref_grain_size_m is not None:
		grain_term = (grain_size_m / ref_grain_size_m) ** grain_size_exponent
	else:
		grain_term = 1.0

	eta = (
		ETA_R
		* grain_term
		* np.exp((E_ACT / R_GAS) * (1.0 / T_K - 1.0 / T_REF_K))
		* np.exp((V_ACT / R_GAS) * (P_Pa / T_K - P_REF_PA / T_REF_K))
		* A_eta
	)

	# ==================================================================
	# Unrelaxed (anharmonic) modulus and compliance
	# ==================================================================
	# Ju is scaled by (1 + delta_poro) for the poroelastic softening of the
	# undrained frame when melt is present. delta_poro ~ 0 for phi << 1.
	mu_U_GPa = MU0_GPA + DMU_DP * P + DMU_DT * (T_K - T_MODULUS_REF_K)
	Ju = (1.0 + delta_poro) / (1.0e9 * mu_U_GPa)

	# ==================================================================
	# Maxwell time and normalised seismic period
	# ==================================================================
	tauM = eta * Ju
	if period_s is None:
		# Heuristic only: wavelength ~ 3*depth, group velocity ~ 4.2 km/s.
		period_s = 3.0 * depth_km / 4.2
	tauS_n = period_s / (2.0 * np.pi * tauM)

	# ==================================================================
	# Peak amplitude A_P(Tn), eq. 13
	# ==================================================================
	#	Tn <  0.91          0.01
	#	0.91 <= Tn < 0.96   0.01 + 0.4*(Tn - 0.91)   ramps 0.01 -> 0.03
	#	0.96 <= Tn < 1      0.03
	#	Tn >= 1             0.03 + beta(phi)
	A_P = np.where(
		Tn < 0.91,
		0.01,
		np.where(
			Tn < 0.96,
			0.01 + 0.4 * (Tn - 0.91),
			np.where(Tn < 1.0, 0.03, 0.03 + beta_phi),
		),
	)

	# ==================================================================
	# Peak width sigma_P(Tn), eq. 14
	# ==================================================================
	#	Tn <  0.92          4
	#	0.92 <= Tn < 1      4 + 37.5*(Tn - 0.92)     ramps 4 -> 7
	#	Tn >= 1             7
	#
	# This broadening is what lets the peak reach the very high normalised
	# frequencies seismic waves actually sample (tauS_n ~ 1e-9), which is
	# why Vs falls off a cliff near the solidus.
	sigma_P = np.where(
		Tn < 0.92,
		4.0,
		np.where(Tn < 1.0, 4.0 + 37.5 * (Tn - 0.92), 7.0),
	)

	# ==================================================================
	# Complex compliance J* = J1 + i*J2, eq. 15
	# ==================================================================
	# J1/Ju is the storage term, i.e. the modulus reduction. Three parts:
	#	1                                     unrelaxed (elastic) part
	#	A_B*tau_n^alpha/alpha                 high-temperature background
	#	sqrt(2pi)/2 * A_P*sigma_P*(1 - erf)   high-frequency peak
	# The erf term is the fraction of the Gaussian peak lying at higher
	# frequencies than the wave, so it is what actually softens the rock.
	#
	# J2/Ju is the loss term, needed for Q and for the exact Vs. Its last
	# term (+ tauS_n) is the viscous (Maxwell) contribution.
	with np.errstate(divide="ignore", invalid="ignore"):
		log_ratio = np.log(TAU_P_N / tauS_n)

		J1_over_Ju = (
			1.0
			+ (A_B * tauS_n ** ALPHA) / ALPHA
			+ (np.sqrt(2.0 * np.pi) / 2.0)
			* A_P
			* sigma_P
			* (1.0 - erf(log_ratio / (np.sqrt(2.0) * sigma_P)))
		)

		J2_over_Ju = (np.pi / 2.0) * (
			A_B * tauS_n ** ALPHA
			+ A_P * np.exp(-(log_ratio ** 2) / (2.0 * sigma_P ** 2))
		) + tauS_n

	# ==================================================================
	# Velocities and attenuation, eq. 16
	# ==================================================================
	# Vs = 1/sqrt(rho*J1) * [(1 + sqrt(1 + (J2/J1)^2))/2]^(-1/2)
	# Dividing by the unrelaxed Vs = 1/sqrt(rho*Ju) leaves a pure ratio.
	J2_J1 = J2_over_Ju / J1_over_Ju
	if exact:
		correction = ((1.0 + np.sqrt(1.0 + J2_J1 ** 2)) / 2.0) ** -0.5
	else:
		correction = 1.0
	v_anelasticity_s = correction / np.sqrt(J1_over_Ju)

	# Shear modulus reduction factor, mu_relaxed / mu_U.
	mu_factor = v_anelasticity_s ** 2

	# Bulk sound speed Vb = sqrt(K/rho). Only the shear modulus relaxes in
	# this model, so K is untouched and Vb is unchanged.
	v_anelasticity_b = np.ones_like(np.asarray(v_anelasticity_s, dtype=float))

	# Vp^2 = (K + 4/3 mu)/rho. Writing K = mu_U*(r^2 - 4/3) with r = Vp/Vs
	# gives an exact ratio with no linearisation needed.
	r2 = vp_vs_ratio ** 2
	v_anelasticity_p = np.sqrt((r2 - 4.0 / 3.0 + (4.0 / 3.0) * mu_factor) / r2)

	if not full_output:
		return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s

	# Qs^-1 = J2/J1 * [(1 + sqrt(1 + (J2/J1)^2))/2]^-1
	Qs_inv = J2_J1 * ((1.0 + np.sqrt(1.0 + J2_J1 ** 2)) / 2.0) ** -1.0
	# With no bulk dissipation, Qp^-1 = (4/3)*(Vs/Vp)^2 * Qs^-1.
	Qp_inv = (4.0 / 3.0) * (1.0 / r2) * Qs_inv

	extras = {
		"Tn": Tn,
		"eta": eta,
		"mu_U_GPa": mu_U_GPa,
		"tauM": tauM,
		"tauS_n": tauS_n,
		"J1_over_Ju": J1_over_Ju,
		"J2_over_Ju": J2_over_Ju,
		"Qs_inv": Qs_inv,
		"Qp_inv": Qp_inv,
	}
	
	if return_extras == True:
		return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s, extras
	else:
		return v_anelasticity_b, v_anelasticity_p, v_anelasticity_s