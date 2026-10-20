import numpy as np

def calculate_hasterok2011_geotherm(SHF,  T_0, max_depth, moho, adiabat=True, BDL_T = 0,kinked = False, thermal_lab = False, thermal_lab_depth = None, thermal_lab_temp = None, **kwargs):

	'''

	A function to calculate the generalized steady-state geotherm for the
	continental lithosphere as taken from Hasterok & Chapman (2011).
	:cite: `hasterok2011`

	All the parameters used in constructing the geotherm are derived from the
	cited study. Generalized continental geotherm can be creating assuming
	26% of the heat generation occurs in the upper crustal layer.

	Input
	-------

	SHF: Surface Heat Flow value in mW/m^2

	BDL_T: is Temperature at the base of the depleted lithosphere in Celsius

	T_0: is Temperature at the surface in Celsius

	max_depth: Maximum depth the geotherm is going to be calculated for.

	moho: The depth of Mohorivicic discontinuity (i.e. crustal thickness)

	kinked: Boolean parameter to determine whether the conductive
	geotherm will be kinked parallel to D-G transition curve after the BDL.


		|  'Array' - A numpy array output

	Output
	-------
	Temperature (Kelvin), Depth (kilometers), Pressure (GPa), index point where the geotherm
	is kinked

	'''

	mechanical_lab_depth = kwargs.pop('mechanical_lab_depth', 120)

	if kinked == False:
		BDL_T = 0
	else:
		BDL_T = BDL_T
		d_g_t = [600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
		d_g_p = [35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60]
		d_g_p = np.array(d_g_p) * 0.1
		d_g_t = np.array(d_g_t)
		d_g_slope = (np.array(d_g_t)[-1] - np.array(d_g_t)[0]) / (d_g_p[-1] - d_g_p[0])

	moho = moho * 1e3 #converting moho depth to meters
	T_0 = T_0 + 273.0 #converting celsius to kelvin
	max_depth = max_depth * 1e3 #converting max depth to meters.
	rho_mantle = kwargs.pop('rho_mantle', 3300.0)
	g = 9806.0
	rho_crust = 2850.0

	A_upper_crust = 1e-3 #mW/m^3
	A_lower_crust = 0.4e-3 #mW/m^3
	heat_prod_mantle = 4e-5 #mW/m^3

	interval = 1000.0 #Depth interval for calculations, in kms
	depth = np.arange(0,(max_depth) + interval, interval) #setting up the depth range in meters

	#Setting up the thermal conductivity parameters.
	k0_upper_1 = 1.496e3
	k0_upper_2 = 2.964e3
	k0_mid_1 = 1.733e3
	k0_mid_2 = 2.717e3
	k0_low_1 = 1.723e3
	k0_low_2 = 2.32e3
	k0_mantle_sp = 2.271e3
	k0_mantle_gt = 2.371e3
	k1_upper_1 = 398.84e3
	k1_upper_2 = -495.29e3
	k1_mid_1 = 194.59e3
	k1_mid_2 = -398.93e3
	k1_low_1 = 219.88e3
	k1_low_2 = -96.88e3
	k1_mantle_sp = 681.12e3
	k1_mantle_gt = 669.4e3
	k2_upper_1 = 4.573e-4
	k2_upper_2 = 0.866e-4
	k2_mid_1 = 2.906e-4
	k2_mid_2 = 0.032e-4
	k2_low_1 = 1.705e-4
	k2_low_2 = -0.981e-4
	k2_mantle_sp = -1.259e-4
	k2_mantle_gt = -1.288e-4
	k3_upper_1 = 0.0950
	k3_upper_2 = 0.0692
	k3_mid_1 = 0.0788
	k3_mid_2 = 0.0652
	k3_low_1 = 0.0520
	k3_low_2 = 0.0384
	k3_mantle_sp = 0.0399
	k3_mantle_gt = 0.0384

	p = np.zeros(len(depth)) #Setting up pressure array
	density_mantle = np.zeros(len(depth)) #setting up density array

	for i in range(1,len(depth)):
		if depth[i] <= moho:
			p[i] = interval * rho_crust * g
			density_mantle[i] = rho_crust
		else:
			p[i] = interval * rho_mantle * g
			density_mantle[i] = rho_mantle

	p = np.cumsum(p) / 1e12 #cumulative addition and converting to GPa

	t_criterion = 844.0 #Kelvin

	def calculate_k_st(k0,k1,k2,k3,Temp,P):

		k = (k0 + (k1 / Temp) + (k2 * Temp**2.0)) * (1 + (k3 * P))

		return k

	def interception(y1,y2):
		#Local function to calculate the interception point of arrays that have
		#used same temperature variation as paramaeter:(y1(T),y2(T))
		idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)

		return idx

	#Setting up parameters as numpy arrays at length of d_z (layers)
	T = np.zeros(len(depth)) #Temperature in Kelvin
	q = np.zeros(len(depth)) #Heat flow
	A_list = np.zeros(len(depth)) #Heat generation
	k = np.zeros(len(depth)) #Thermal conductivity
	#Setting up the first layers
	T[0] = T_0 #Temperature at surface
	q[0] = SHF #Heat flow at surface
	k[0] = (k0_upper_1 + (k1_upper_1 / T[0]) + (k2_upper_1 * T[0]**2.0)) #Thermal
	#conductivity at surface, calculated via Hasterok (2010, PhD thesis).

	#Logical parameter that will be changed when sp-garnet transition
	#going to be intercepted.
	transition_sp_search = True

	hr = 16 * 1e3
	A_upper_crust = (SHF*0.26) / 16000.0 #26% of heat generation happens in first 16 km

	#Setting up for loop that will iterate for number of layers.
	for i in range(0,len(depth)):
		#If conditionals setting the heat production parameters based on given input.

		if depth[i] <= hr:
			A_list[i] = A_upper_crust
		elif (depth[i] > hr) and (depth[i] <= moho):
			A_list[i] = A_lower_crust
		elif (depth[i] > moho) and (depth[i] <= max_depth):
			A_list[i] = heat_prod_mantle

		if i != 0:

			#Solving the temperature with formula (?)
			T[i] = T[i-1] + ((q[i-1] * interval) / k[i-1]) - ((A_list[i-1] * interval**2.0) / (2.0 * k[i-1]))
			T_avg = (T[i-1] + T[i]) / 2.0 #Averaging T

			p_spinel = 1.4209 + np.exp((3.9073 * 1e-3 * T) - 6.8041)  #Spinel-Garnet transition line calculated with T
			#Taken from Hasterok & Chapman (2011)
			p_spinel = p_spinel[np.nonzero(T)] #Getting rid of zeros


			#Calculating T-dependent thermal conductivity
			#If conditionals for different empirical parameters for critical temperature (Hasterok, 2010)

			if T[i] <= t_criterion:
				if depth[i] <= hr:
					k[i] = calculate_k_st(k0_upper_1,k1_upper_1,k2_upper_1,k3_upper_1,T_avg,p[i])
				elif (depth[i] > hr) and (depth[i] <= moho):
					k[i] = calculate_k_st(k0_low_1,k1_low_1,k2_low_1,k3_low_1,T_avg,p[i])
				else:
					#In the mantle searching for sp-gt transition
					if transition_sp_search == True:
						idx_sp = interception(p_spinel,p[:len(p_spinel)])
						if idx_sp != None:
							#If producted geotherm intercepts P-T line of sp-gt transition
							#change logical parameter to False for changing to gt.
							transition_sp_search = False
							depth_spinel = depth[i] #Depth of transition
							p_spinel_trans = p[i]
							index_spinel = i
						k[i] = calculate_k_st(k0_mantle_sp,k1_mantle_sp,k2_mantle_sp,k3_mantle_sp,T_avg,p[i])
					else:
						k[i] = calculate_k_st(k0_mantle_gt,k1_mantle_gt,k2_mantle_gt,k3_mantle_gt,T_avg,p[i])
			else:
				#Same if tree for temperatures bigger than critical temperature : (t_criterion)
				if depth[i] <= hr:
					k[i] = calculate_k_st(k0_upper_2,k1_upper_2,k2_upper_2,k3_upper_2,T_avg,p[i])
				elif (depth[i] > hr) and (depth[i] <= moho):
					k[i] = calculate_k_st(k0_low_2,k1_low_2,k2_low_2,k3_low_2,T_avg,p[i])
				else:
					if transition_sp_search == True:
						idx_sp = interception(p_spinel,p[:len(p_spinel)])
						if idx_sp != None:
							transition_sp_search = False
							depth_spinel = depth[i]
							p_spinel_trans = p[i]
						k[i] = calculate_k_st(k0_mantle_sp,k1_mantle_sp,k2_mantle_sp,k3_mantle_sp,T_avg,p[i])
					else:
						k[i] = calculate_k_st(k0_mantle_gt,k1_mantle_gt,k2_mantle_gt,k3_mantle_gt,T_avg,p[i])

			q[i] = q[i-1] - (A_list[i-1] * interval)
			if depth[i] == moho:
				crust_production = SHF - q[i]
				moho_heat_flow = q[i]

	if kinked == True:

		idx_geotherm_nearest = (np.abs(T-BDL_T)).argmin()
		t_geotherm_inflict = T[idx_geotherm_nearest]
		p_geotherm_inflict = p[idx_geotherm_nearest]
		depth_geotherm_inflict = depth[idx_geotherm_nearest]

		kinked_geotherm_add = t_geotherm_inflict - (d_g_slope * p_geotherm_inflict)
		kinked_geotherm = kinked_geotherm_add + (d_g_slope * p)

		T[idx_geotherm_nearest:] = kinked_geotherm[idx_geotherm_nearest:]

		if adiabat == True:
			adiabat = False
			print('Both "kinked" and "adiabat" cannot be True. Turning adiabat==False.')

	else:

		idx_geotherm_nearest = 0

	if adiabat == True:
		
		if thermal_lab and thermal_lab_temp is not None:
			# Find where conductive geotherm reaches the LAB temperature
			thermal_lab_temp_K = thermal_lab_temp + 273.15
			idx_LAB = np.argmin(np.abs(T - thermal_lab_temp_K))
			
			# If conductive geotherm never reaches LAB temp, use the last point
			if T[idx_LAB] < thermal_lab_temp_K * 0.95:
				idx_LAB = len(T) - 1
			
			# Build adiabat starting from LAB temperature
			grad_adiabat = T_Katsura_2022_Adiabat(p, out_type='grad')
			T[idx_LAB] = thermal_lab_temp_K
			for idx_rc in range(idx_LAB + 1, len(T)):
				T[idx_rc] = T[idx_rc - 1] + (grad_adiabat[idx_rc] * (depth[idx_rc] - depth[idx_rc - 1]) * 1e-3)


		elif thermal_lab and thermal_lab_depth is not None:
			
			idx_LAB = np.argmin(np.abs(depth/1e3 - thermal_lab_depth))
			
			T_at_LAB = T[idx_LAB]
			
			grad_adiabat = T_Katsura_2022_Adiabat(p, out_type='grad')
			for idx_rc in range(idx_LAB + 1, len(T)):
				T[idx_rc] = T[idx_rc-1] + (grad_adiabat[idx_rc] * (depth[idx_rc] - depth[idx_rc-1]) * 1e-3)
		
		else:
			# Original behaviour: find intersection with Katsura adiabat
			T_C_Adiabat, T_K_Adiabat = T_Katsura_2022_Adiabat(p)
			idx_LAB = np.argwhere(np.diff(np.sign(T - T_K_Adiabat)) != 0)
			try:
				T[idx_LAB[0][0]:] = T_K_Adiabat[idx_LAB[0][0]:]
			except IndexError:
				pass

	if kinked == False:
		if adiabat == True:
			return T, depth/1e3, p, idx_LAB
		else:
			return T, depth/1e3, p
	else:
		if adiabat == True:
			return T, depth/1e3, p, idx_LAB, idx_geotherm_nearest
		else:
			return T, depth/1e3, p, idx_geotherm_nearest
	
def T_Katsura_2022_Adiabat(P_input, out_type = 'T'):

	'''
	A function that calculates the mantle adiabat temperature for given pressure
	range. Taken from Katsura (2022).

	###Parameters###
	P_input: Pressure in GPa
	out_type: 'T' for temperature, 'grad' for dT/dz gradient in K/km

	'''

	Depth = [50,70,90,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,
	400,410,410,420,440,460,480,500]
	P = [1.5,2.1,2.8,3.1,3.8,4.4,5.1,5.8,6.4,7.1,7.8,8.5,9.2,9.9,10.6,11.2,
	11.9,12.6,13.4,13.7,13.7,14.1,14.9,15.6,16.4,17.1]
	T = [1646,1657,1667,1672,1682,1691,1700,1709,1718,1726,1735,1743,1751,
	1759,1766,1774,1781,1788,1796,1799,1860,1863,1871,1878,1885,1892]
	grad = [0.54,0.53,0.51,0.5,0.48,0.47,0.46,0.44,0.43,0.42,0.41,0.41,
	0.4,0.39,0.38,0.37,0.36,0.36,0.37,0.36,0.36,0.36,0.36,0.36,0.35,0.35]

	Depth = np.array(Depth)
	T_C = np.array(T) - 273.15

	if out_type == 'T':
		T_C_out = np.interp(P_input ,P ,T_C)
		T_K_out = T_C_out + 273.15
		
		return T_C_out, T_K_out
		
	elif out_type == 'grad':
	
		# dT/dz in K/km from the tabulated values
		dT = np.diff(T)
		dz = np.diff(Depth)
		grad = dT / dz  # K/km at midpoints

		# Pressures at midpoints
		P_mid = (np.array(P[:-1]) + np.array(P[1:])) / 2.0

		# Interpolate gradient onto input pressures
		grad_out = np.interp(P_input, P_mid, grad)

		return grad_out

#k = (k0 + k1/T + k2 T^2)(1 + k3 P), T in K, P in GPa, k in mW/m/K
K_MANTLE_SP = (2.271e3, 681.12e3, -1.259e-4, 0.0399)
K_MANTLE_GT = (2.371e3, 669.4e3, -1.288e-4, 0.0384)

A_MANTLE = 4e-5      #mW/m^3
RHO_CRUST = 2850.0
RHO_MANTLE = 3300.0
G = 9806.0
KELVIN = 273.0


def _k_mantle(T_K, P_GPa):
	"""Hasterok (2010) mantle conductivity, mW/m/K, with the sp-gt switch."""

	T_K = np.asarray(T_K, dtype = float)
	P = np.asarray(P_GPa, dtype = float)

	def f(c):
		return (c[0] + (c[1] / T_K) + (c[2] * T_K ** 2.0)) * (1.0 + c[3] * P)

	p_transition = 1.4209 + np.exp((3.9073e-3 * T_K) - 6.8041)

	return np.where(P >= p_transition, f(K_MANTLE_GT), f(K_MANTLE_SP))


def _pressure(depth, moho):
	"""GPa. Matches the cumulative sum in the Hasterok function."""

	z = np.asarray(depth, dtype = float)

	return (RHO_CRUST * np.minimum(z, moho)
		+ RHO_MANTLE * np.maximum(z - moho, 0.0)) * G * 1e-9


def calculate_conductive_geotherm(T_moho, q_moho, moho = 40.0, T_0 = 0.0,
	max_depth = 300.0, dz = 1.0, A_mantle = A_MANTLE, k_scale = 1.0,
	thermal_lab = True, thermal_lab_temp = 1300.0, adiabat_gradient = 0.4):
	"""
	Steady state conductive geotherm anchored at the Moho. Stage one.

	Solves dT/dz = q/k and dq/dz = -A_mantle downward from T_moho with the
	Hasterok mantle conductivity, switches to an adiabat at the thermal LAB.

	The crust is a straight line to the anchor and carries no information.

	Parameters
	----------
	T_moho : degrees C
	q_moho : mW/m^2
	moho, max_depth, dz : km
	A_mantle : mW/m^3
	k_scale : multiplier on the mantle conductivity. Very nearly degenerate
		with 1/q_moho. Only A_mantle separates them, and weakly: refitting
		T_moho and q_moho to hold the 50 and 120 km temperatures fixed, a
		k_scale of 0.85 to 1.15 moves the middle of the column by at most
		1.4 K, and even a factor of two moves it 8 K. Against a 100 K bow
		that is nothing. Use it for a one-off sensitivity test, never as a
		free parameter.
	thermal_lab, thermal_lab_temp, adiabat_gradient : the LAB switch

	Returns
	-------
	T : KELVIN
	depth : km
	P : GPa
	idx_lab : index of the thermal LAB, or -1

	Notes
	-----
	This family has TWO real degrees of freedom, level and slope, and no
	useful curvature freedom. That is the point: it is a stiff model, so its
	failure is informative. Neither candidate third parameter helps. k_scale
	buys 1.4 K of curvature over its plausible range, see above. A_mantle
	cannot do it either, because a layer pinned at both ends bows by
	A L^2 / 8k, so 100 K over 70 km needs A near 0.49 microW/m^3, twelve
	times depleted peridotite.

	Unlike the gradient geotherm, this one does return a LAB, so lab_depth
	can be fitted in stage one.
	"""

	if q_moho <= 0.0:
		raise ValueError('q_moho must be positive.')
	if moho >= max_depth:
		raise ValueError('moho must be shallower than max_depth.')

	depth = np.arange(0.0, max_depth + 0.5 * dz, dz)
	P = _pressure(depth, moho)
	T = np.zeros(len(depth), dtype = float)

	T_moho_K = T_moho + KELVIN
	T_0_K = T_0 + KELVIN

	in_crust = depth <= moho
	T[in_crust] = T_0_K + (T_moho_K - T_0_K) * depth[in_crust] / moho

	i_start = int(np.argmax(depth > moho))
	if i_start == 0:
		return T, depth, P, -1

	T_cur = float(T_moho_K)
	q = float(q_moho)
	z_cur = float(moho)

	idx_lab = -1
	T_lab_K = thermal_lab_temp + 273.15

	for i in range(i_start, len(depth)):

		h = depth[i] - z_cur

		#q in mW/m^2 over k in mW/m/K gives K/m, so scale by 1e3 for K/km.
		#A in mW/m^3 times km times 1e3 gives mW/m^2.
		k1 = k_scale * float(_k_mantle(T_cur, P[i - 1]))
		T_mid = T_cur + (q / k1 * 1e3) * 0.5 * h
		k2 = k_scale * float(_k_mantle(T_mid, 0.5 * (P[i - 1] + P[i])))

		T_next = T_cur + ((q - A_mantle * 0.5 * h * 1e3) / k2 * 1e3) * h
		q_next = q - A_mantle * h * 1e3

		if thermal_lab and T_next >= T_lab_K:
			#land the LAB on the isotherm, then follow the adiabat
			frac = (T_lab_K - T_cur) / (T_next - T_cur) if T_next > T_cur else 0.0
			z_lab = z_cur + frac * h
			idx_lab = i
			T[i:] = T_lab_K + adiabat_gradient * (depth[i:] - z_lab)
			break

		T[i] = T_next
		T_cur, q, z_cur = T_next, q_next, depth[i]

		if q <= 0.0:
			T[i:] = T_cur
			break

	return T, depth, P, idx_lab


def calculate_gradient_geotherm(T_moho, gradients, gradient_depths,
	moho = 40.0, T_0 = 0.0, max_depth = 300.0, dz = 1.0,
	interp_gradients = True):
	"""
	Geotherm from explicit interval gradients below the Moho.

	Interval i runs from gradient_depths[i-1] to gradient_depths[i], with the
	Moho as the lower edge of the first. Positive gradients make the profile
	monotonic automatically, so nothing is ever rejected for an inverted
	geotherm.

	No LAB and no adiabat are imposed. With the gradients bounded down to
	0.3 K/km the adiabat is reachable as a special case, and the LAB comes
	out afterwards as the shallowest depth where the sampled gradient falls
	below about 1 K/km. Imposing a cap instead would make every gradient
	below it inert, with a flat posterior.

	Parameters
	----------
	T_moho : degrees C
	gradients : K/km, one per interval
	gradient_depths : km, strictly increasing, all deeper than moho
	moho : km
	T_0 : degrees C at the surface, for the placeholder crustal line only
	max_depth : km
	dz : km
	interp_gradients : True for a smooth profile, False for kinks at nodes

	Returns
	-------
	T : KELVIN
	depth : km
	P : GPa
	idx_lab : index where the gradient first falls below 1 K/km, else -1
	"""

	grad = np.asarray(gradients, dtype = float)
	gd = np.asarray(gradient_depths, dtype = float)

	if len(grad) != len(gd):
		raise ValueError('gradients and gradient_depths must match in length, '
			'got %d and %d.' % (len(grad), len(gd)))
	if np.any(np.diff(gd) <= 0.0):
		raise ValueError('gradient_depths must be strictly increasing.')
	if gd[0] <= moho:
		raise ValueError('every gradient depth must be deeper than the Moho.')
	if np.any(grad <= 0.0):
		raise ValueError('gradients must be positive.')

	depth = np.arange(0.0, max_depth + 0.5 * dz, dz)
	P = _pressure(depth, moho)
	T = np.zeros(len(depth), dtype = float)

	T_moho_K = T_moho + KELVIN
	T_0_K = T_0 + KELVIN

	#placeholder crustal line, so the array starts at the surface
	in_crust = depth <= moho
	T[in_crust] = T_0_K + (T_moho_K - T_0_K) * depth[in_crust] / moho

	#gradient as a function of depth
	edges = np.concatenate(([moho], gd))
	if interp_gradients:
		mids = 0.5 * (edges[:-1] + edges[1:])
		g_of_z = np.interp(depth, mids, grad, left = grad[0], right = grad[-1])
	else:
		idx = np.clip(np.searchsorted(edges, depth, side = 'right') - 1,
			0, len(grad) - 1)
		g_of_z = grad[idx]

	#cumulative integral from the anchor
	i_start = int(np.argmax(depth > moho))
	if i_start == 0:
		return T, depth, P, -1

	T[i_start:] = T_moho_K + np.cumsum(
		0.5 * (g_of_z[i_start - 1:-1] + g_of_z[i_start:])
		* np.diff(depth[i_start - 1:]))

	#derived LAB
	idx_lab = -1
	below = np.where((depth > moho) & (g_of_z < 1.0))[0]
	if len(below) > 0:
		idx_lab = int(below[0])

	return T, depth, P, idx_lab


def conductive_gradients_at(depths, T_moho, q_moho, moho = 40.0,
	A_mantle = A_MANTLE, dz = 1.0):
	"""
	The interval gradients a conductive solution gives, for a prior mean.

	Integrates dT/dz = q/k and dq/dz = -A_mantle downward from the anchor
	with the Hasterok mantle conductivity, then averages over each interval
	of [moho, depths[0], depths[1], ...].

	Returns
	-------
	gradients : K/km, len(depths) values
	"""

	d = np.asarray(depths, dtype = float)

	#LAB off, so the conductive branch runs the whole way and the gradients
	#are the pure conductive ones at every node
	T, depth, _P, _i = calculate_conductive_geotherm(T_moho, q_moho,
		moho = moho, max_depth = float(d[-1]) + 10.0, dz = dz,
		A_mantle = A_mantle, thermal_lab = False)

	edges = np.concatenate(([moho], d))
	T_edges = np.interp(edges, depth, T)

	return np.diff(T_edges) / np.diff(edges)