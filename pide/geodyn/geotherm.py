import numpy as np

def calculate_hasterok2011_geotherm(SHF,  T_0, max_depth, moho, adiabat=True, BDL_T = 0,kinked = False, **kwargs):

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

	g = 9806.0
	rho_crust = 2850.0
	rho_mantle = 3340.0
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
	
def T_Katsura_2022_Adiabat(P_input):

	'''
	A function that calculates the mantle adiabat temperature for given pressure
	range. Taken from Katsura (2022).

	###Parameters###
	P_input: Pressure in GPa

	'''

	Depth = [50,70,90,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,
	400,410,410,420,440,460,480,500]
	P = [1.5,2.1,2.8,3.1,3.8,4.4,5.1,5.8,6.4,7.1,7.8,8.5,9.2,9.9,10.6,11.2,
	11.9,12.6,13.4,13.7,13.7,14.1,14.9,15.6,16.4,17.1]
	T = [1646,1657,1667,1672,1682,1691,1700,1709,1718,1726,1735,1743,1751,
	1759,1766,1774,1781,1788,1796,1799,1860,1863,1871,1878,1885,1892]

	Depth = np.array(Depth)
	T_C = np.array(T) - 273.15

	T_C_out = np.interp(P_input ,P ,T_C)
	T_K_out = T_C_out + 273.15

	return T_C_out, T_K_out