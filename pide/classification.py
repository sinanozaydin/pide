import numpy as np
import pide
from pide.material import Material
from pide.geodyn.mantlemelting import katz_2003
from pide.inversion import conductivity_solver_single_param
import multiprocessing
import os
from functools import partial
	
def _katz2003_wet_parallel(index, max_water, d_per_melt, cpx_frac , P, T, f_dummy):
	
	"A function calling katz2003 peridotite melting scheme for detecting if melt can exist in the given node [index] in a parallel manner."

	bw_list = (max_water[index] * (f_dummy + ((1-f_dummy) * d_per_melt[index]))) / d_per_melt[index] #creating a bulk water list from maximum solubility of peridotite and varying melt content
	
	for j in range(len(bw_list)):
		F_calced = katz_2003.F_wet(P = P[index], T = T[index]-273.15, X = bw_list[j]/1e4, D=d_per_melt[index], M = cpx_frac[index]) #calculating if melt exist with the given P,T,bulk water
		if np.mean(F_calced) > 0.0:
			break
		else:
			pass
	return np.mean(F_calced) > 0.0

def mantle_classification(cond, T, P, Vp = None, Vs = None, material = None,num_cpu = 1,**kwargs):

	import ipdb
	import time
	
	print('ENTERED')
	if material is None:
		material = Material(composition={"ol":0.6, "opx":0.3, "cpx": 0.1, "garnet":0.0}) 
	if (Vp is not None) and (Vs is not None):
		if len(cond) != len(T) != len(P):
			raise ValueError('The length of conductivity, pressure, temperature arrays are not equal!')
	else:
		if len(cond) != len(T) != len(P) != len(Vs) != len(Vp):
			raise ValueError('The length of conductivity, pressure, temperature, Vs and Vp arrays are not equal!')
			
	cond_uncert = kwargs.pop("cond_uncert",np.ones(len(T)) * 0.1)
	vp_uncert = kwargs.pop("vp_uncert",np.ones(len(T)) * 0.1)
	vs_uncert = kwargs.pop("vs_uncert",np.ones(len(T)) * 0.1)
	depleted_water_cutoff = kwargs.pop("depleted_water_cutoff", 30)
	
	T = np.array(T)
	P = np.array(P)
	cond = np.array(cond)
	if Vp is not None:
		Vp = np.array(Vp)
	if Vs is not None:
		Vs = np.array(Vs)
	
	p_obj = pide.pide()
	p_obj.set_temperature(T)
	p_obj.set_pressure(P)
	#setting the composition ol-opx-cpx-garnet
	p_obj.set_composition_solid_mineral(ol = material.composition["ol"], opx = material.composition["opx"],
	cpx = material.composition["cpx"], garnet = material.composition["garnet"])
	#setting the mantle water solubility options
	p_obj.set_mantle_water_solubility(ol = material.mantle_water_solubility["ol"], opx = material.mantle_water_solubility["opx"],
	cpx = material.mantle_water_solubility["cpx"], garnet = material.mantle_water_solubility["garnet"])
	#setting the mantle water partitionings
	p_obj.set_mantle_water_partitions(opx_ol = material.mantle_water_part["opx_ol"], cpx_ol = material.mantle_water_part["cpx_ol"],
	garnet_ol = material.mantle_water_part["garnet_ol"], ol_melt = material.mantle_water_part["ol_melt"],
	opx_melt = material.mantle_water_part["opx_melt"], cpx_melt = material.mantle_water_part["cpx_melt"],
	garnet_melt = material.mantle_water_part["garnet_melt"])
	
	p_obj.set_solid_phs_mix_method(1) #H-S lower bound
	if material.mantle_water_solubility["ol"] == 4:
		p_obj.set_parameter('ti_ol', 0.01)
	
	mask_unclassified = np.array([True] * len(T))
	mantle_class = np.zeros(len(T))
	melt_bool_dry = np.array([False] * len(T))
	melt_bool_wet = np.array([False] * len(T))
	
	#Step 1: 
	print("STEP 1 Started")
	
	#Checking if is there any melting without any addition of volatiles after Katz et al. (2003)
	dry_melt_frac = katz_2003.F_dry(P=P, T=T-273.15, M=p_obj.cpx_frac)
	melt_bool_dry = ~np.isnan(dry_melt_frac)
	
	
	#calculating max water content to check whether it will depress the solidus enough to enter the melting domain
	max_water = p_obj.calculate_bulk_mantle_water_solubility()
	d_per_melt = (p_obj.ol_frac * p_obj.d_melt_ol) +\
					(p_obj.opx_frac * p_obj.d_melt_opx) +\
					(p_obj.cpx_frac * p_obj.d_melt_cpx) +\
					(p_obj.garnet_frac * p_obj.d_melt_garnet)
	
	#creating melt fraction form 0.5 to 0, to calculate dummy. Reverse direction to save time in processing.
	f_dummy = np.arange(0.5,0,-0.01)

	#Calling katz2003 for peridotite melting detection.
	if num_cpu == 1:
		for i in range(0,len(max_water)):
			bw_list = (max_water[i] * (f_dummy + ((1-f_dummy) * d_per_melt[i]))) / d_per_melt[i] #creating a bulk water list from maximum solubility of peridotite and varying melt content
			for j in range(len(bw_list)):
				F_calced = katz_2003.F_wet(P = P[i], T = T[i]-273.15, X = bw_list[j]/1e4, D=d_per_melt[i], M = p_obj.cpx_frac[i]) #calculating if melt exist with the given P,T,bulk water
				if np.mean(F_calced) > 0.0:
					melt_bool_wet[i] = True
					break
				else:
					pass
	else:
		index_list = np.array(list(range(0,len(max_water))))
		with multiprocessing.Pool(processes=num_cpu) as pool:
			process_item_partial = partial(_katz2003_wet_parallel, max_water = max_water, 
			d_per_melt = d_per_melt, cpx_frac = p_obj.cpx_frac , P = P, T = T, f_dummy = f_dummy)
			
			c_sol = pool.map(process_item_partial, index_list)
		melt_bool_wet = np.array(c_sol)
	

	
	#assigning all possible melting areas from temperature only.
	melt_bool = melt_bool_dry + melt_bool_wet
	print('STEP 1 Finished...')
	#Step 2:
	#Binning every bit where observed conductivity is below the conductivity of the depleted mantle.
	cond_dry_mantle = p_obj.calculate_conductivity()
	diff_mantle = cond - cond_dry_mantle

	crit_0, = np.where(diff_mantle <= 0.0)

	#assigning to class 1 - Dry/Depleted Peridotite
	mantle_class[crit_0] = 1
	mask_unclassified[mantle_class != 0] = False
	print("STEP 2 Finished...")
	#looking from the difference, removing mantle class = 1 from the possbiel melt search areas:
	melt_bool[mantle_class == 1] = False
	
	#STEP 3:
	#Performing possible existence of melt
	pass
	print("STEP 3 Finished...")
	#STEP_4:
	#Defining the unclassified nodes
	T_step_4 = T[mask_unclassified].copy()
	P_step_4 = P[mask_unclassified].copy()
	
	p_obj = pide.pide()
	p_obj.set_temperature(T_step_4)
	p_obj.set_pressure(P_step_4)
	
	#setting the composition ol-opx-cpx-garnet
	p_obj.set_composition_solid_mineral(ol = material.composition["ol"], opx = material.composition["opx"],
	cpx = material.composition["cpx"], garnet = material.composition["garnet"])
	#setting the mantle water solubility options
	p_obj.set_mantle_water_solubility(ol = material.mantle_water_solubility["ol"], opx = material.mantle_water_solubility["opx"],
	cpx = material.mantle_water_solubility["cpx"], garnet = material.mantle_water_solubility["garnet"])
	#setting the mantle water partitionings
	p_obj.set_mantle_water_partitions(opx_ol = material.mantle_water_part["opx_ol"], cpx_ol = material.mantle_water_part["cpx_ol"],
	garnet_ol = material.mantle_water_part["garnet_ol"], ol_melt = material.mantle_water_part["ol_melt"],
	opx_melt = material.mantle_water_part["opx_melt"], cpx_melt = material.mantle_water_part["cpx_melt"],
	garnet_melt = material.mantle_water_part["garnet_melt"])
	
	p_obj.set_solid_phs_mix_method(1) #H-S lower bound
	print('STEP 4 Inversion Started...')
	water_solv = conductivity_solver_single_param(object = p_obj, cond_list = cond[mask_unclassified],param_name = 'bulk_water',
	upper_limit_list=max_water[mask_unclassified], lower_limit_list=np.zeros(len(T_step_4)),acceptence_threshold=0.1,num_cpu = num_cpu)

	
	import ipdb
	ipdb.set_trace()