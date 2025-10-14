import numpy as np
import pide
from pide.material import Material
from pide.geodyn.mantlemelting import katz_2003

def _crit_1(cond_difference, T, P):
	
	perid_class = True
	
	return perid_class

def mantle_classification(cond, T, P, Vp = None, Vs = None, material = None):
	
	if material is None:
		material = Material(composition={"ol":0.6, "opx":0.3, "cpx": 0.1, "garnet":0.0}) 
	if (Vp is not None) and (Vs is not None):
		if len(cond) != len(T) != len(P):
			raise ValueError('The length of conductivity, pressure, temperature arrays are not equal!')
	else:
		if len(cond) != len(T) != len(P) != len(Vs) != len(Vp):
			raise ValueError('The length of conductivity, pressure, temperature, Vs and Vp arrays are not equal!')
			
	p_obj = pide.pide()
	p_obj.set_temperature(T)
	p_obj.set_pressure(P)
	p_obj.set_composition_solid_mineral(ol = material.composition["ol"], opx = material.composition["opx"],
	cpx = material.composition["cpx"], garnet = material.composition["garnet"])
	p_obj.set_solid_phs_mix_method(1)
	
	mask_unclassified = np.array([True] * len(T))
	mantle_class = np.zeros(len(T))
	
	#Crit_0
	cond_dry_mantle = p_obj.calculate_conductivity()
	diff_mantle = cond - cond_dry_mantle

	crit_0, = np.where(diff_mantle <= 0.0)

	#assigning to class 1 - Dry/Depleted Peridotite
	mantle_class[crit_0] = 1
	mask_unclassified[mantle_class != 0] = False
	
	
	
	

	
	import ipdb
	ipdb.set_trace()