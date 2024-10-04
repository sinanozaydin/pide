#!/usr/bin/env python3

import numpy as np

from .pide import pide
from .geodyn.material_process import return_material_bool

#importing the function
from .geodyn.deform_cond import plastic_strain_2_conductivity
from .utils.utils import text_color, sort_through_external_list, check_type

def run_model(index_list, material, pide_object, t_array, p_array, melt_array, type = 'conductivity'):

	"""
	A global function to run conductivity model of a Model object. Designed to be global def to run parallel
	with multiprocessing library. The users are not encouraged to use this function directly. 
	
	Input:
	array: index_list - list of indexes associated with material coordinates 
	object: material - material object associated with the indexes.
	object: pide_object - pide object used in calculation.
	array: t_array - temperature array || in K.
	array: p_array - pressure array || in GPa.
	array: melt_array - array of melt content || in fraction.
	str: type - type of calculation. || 'conductivity' or 'seismic' or 'both'.
	
	Output:
	Depending on the 'type':
	array: conductivity || in S/m.
	array: seismic velocity || in km/s.
	list of arrays: [conductivity, seismic_velocity]
	
	Example:
	cond = run_model(index_list= sliced_material_idx, material = material_list_holder[l][i], pide_object = mat_pide_obj,
							t_array = self.T, p_array=self.P, melt_array=self.melt_frac)
	"""

	#global function to run conductivity model. designed to be global def to run parallel with multiprocessing
	
	#setting temperatures at the pide_object
	pide_object.set_temperature(t_array[index_list])
	pide_object.set_pressure(p_array[index_list])
	
	pide_object.set_o2_buffer(material.o2_buffer)
	pide_object.set_solid_phs_mix_method(material.solid_phase_mixing_idx)
	pide_object.set_solid_melt_fluid_mix_method(material.melt_fluid_phase_mixing_idx)
	
	if melt_array[index_list].any() > 0.0:
		material.melt_fluid_incorporation_method == 'field' #if melt exists it overwrites the field value.
		
	if material.melt_fluid_incorporation_method == 'value':

		if material.melt_or_fluid == 'melt':
			pide_object.set_melt_or_fluid_mode('melt')
		else:
			pide_object.set_melt_or_fluid_mode('fluid')

		pide_object.set_melt_fluid_frac(material.melt_fluid_frac)

	elif material.melt_fluid_incorporation_method == 'field':

		pide_object.set_melt_or_fluid_mode('melt') #only melt field can be taken from the area
		pide_object.set_melt_fluid_frac(melt_array[index_list])

	else:
		pass
		
	pide_object.set_watercalib(ol = material.water_calib['ol'], px_gt = material.water_calib['px_gt'], feldspar = material.water_calib['feldspar'])
	pide_object.set_o2_buffer(o2_buffer = material.o2_buffer)
	pide_object.set_solid_phase_method(material.calculation_type)
	pide_object.set_alopx(material.al_opx)

	#adjusting material parameters for the pide_object
	if material.calculation_type == 'mineral':
	
		pide_object.set_composition_solid_mineral(ol = material.composition['ol'],
		opx = material.composition['opx'],
		cpx = material.composition['cpx'],
		garnet = material.composition['garnet'],
		mica = material.composition['mica'],
		amp = material.composition['amp'],
		quartz = material.composition['quartz'],
		plag = material.composition['plag'],
		kfelds = material.composition['kfelds'],
		sulphide = material.composition['sulphide'],
		graphite = material.composition['graphite'],
		mixture = material.composition['mixture'],
		sp = material.composition['sp'],
		wds = material.composition['wds'],
		rwd = material.composition['rwd'],
		perov = material.composition['perov'],
		other = material.composition['other'])
		
		pide_object.set_param1_mineral(ol = material.param1['ol'],
		opx = material.param1['opx'],
		cpx = material.param1['cpx'],
		garnet = material.param1['garnet'],
		mica = material.param1['mica'],
		amp = material.param1['amp'],
		quartz = material.param1['quartz'],
		plag = material.param1['plag'],
		kfelds = material.param1['kfelds'],
		sulphide = material.param1['sulphide'],
		graphite = material.param1['graphite'],
		mixture = material.param1['mixture'],
		sp = material.param1['sp'],
		wds = material.param1['wds'],
		rwd = material.param1['rwd'],
		perov = material.param1['perov'],
		other = material.param1['other'])
		
		if material.solid_phase_mixing_idx == 0:
		
			pide_object.set_phase_interconnectivities(ol = material.interconnectivities['ol'],
			opx = material.interconnectivities['opx'],
			cpx = material.interconnectivities['cpx'],
			garnet = material.interconnectivities['garnet'],
			mica = material.interconnectivities['mica'],
			amp = material.interconnectivities['amp'],
			quartz = material.interconnectivities['quartz'],
			plag = material.interconnectivities['plag'],
			kfelds = material.interconnectivities['kfelds'],
			sulphide = material.interconnectivities['sulphide'],
			graphite = material.interconnectivities['graphite'],
			mixture = material.interconnectivities['mixture'],
			sp = material.interconnectivities['sp'],
			wds = material.interconnectivities['wds'],
			rwd = material.interconnectivities['rwd'],
			perov = material.interconnectivities['perov'],
			other = material.interconnectivities['other'],
			)
			
		if material.water_distr == False:
		
			pide_object.set_mineral_water(ol = material.water['ol'],
			opx = material.water['opx'],
			cpx = material.water['cpx'],
			garnet = material.water['garnet'],
			mica = material.water['mica'],
			amp = material.water['amp'],
			quartz = material.water['quartz'],
			plag = material.water['plag'],
			kfelds = material.water['kfelds'],
			sulphide = material.water['sulphide'],
			graphite = material.water['graphite'],
			mixture = material.water['mixture'],
			sp = material.water['sp'],
			wds = material.water['wds'],
			rwd = material.water['rwd'],
			perov = material.water['perov'],
			other = material.water['other'])
			
		else:
		
			pide_object.set_bulk_water(material.water['bulk'])
			pide_object.set_mantle_water_partitions(opx_ol = material.mantle_water_part['opx_ol'],
			cpx_ol = material.mantle_water_part['cpx_ol'],
			garnet_ol = material.mantle_water_part['garnet_ol'],
			ol_melt = material.mantle_water_part['ol_melt'],
			opx_melt = material.mantle_water_part['opx_melt'],
			cpx_melt = material.mantle_water_part['cpx_melt'],
			garnet_melt = material.mantle_water_part['garnet_melt'])
			pide_object.mantle_water_distribute(method = 'array')
			
		pide_object.set_mineral_conductivity_choice(ol = material.el_cond_selections['ol'],
			opx = material.el_cond_selections['opx'],
			cpx = material.el_cond_selections['cpx'],
			garnet = material.el_cond_selections['garnet'],
			mica = material.el_cond_selections['mica'],
			amp = material.el_cond_selections['amp'],
			quartz = material.el_cond_selections['quartz'],
			plag = material.el_cond_selections['plag'],
			kfelds = material.el_cond_selections['kfelds'],
			sulphide = material.el_cond_selections['sulphide'],
			graphite = material.el_cond_selections['graphite'],
			mixture = material.el_cond_selections['mixture'],
			sp = material.el_cond_selections['sp'],
			wds = material.el_cond_selections['wds'],
			rwd = material.el_cond_selections['rwd'],
			perov = material.el_cond_selections['perov'],
			other = material.el_cond_selections['other'])
									
	elif material.calculation_type == 'rock':
	
		pide_object.set_composition_solid_rock(granite = material.composition['granite'],
		granulite = material.composition['granulite'],
		sandstone = material.composition['sandstone'],
		gneiss = material.composition['gneiss'],
		amphibolite = material.composition['amphibolite'],
		basalt = material.composition['basalt'],
		mud = material.composition['mud'],
		gabbro = material.composition['gabbro'],
		other_rock = material.composition['other_rock'])
		
		pide_object.set_param1_rock(granite = material.param1['granite'],
		granulite = material.param1['granulite'],
		sandstone = material.param1['sandstone'],
		gneiss = material.param1['gneiss'],
		amphibolite = material.param1['amphibolite'],
		basalt = material.param1['basalt'],
		mud = material.param1['mud'],
		gabbro = material.param1['gabbro'],
		other_rock = material.param1['other_rock'])
		
		if material.solid_phase_mixing_idx == 0:
			
			pide_object.set_phase_interconnectivities(granite = material.interconnectivities['granite'],
			granulite = material.interconnectivities['granulite'],
			sandstone = material.interconnectivities['sandstone'],
			gneiss = material.interconnectivities['gneiss'],
			amphibolite = material.interconnectivities['amphibolite'],
			basalt = material.interconnectivities['basalt'],
			mud = material.interconnectivities['mud'],
			gabbro = material.interconnectivities['gabbro'],
			other_rock = material.interconnectivities['other_rock'])
			
		if material.water_distr == False:
		
			pide_object.set_rock_water(granite = material.water['granite'],
			granulite = material.water['granulite'],
			sandstone = material.water['sandstone'],
			gneiss = material.water['gneiss'],
			amphibolite = material.water['amphibolite'],
			basalt = material.water['basalt'],
			mud = material.water['mud'],
			gabbro = material.water['gabbro'],
			other_rock = material.water['other_rock'])
			
		pide_object.set_rock_conductivity_choice(granite = material.el_cond_selections['granite'],
			granulite = material.el_cond_selections['granulite'],
			sandstone = material.el_cond_selections['sandstone'],
			gneiss = material.el_cond_selections['gneiss'],
			amphibolite = material.el_cond_selections['amphibolite'],
			basalt = material.el_cond_selections['basalt'],
			mud = material.el_cond_selections['mud'],
			gabbro = material.el_cond_selections['gabbro'],
			other_rock = material.el_cond_selections['other_rock'])
	
	if material.melt_fluid_cond_selection != None:
		if material.melt_or_fluid == 'melt':
			
			pide_object.set_melt_fluid_conductivity_choice(melt = material.melt_fluid_cond_selection)
			if material.melt_fluid_phase_mixing_idx == 0:
				pide_object.set_melt_fluid_interconnectivity(material.melt_fluid_m)
		elif material.melt_or_fluid == 'fluid':
			pide_object.set_melt_fluid_conductivity_choice(fluid = material.melt_fluid_cond_selection)
			pide_object.set_fluid_properties(salinity = material.fluid_salinity)
			if material.melt_fluid_phase_mixing_idx == 0:
				pide_object.set_melt_fluid_interconnectivity(material.melt_fluid_m)
	
	if type == 'conductivity':
		c = pide_object.calculate_conductivity(method = 'array')
		return c
	elif type == 'seismic':
		c = pide_object.calculate_seismic_velocities(method = 'array')
		return c
	elif type == 'both':
		c = pide_object.calculate_conductivity(method = 'array')
		s = pide_object.calculate_seismic_velocities(method = 'array')
		return c, s
		
	else:
		raise NameError('The type for "run_model" function entered wrongly...It has to be either "conductivity", "seismic" or "both".')
	
def run_deform2cond(index_number,p_strain, background_cond, max_cond, low_deformation_threshold, high_deformation_threshold, function_method,
	conductivity_decay_factor, conductivity_decay_factor_2, strain_decay_factor, model_type):

	"""
	A function to run deformation related conductivity.

	Inputs:
	index_number: array of index numbers.
	p_strain: Plastic strain field  - np.ndarray
	background_cond: lower conductivity to associated with low plastic strain - np.ndarray
	max_cond: highest conductivity to associate with high plastic strain - np.ndarray
	low_deformation_threshold: Threshold point in plastic strain to start the linking function - float
	high_deformation_threshold: Threshold point in plastic strain to end the linking function - float
	function_method: The type of function to fit: 'exponential', 'linear, 'logarithmic'.
	conductivity_decay_factor: How much conductivity mid-point is moved from the middle point: -1 to 1
	strain_decay_factor:
	model_type:

	To understand how this function works. Try to use the notebook deformation_related_conductivity notebook in examples.

	"""
	if model_type == 'underworld_2d':
		#Global function to link deformation and conductivity.
		if (background_cond[index_number[0],index_number[1]] == np.nan) or (max_cond[index_number[0],index_number[1]] == np.nan): 
		
			c = np.nan
			rms = np.nan
		
		else:
		
			if p_strain[index_number[0],index_number[1]] <= low_deformation_threshold:
				
				c = background_cond[index_number[0],index_number[1]]
				rms = np.nan
				
			elif p_strain[index_number[0],index_number[1]] >= high_deformation_threshold:
				
				c = max_cond[index_number[0],index_number[1]]
				rms = np.nan
				
			else:
				
				c, rms = plastic_strain_2_conductivity(strain = p_strain[index_number[0],index_number[1]],low_cond = background_cond[index_number[0],index_number[1]],
					high_cond=max_cond[index_number[0],index_number[1]],low_strain=low_deformation_threshold, high_strain=high_deformation_threshold,
					function_method = function_method, conductivity_decay_factor = conductivity_decay_factor, 
					conductivity_decay_factor_2 = conductivity_decay_factor_2,strain_decay_factor = strain_decay_factor, return_all_params = True)
	
	elif model_type == 'underworld_3d':
	
		if (background_cond[index_number] == np.nan) or (max_cond[index_number] == np.nan): 
		
			c = np.nan
			rms = np.nan
		
		else:
		
			if p_strain[index_number] <= low_deformation_threshold:
				
				c = background_cond[index_number]
				rms = 0.0
				
			elif p_strain[index_number] >= high_deformation_threshold:

				c = max_cond[index_number]
				rms = 0.0
				
			else:
				
				c, rms = plastic_strain_2_conductivity(strain = p_strain[index_number],low_cond = background_cond[index_number],
					high_cond=max_cond[index_number],low_strain=low_deformation_threshold, high_strain=high_deformation_threshold,
					function_method = function_method, conductivity_decay_factor = conductivity_decay_factor,
					conductivity_decay_factor_2 = conductivity_decay_factor_2, strain_decay_factor = strain_decay_factor, return_all_params = True, index_plot = index_number)
				
	else:
	
		raise ValueError('The model type can only be underworld_2d or underworld_3d.')

	return c, rms

class Model(object):

	def __init__(self, material_list, material_array = None, T = None, P = None, depth = None, model_type = 'underworld_2d', material_list_2 = None,
	melt = None, p_strain = None, strain_rate = None, material_node_skip_rate_list = None):
		
		"""
		Model object: This is an object used to append materials in a 3D or 2D space, then perform calculations in batch.
		
		Inputs:
		array: material_list - list of pide.Material objects that will be associated with the Model object.
		array: material_array - list of indexes associated with the pide.Material objects.
		array: T - temperature array || in K.
		array: P - pressure array || in GPa.
		array: depth - depth array || in km.
		str: model_type - type of the model to used in calculations. || 'underworld_2d' or 'underworld_3d' or any other string.
																		Only underworld_2d has special indexing conditions.
																		Whatever else is entered for model_type will be the
																		same methodology of single array indexing of underworld_3d.
																		
																		underworld_2d - Special arrays for every parameter (e.g., T,P) that
																		should not be used for other purposes.
																		anything else - 1D numpy arrays.
																		
		array: material_list_2 - list of pide.Material objects that can be used to calculate a connection between the first material list.
		array: melt - melt array || in fraction.
		array: p_strain - plastic strain || in plastic strain rate.
		array: material_node_skip_rate_list -  an array to identify skipping rate for calculation for each pide.Material object listed in
																		material_list. 
		
		
		Methods:
		----------------------------------------------- -----------------------------------------
		Methods                             			Description
		----------------------------------------------- -----------------------------------------

		calculate_model					                Calculates the conductivity or seismic
														velocities of the model object that is setup.
														
		calculate_geothermal_block                      Calculates the conductivity or seismic
														velocities of the model object with
														material objects.
		
		calculate_deformation_related_conductivity		Calculates the deformation related conductivity
														withe the given information loaded onto material
														objects.								

		
		"""

		self.material_list = material_list
		self.material_list_2 = material_list_2
		self.material_array = material_array
		self.T = T
		self.P = P
		self.depth = depth
		self.model_type = model_type
		self.melt_frac = melt
		self.p_strain = p_strain
		self.strain_rate = strain_rate
		self.material_node_skip_rate_list = material_node_skip_rate_list

		if self.T is not None:
			if len(self.T) != len(self.P):
				raise ValueError(text_color.RED + f'The length of T and P arrays do not match!' + text_color.END)
		
	def calculate_model(self,type = 'conductivity', num_cpu = 1):

		"""
		Calculates conductivity for the given 2 or 3D model with given material_list [Material object] and associated material
		array, where indexes of materials are stored. These index arrays has to be same length with T, P and all other comp-
		positional arrays. These materials and indexes can be imported form a thermomechanical model, a geological model or
		anything in particular as long as they are entered in the right format.
		
		Inputs:
		str: type - type of the calculation || 'conductivity' or 'seismic' or 'both'.
		int: num_cpu - number of cpus used in the calculation.
		
		Output:
		Depending on the 'type':
		array: conductivity || in S/m.
		array: seismic velocity || in km/s.
		list of arrays: [conductivity, seismic_velocity]

		"""
		initial_cpu = int(num_cpu)

		if num_cpu > 1:
		
			import multiprocessing
			import os
			from functools import partial
			
			max_num_cores = os.cpu_count()
			
			if num_cpu > max_num_cores:
				raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
		
		acceptable_types = ['conductivity','seismic','both']
		
		if type in acceptable_types:
			pass
		else:
			raise NameError('Invalid entry for type entry.')
			
		cond = np.zeros_like(self.T)
		v_p = np.zeros_like(self.T)
		v_s = np.zeros_like(self.T)
		
		material_list_holder = [self.material_list]
					
		for l in range(0,len(material_list_holder)):

			print(text_color.RED + 'Initiating calculation for the materials appended to the model.' + text_color.END)
			print('##############################################################')
			for i in range(0,len(material_list_holder[l])):
	
				#determining the material indexes from the material array
				
				if self.material_node_skip_rate_list != None:
					if self.material_node_skip_rate_list[i] != None:
						mat_skip = self.material_node_skip_rate_list[i]
					else:
						mat_skip = None
				else:
					mat_skip = None
				
				#getting the relevant indexes with information given for it
				
				material_idx = return_material_bool(material_index = material_list_holder[l][i].material_index,
				model_array = self.material_array, material_skip = mat_skip, model_type = self.model_type)		
				
				#setting up the object for the material
				mat_pide_obj = pide()
				
				if material_list_holder[l][i].calculation_type != 'value':
					mat_pide_obj.set_solid_phase_method(material_list_holder[l][i].calculation_type)
				else:
					num_cpu = 1 #defaulting num_cpu for 1 since it is not needed for value-method
				
				#Slicing the array for parallel calculation
				if num_cpu > 1:
					
					if len(material_idx) == 2: #if clause for 2D underworld model, yes the number is 3 for 2D and 1 for 3D. You do not need to confuse!
						#condition to check if array is too small to parallelize for the material num_cpu*num_cpu 
						
						if len(material_idx[0]) > (num_cpu*num_cpu):
							size_arrays = len(material_idx[0]) // num_cpu
							sliced_material_idx = []
							#adjusting the material_index_array
							for idx in range(0, len(material_idx[0]), size_arrays):
								if idx >= (size_arrays*num_cpu):
									sliced_material_idx.append(tuple((material_idx[0][idx:len(material_idx[0])], material_idx[1][idx:len(material_idx[1])])))
								else:
									sliced_material_idx.append(tuple((material_idx[0][idx:idx+size_arrays], material_idx[1][idx:idx+size_arrays])))
						else:
							#revert back to the single cpu if the material_idx_list is not long enough
							num_cpu = 1
							sliced_material_idx = material_idx
					else: #if clause for 3D underworld model
						if len(material_idx) > (num_cpu*num_cpu):
							size_arrays = len(material_idx) // num_cpu
							sliced_material_idx = []
							for idx in range(0, len(material_idx), size_arrays):
								if idx >= (size_arrays*num_cpu):
									sliced_material_idx.append(material_idx[idx:len(material_idx)])
								else:
									sliced_material_idx.append(material_idx[idx:idx+size_arrays])
						else:	
							#revert back to the single_cpu if the material_idx is not long enough
							num_cpu = 1
							sliced_material_idx = material_idx
						
				#if not parallel material_idx stays the same for both 2-D and 3-D cases.
				elif num_cpu == 1:
					
					sliced_material_idx = material_idx

				if material_list_holder[l][i].calculation_type == 'value':
					#calculation not necessary for value method so automatically not parallel and indexed into sliced_material_idx
					if type == 'conductivity':
						cond[sliced_material_idx] = 1.0 / material_list_holder[l][i].resistivity_medium
					elif type == 'seismic':
						v_p[sliced_material_idx] = material_list_holder[l][i].vp_medium
						v_s[sliced_material_idx] = material_list_holder[l][i].vs_medium
				else:
				
					if num_cpu == 1:
						
						if type == 'conductivity':
							cond[sliced_material_idx] = run_model(index_list= sliced_material_idx, material = material_list_holder[l][i], pide_object = mat_pide_obj,
							t_array = self.T, p_array=self.P, melt_array=self.melt_frac)
						elif type == 'seismic':
							v_bulk, v_p[sliced_material_idx], v_s[sliced_material_idx] = run_model(index_list= sliced_material_idx, material = material_list_holder[l][i], pide_object = mat_pide_obj,
							t_array = self.T, p_array=self.P, melt_array=self.melt_frac, type='seismic') 
						
					else:
						#solving for parallel with multiprocessing
						with multiprocessing.Pool(processes=num_cpu) as pool:
							
							process_item_partial = partial(run_model, material =  material_list_holder[l][i], pide_object = mat_pide_obj, t_array = self.T,
							p_array = self.P, melt_array = self.melt_frac, type = type)
							
							c = pool.map(process_item_partial, sliced_material_idx)
						
						#assigning to the global cond list
						for idx in range(0,len(sliced_material_idx)):
							
							if type == 'conductivity':
							
								cond[sliced_material_idx[idx]] = c[idx]
							
							elif type == 'seismic':
								
								v_p[sliced_material_idx[idx]] = c[idx][1]
								v_s[sliced_material_idx[idx]] = c[idx][2]
				
				print(f'The conductivity for the material  {material_list_holder[l][i].name}  is calculated.')
				num_cpu = initial_cpu
						
			#converting all zero vals in the cond to None values
			if type == 'conductivity':
				cond[cond == 0.0] = np.nan
			elif type == 'seismic':
				v_p[v_p == 0.0] = np.nan
				v_s[v_s == 0.0] = np.nan
				
			print('##############################################################')
				
		if type == 'conductivity':
			return cond
		elif type == 'seismic':
			return v_p, v_s

	def calculate_geothermal_block(self, type = 'conductivity'):
	
		"""
		A function to calculated conductivity and seismic velocities built over material arrays specifically designed for calculations
		along a geotherm.
		
		Input:
		str: type - type of the calculation || 'conductivity' or 'seismic' or 'both'.
		
		Output:
		Depending on the 'type':
		array: conductivity || in S/m.
		array: seismic velocity || in km/s.
		list of arrays: [conductivity, seismic_velocity]

		"""

		if self.depth is None:
			raise ValueError(text_color.RED + 'The depth array of the geotherm has to be entered...')
				
		if self.melt_frac is None:
			self.melt_frac = np.zeros(len(self.depth))
		
		#Getting into layer
		top_bottom_list = []
		for item in self.material_list:
			top_bottom_list.append([item.top,item.bottom])
		top_bottom_list = np.array(top_bottom_list)

		if np.any(top_bottom_list == None):
			raise ValueError(text_color.RED + 'There is None encountered in the material top-bottom values. You have to enter top bottom attributes for each material object.')
		
		#Order the material dependent on their 
		self.material_list = sort_through_external_list(top_bottom_list[:,0],self.material_list)
		layer_end_list = []
		for i in range(1,len(self.material_list)):
			if self.material_list[i].top < self.material_list[i-1].bottom:
				raise ValueError(text_color.RED + f'The depth of the "top" value is found to be over the "bottom value of the consequent layer at layer no:' +  text_color.GREEN + f' {str(i)}' + text_color.END)
			else:
				if i == len(self.material_list):
					layer_end_list.append(len(self.depth))
				else:
					idx_match = (np.abs(self.depth-self.material_list[i-1].bottom)).argmin()
					layer_end_list.append(idx_match)

		layer_end_list.append(len(self.depth))

		sliced_material_idx = [list(range(0,layer_end_list[0]))]
		for idx in range(1,len(layer_end_list)):
			sliced_material_idx.append(list(range(layer_end_list[idx-1],layer_end_list[idx])))

		mat_pide_obj = pide()

		if type == 'conductivity':

			cond = np.zeros(len(self.depth))

			for layer_idx in range(0,len(sliced_material_idx)):
				
				cond[sliced_material_idx[layer_idx]] = run_model(index_list= sliced_material_idx[layer_idx], material = self.material_list[layer_idx],
							pide_object = mat_pide_obj,	t_array = self.T, p_array=self.P, melt_array=self.melt_frac)
			return cond
		elif type == 'seismic':

			vel_bulk = np.zeros(len(self.depth))
			vel_p =  np.zeros(len(self.depth))
			vel_s = np.zeros(len(self.depth))

			for layer_idx in range(0,len(sliced_material_idx)):
				
				vel_bulk[sliced_material_idx[layer_idx]], vel_p[sliced_material_idx[layer_idx]],vel_s[sliced_material_idx[layer_idx]] = \
					run_model(index_list = sliced_material_idx[layer_idx], material = self.material_list[layer_idx],
							pide_object = mat_pide_obj,	t_array = self.T, p_array=self.P, melt_array=self.melt_frac,
							type = 'seismic')
				
			return vel_bulk,vel_p,vel_s
		
		else:
			raise ValueError(text_color.RED + f'The type is entered wrongly. The available type string inputs are: "conductivity"' + text_color.END)
	
	def calculate_deformation_related_conductivity(self, cond_min, cond_max, method = 'plastic_strain',
		low_deformation_threshold = 1e-2, high_deformation_threshold = 100, num_cpu = 1):
		
		"""
		A function to calculate deformation related conductivity via given method parameter, minimum and maximum conductivity arrays.
		The structure of the input arrays will depend on whether the model_type is entered as 'underworld_2d' or not. 
		
		Input:
		array: cond_min - array of conductivity values for minimum conductivities. || in S/m
		array: cond_max - array of conductivity values for maximum conductivities. || in S/m
		str: method - method of linking two conductivites.
		float: low_deformation_threshold - low plastic strain treshold where one can assume the material can be fully associated with cond_min
		float: high_deformation_threshold - high plastic strain treshold where one can assume the material can be fully associated with cond_max
		num_cpu: number of cpus to perform the calculation
		
		Output:
		[deform_cond, rms]
		array: deform_cond - calculated deformation related conductivity || in S/m
		array: rms - local rms misfit array from the fitting of the conductivity-strain relationship curves.
		"""
		
		if num_cpu > 1:
		
			import multiprocessing
			import os
			from functools import partial
			
			max_num_cores = os.cpu_count()
			
			if num_cpu > max_num_cores:
				raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')			
	
		deform_cond = np.zeros_like(self.T)
		rms = np.zeros_like(self.T)
		
		if method == 'plastic_strain':
						
			for i in range(0,len(self.material_list)):
			
				if self.material_node_skip_rate_list != None:
					if self.material_node_skip_rate_list[i] != None:
						mat_skip = self.material_node_skip_rate_list[i]
					else:
						mat_skip = None
				else:
					mat_skip = None
				
				#getting material index for each material
				material_idx = return_material_bool(material_index = self.material_list[i].material_index, model_array = self.material_array,
				material_skip = mat_skip, model_type=self.model_type)

				#turning material index list to be useable format for the np.ndarray fields

				if self.model_type == "underworld_2d":
					material_idx_list = [[material_idx[0][idx],material_idx[1][idx]] for idx in range(0,len(material_idx[0]))]
				else:
					material_idx_list = material_idx

				if num_cpu > 1:
				
					#multiprocessing loop for each material
					with multiprocessing.Pool(processes=num_cpu) as pool:
						
						process_item_partial = partial(run_deform2cond, p_strain = self.p_strain, background_cond = cond_min,
						max_cond = cond_max, low_deformation_threshold = low_deformation_threshold,
						high_deformation_threshold = high_deformation_threshold,
						function_method = self.material_list[i].deformation_dict['function_method'],
						conductivity_decay_factor = self.material_list[i].deformation_dict['conductivity_decay_factor'],
						conductivity_decay_factor_2 = self.material_list[i].deformation_dict['conductivity_decay_factor_2'],
						strain_decay_factor = self.material_list[i].deformation_dict['strain_decay_factor'], model_type = self.model_type)
						
						c = pool.map(process_item_partial, material_idx_list)
					
					deform_cond[material_idx] = [x[0] for x in c]
					rms[material_idx] = [x[1] for x in c]
					
				else:
					
					for idx_ in range(0,len(material_idx_list)):
					
						condd, rmss = run_deform2cond(material_idx_list[idx_],p_strain = self.p_strain, background_cond = cond_min,
						max_cond = cond_max, low_deformation_threshold = low_deformation_threshold,
						high_deformation_threshold = high_deformation_threshold,
						function_method = self.material_list[i].deformation_dict['function_method'],
						conductivity_decay_factor = self.material_list[i].deformation_dict['conductivity_decay_factor'],
						conductivity_decay_factor_2 = self.material_list[i].deformation_dict['conductivity_decay_factor_2'],
						strain_decay_factor = self.material_list[i].deformation_dict['strain_decay_factor'], model_type = self.model_type)
						try:
							if self.model_type == "underworld_3d":
								deform_cond[material_idx_list[idx_]] = condd
								rms[material_idx_list[idx_]] = rmss
							elif self.model_type == "underworld_2d":
								deform_cond[material_idx] = condd
								rms[material_idx] = rmss
							else:
								deform_cond[material_idx_list[idx_]] = condd
								rms[material_idx_list[idx_]] = rmss
						except IndexError:
							import ipdb
							ipdb.set_trace()
										
				print(f'The deformation related conductivity for the material  {self.material_list[i].name}  is calculated.')
			
			#converting all zero vals in the cond to None values
			deform_cond[deform_cond == 0.0] = np.nan
			
		return deform_cond, rms
		




