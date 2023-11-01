#!/usr/bin/env python3

import numpy as np

import SEL
from geodyn.material_process import return_material_bool

#importing the function
from sel_src.deformation.deform_cond import plastic_strain_2_conductivity

def run_deform2cond(index_number,p_strain, background_cond, max_cond, low_deformation_threshold, high_deformation_threshold, function_method, conductivity_decay_factor,strain_decay_factor):
	
	#Global function to link deformation and conductivity.
	if (background_cond[index_number[0],index_number[1],index_number[2]] == np.nan) or (max_cond[index_number[0],index_number[1],index_number[2]] == np.nan): 
	
		c = np.nan
	
	else:
	
		if p_strain[index_number[0],index_number[1],index_number[2]] <= low_deformation_threshold:
			
			c = background_cond[index_number[0],index_number[1],index_number[2]]
			
		elif p_strain[index_number[0],index_number[1],index_number[2]] >= high_deformation_threshold:
			
			c = max_cond[index_number[0],index_number[1],index_number[2]]
			
		else:
			
			c = plastic_strain_2_conductivity(strain = p_strain[index_number[0],index_number[1],index_number[2]],low_cond = background_cond[index_number[0],index_number[1],index_number[2]],
				high_cond=max_cond[index_number[0],index_number[1],index_number[2]],low_strain=low_deformation_threshold, high_strain=high_deformation_threshold,
				function_method = function_method, conductivity_decay_factor = conductivity_decay_factor, strain_decay_factor = strain_decay_factor)
	
	
	return c

class Model(object):

	def __init__(self, material_list, material_array, T, P, model_type = 'underworld', material_list_2 = None, melt = None, p_strain = None, strain_rate = None, material_node_skip_rate_list = None):

		self.material_list = material_list
		self.material_list_2 = material_list_2
		self.material_array = material_array
		self.T = T
		self.P = P
		self.model_type = model_type
		self.melt_frac = melt
		self.p_strain = p_strain
		self.strain_rate = strain_rate
		self.material_node_skip_rate_list = material_node_skip_rate_list

	def calculate_conductivity(self,type = 'background'):

		cond = np.zeros_like(self.T)
		
		if type == 'background':
			material_list_holder = [self.material_list]
		elif type == 'maximum':
			material_list_holder = [self.material_list_2]
		else:
			material_list_holder = [self.material_list]
			if self.material_list_2 != None:
				material_list_holder.append(self.material_list_2)
		
		cond_list = []
			
		for l in range(0,len(material_list_holder)):
		
			for i in range(0,len(self.material_list)):
	
				#determining the material indexes from the material array
				
				if self.material_node_skip_rate_list != None:
					if self.material_node_skip_rate_list[i] != None:
						mat_skip = self.material_node_skip_rate_list[i]
					else:
						mat_skip = None
				else:
					mat_skip = None
				
				#getting the relevant indexes with information given for it
					
				material_idx = return_material_bool(material_index = material_list_holder[l][i].material_index, model_array = self.material_array, material_skip = mat_skip)		
	
				#getting only the relevant arrays for calculation
				t_relevant = self.T[material_idx]
				p_relevant = self.P[material_idx]
				melt_relevant = self.melt_frac[material_idx]
				
				mat_sel_obj = SEL.SEL()
				mat_sel_obj.set_temperature(t_relevant)
				mat_sel_obj.set_pressure(p_relevant)
				if material_list_holder[l][i].calculation_type != 'value':
					mat_sel_obj.set_solid_phase_method(material_list_holder[l][i].calculation_type)
				
				mat_sel_obj.set_o2_buffer(material_list_holder[l][i].o2_buffer)
				mat_sel_obj.set_melt_or_fluid_mode('melt')
				
				mat_sel_obj.set_melt_fluid_frac(melt_relevant)
				
				if material_list_holder[l][i].calculation_type == 'mineral':
				
					mat_sel_obj.set_composition_solid_mineral(ol = material_list_holder[l][i].composition['ol'],
					opx = material_list_holder[l][i].composition['opx'],
					cpx = material_list_holder[l][i].composition['cpx'],
					garnet = material_list_holder[l][i].composition['garnet'],
					mica = material_list_holder[l][i].composition['mica'],
					amp = material_list_holder[l][i].composition['amp'],
					quartz = material_list_holder[l][i].composition['quartz'],
					plag = material_list_holder[l][i].composition['plag'],
					kfelds = material_list_holder[l][i].composition['kfelds'],
					sulphide = material_list_holder[l][i].composition['sulphide'],
					graphite = material_list_holder[l][i].composition['graphite'],
					mixture = material_list_holder[l][i].composition['mixture'],
					sp = material_list_holder[l][i].composition['sp'],
					wds = material_list_holder[l][i].composition['wds'],
					rwd = material_list_holder[l][i].composition['rwd'],
					perov = material_list_holder[l][i].composition['perov'],
					other = material_list_holder[l][i].composition['other'])
					
					if material_list_holder[l][i].phase_mixing_idx == 0:
					
						mat_sel_obj.set_phase_interconnectivities(ol = material_list_holder[l][i].interconnectivities['ol'],
						opx = material_list_holder[l][i].interconnectivities['opx'],
						cpx = material_list_holder[l][i].interconnectivities['cpx'],
						garnet = material_list_holder[l][i].interconnectivities['garnet'],
						mica = material_list_holder[l][i].interconnectivities['mica'],
						amp = material_list_holder[l][i].interconnectivities['amp'],
						quartz = material_list_holder[l][i].interconnectivities['quartz'],
						plag = material_list_holder[l][i].interconnectivities['plag'],
						kfelds = material_list_holder[l][i].interconnectivities['kfelds'],
						sulphide = material_list_holder[l][i].interconnectivities['sulphide'],
						graphite = material_list_holder[l][i].interconnectivities['graphite'],
						mixture = material_list_holder[l][i].interconnectivities['mixture'],
						sp = material_list_holder[l][i].interconnectivities['sp'],
						wds = material_list_holder[l][i].interconnectivities['wds'],
						rwd = material_list_holder[l][i].interconnectivities['rwd'],
						perov = material_list_holder[l][i].interconnectivities['perov'],
						other = material_list_holder[l][i].interconnectivities['other'])
						
					if material_list_holder[l][i].water_distr == False:
					
						mat_sel_obj.set_mineral_water(ol = material_list_holder[l][i].water['ol'],
						opx = material_list_holder[l][i].water['opx'],
						cpx = material_list_holder[l][i].water['cpx'],
						garnet = material_list_holder[l][i].water['garnet'],
						mica = material_list_holder[l][i].water['mica'],
						amp = material_list_holder[l][i].water['amp'],
						quartz = material_list_holder[l][i].water['quartz'],
						plag = material_list_holder[l][i].water['plag'],
						kfelds = material_list_holder[l][i].water['kfelds'],
						sulphide = material_list_holder[l][i].water['sulphide'],
						graphite = material_list_holder[l][i].water['graphite'],
						mixture = material_list_holder[l][i].water['mixture'],
						sp = material_list_holder[l][i].water['sp'],
						wds = material_list_holder[l][i].water['wds'],
						rwd = material_list_holder[l][i].water['rwd'],
						perov = material_list_holder[l][i].water['perov'],
						other = material_list_holder[l][i].water['other'])
						
					else:
					
						mat_sel_obj.set_bulk_water(material_list_holder[l][i].water['bulk'])
						mat_sel_obj.set_mantle_water_partitions(opx_ol = material_list_holder[l][i].mantle_water_part['opx_ol'],
						cpx_ol = material_list_holder[l][i].mantle_water_part['cpx_ol'],
						garnet_ol = material_list_holder[l][i].mantle_water_part['garnet_ol'],
						ol_melt = material_list_holder[l][i].mantle_water_part['ol_melt'],
						opx_melt = material_list_holder[l][i].mantle_water_part['opx_melt'],
						cpx_melt = material_list_holder[l][i].mantle_water_part['cpx_melt'],
						garnet_melt = material_list_holder[l][i].mantle_water_part['garnet_melt'])
						mat_sel_obj.mantle_water_distribute(method = 'array')
						
					mat_sel_obj.set_mineral_conductivity_choice(ol = material_list_holder[l][i].el_cond_selections['ol'],
						opx = material_list_holder[l][i].el_cond_selections['opx'],
						cpx = material_list_holder[l][i].el_cond_selections['cpx'],
						garnet = material_list_holder[l][i].el_cond_selections['garnet'],
						mica = material_list_holder[l][i].el_cond_selections['mica'],
						amp = material_list_holder[l][i].el_cond_selections['amp'],
						quartz = material_list_holder[l][i].el_cond_selections['quartz'],
						plag = material_list_holder[l][i].el_cond_selections['plag'],
						kfelds = material_list_holder[l][i].el_cond_selections['kfelds'],
						sulphide = material_list_holder[l][i].el_cond_selections['sulphide'],
						graphite = material_list_holder[l][i].el_cond_selections['graphite'],
						mixture = material_list_holder[l][i].el_cond_selections['mixture'],
						sp = material_list_holder[l][i].el_cond_selections['sp'],
						wds = material_list_holder[l][i].el_cond_selections['wds'],
						rwd = material_list_holder[l][i].el_cond_selections['rwd'],
						perov = material_list_holder[l][i].el_cond_selections['perov'],
						other = material_list_holder[l][i].el_cond_selections['other'])
												
				elif material_list_holder[l][i].calculation_type == 'rock':
				
					mat_sel_obj.set_composition_solid_rock(granite = material_list_holder[l][i].composition['granite'],
					granulite = material_list_holder[l][i].composition['granulite'],
					sandstone = material_list_holder[l][i].composition['sandstone'],
					gneiss = material_list_holder[l][i].composition['gneiss'],
					amphibolite = material_list_holder[l][i].composition['amphibolite'],
					basalt = material_list_holder[l][i].composition['basalt'],
					mud = material_list_holder[l][i].composition['mud'],
					gabbro = material_list_holder[l][i].composition['gabbro'],
					other_rock = material_list_holder[l][i].composition['other_rock'])
					
					if material_list_holder[l][i].phase_mixing_idx == 0:
						
						mat_sel_obj.set_phase_interconnectivities(granite = material_list_holder[l][i].interconnectivities['granite'],
						granulite = material_list_holder[l][i].interconnectivities['granulite'],
						sandstone = material_list_holder[l][i].interconnectivities['sandstone'],
						gneiss = material_list_holder[l][i].interconnectivities['gneiss'],
						amphibolite = material_list_holder[l][i].interconnectivities['amphibolite'],
						basalt = material_list_holder[l][i].interconnectivities['basalt'],
						mud = material_list_holder[l][i].interconnectivities['mud'],
						gabbro = material_list_holder[l][i].interconnectivities['gabbro'],
						other_rock = material_list_holder[l][i].interconnectivities['other_rock'])
						
					if material_list_holder[l][i].water_distr == False:
					
						mat_sel_obj.set_rock_water(granite = material_list_holder[l][i].water['granite'],
						granulite = material_list_holder[l][i].water['granulite'],
						sandstone = material_list_holder[l][i].water['sandstone'],
						gneiss = material_list_holder[l][i].water['gneiss'],
						amphibolite = material_list_holder[l][i].water['amphibolite'],
						basalt = material_list_holder[l][i].water['basalt'],
						mud = material_list_holder[l][i].water['mud'],
						gabbro = material_list_holder[l][i].water['gabbro'],
						other_rock = material_list_holder[l][i].water['other_rock'])
						
					mat_sel_obj.set_rock_conductivity_choice(granite = material_list_holder[l][i].el_cond_selections['granite'],
						granulite = material_list_holder[l][i].el_cond_selections['granulite'],
						sandstone = material_list_holder[l][i].el_cond_selections['sandstone'],
						gneiss = material_list_holder[l][i].el_cond_selections['gneiss'],
						amphibolite = material_list_holder[l][i].el_cond_selections['amphibolite'],
						basalt = material_list_holder[l][i].el_cond_selections['basalt'],
						mud = material_list_holder[l][i].el_cond_selections['mud'],
						gabbro = material_list_holder[l][i].el_cond_selections['gabbro'],
						other_rock = material_list_holder[l][i].el_cond_selections['other_rock'])
												
				elif material_list_holder[l][i].calculation_type == 'value':
				
					self.cond_backgr = 1.0 / material_list_holder[l][i].resistivity_medium
								
				if material_list_holder[l][i].calculation_type == 'value':
					cond[material_idx] = self.cond_backgr
				else:
					cond[material_idx] = mat_sel_obj.calculate_conductivity(method = 'array')
					
				print('The conductivity for the material ' + material_list_holder[l][i].name + ' is calculated.')
				
			#converting all zero vals in the cond to None values
			cond[cond == 0.0] = np.nan
			
			cond_list.append(cond)
		
		if type == 'background':
			self.background_cond = cond_list[0]
			return cond_list[0]
		elif type == 'maximum':
			self.maximum_cond = cond_list[0]
			return cond_list[0]
		else:
			self.background_cond = cond_list[0]
			self.maximum_cond = cond_list[1]
			return cond_list[0],cond_list[1]
	
	def calculate_deformation_related_conductivity(self, method = 'plastic_strain', function_method = 'linear',
		low_deformation_threshold = 1e-2, high_deformation_threshold = 100, num_cpu = 1):
	
		try:
			self.background_cond
			self.maximum_cond
			
		except NameError:
			raise NameError('Background and maximum conductivities has to be assigned first to the model object to use deformation related conductivity function.')
	
		if num_cpu > 1:
			import multiprocessing
			import os
			from functools import partial
			
			max_num_cores = os.cpu_count()
			
			if num_cpu > max_num_cores:
				raise ValueError('There are not enough cpus in the machine to run this action with ' + str(num_cpu) + ' cores.')
	
		deform_cond = np.zeros_like(self.T)
		
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
				material_idx = return_material_bool(material_index = self.material_list[i].material_index, model_array = self.material_array, material_skip = mat_skip)
				#turning material index list to be useable format for the np.ndarray fields
				material_idx_list = [[material_idx[0][idx],material_idx[1][idx],material_idx[2][idx]] for idx in range(0,len(material_idx[0]))]
				
				#multiprocessing loop for each material
				with multiprocessing.Pool(processes=num_cpu) as pool:
					
					process_item_partial = partial(run_deform2cond, p_strain = self.p_strain, background_cond = self.background_cond,
					max_cond = self.maximum_cond, low_deformation_threshold = low_deformation_threshold,
					high_deformation_threshold = high_deformation_threshold, function_method = function_method,
					conductivity_decay_factor = self.material_list[i].deformation_dict['conductivity_decay_factor'],
					strain_decay_factor = self.material_list[i].deformation_dict['strain_decay_factor'])
					
					c = pool.map(process_item_partial, material_idx_list)
									
				deform_cond[material_idx] = c
				
				print('The deformation related conductivity for the material ' + self.material_list[i].name + ' is calculated.')
			
			#converting all zero vals in the cond to None values
			deform_cond[deform_cond == 0.0] = np.nan
			
		return deform_cond
		




