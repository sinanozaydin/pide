#!/usr/bin/env python3

import numpy as np

import SEL
from geodyn.material_process import return_material_bool

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
				
			#converting all zero vals in the cond tuple to None values
			cond = tuple(np.where(array == 0, np.nan, array) for array in cond)
			cond_list.append(cond)
		
		if type == 'background':
			return cond_list[0]
		elif type == 'maximum':
			return cond_list[0]
		else:
			return cond_list[0],cond_list[1]
			




