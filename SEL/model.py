#!/usr/bin/env python3

import numpy as np

import SEL
from geodyn.material_process import return_material_bool

class Model(object):

	def __init__(self, material_list, material_array, T, P, melt = None, p_strain = None, strain_rate = None):

		self.material_list = material_list
		self.material_array = material_array
		self.T = T
		self.P = P
		self.melt_frac = melt
		self.p_strain = p_strain
		self.strain_rate = strain_rate

	def calculate_background_conductivity(self):

		cond = np.zeros_like(self.T)

		for i in range(0,len(self.material_list)):

			#determining the material indexes from the material array
			material_idx = return_material_bool(self.material_list[i].material_index, self.material_array)

			#getting only the relevant arrays for calculation
			t_relevant = self.T[material_idx]
			p_relevant = self.P[material_idx]
			melt_relevant = self.melt_frac[material_idx]

			mat_sel_obj = SEL.SEL()
			mat_sel_obj.set_temperature(t_relevant)
			mat_sel_obj.set_pressure(p_relevant)

			if self.material_list[i].calculation_type == 'mineral':
				mat_sel_obj.set_composition_solid_mineral(ol = self.material_list[i].composition['ol'],
				opx = self.material_list[i].composition['opx'],
				cpx = self.material_list[i].composition['cpx'],
				garnet = self.material_list[i].composition['garnet'],
				mica = self.material_list[i].composition['mica'],
				amp = self.material_list[i].composition['amp'],
				quartz = self.material_list[i].composition['quartz'],
				plag = self.material_list[i].composition['plag'],
				kfelds = self.material_list[i].composition['kfelds'],
				sulphide = self.material_list[i].composition['sulphide'],
				graphite = self.material_list[i].composition['graphite'],
				mixture = self.material_list[i].composition['mixture'],
				sp = self.material_list[i].composition['sp'],
				wds = self.material_list[i].composition['wds'],
				rwd = self.material_list[i].composition['rwd'],
				perov = self.material_list[i].composition['perov'],
				other = self.material_list[i].composition['other'])
				
				if self.material_list[i].phase_mixing_idx == 0:
					mat_sel_obj.set_phase_interconnectivities(ol = self.material_list[i].interconnectivities['ol'],
					opx = self.material_list[i].interconnectivities['opx'],
					cpx = self.material_list[i].interconnectivities['cpx'],
					garnet = self.material_list[i].interconnectivities['garnet'],
					mica = self.material_list[i].interconnectivities['mica'],
					amp = self.material_list[i].interconnectivities['amp'],
					quartz = self.material_list[i].interconnectivities['quartz'],
					plag = self.material_list[i].interconnectivities['plag'],
					kfelds = self.material_list[i].interconnectivities['kfelds'],
					sulphide = self.material_list[i].interconnectivities['sulphide'],
					graphite = self.material_list[i].interconnectivities['graphite'],
					mixture = self.material_list[i].interconnectivities['mixture'],
					sp = self.material_list[i].interconnectivities['sp'],
					wds = self.material_list[i].interconnectivities['wds'],
					rwd = self.material_list[i].interconnectivities['rwd'],
					perov = self.material_list[i].interconnectivities['perov'],
					other = self.material_list[i].interconnectivities['other'])
					
				if self.material_list[i].water_distr == False:
					mat_sel_obj.set_mineral_water(ol = self.material_list[i].water['ol'],
					opx = self.material_list[i].water['opx'],
					cpx = self.material_list[i].water['cpx'],
					garnet = self.material_list[i].water['garnet'],
					mica = self.material_list[i].water['mica'],
					amp = self.material_list[i].water['amp'],
					quartz = self.material_list[i].water['quartz'],
					plag = self.material_list[i].water['plag'],
					kfelds = self.material_list[i].water['kfelds'],
					sulphide = self.material_list[i].water['sulphide'],
					graphite = self.material_list[i].water['graphite'],
					mixture = self.material_list[i].water['mixture'],
					sp = self.material_list[i].water['sp'],
					wds = self.material_list[i].water['wds'],
					rwd = self.material_list[i].water['rwd'],
					perov = self.material_list[i].water['perov'],
					other = self.material_list[i].water['other'])
					
				else:
					mat_sel_obj.set_bulk_water(self.material_list[i].water['bulk'])
					mat_sel_obj.mantle_water_distribute(method = 'array')
					
				mat_sel_obj.set_mineral_conductivity_choice(ol = self.material_list[i].el_cond_selections['ol'],
					opx = self.material_list[i].el_cond_selections['opx'],
					cpx = self.material_list[i].el_cond_selections['cpx'],
					garnet = self.material_list[i].el_cond_selections['garnet'],
					mica = self.material_list[i].el_cond_selections['mica'],
					amp = self.material_list[i].el_cond_selections['amp'],
					quartz = self.material_list[i].el_cond_selections['quartz'],
					plag = self.material_list[i].el_cond_selections['plag'],
					kfelds = self.material_list[i].el_cond_selections['kfelds'],
					sulphide = self.material_list[i].el_cond_selections['sulphide'],
					graphite = self.material_list[i].el_cond_selections['graphite'],
					mixture = self.material_list[i].el_cond_selections['mixture'],
					sp = self.material_list[i].el_cond_selections['sp'],
					wds = self.material_list[i].el_cond_selections['wds'],
					rwd = self.material_list[i].el_cond_selections['rwd'],
					perov = self.material_list[i].el_cond_selections['perov'],
					other = self.material_list[i].el_cond_selections['other'])
					
				cond = mat_sel_obj.calculate_conductivity(method = 'array')
			




