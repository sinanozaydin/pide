#!/usr/bin/env python3

import numpy as np

def brute_water_solver(sol_index, water_search_start, water_search_end, water_search_increment):

	self.water_end = False
	water_search = np.arange(water_search_start, water_search_end[sol_index],water_search_increment)
	restart = True
	while restart:

		self.residual_list = []
		res_check = []
		restart = False
		for j in range(0,len(water_search)):
			self.h2o[sol_index] = water_search[j]
			self.calculate_water(method = 'index', idx = sol_index)
			self.calculate_ol_conductivity(method = 'index', sol_idx = sol_index)
			self.calculate_opx_conductivity(method = 'index', sol_idx = sol_index)
			self.calculate_cpx_conductivity(method = 'index', sol_idx = sol_index)
			self.calculate_gt_conductivity(method = 'index', sol_idx = sol_index)
			if self.pl_method == '0':
				self.calculate_pl_conductivity(method = 'index', sol_idx = sol_index)
			if self.amp_method == '0':
				self.calculate_amp_conductivity(method = 'index', sol_idx = sol_index)
			self.calculate_sp_chr_conductivity(method = 'index', sol_idx = sol_index)
			if self.melt_method == '0':
				self.calculate_melt_conductivity(method = 'index', sol_idx = sol_index)
			self.phase_mixing_function(method = self.phs_mix_method, melt_method = self.phs_melt_mix_method, indexing_method = 'index', sol_idx = sol_index)

			residual = (1.0/self.bulk_cond[sol_index] - res_to_invert[sol_index])

			if j == 0:
				if residual < 0.0:
					if water_search_start != 0.0:
						water_search_start = water_search_start - (water_search_start * 0.1)
						if water_search_start <= 20:
							water_search_start = 0
						water_search = np.arange(water_search_start,water_search_end[sol_index],water_search_increment)
						restart = True
						break

			self.residual_list.append(abs(residual))
			res_check.append(residual)

			chek = 0
			if len(res_check) > 5:
				for k in range(-1,-4,-1):
					if res_check[k] < 0.0:
						chek = chek + 1
				if chek > 2:
					water_search_increment = water_search_increment / 2.0
					if water_search[j] >= water_search_increment:
						water_search_start = water_search[j-5]
					else:
						water_search_start = 0.0

					water_search = np.arange(water_search_start,water_search_end[sol_index],water_search_increment)
					if water_search_increment <= 0.5:
						restart = False
						break
					else:
						restart = True
						break
	if len(self.residual_list) == 0:
		self.h2o[sol_index] = 0.0
		self.residual_to_write = -1
	else:
		min_idx = self.residual_list.index(min(self.residual_list))
		self.h2o[sol_index] = water_search[min_idx]
		self.residual_to_write = min(self.residual_list)

	
