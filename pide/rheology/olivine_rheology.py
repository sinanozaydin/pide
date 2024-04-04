#!/usr/bin/env python3

import numpy as np


class olivine_rheology(object):

	def __init__(self, T, P, water = 0, xFe = 0.1, difffusion_model = None, dislocation_model = None, GBS_model = None,):
	
		self.R_const = 8.3144621
		self.Water_cutoff_rate = 5.0 #ppm
		self.Water_cutoff_rate_hsi = self.Water_cutoff_rate * (4.39e4) / 2695.0 
		self.fugacity_calculated = False
		self.T = T
		self.P = P
		
		if isinstance(self.T, (np.ndarray)):
			if isinstance(self.P, (np.ndarray)):
				if len(self.T) != len(self.P):
					raise ValueError('The length of the T and P arrays are not the same.')
		
		if isinstance(xFe, (int,float)):
			if isinstance(self.T, (np.ndarray, list)):
				self.xFe = xFe * np.ones(len(self.T))
				self.water = water * np.ones(len(self.T))
			else:
				self.xFe = xFe
				self.water = water
				
		else:
			self.xFe = xFe
			self.water = water
		
	def calculate_effective_viscosity(stress, strain_diff = 0.0, strain_disl = 0.0, strain_GBS = 0.0,other_strain_list = None):

		#calculates the effective viscosity calculated from given stress and calculated strain rates
		#stress in MPa
				
		total_strain = [strain_diff, strain_disl, strain_GBS]
		
		if other_strain_list != None:
		
			if isinstance(other_strain_list, (list,np.ndarray)):
			
				for item in other_strain_list:
				
					total_strain.append(total_strain)
					
			else:
			
				raise TypeError('The type for other_strain_list is not supported. It has to be list or a numpy.ndarray.')
		
		eff_visc = (stress * 1e6) / (2 * (sum(total_strain))) #pa s-1
		
		return eff_visc
		
	def convert_water_to_fh2o(self, fugacity_model = 'Zhao2004'):
		
		#Water in ppm
		#P in GPa
		#T in K
				
		water = self.water * 16.289430602953274 #converting to H/Si^6
		
		if fugacity_model == 'Zhao2004':
			A_fugacity = 90
			E_fugacity = 50000.0
			dV_fugacity = 10e-6
			alpha_fugacity = 97000.0
		elif fugacity_model == 'Kohlstedt1996':
			A_fugacity = 1.1
			E_fugacity = 0
			dV_fugacity = 10.6e-6
			alpha_fugacity = 0
		else:
			raise ValueError('There is no such reference model for fugacity: ' + str(fugacity_model))
	
		
		f_h2o = water  / (A_fugacity * np.exp(- (E_fugacity + (dV_fugacity * 1e9 * self.P)) / (self.R_const * self.T)) * np.exp((alpha_fugacity * self.xFe) / (self.R_const * self.T)))
		
		self.fugacity_calculated = True
		
		return f_h2o #out in MPa
		
	def Hirth_Kohlstedt_2003_diff_fugacity(self, gr_sz, stress, melt, fugacity_model = 'Zhao2004', calibration_model = 'Paterson1982'):
	
		#P in GPa
		#T in Kelvin
		#Water in ppm
		#stress in MPa
		#melt in fraction
		#xFe in fraction
		#fugacity_model Zhao2004 or Kohlstedth1996
			
		gr_sz = gr_sz * 1e3 #mm to micron
		
		if self.fugacity_calculated == False:
			self.fh2o = self.convert_water_to_fh2o(water = self.water, fugacity_model = fugacity_model)
		
		strain_array = np.zeros(len(self.T))
	
		n = 1
		p = 3
		alpha = 30.0
		
		for i in range(0,len(self.T)):
		
			if self.water[i] >= self.Water_cutoff_rate:
			
				A = 2.5e7
				r = 1
				dV = 15e-6
				E = 375000.0
				if calibration_model == 'Paterson1982':
					water_corr = 1.0
				elif calibration_model == 'Bell2003':
					water_corr = 3.0**r
				elif calibration_model == 'Withers2012':
					water_corr = 1.9**r
				else:
					raise ValueError('There is no such reference model for water calibration: ' + str(calibration_model))
				
			else:
			
				A = 1.5e9
				r = 0.0
				dV = 5e-6 #the middle of the range
				E = 375000.0
				water_corr = 1.0

			strain = A * (stress**n) * (gr_sz**-p) * ((self.fh2o[i]**r)) * np.exp(melt * alpha) * np.exp(-(E + (self.P[i]*1e9*dV)) / (self.R_const * self.T[i]))
			
			strain_array[i] = strain * water_corr
			
		return strain_array
	
	def Hirth_Kohlstedt_2003_dislocation_fugacity(self, stress, melt, fugacity_model = 'Zhao2004', calibration_model = 'Paterson1982'):
	
		n = 3.5
		p = 0
		alpha = 37.5
		
		strain_array = np.zeros(len(self.T))

		if self.fugacity_calculated == False:
			self.fh2o = self.convert_water_to_fh2o(water = self.water, fugacity_model = fugacity_model)
		
		for i in range(0,len(self.T)):
		
			if self.water[i] >= self.Water_cutoff_rate:

				A = 1600.0
				r = 1.2
				dV = 22e-6
				E = 520000.0
				if calibration_model == 'Paterson1982':
					water_corr = 1.0
				elif calibration_model == 'Bell2003':
					water_corr = 3.0**r
				elif calibration_model == 'Withers2012':
					water_corr = 1.9**r
				else:
					raise ValueError('There is no such reference model for water calibration: ' + str(calibration_model))
				
			else:
			
				A = 1.1e5
				r = 0.0
				dV = 15e-6 #the middle of the range
				E = 530000.0
				water_corr = 1.0

			strain = A * (stress**n) * ((self.fh2o[i]**r)) * np.exp(melt * alpha) * np.exp(-(E + (self.P*1e9*dV)) / (self.R_const * self.T[i]))

			strain_array[i] = strain * water_corr

		return strain_array
		
	def Ohuchi_et_al_2014_GBS(self, gr_sz, stress, fugacity_model = 'Zhao2004',calibration_model = 'Paterson1982'):
	
		gr_sz = gr_sz * 1e-3 #to m
		
		if self.fugacity_calculated == False:
			self.fh2o = self.convert_water_to_fh2o(water = self.water, fugacity_model = fugacity_model)
	
		A = 10**-4.89
		n = 3
		p = 1
		r = 1.25
		E = 423000.0
		dV = 17.6e-6
		if calibration_model == 'Paterson1982':
			water_corr = 1.0
		elif calibration_model == 'Bell2003':
			water_corr = 3.0**r
		elif calibration_model == 'Withers2012':
			water_corr = 1.9**r
		else:
			raise ValueError('There is no such reference model for water calibration: ' + str(calibration_model))
	
		strain = A * ((stress**n) / (gr_sz)) * (self.fh2o**r) * np.exp(-(E + (self.P*1e9*dV)) / (self.R_const * self.T))
		
		strain = strain * water_corr
		
		return strain