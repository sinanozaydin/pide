#!/usr/bin/env python3

import SEL

class Material(object):

	def __init__(self, name = "Unnamed", mineral_or_rock = 'mineral', composition = {'ol':1}, interconnectivities = {'ol':1}, el_cond_selections = {'ol':0}, water_distr = False,
	water = {'ol':30}, xfe = {'ol':0.9}, phase_mixing_idx = 0, **kwargs):
	
		self.mineral_list = ['ol','opx','cpx','garnet','mica','amp','quartz','plag','kfelds','sulphide','graphite','mixture','sp','wds','rwd','perov','other','bulk']
	
		self.name = name
		self.mineral_or_rock = mineral_or_rock
		
		self._composition = None
		self.composition = composition
		
		self._interconnectivities = None
		self.interconnectivities = interconnectivities
		
		self._el_cond_selections = None
		self.el_cond_selections = el_cond_selections
		
		self.water_distr = water_distr
		
		self._water = None
		self.water = water
		
		self._xfe = None
		self.xfe = xfe
		
		self._phase_mixing_idx = None
		self.phase_mixing_idx = phase_mixing_idx
				
		self._mantle_water_part = None
		self.mantle_water_part = kwargs.pop('mantle_water_part', {'ol':0,'opx':0,'cpx':0,'garnet':0})
		
		self.mantle_water_sol_ref = kwargs.pop('mantle_water_sol_ref', 'ol')
		
		self.surpassing_resistivity = kwargs.pop('surpassing_resistivity',100)
			
	def check_mineral_vals(self,value,type):
		
		for item in value:
			if (item in self.mineral_list) == False:
				raise ValueError('The mineral ' + item + ' is wrongly defined in the composition dictionary. The possible mineral names are:' + str(self.mineral_list))
			
		for item in self.mineral_list:
			if item not in value:
				if type == 'comp':
					value[item] = 0
				elif type == 'archie':
					value[item] = 8.0
					
		return value
		
	#attributes listing here
	@property
	def composition(self):
		return self._composition
		
	@composition.setter
	def composition(self, value):
		self._composition = self.check_mineral_vals(value=value,type = 'comp')
		
	@property
	def interconnectivities(self):
		return self._interconnectivities
		
	@interconnectivities.setter
	def interconnectivities(self, value):
		self._interconnectivities = self.check_mineral_vals(value=value,type = 'archie')
		
	@property
	def el_cond_selections(self):
		return self._el_cond_selections
		
	@el_cond_selections.setter
	def el_cond_selections(self, value):
		self._el_cond_selections = self.check_mineral_vals(value=value,type = 'comp')
		
	@property
	def xfe(self):
		return self._xfe
		
	@xfe.setter
	def xfe(self, value):
		self._xfe = self.check_mineral_vals(value=value,type = 'comp')
		
	@property
	def solid_phase_mixing_idx(self):
		return self._solid_phase_mixing_idx
		
	@solid_phase_mixing_idx.setter
	def solid_phase_mixing_idx(self, value):
		self._solid_phase_mixing_idx = value
		

		
		