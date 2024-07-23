#!/usr/bin/env python3

from .utils.utils import check_type, text_color
from pide import pide

class Material(object):

	def __init__(self, name = "Unnamed", material_index = None, calculation_type = 'mineral', composition = None, melt_fluid_frac = 0.0,
	interconnectivities = None, param1 = None, el_cond_selections = None, melt_fluid_incorporation_method = 'Field', melt_or_fluid = 'melt', melt_fluid_m = 8.0,
	melt_properties = None,	melt_fluid_cond_selection = None, water_distr = False, water = None, xfe = None, solid_phase_mixing_idx = 0, melt_fluid_phase_mixing_idx = 0,
	deformation_dict = None, top = None, bottom = None, **kwargs):
	
		"""
		This is the core object to create a material object for pide.
		
		---------------------------------- -----------------------------------------
		Methods                             Description
		---------------------------------- -------------------------------------------
		
		
		
		"""
	
		self.mineral_list = ['ol','opx','cpx','garnet','mica','amp','quartz','plag','kfelds','sulphide','graphite','mixture','sp','wds','rwd','perov','other','bulk']
		self.rock_list = ['granite', 'granulite', 'sandstone', 'gneiss', 'amphibolite', 'basalt', 'mud', 'gabbro', 'other_rock']
		self.name = name
		self.material_index = material_index
		self.calculation_type = calculation_type
		self.melt_or_fluid = melt_or_fluid
		self.melt_fluid_frac = melt_fluid_frac
		self.melt_fluid_incorporation_method = melt_fluid_incorporation_method
		
		if composition == None:
			if self.calculation_type == 'rock':
				composition = {'granite':1}
			else:
				composition = {'ol':1}
		self._composition = None
		self.composition = composition
		
		#magnetotellurics
		if interconnectivities == None:
			if self.calculation_type == 'rock':
				interconnectivities = {'granite':1}
			else:
				interconnectivities = {'ol':1}
		self._interconnectivities = None
		self.interconnectivities = interconnectivities
		
		if el_cond_selections == None:
			if self.calculation_type == 'rock':
				el_cond_selections = {'granite':0}
			else:
				el_cond_selections = {'ol':0}
		self._el_cond_selections = None
		self.el_cond_selections = el_cond_selections
		
		if melt_fluid_cond_selection == None:
			if self.melt_or_fluid == 'melt':
				melt_fluid_cond_selection = 0
			else:
				melt_fluid_cond_selection = 0
		self._melt_fluid_cond_selection = None
		self.melt_fluid_cond_selection = melt_fluid_cond_selection
		
		self.water_distr = water_distr
		
		if water == None:
			if self.calculation_type == 'rock':
				water = {'granite':0}
			else:
				water = {'ol':0}
		self._water = None
		self.water = water
		
		if param1 == None:
			if self.calculation_type == 'rock':
				param1 = {'granite':0}
			else:
				param1 = {'ol':0}
				
		self._param1 = None
		self.param1 = param1

		self._top = None
		self.top = top

		self._bottom = None
		self.bottom = bottom
						
		if xfe == None:
			if self.calculation_type == 'rock':
				xfe = {'granite':0.1}
			else:
				xfe = {'ol':0.1}
		self._xfe = None
		self.xfe = xfe
		
		self._solid_phase_mixing_idx = None
		self.solid_phase_mixing_idx = solid_phase_mixing_idx
		
		self._melt_fluid_phase_mixing_idx = None
		self.melt_fluid_phase_mixing_idx = melt_fluid_phase_mixing_idx
		
		self._melt_fluid_m = None,
		self.melt_fluid_m = melt_fluid_m
		
		if melt_properties == None:
			melt_properties = {'water': 20, 'co2': 0, 'na2o': 0, 'k2o':0}
		self._melt_properties = None
		self.melt_properties = melt_properties
		
		if deformation_dict == None:
			deformation_dict = {'function_method':'linear','conductivity_decay_factor':0, 'strain_decay_factor':0, 'strain_percolation_threshold': None}
		
		self._deformation_dict = None
		self.deformation_dict = deformation_dict
				
		self._mantle_water_part = None
		self.mantle_water_part = kwargs.pop('mantle_water_part', {'opx_ol':0,'cpx_ol':0,'garnet_ol':0, 'ol_melt':0, 'opx_melt':0, 'cpx_melt':0,'garnet_melt':0})
		
		self.mantle_water_sol_ref = kwargs.pop('mantle_water_sol_ref', 'ol')
		
		self.resistivity_medium = kwargs.pop('resistivity_medium', None)
		self.vp_medium = kwargs.pop('vp_medium', None)
		self.vs_medium = kwargs.pop('vs_medium', None)
		
		self.water_calib = kwargs.pop('water_calib', {'ol':3,'px_gt':2,'feldspar':2})
		
		self.o2_buffer = kwargs.pop('o2_buffer', 0)
		
		self.linked_material_index = kwargs.pop('linked_material_index', None)
		
		self.fluid_salinity = kwargs.pop('fluid_salinity', 0.0)
		
		self.al_opx = kwargs.pop('al_opx', 0.0)				
		
		if (self.calculation_type == 'value') and (self.resistivity_medium == None):
		
			raise AttributeError('Calculation type is selected as value. You have to set resistivity medium as a floating number in Ohm meters.')
		
		#magnetics
		self.magnetic_susceptibility = kwargs.pop('magnetic_susceptivility', 1e-6)
		
	def calculate_conductivity(self, T, P, melt = None):
	
		from .model import run_model
				
		if check_type(T) == 'array':
			if len(T) != len(P):
			
				raise IndexError('The length of the temperature and pressure arrays do not match')
				
		if melt is None:
			import numpy as np
			melt = np.zeros(len(T))
		
		p_obj_material = pide()
		
		cond = run_model(index_list=list(range(0,len(T))), material = self, pide_object=p_obj_material, 
		t_array=T, p_array=P, melt_array=melt)
		
		return cond
		
	def calculate_seismic_velocity(self, T, P, melt = None):
	
		from .model import run_model
		
		if check_type(T) == 'array':
			if len(T) != len(P):
			
				raise IndexError('The length of the temperature and pressure arrays do not match')
				
		if melt is None:
			import numpy as np
			melt = np.zeros(len(T))
				
		p_obj_material = pide()
		
		v_bulk, v_p, v_s = run_model(index_list=list(range(0,len(T))), material = self, pide_object=p_obj_material, 
		t_array=T, p_array=P, melt_array=melt, type = 'seismic')
		
		return v_bulk, v_p, v_s
			
	def check_vals(self,value,type):
		
		for item in value:
			
			if self.calculation_type == 'mineral':
				if (item in self.mineral_list) == False:
					raise ValueError('The mineral ' + item + ' is wrongly defined in the composition dictionary. The possible mineral names are:' + str(self.mineral_list))
			elif self.calculation_type == 'rock':
				if (item in self.rock_list) == False:
					raise ValueError('The rock ' + item + ' is wrongly defined in the composition dictionary. The possible rock names are:' + str(self.rock_list))
			elif self.calculation_type == 'value':
				pass
			else:
				raise ValueError('The calculation type is wrongly defined. It has to be one of those three: 1.mineral, 2.rock, 3.value.')
		
		if self.calculation_type == 'mineral':
			list2check = self.mineral_list
		elif self.calculation_type == 'rock':
			list2check = self.rock_list
		else:
			list2check = None
		
		if list2check != None:
		
			for item in list2check:
				if item not in value:
					if type == 'comp':
						value[item] = 0
					elif type == 'archie':
						value[item] = 8.0
		else:
			value = None
			
		return value
		
	def check_property(self, value, list_vals):
	
		for item in list_vals:
			if item not in value:
				value[item] = 0.0
				
		return value
		
	def set_parameter(self, param_name, value):
	
		setattr(self, param_name, value)
		
	def copy_attributes(self, dest_object):
	
		if type(dest_object) == list:
			for idx_list in range(0,len(dest_object)):
				for attr_name in vars(self):
					setattr(dest_object[idx_list], attr_name, getattr(self, attr_name))
		else:
		
			for attr_name in vars(self):
				setattr(dest_object, attr_name, getattr(self, attr_name))
			
	
		
	#attributes listing here
	@property
	def composition(self):
		return self._composition
		
	@composition.setter
	def composition(self, value):
		self._composition = self.check_vals(value=value,type = 'comp')
		
	@property
	def interconnectivities(self):
		return self._interconnectivities
		
	@interconnectivities.setter
	def interconnectivities(self, value):
		self._interconnectivities = self.check_vals(value=value,type = 'archie')
		
	@property
	def param1(self):
		return self._param1
		
	@param1.setter
	def param1(self, value):
		self._param1 = self.check_vals(value=value,type = 'comp')
				
	@property
	def water(self):
		return self._water
		
	@water.setter
	def water(self, value):
		self._water = self.check_vals(value=value,type = 'comp')
		
	@property
	def melt_properties(self):
		return self._melt_properties
		
	@melt_properties.setter
	def melt_properties(self, value):
		self._melt_properties = self.check_property(value=value, list_vals = ['water','co2','na2o','k2o'])
		
	@property
	def el_cond_selections(self):
		return self._el_cond_selections
		
	@el_cond_selections.setter
	def el_cond_selections(self, value):
		self._el_cond_selections = self.check_vals(value=value,type = 'comp')
		
	@property
	def melt_fluid_cond_selection(self):
		return self._melt_fluid_cond_selection
		
	@melt_fluid_cond_selection.setter
	def melt_fluid_cond_selection(self, value):
		self._melt_fluid_cond_selection = value
		
	@property
	def xfe(self):
		return self._xfe
		
	@xfe.setter
	def xfe(self, value):
		self._xfe = self.check_vals(value=value,type = 'comp')
		
	@property
	def solid_phase_mixing_idx(self):
		return self._solid_phase_mixing_idx
		
	@solid_phase_mixing_idx.setter
	def solid_phase_mixing_idx(self, value):
		self._solid_phase_mixing_idx = value
		
	@property
	def melt_fluid_phase_mixing_idx(self):
		return self._melt_fluid_phase_mixing_idx
		
	@melt_fluid_phase_mixing_idx.setter
	def melt_fluid_phase_mixing_idx(self, value):
		self._melt_fluid_phase_mixing_idx = value
		
	@property
	def melt_fluid_m(self):
		return self._melt_fluid_m
		
	@melt_fluid_m.setter
	def melt_fluid_m(self, value):
		self._melt_fluid_m = value
		
	@property
	def deformation_dict(self):
		return self._deformation_dict
		
	@deformation_dict.setter
	def deformation_dict(self, value):
		self._deformation_dict = value
		
	@property
	def melt_fluid_frac(self):
		return self._melt_fluid_frac
		
	@melt_fluid_frac.setter
	def melt_fluid_frac(self, value):
		self._melt_fluid_frac = value
		
	@property
	def water_calib(self):
		return self._water_calib
		
	@water_calib.setter
	def water_calib(self, value):
		self._water_calib = value
		
	@property
	def o2_buffer(self):
		return self._o2_buffer
		
	@o2_buffer.setter
	def o2_buffer(self, value):
		self._o2_buffer = value
		
	@property
	def magnetic_susceptibility(self):
		return self._magnetic_susceptibility
		
	@magnetic_susceptibility.setter
	def magnetic_susceptibility(self, value):
		self._magnetic_susceptibility = value

	@property
	def top(self):
		return self._top
		
	@top.setter
	def top(self, value):
		self._top = value

	@property
	def bottom(self):
		return self._bottom
		
	@bottom.setter
	def bottom(self, value):
		self._bottom = value