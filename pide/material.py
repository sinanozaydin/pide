#!/usr/bin/env python3

from .utils.utils import check_type, text_color
from pide import pide

class Material(object):

	def __init__(self, name = "Unnamed", material_index = None, calculation_type = 'mineral', composition = None, melt_fluid_frac = 0.0,
	interconnectivities = None, param1 = None, el_cond_selections = None, melt_fluid_incorporation_method = 'Field', melt_or_fluid = 'melt', melt_fluid_m = 8.0,
	melt_properties = None, fluid_properties = None, melt_fluid_cond_selection = None, water_distr = False, water = None, xfe = None, solid_phase_mixing_idx = 0, melt_fluid_phase_mixing_idx = 0,
	deformation_dict = None, top = None, bottom = None, **kwargs):
	
		"""
		Material object: an python object to define properties of a certain material. 
		
		Input:
		str: name - Name of the material
		
		int: material_index - Index of the material for calculations within a pide.Model object.
		
		str: calculation_type - Type of calculation for a solid phases || 'mineral' or 'rock'.
		
		dict: composition -  dictionary entered in mineral_name:fraction or rock_name:fraction fashion.
											ol,opx,cpx... granite,granulite,sandstone...
											
		array: melt_fluid_frac - melt/fluid content array || in fraction.
		
		dict: interconnectivities - dictionary entered in mineral_name:cementation_exponent(m||float) fashion.
		
		dict: param1 - dictionary entered in mineral_name:param1 fashion.
		
		dict: el_cond_selections - electrical conductivirty model dictionary entered in mineral_name:el_cond_selection(int)
														or rock_name:el_cond_selection(int) fashion.
														
		str: melt_or_fluid - liquid phase type for the given material. || 'melt' or 'fluid.
		
		float: melt_fluid_m - interconnectivity of liquid phase used in Archie-type mixing models.
		
		dict: melt_properties - properties of melt entered in property_name:value fashion. 
											      'water','co2','na2o','k2o'
											      
		dict: fluid_properties - properties of fluid entered in property_name:value fashion.
												   'salinity'
											      
		int: melt_fluid_cond_selection - electrical conductivity model chosen for melt/fluid phase.
		
		bool: water_distr - Boolean to automatically distribute water among the phases or not.
		
		dict: water - dictionary to enter water content of the solid phases or bulk water content in mineral_name:water_content
																rock_name:water_content or bulk:water_content fashion.
		
		dict: xfe - dictionary to enter xFe value of the solid phases in mineral_name:xFe or rock_name:xFe fashion.
		
		int: solid_phase_mixing_idx - Solid phase mixing method selection.
		
		int: melt_fluid_phase_mixing_idx - solid/melt phase mixing method selection.
		
		dict: deformation_dict - dictionary to perform deform2cond function for the given material.
												'function_method', 'conductivity_decay_factor', 
												'conductivity_decay_factor_2' , 'strain_decay_factor'
		
		float: top - top of the material if geotherm calculation method will be used in pide.Model object. || in km.
		
		float: bottom - bottom of the material if geotherm calculation method will be used in pide.Model object. || in km.
		
		dict: mantle_water_part - mantle water partitioning dictionart entered in mineral_name:index fashion.
		
		float: resistivity_medium - resistivity of the material entered as a value instead of empirical calculation. || in ohm.m
		
		float: vp_medium - P-wave velocity value of the material entered as a value instead of empirical calculation. || in km/s
		
		float: vs_medium - S-wave velocity value of the material entered as a value instead of empirical calculation. || in km/s

		dict: water_calib - dictionary to enter water calibration methods for water calculations entered in material:mode fashion.
												'ol', 'px_gt', 'plag'
		int: o2_buffer - selection of o2 buffer with denoted integer value.
		
		array: al_opx - Al in orthopyroxene in array form || in w.t.%
		
		float: magnetic_susceptibility - magnetic susceptibility of the material || in Am^-1
		
		
		---------------------------------- -----------------------------------------
		Methods                             Description
		---------------------------------- -------------------------------------------
		calculate_conductivity              method to calculate conductivity of pide.Material object.
		
		calculate_seismic_velocity          method to calculate seismic velocity of pide.Material object.
		
		set_parameter                       method to set extra parameters that are not defined in default parameters.
		
		copy_attributes                     method to copy attributes to another pide.Material object.
	
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
		
		if fluid_properties == None:
			fluid_properties = {'salinity': 0.1}
		self._fluid_properties = None
		self.fluid_properties = fluid_properties
		
		if deformation_dict == None:
			deformation_dict = {'function_method':'linear','conductivity_decay_factor':0, 'conductivity_decay_factor_2':0,
			'strain_decay_factor':0, 'strain_percolation_threshold': None}
		
		self._deformation_dict = None
		self.deformation_dict = deformation_dict
				
		self._mantle_water_part = None
		self.mantle_water_part = kwargs.pop('mantle_water_part', {'opx_ol':0,'cpx_ol':0,'garnet_ol':0, 'ol_melt':0, 'opx_melt':0, 'cpx_melt':0,'garnet_melt':0})
				
		self.resistivity_medium = kwargs.pop('resistivity_medium', None)
		self.vp_medium = kwargs.pop('vp_medium', None)
		self.vs_medium = kwargs.pop('vs_medium', None)
		
		self.water_calib = kwargs.pop('water_calib', {'ol':3,'px_gt':2,'feldspar':2})
		
		self.o2_buffer = kwargs.pop('o2_buffer', 0)
						
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
	def fluid_properties(self):
		return self._fluid_properties
		
	@fluid_properties.setter
	def fluid_properties(self, value):
		self._fluid_properties = self.check_property(value=value, list_vals = ['salinity'])
		
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