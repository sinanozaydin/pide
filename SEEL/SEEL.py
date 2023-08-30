#!/usr/bin/env python3

import os

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , 'seel_src')

import sys, csv, platform, warnings
import numpy as np
import iapws

#Importing external functions

sys.path.append(core_path_ext)

#importing odd melt/fluid functions
from seel_src.cond_models.melt_odd import * 
from seel_src.cond_models.fluids_odd import * 
#importing odd rock functions
from seel_src.cond_models.rocks.granite_odd import * 
from seel_src.cond_models.rocks.granulite_odd import *
from seel_src.cond_models.rocks.sandstone_odd import *
from seel_src.cond_models.rocks.gneiss_odd import *
from seel_src.cond_models.rocks.amphibolite_odd import *
from seel_src.cond_models.rocks.basalt_odd import *
from seel_src.cond_models.rocks.mud_odd import *
from seel_src.cond_models.rocks.gabbro_odd import *
from seel_src.cond_models.rocks.other_rocks_odd import *
#importing odd mineral functions
from seel_src.cond_models.minerals.quartz_odd import *
from seel_src.cond_models.minerals.plag_odd import *
from seel_src.cond_models.minerals.amp_odd import *
from seel_src.cond_models.minerals.kfelds_odd import *
from seel_src.cond_models.minerals.opx_odd import *
from seel_src.cond_models.minerals.cpx_odd import *
from seel_src.cond_models.minerals.mica_odd import *
from seel_src.cond_models.minerals.garnet_odd import *
from seel_src.cond_models.minerals.ol_odd import *
from seel_src.cond_models.minerals.mixtures_odd import *
from seel_src.cond_models.minerals.other_odd import *

warnings.filterwarnings("ignore", category=RuntimeWarning) #ignoring many RuntimeWarning printouts that are useless

#Version 0.1, June. 2023.
#SEEL - (S)ynthetic (E)lectrical (E)arth (L)ibrary
#Program written by Sinan Ozaydin (Macquarie University, School of Natural Sciences
#sciences, Australia).

#Indentation method: hard tabs ('\t')

#Works with Python3
#Required libraries: numpy,matplotlib,PyQt5
#optional libraries: pyperclip

class SEEL(object):
	
	def __init__(self, core_path = core_path_ext):

		self.core_path = core_path
				
		self.home()
		
	def home(self):
		
		#Setting up initial variables.

		SEEL.loaded_file = False
		self.cond_calculated = False

		self.init_params = self.read_csv(filename = os.path.join(self.core_path,'init_param.csv'),delim = ',') #loading the blueprint parameter file.
		
		self.read_cond_models()
		self.read_params()
		
		#setting up default values for the SEEL object
		self.set_temperature(900.0) #in Kelvin
		self.set_pressure(1.0) #in GPa
		self.set_mineral_conductivity_choice()
		self.set_rock_conductivity_choice()
		self.set_mineral_water()
		self.set_watercalib()
		self.set_o2_buffer()
		self.set_param1_mineral()
		self.set_param2_mineral()
		self.set_param1_rock()
		self.set_param2_rock()
		self.set_melt_or_fluid_mode(mode = 1) #default choice is melt - 1
		self.set_solid_phase_method(mode = 2) #default choice is mineral - 2
		self.set_melt_fluid_conductivity_choice()
		self.set_melt_fluid_frac()
		self.set_melt_properties()
		self.set_fluid_properties()
		self.set_phase_interconnectivities()
		
	def read_csv(self,filename,delim):

		#Simple function for reading csv files and give out filtered output for given delimiter (delim)

		file_obj = open(filename,'rt',encoding = "utf8") #Creating file object
		file_csv = csv.reader(file_obj,delimiter = delim) #Reading the file object with csv module, delimiter assigned to ','
		data = [] #Creating empty array to append data

		#Appending data from csb object
		for row in file_csv:
			data.append(row)

		#Filtering data for None elements read.
		for j in range(0,len(data)):
			data[j] = list(filter(None,data[j]))
		data = list(filter(None,data))

		return data

	def read_cond_models(self):

		#A function that reads conductivity model files and get the data.

		self.fluid_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'fluids.csv'),delim = ',') 
		self.melt_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'melt.csv'),delim = ',')

		#reading rocks
		self.granite_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'granite.csv'),delim = ',')
		self.granulite_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'granulite.csv'),delim = ',')
		self.sandstone_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'sandstone.csv'),delim = ',')
		self.gneiss_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'gneiss.csv'),delim = ',')
		self.amphibolite_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'amphibolite.csv'),delim = ',')
		self.basalt_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'basalt.csv'),delim = ',')
		self.mud_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'mud.csv'),delim = ',')
		self.gabbro_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'gabbro.csv'),delim = ',')
		self.other_rock_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'rocks', 'other_rock.csv'),delim = ',')

		#reading minerals
		self.quartz_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'quartz.csv'),delim = ',')
		self.plag_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'plag.csv'),delim = ',')
		self.amp_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'amp.csv'),delim = ',')
		self.kfelds_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'kfelds.csv'),delim = ',')
		self.opx_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'opx.csv'),delim = ',')
		self.cpx_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'cpx.csv'),delim = ',')
		self.mica_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'mica.csv'),delim = ',')
		self.garnet_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'garnet.csv'),delim = ',')
		self.sulphides_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'sulphides.csv'),delim = ',')
		self.graphite_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'graphite.csv'),delim = ',')
		self.ol_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'ol.csv'),delim = ',')
		self.mixture_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'mixtures.csv'),delim = ',')
		self.other_cond_data = self.read_csv(os.path.join(self.core_path, 'cond_models' , 'minerals', 'other.csv'),delim = ',')
		

		self.cond_data_array = [self.fluid_cond_data, self.melt_cond_data, self.granite_cond_data, self.granulite_cond_data,
			  self.sandstone_cond_data, self.gneiss_cond_data, self.amphibolite_cond_data, self.basalt_cond_data, self.mud_cond_data,
			   self.gabbro_cond_data, self.other_rock_cond_data, self.quartz_cond_data, self.plag_cond_data,
			  self.amp_cond_data, self.kfelds_cond_data, self.opx_cond_data, self.cpx_cond_data, self.mica_cond_data,
			  self.garnet_cond_data, self.sulphides_cond_data, self.graphite_cond_data, self.ol_cond_data, self.mixture_cond_data,
			  self.other_cond_data]

		len_fluid = len(self.fluid_cond_data) - 1 
		len_melt = len(self.melt_cond_data) - 1

		self.fluid_num = 2

		len_granite = len(self.granite_cond_data) - 1
		len_granulite = len(self.granulite_cond_data) - 1
		len_sandstone = len(self.sandstone_cond_data) - 1
		len_gneiss = len(self.gneiss_cond_data) - 1
		len_amphibolite = len(self.amphibolite_cond_data) - 1
		len_basalt = len(self.basalt_cond_data) - 1
		len_mud = len(self.mud_cond_data) - 1
		len_gabbro = len(self.gabbro_cond_data) - 1
		len_other_rock = len(self.other_rock_cond_data) - 1
		
		self.rock_num = 9

		len_quartz = len(self.quartz_cond_data) - 1
		len_plag = len(self.plag_cond_data) - 1
		len_amp = len(self.amp_cond_data) - 1
		len_kfelds = len(self.kfelds_cond_data) - 1
		len_opx = len(self.opx_cond_data) - 1
		len_cpx = len(self.cpx_cond_data) - 1
		len_mica = len(self.mica_cond_data) - 1
		len_garnet = len(self.garnet_cond_data) - 1
		len_sulphides = len(self.sulphides_cond_data) - 1
		len_graphite = len(self.graphite_cond_data) - 1
		len_ol = len(self.ol_cond_data) - 1
		len_mixture = len(self.mixture_cond_data) - 1
		len_other = len(self.other_cond_data) - 1

		self.mineral_num = 15

		#Creating empty arrays for appending new data.
		SEEL.name = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		SEEL.type = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		SEEL.t_min = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		SEEL.t_max = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.p_min = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.p_max = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.w_calib = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.mg_cond = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_i =  [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_i_err =  [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_i =  [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_i_err =  [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_pol = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_pol_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_pol = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_pol_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_p = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.sigma_p_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_p = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.h_p_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.r = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.r_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.alpha_p = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.alpha_p_err = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.wtype = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]
		self.dens_mat = [[None] * len_fluid, [None] * len_melt, [None] * len_granite, [None] * len_granulite, [None] * len_sandstone, [None] * len_gneiss,
		   [None] * len_amphibolite, [None] * len_basalt, [None] * len_mud, [None] * len_gabbro, [None] * len_other_rock, [None] * len_quartz,
			[None] * len_plag, [None] * len_amp, [None] * len_kfelds, [None] * len_opx, [None] * len_cpx, [None] * len_mica,
			[None] * len_garnet, [None] * len_sulphides, [None] * len_graphite, [None] * len_ol, [None] * len_mixture, [None] * len_other]

		#Filling up the arrays.
		for i in range(0,len(SEEL.type)):
			count = 1
			for j in range(0,len(SEEL.type[i])):
				SEEL.name[i][count-1] = self.cond_data_array[i][count][0]
				SEEL.type[i][count-1] = self.cond_data_array[i][count][1]
				SEEL.t_min[i][count-1] = float(self.cond_data_array[i][count][2])
				SEEL.t_max[i][count-1] = float(self.cond_data_array[i][count][3])
				self.p_min[i][count-1] = float(self.cond_data_array[i][count][4])
				self.p_max[i][count-1] = float(self.cond_data_array[i][count][5])
				self.w_calib[i][count-1] = int(self.cond_data_array[i][count][6])
				self.mg_cond[i][count-1] = float(self.cond_data_array[i][count][7])
				self.sigma_i[i][count-1] = float(self.cond_data_array[i][count][8])
				self.sigma_i_err[i][count-1] = float(self.cond_data_array[i][count][9])
				self.h_i[i][count-1] = float(self.cond_data_array[i][count][10])
				self.h_i_err[i][count-1] = float(self.cond_data_array[i][count][11])
				self.sigma_pol[i][count-1] = float(self.cond_data_array[i][count][12])
				self.sigma_pol_err[i][count-1] = float(self.cond_data_array[i][count][13])
				self.h_pol[i][count-1] = float(self.cond_data_array[i][count][14])
				self.h_pol_err[i][count-1] = float(self.cond_data_array[i][count][15])
				self.sigma_p[i][count-1] = float(self.cond_data_array[i][count][16])
				self.sigma_p_err[i][count-1] = float(self.cond_data_array[i][count][17])
				self.h_p[i][count-1] = float(self.cond_data_array[i][count][18])
				self.h_p_err[i][count-1] = float(self.cond_data_array[i][count][19])
				self.r[i][count-1] = float(self.cond_data_array[i][count][20])
				self.r_err[i][count-1] = float(self.cond_data_array[i][count][21])
				self.alpha_p[i][count-1] = float(self.cond_data_array[i][count][22])
				self.alpha_p_err[i][count-1] = float(self.cond_data_array[i][count][23])
				self.wtype[i][count-1] = int(self.cond_data_array[i][count][24])
				self.dens_mat[i][count-1] = float(self.cond_data_array[i][count][25])
			
				count += 1

	def read_params(self):

		#READING THE PARAMETERS IN PARAMS.CSV WHICH ARE GENERAL PHYSICAL CONSTANTS
		#AND PROPERTIES OF MATERIALS

		params_dat = self.read_csv(os.path.join(self.core_path, 'params.csv'), delim = ',')

		self.g = float(params_dat[0][1]) # in kg/
		self.R = float(params_dat[1][1])
		self.avog = float(params_dat[2][1])
		self.boltz = float(params_dat[3][1])
		self.el_q = float(params_dat[4][1])
		SEEL.spreadsheet = str(params_dat[5][1])
		self.mu = 4.0 * np.pi * 10**(-7)
		
	def set_composition_solid_mineral(self, **kwargs):
	
		#Enter composition in fraction 0.6 == 60% volumetric percentage
	
		self.ol_frac = np.array(kwargs.pop('ol_frac', 0.0))
		self.opx_frac = np.array(kwargs.pop('opx_frac', 0.0))
		self.cpx_frac = np.array(kwargs.pop('cpx_frac', 0.0))
		self.garnet_frac = np.array(kwargs.pop('garnet_frac', 0.0))
		self.mica_frac = np.array(kwargs.pop('mica_frac', 0.0))
		self.amp_frac = np.array(kwargs.pop('amp_frac', 0.0))
		self.quartz_frac = np.array(kwargs.pop('quartz_frac', 0.0))
		self.plag_frac = np.array(kwargs.pop('plag_frac', 0.0))
		self.kfelds_frac = np.array(kwargs.pop('kfelds_frac', 0.0))
		self.graphite_frac = np.array(kwargs.pop('graphite_frac', 0.0))
		self.sulphide_frac = np.array(kwargs.pop('sulphide_frac', 0.0))
		self.mixture_frac = np.array(kwargs.pop('mixture_frac', 0.0))
		self.other_frac = np.array(kwargs.pop('other_frac', 0.0))
		
		bool_composition = self.check_composition(method = 'mineral')

		if bool_composition == False:
		
			raise ValueError('Some of the values entered in mineral composition do not add up to 1.')
			
	
	def set_composition_solid_rock(self, **kwargs):
	
		#Enter composition in fraction 0.6 == 60% volumetric percentage
	
		self.granite_frac = np.array(kwargs.pop('granite_frac', 0.0))
		self.granulite_frac = np.array(kwargs.pop('granulite_frac', 0.0))
		self.sandstone_frac = np.array(kwargs.pop('sandstone_frac', 0.0))
		self.gneiss_frac = np.array(kwargs.pop('gneiss_frac', 0.0))
		self.amphibolite_frac = np.array(kwargs.pop('amphibolite_frac', 0.0))
		self.basalt_frac = np.array(kwargs.pop('basalt_frac', 0.0))
		self.mud_frac = np.array(kwargs.pop('mud_frac', 0.0))
		self.gabbro_frac = np.array(kwargs.pop('gabbro_frac', 0.0))
		self.other_rock_frac = np.array(kwargs.pop('other_rock_frac', 0.0))
	
		
		bool_composition = self.check_composition(method = 'rock')

		if bool_composition == False:
		
			raise ValueError('Some of the values entered in mineral composition do not add up to 1.')
			
			
	def set_temperature(self,T):
	
		try: 
			len(T)
			self.T = T
		except TypeError:
			self.T = np.array(T)
		
	def set_pressure(self,P):
		
		try: 
			len(P)
			self.p = P
			t_check = self.check_p_n_T()
			if t_check == False:
				raise ValueError('The arrays of pressure and temperature are not the same...')
				sys.exit()
		except TypeError:
			try:
				self.p = np.ones(len(self.T)) * P
			except TypeError:
				self.p = np.ones(1) * P
			
	def check_p_n_T(self):
	
		if len(self.T) != len(self.p):
			T_check = False
		else:
			T_check = True
			
		return T_check
		
	def set_watercalib(self,**kwargs):
	
		SEEL.ol_calib = kwargs.pop('ol', 3)
		SEEL.px_gt_calib = kwargs.pop('px-gt', 2)
		SEEL.feldspar_calib = kwargs.pop('feldspar', 2)
		
	def set_o2_buffer(self,**kwargs):
		
		SEEL.o2_buffer = kwargs.pop('o2_buffer', 0)

	def check_composition(self, method = None):

		continue_adjusting = True

		if method == 'rock':

			tot = self.granite_frac + self.granulite_frac + self.sandstone_frac +\
			self.gneiss_frac + self.amphibolite_frac + self.basalt_frac + self.mud_frac +\
				 self.gabbro_frac + self.other_rock_frac
			
			if any(item <= 0.99 for item in tot) == True:
				continue_adjusting = False
			
			if any(item >= 1.01 for item in tot) == True:

				continue_adjusting = False
				
		elif method == 'mineral':
			
			tot = self.quartz_frac + self.plag_frac + self.amp_frac + self.kfelds_frac +\
			self.opx_frac + self.cpx_frac + self.mica_frac + self.garnet_frac + self.sulphide_frac + self.graphite_frac +\
			self.ol_frac + self.mixture_frac + self.other_frac

			if any(item <= 0.99 for item in tot) == True:
			
				continue_adjusting = False
			
			if any(item >= 1.01 for item in tot) == True:

				continue_adjusting = False

		return continue_adjusting
	
	def match_mineral_index(self, mineral_name):
	
		if (mineral_name == 'ol') or (mineral_name == 'olivine'):
			min_index = 21
		elif (mineral_name == 'opx') or (mineral_name == 'orthopyroxene'):
			min_index = 15
		elif (mineral_name == 'cpx') or (mineral_name == 'clinopyroxene'):
			min_index = 16
		elif (mineral_name == 'garnet') or (mineral_name == 'gt'):
			min_index = 18
		elif (mineral_name == 'mica') or (mineral_name == 'Mica'):
			min_index = 17
		elif (mineral_name == 'amp') or (mineral_name == 'amphibole'):
			min_index = 13
		elif (mineral_name == 'quartz') or (mineral_name == 'qtz'):
			min_index = 11
		elif (mineral_name == 'plag') or (mineral_name == 'plagioclase'):
			min_index = 12
		elif (mineral_name == 'kfelds') or (mineral_name == 'kfeldspar'):
			min_index = 14
		elif (mineral_name == 'sulphide') or (mineral_name == 'Sulphide'):
			min_index = 19
		elif (mineral_name == 'graphite') or (mineral_name == 'Graphite'):
			min_index = 20
		elif (mineral_name == 'mixture') or (mineral_name == 'mixtures'):
			min_index = 22
		elif (mineral_name == 'other') or (mineral_name == 'other'):
			min_index = 23
			
		else:
		
			raise ValueError('There is no such a mineral specifier called :' + mineral_name)
			
		return min_index
		
	def list_mineral_econd_models(self, mineral_name):
		
		if (mineral_name == 'ol') or (mineral_name == 'olivine'):
			min_index = 21
		elif (mineral_name == 'opx') or (mineral_name == 'orthopyroxene'):
			min_index = 15
		elif (mineral_name == 'cpx') or (mineral_name == 'clinopyroxene'):
			min_index = 16
		elif (mineral_name == 'garnet') or (mineral_name == 'gt'):
			min_index = 18
		elif (mineral_name == 'mica') or (mineral_name == 'Mica'):
			min_index = 17
		elif (mineral_name == 'amp') or (mineral_name == 'amphibole'):
			min_index = 13
		elif (mineral_name == 'quartz') or (mineral_name == 'qtz'):
			min_index = 11
		elif (mineral_name == 'plag') or (mineral_name == 'plagioclase'):
			min_index = 12
		elif (mineral_name == 'kfelds') or (mineral_name == 'kfeldspar'):
			min_index = 14
		elif (mineral_name == 'sulphide') or (mineral_name == 'Sulphide'):
			min_index = 19
		elif (mineral_name == 'graphite') or (mineral_name == 'Graphite'):
			min_index = 20
		elif (mineral_name == 'mixture') or (mineral_name == 'mixtures'):
			min_index = 22
		elif (mineral_name == 'other') or (mineral_name == 'other'):
			min_index = 23
			
		def print_lists(min_idx):
		
			for i in range(0,len(self.name[min_idx])):
				print(str(i) + '.  ' + self.name[min_idx][i])
			
		print_lists(min_idx = min_index)
		
		return self.name[min_index]
		
	def list_rock_econd_models(self, rock_name):
		
		if (rock_name == 'granite'):
			rock_idx = 2
		elif (rock_name == 'granulite'):
			rock_idx = 3
		elif (rock_name == 'sandstone'):
			rock_idx = 4
		elif (rock_name == 'gneiss'):
			rock_idx = 5
		elif (rock_name == 'amphibolite'):
			rock_idx = 6
		elif (rock_name == 'basalt'):
			rock_idx = 7
		elif (rock_name == 'mud'):
			rock_idx = 8
		elif (rock_name == 'gabbro'):
			rock_idx = 9
		elif (rock_name == 'other_rock'):
			rock_idx = 10
			
		def print_lists(rock_idx):
		
			for i in range(0,len(self.name[rock_idx])):
				print(str(i) + '.  ' + self.name[rock_idx][i])
			
		print_lists(rock_idx = rock_idx)
		
		return self.name[rock_idx]
		
	def list_melt_econd_models(self):
	
		for i in range(0,len(self.name[1])):
			print(str(i) + '.  ' + self.name[1][i])
			
		return self.name[1]
		
	def list_fluid_econd_models(self):
	
		for i in range(0,len(self.name[0])):
			print(str(i) + '.  ' + self.name[0][i])
			
		return self.name[0]
		
	def set_melt_fluid_conductivity_choice(self,**kwargs):
	
		SEEL.melt_cond_selection = kwargs.pop('melt', 0)
		SEEL.fluid_cond_selection = kwargs.pop('fluid', 0)
		
	def set_mineral_conductivity_choice(self,**kwargs):
	
		SEEL.ol_cond_selection = kwargs.pop('ol', 0)
		SEEL.opx_cond_selection = kwargs.pop('opx', 0)
		SEEL.cpx_cond_selection = kwargs.pop('cpx', 0)
		SEEL.garnet_cond_selection = kwargs.pop('garnet', 0)
		SEEL.mica_cond_selection = kwargs.pop('mica', 0)
		SEEL.amp_cond_selection = kwargs.pop('amp', 0)
		SEEL.quartz_cond_selection = kwargs.pop('quartz', 0)
		SEEL.plag_cond_selection = kwargs.pop('plag', 0)
		SEEL.kfelds_cond_selection = kwargs.pop('kfelds', 0)
		SEEL.sulphide_cond_selection = kwargs.pop('sulphide', 0)
		SEEL.graphite_cond_selection = kwargs.pop('graphite', 0)
		SEEL.mixture_cond_selection = kwargs.pop('mixture', 0)
		SEEL.other_cond_selection = kwargs.pop('other', 0)
		
		SEEL.minerals_cond_selections = [SEEL.quartz_cond_selection, SEEL.plag_cond_selection, SEEL.amp_cond_selection, SEEL.kfelds_cond_selection, SEEL.opx_cond_selection,
				   SEEL.cpx_cond_selection, SEEL.mica_cond_selection, SEEL.garnet_cond_selection, SEEL.sulphide_cond_selection,
				   SEEL.graphite_cond_selection, SEEL.ol_cond_selection, SEEL.mixture_cond_selection, SEEL.other_cond_selection]
				   
	
	def set_rock_conductivity_choice(self,**kwargs):
	
		SEEL.granite_cond_selection = kwargs.pop('granite', 0)
		SEEL.granulite_cond_selection = kwargs.pop('granulite', 0)
		SEEL.sandstone_cond_selection = kwargs.pop('sandstone', 0)
		SEEL.gneiss_cond_selection = kwargs.pop('gneiss', 0)
		SEEL.amphibolite_cond_selection = kwargs.pop('amphibolite', 0)
		SEEL.basalt_cond_selection = kwargs.pop('basalt', 0)
		SEEL.mud_cond_selection = kwargs.pop('mud', 0)
		SEEL.gabbro_cond_selection = kwargs.pop('gabbro', 0)
		SEEL.other_rock_cond_selection = kwargs.pop('other_rock', 0)
		
		SEEL.rock_cond_selections = [SEEL.granite_cond_selection, SEEL.granulite_cond_selection, SEEL.sandstone_cond_selection, SEEL.gneiss_cond_selection,
				   SEEL.amphibolite_cond_selection, SEEL.basalt_cond_selection, SEEL.mud_cond_selection, SEEL.gabbro_cond_selection, SEEL.other_rock_cond_selection]
				   
	def set_mineral_water(self, **kwargs):
	
		SEEL.ol_water = np.array(kwargs.pop('ol', 0))
		SEEL.opx_water = np.array(kwargs.pop('opx', 0))
		SEEL.cpx_water = np.array(kwargs.pop('cpx', 0))
		SEEL.garnet_water = np.array(kwargs.pop('garnet', 0))
		SEEL.mica_water = np.array(kwargs.pop('mica', 0))
		SEEL.amp_water = np.array(kwargs.pop('amp', 0))
		SEEL.quartz_water = np.array(kwargs.pop('quartz', 0))
		SEEL.plag_water = np.array(kwargs.pop('plag', 0))
		SEEL.kfelds_water = np.array(kwargs.pop('kfelds', 0))
		SEEL.sulphide_water = np.array(kwargs.pop('sulphide', 0))
		SEEL.graphite_water = np.array(kwargs.pop('graphite', 0))
		SEEL.mixture_water = np.array(kwargs.pop('mixture', 0))
		SEEL.other_water = np.array(kwargs.pop('other', 0))

		SEEL.mineral_water_list = [SEEL.quartz_water, SEEL.plag_water, SEEL.amp_water, SEEL.kfelds_water,
			 SEEL.opx_water, SEEL.cpx_water, SEEL.mica_water, SEEL.garnet_water, SEEL.sulphide_water,
				   SEEL.graphite_water, SEEL.ol_water, SEEL.mixture_water, SEEL.other_water]
				   
	def set_rock_water(self, **kwargs):
	
		SEEL.granite_water = np.array(kwargs.pop('granite', 0))
		SEEL.granulite_water = np.array(kwargs.pop('granulite', 0))
		SEEL.sandstone_water = np.array(kwargs.pop('sandstone', 0))
		SEEL.gneiss_water = np.array(kwargs.pop('gneiss', 0))
		SEEL.amphibolite_water = np.array(kwargs.pop('amphibolite', 0))
		SEEL.basalt_water = np.array(kwargs.pop('basalt', 0))
		SEEL.mud_water = np.array(kwargs.pop('mud', 0))
		SEEL.gabbro_water = np.array(kwargs.pop('gabbro', 0))
		SEEL.other_rock_water = np.array(kwargs.pop('other_rock', 0))
		
		SEEL.rock_water_list = [SEEL.granite_water, SEEL.granulite_water,
			SEEL.sandstone_water, SEEL.gneiss_water, SEEL.amphibolite_water, SEEL.basalt_water,
			SEEL.mud_water, SEEL.gabbro_water, SEEL.other_rock_water]
			
	def set_param1_mineral(self, **kwargs):
	
		SEEL.ol_param1 = np.array(kwargs.pop('ol', 0))
		SEEL.opx_param1 = np.array(kwargs.pop('opx', 0))
		SEEL.cpx_param1 = np.array(kwargs.pop('cpx', 0))
		SEEL.garnet_param1 = np.array(kwargs.pop('garnet', 0))
		SEEL.mica_param1 = np.array(kwargs.pop('mica', 0))
		SEEL.amp_param1 = np.array(kwargs.pop('amp', 0))
		SEEL.quartz_param1 = np.array(kwargs.pop('quartz', 0))
		SEEL.plag_param1 = np.array(kwargs.pop('plag', 0))
		SEEL.kfelds_param1 = np.array(kwargs.pop('kfelds', 0))
		SEEL.sulphide_param1 = np.array(kwargs.pop('sulphide', 0))
		SEEL.graphite_param1 = np.array(kwargs.pop('graphite', 0))
		SEEL.mixture_param1 = np.array(kwargs.pop('mixture', 0))
		SEEL.other_param1 = np.array(kwargs.pop('other', 0))
		
		SEEL.param1_mineral_list = [SEEL.quartz_param1, SEEL.plag_param1, SEEL.amp_param1, SEEL.kfelds_param1,
			 SEEL.opx_param1, SEEL.cpx_param1, SEEL.mica_param1, SEEL.garnet_param1, SEEL.sulphide_param1,
				   SEEL.graphite_param1, SEEL.ol_param1, SEEL.mixture_param1, SEEL.other_param1]
				   
	def set_param1_rock(self, **kwargs):
	
		SEEL.granite_param1 = np.array(kwargs.pop('granite', 0))
		SEEL.granulite_param1 = np.array(kwargs.pop('granulite', 0))
		SEEL.sandstone_param1 = np.array(kwargs.pop('sandstone', 0))
		SEEL.gneiss_param1 = np.array(kwargs.pop('gneiss', 0))
		SEEL.amphibolite_param1 = np.array(kwargs.pop('amphibolite', 0))
		SEEL.basalt_param1 = np.array(kwargs.pop('basalt', 0))
		SEEL.mud_param1 = np.array(kwargs.pop('mud', 0))
		SEEL.gabbro_param1 = np.array(kwargs.pop('gabbro', 0))
		SEEL.other_rock_param1 = np.array(kwargs.pop('other_rock', 0))
		
		SEEL.param1_rock_list = [SEEL.granite_param1, SEEL.granulite_param1,
			SEEL.sandstone_param1, SEEL.gneiss_param1, SEEL.amphibolite_param1, SEEL.basalt_param1,
			SEEL.mud_param1, SEEL.gabbro_param1, SEEL.other_rock_param1]
			
	def set_param2_mineral(self, **kwargs):
	
		SEEL.ol_param2 = np.array(kwargs.pop('ol', 0))
		SEEL.opx_param2 = np.array(kwargs.pop('opx', 0))
		SEEL.cpx_param2 = np.array(kwargs.pop('cpx', 0))
		SEEL.garnet_param2 = np.array(kwargs.pop('garnet', 0))
		SEEL.mica_param2 = np.array(kwargs.pop('mica', 0))
		SEEL.amp_param2 = np.array(kwargs.pop('amp', 0))
		SEEL.quartz_param2 = np.array(kwargs.pop('quartz', 0))
		SEEL.plag_param2 = np.array(kwargs.pop('plag', 0))
		SEEL.kfelds_param2 = np.array(kwargs.pop('kfelds', 0))
		SEEL.sulphide_param2 = np.array(kwargs.pop('sulphide', 0))
		SEEL.graphite_param2 = np.array(kwargs.pop('graphite', 0))
		SEEL.mixture_param2 = np.array(kwargs.pop('mixture', 0))
		SEEL.other_param2 = np.array(kwargs.pop('other', 0))
		
		SEEL.param2_mineral_list = [SEEL.quartz_param2, SEEL.plag_param2, SEEL.amp_param2, SEEL.kfelds_param2,
			SEEL.opx_param2, SEEL.cpx_param2, SEEL.mica_param2, SEEL.garnet_param2, SEEL.sulphide_param2,
				   SEEL.graphite_param2, SEEL.ol_param2, SEEL.mixture_param2, SEEL.other_param2]
				   
	def set_param2_rock(self, **kwargs):
	
		SEEL.granite_param2 = np.array(kwargs.pop('granite', 0))
		SEEL.granulite_param2 = np.array(kwargs.pop('granulite', 0))
		SEEL.sandstone_param2 = np.array(kwargs.pop('sandstone', 0))
		SEEL.gneiss_param2 = np.array(kwargs.pop('gneiss', 0))
		SEEL.amphibolite_param2 = np.array(kwargs.pop('amphibolite', 0))
		SEEL.basalt_param2 = np.array(kwargs.pop('basalt', 0))
		SEEL.mud_param2 = np.array(kwargs.pop('mud', 0))
		SEEL.gabbro_param2 = np.array(kwargs.pop('gabbro', 0))
		SEEL.other_rock_param2 = np.array(kwargs.pop('other_rock', 0))
		
		SEEL.param2_rock_list = [SEEL.granite_param2, SEEL.granulite_param2,
			SEEL.sandstone_param2, SEEL.gneiss_param2, SEEL.amphibolite_param2, SEEL.basalt_param2,
			SEEL.mud_param2, SEEL.gabbro_param2, SEEL.other_rock_param2]
			
	def set_melt_fluid_frac(self, **kwargs):
	
		self.melt_fluid_mass_frac = np.array(kwargs.pop('frac', 0))
		
	def set_melt_or_fluid_mode(self,mode):
	
		SEEL.fluid_or_melt_method = mode
		
	def set_solid_phase_method(self,mode):
	
		SEEL.solid_phase_method = mode
			
	def set_melt_properties(self, **kwargs):
	
		self.co2_melt = np.array(kwargs.pop('co2', 0)) #in ppm
		self.h2o_melt = np.array(kwargs.pop('water', 0)) #in ppm
		self.na2o_melt = np.array(kwargs.pop('na2o', 0)) #in wt
		self.k2o_melt = np.array(kwargs.pop('k2o', 0)) #in wt
				
	def set_fluid_properties(self, **kwargs):
	
		self.salinity_fluid = np.array(kwargs.pop('salinity', 0))
		
	def set_phase_interconnectivities(self,**kwargs):
	
		SEEL.ol_m = np.array(kwargs.pop('ol', 4))
		SEEL.opx_m = np.array(kwargs.pop('opx', 4))
		SEEL.cpx_m = np.array(kwargs.pop('cpx', 4))
		SEEL.garnet_m = np.array(kwargs.pop('garnet', 4))
		SEEL.mica_m = np.array(kwargs.pop('mica', 4))
		SEEL.amp_m = np.array(kwargs.pop('amp', 4))
		SEEL.quartz_m = np.array(kwargs.pop('quartz', 4))
		SEEL.plag_m = np.array(kwargs.pop('plag', 4))
		SEEL.kfelds_m = np.array(kwargs.pop('kfelds', 4))
		SEEL.sulphide_m = np.array(kwargs.pop('sulphide', 4))
		SEEL.graphite_m = np.array(kwargs.pop('graphite', 4))
		SEEL.mixture_m = np.array(kwargs.pop('mixture', 4))
		SEEL.other_m = np.array(kwargs.pop('other', 4))
		
		SEEL.granite_m = np.array(kwargs.pop('granite', 4))
		SEEL.granulite_m = np.array(kwargs.pop('granulite', 4))
		SEEL.sandstone_m = np.array(kwargs.pop('sandstone', 4))
		SEEL.gneiss_m = np.array(kwargs.pop('gneiss', 4))
		SEEL.amphibolite_m = np.array(kwargs.pop('amphibolite', 4))
		SEEL.basalt_m = np.array(kwargs.pop('basalt', 4))
		SEEL.mud_m = np.array(kwargs.pop('mud', 4))
		SEEL.gabbro_m = np.array(kwargs.pop('gabbro', 4))
		SEEL.other_rock_m = np.array(kwargs.pop('other_rock', 4))

	def calculate_arrhenian_single(self, T, sigma, E, r, alpha, water):

		if (sigma == 0.0) and (E == 0.0):
			cond = 0.0
		else:
			cond = (10.0**sigma) * (water**r) * np.exp(-(E + (alpha * water)**(1.0/3.0)) / (self.R * T))

		return cond
	
	def calculate_fluids_conductivity(self, method, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond_fluids = np.zeros(len(self.T))

		if SEEL.type[0][SEEL.fluid_cond_selection] == '0':

			cond_fluids[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[0][SEEL.fluid_cond_selection],
								   E = self.h_i[0][SEEL.fluid_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[0][SEEL.fluid_cond_selection],
								   E = self.h_pol[0][SEEL.fluid_cond_selection],r = 0, alpha = 0, water = 0)
			
		elif SEEL.type[0][SEEL.fluid_cond_selection] == '1':

			cond_fluids[idx_node] =  self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[0][SEEL.fluid_cond_selection],
								   E = self.h_i[0][SEEL.fluid_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[0][SEEL.fluid_cond_selection],
								   E = self.h_pol[0][SEEL.fluid_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[0][SEEL.fluid_cond_selection],
								   E = self.h_p[0][SEEL.fluid_cond_selection],r = 0, alpha = 0, water = 0)
			
		elif SEEL.type[0][SEEL.fluid_cond_selection] == '3':

			if ('*' in SEEL.name[0][SEEL.fluid_cond_selection]) == True:

				fluids_odd_function = SEEL.name[0][SEEL.fluid_cond_selection].replace('*','')

			else:

				fluids_odd_function = SEEL.name[0][SEEL.fluid_cond_selection]

			cond_fluids[idx_node] = eval(fluids_odd_function + '(T = self.T[idx_node], P = self.p[idx_node], salinity = self.salinity_fluid[idx_node], method = method)')
	
		return cond_fluids

	def calculate_melt_conductivity(self, method, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond_melt = np.zeros(len(self.T))
		
		if ("Wet" in self.name[1][SEEL.melt_cond_selection]) == True:
					
			if self.wtype[1][SEEL.melt_cond_selection] == 0:
				water_corr_factor = 1e4 #converting to wt % if the model requires
			else:
				water_corr_factor = 1.0
			
		else:
			
			water_corr_factor = 1.0

		if SEEL.type[1][SEEL.melt_cond_selection] == '0':

			cond_melt[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[1][SEEL.melt_cond_selection],
								   E = self.h_i[1][SEEL.melt_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[1][SEEL.melt_cond_selection],
								   E = self.h_pol[1][SEEL.melt_cond_selection],r = 0, alpha = 0, water = 0)
			
		elif SEEL.type[1][SEEL.melt_cond_selection] == '1':

			cond_melt[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[1][SEEL.melt_cond_selection],
								   E = self.h_i[1][SEEL.melt_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[1][SEEL.melt_cond_selection],
								   E = self.h_pol[1][SEEL.melt_cond_selection],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[1][SEEL.melt_cond_selection],
								   E = self.h_p[1][SEEL.melt_cond_selection], r = self.r[1][SEEL.melt_cond_selection], alpha = self.alpha_p[1][SEEL.melt_cond_selection],
								   water = self.h2o_melt/water_corr_factor)
			
		elif SEEL.type[1][SEEL.melt_cond_selection] == '3':

			if ('*' in SEEL.name[1][SEEL.melt_cond_selection]) == True:

				melt_odd_function = SEEL.name[1][SEEL.melt_cond_selection].replace('*','')

			else:

				melt_odd_function = SEEL.name[1][SEEL.melt_cond_selection]

			cond_melt[idx_node] = eval(melt_odd_function + '(T = self.T[idx_node], P = self.p[idx_node], Melt_H2O = self.h2o_melt[idx_node]/water_corr_factor,' +
			'Melt_CO2 = self.co2_melt, Melt_Na2O = self.na2o_melt[idx_node], Melt_K2O = self.k2o_melt[idx_node], method = method)')
		
		return cond_melt

	def calculate_rock_conductivity(self, method, rock_idx = None, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond = np.zeros(len(self.T))

		rock_sub_idx = rock_idx - self.fluid_num
		
		if ("Wet" in self.name[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]]) == True:
					
			if self.wtype[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]] == 0:
				water_corr_factor = water_corr_factor * 1e4 #converting to wt % if the model requires
			else:
				water_corr_factor = 1.0
			
		else:
			
			water_corr_factor = 1.0

		if SEEL.type[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]] == '0':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   E = self.h_i[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   E = self.h_pol[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],r = 0, alpha = 0, water = 0)
			
		elif SEEL.type[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]] == '1':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   E = self.h_i[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   E = self.h_pol[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   E = self.h_p[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]], r = self.r[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]],
								   alpha = self.alpha_p[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]], water = SEEL.rock_water_list[rock_sub_idx][idx_node] / water_corr_factor)
			
		elif SEEL.type[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]] == '3':

			if ('*' in SEEL.name[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]]) == True:

				odd_function = SEEL.name[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]].replace('*','')

			else:

				odd_function = SEEL.name[rock_idx][SEEL.rock_cond_selections[rock_sub_idx]]

			if ('fo2' in odd_function) == True:
				cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], water = SEEL.rock_water_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node] / water_corr_factor, param1 = SEEL.param1_rock_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node], param2 = SEEL.param2_rock_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node], fo2 = self.calculate_fugacity(SEEL.o2_buffer),fo2_ref = self.calculate_fugacity(3), method = method)')
			else:
				cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], water = SEEL.rock_water_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node] / water_corr_factor, param1 = SEEL.param1_rock_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node], param2 = SEEL.param2_rock_list[SEEL.rock_cond_selections[rock_sub_idx]][idx_node], method = method)')
		
		return cond
	
	def calculate_mineral_conductivity(self, method, min_idx = None, **kwargs):
	
		sol_idx = kwargs.pop('sol_idx', 0)
	
		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond = np.zeros(len(self.T))

		min_sub_idx = min_idx - self.fluid_num - self.rock_num
		
		if ("Wet" in self.name[min_idx][SEEL.minerals_cond_selections[min_sub_idx]]) == True:
		
			water_corr_factor = self.Water_correction(min_idx = min_idx)
			
			if self.wtype[min_idx][SEEL.minerals_cond_selections[min_sub_idx]] == 0:
				water_corr_factor = water_corr_factor * 1e4 #converting to wt % if the model requires

			else:
				pass
			
		else:
			
			water_corr_factor = 1.0

		
		if SEEL.type[min_idx][SEEL.minerals_cond_selections[min_sub_idx]] == '0':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   E = self.h_i[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   E = self.h_pol[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],r = 0, alpha = 0, water = 0)
			
		elif SEEL.type[min_idx][SEEL.minerals_cond_selections[min_sub_idx]] == '1':
		

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_i[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   E = self.h_i[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   E = self.h_pol[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],r = 0, alpha = 0, water = 0) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   E = self.h_p[min_idx][SEEL.minerals_cond_selections[min_sub_idx]], r = self.r[min_idx][SEEL.minerals_cond_selections[min_sub_idx]],
								   alpha = self.alpha_p[min_idx][SEEL.minerals_cond_selections[min_sub_idx]], water = SEEL.mineral_water_list[min_sub_idx][idx_node] / water_corr_factor)
			
		elif SEEL.type[min_idx][SEEL.minerals_cond_selections[min_sub_idx]] == '3':

			if ('*' in SEEL.name[min_idx][SEEL.minerals_cond_selections[min_sub_idx]]) == True:

				odd_function = SEEL.name[min_idx][SEEL.minerals_cond_selections[min_sub_idx]].replace('*','')

			else:

				odd_function = SEEL.name[min_idx][SEEL.minerals_cond_selections[min_sub_idx]]

			if ('fo2' in odd_function) == True:

				cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], water = SEEL.mineral_water_list[min_sub_idx][idx_node] / water_corr_factor, param1 = SEEL.param1_mineral_list[min_sub_idx][idx_node], param2 = SEEL.param2_mineral_list[min_sub_idx][idx_node], fo2 = self.calculate_fugacity(SEEL.o2_buffer),fo2_ref = self.calculate_fugacity(3), method = method)')

			else:

				cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], water = SEEL.mineral_water_list[min_sub_idx][idx_node] / water_corr_factor, param1 = SEEL.param1_mineral_list[min_sub_idx][idx_node], param2 = SEEL.param2_mineral_list[min_sub_idx][idx_node], fo2 = None, fo2_ref = None, method = method)')

		return cond
		
	def Water_correction(self, min_idx = None):

		#A function that corrects the water content to desired calibration. Numbers are taken from Demouchy and Bolfan-Casanova (2016, Lithos) for olivine, pyroxene and garnet.
		#For feldspars, the correction number is taken from Mosenfelder et al. (2015, Am. Min.)

		if min_idx == 21:
		#olivine

			calib_object = self.w_calib[min_idx][SEEL.ol_cond_selection]
			calib_object_2 = SEEL.ol_calib

			if calib_object_2 == 0:

				if calib_object == 0:
					
					CORR_Factor = self.pat2with #Paterson to Withers

				elif calib_object == 1:

					CORR_Factor = self.bell2with #Bell to Withers from Demouchy and Bolfan Casanova

				else:

					CORR_Factor= 1.0 #Withers to Withers

			elif calib_object_2 == 1:

				if calib_object == 0:

					CORR_Factor = self.pat2bell #Paterson to Bell

				elif calib_object == 1:

					CORR_Factor = 1.0 #Bell to Bell

				else:

					CORR_Factor = self.with2bell #Withers to Bell


			elif calib_object_2 == 2:

				if calib_object == 0:


					CORR_Factor = 1.0 #Paterson to paterson

				elif calib_object == 1:

					CORR_Factor = self.bell2path #Bell to Paterson

				else:

					CORR_Factor = self.with2pat #Withers to Paterson
					
			else:
			
				CORR_Factor = 1.0

		elif (min_idx == 15) or (min_idx == 16) or (min_idx == 18): #opx, cpx and garnet

			calib_object = self.w_calib[min_idx][SEEL.minerals_cond_selections[min_idx - self.rock_num - self.fluid_num]]
			calib_object_2 = SEEL.px_gt_calib

			if calib_object_2 == 0:

				if calib_object == 0:

					CORR_Factor = self.pat2bell95

				else:

					CORR_Factor = 1.0

			elif calib_object_2 == 1:

				if calib_object == 1:

					CORR_Factor = self.bell952pat

				else:

					CORR_Factor = 1.0
			
			else:
			
				CORR_Factor = 1.0

		elif (min_idx == 12) or (min_idx == 14): #Plagioclase and k-feldspar
			
			calib_object = self.w_calib[min_idx][SEEL.minerals_cond_selections[min_idx - self.rock_num - self.fluid_num]]
			calib_object_2 = SEEL.feldspar_calib
				
			if calib_object_2 == 0:

				if calib_object == 0:

					CORR_Factor = self.john2mosen

				else:

					CORR_Factor = 1.0

			elif calib_object_2 == 1:

				if calib_object == 1:

					CORR_Factor = self.mosen2john

				else:

					CORR_Factor = 1.0
					
			else:
			
				CORR_Factor = 1.0

					
		else:
		
			CORR_Factor = 1.0
			
		
		return CORR_Factor
	
	def phase_mixing_function(self, method = None, melt_method = None, indexing_method = None, sol_idx = None):

		self.bulk_cond = np.zeros(len(self.T)) #setting up an empty bulk conductivity array for all methods
		self.dens_melt_fluid = np.zeros(len(self.T))

		if indexing_method == 'array':
			idx_node = None
		elif indexing_method == 'index':
			idx_node = sol_idx
			
		if method == 0:

			#Calculating phase exponent of the abundant mineral to make connectedness equal to unity.
			#From Glover (2010, Geophysics), analytic solution.

			#creating search limits for different indexing methods.
			if indexing_method == 'array':
				start_idx = 0
				end_idx = len(self.T)
			elif indexing_method == 'index':
				start_idx = sol_idx
				end_idx = sol_idx + 1
				
			for i in range(start_idx,end_idx):
			
				if SEEL.solid_phase_method == 1:
					phase_list = [self.granite_frac[i],self.granulite_frac[i],self.sandstone_frac[i],
					self.gneiss_frac[i], self.amphibolite_frac[i], self.basalt_frac[i], self.mud_frac[i],
					 self.gabbro_frac[i], self.other_rock_frac[i]]
					m_list = [SEEL.granite_m[i],SEEL.granulite_m[i],SEEL.sandstone_m[i],
					SEEL.gneiss_m[i], SEEL.amphibolite_m[i], SEEL.basalt_m[i], self.mud_m[i],
					 self.gabbro_m[i], SEEL.other_rock_m[i]]
				elif SEEL.solid_phase_method == 2:
					phase_list = [self.quartz_frac[i], self.plag_frac[i], self.amp_frac[i], self.kfelds_frac[i],
					self.opx_frac[i], self.cpx_frac[0], self.mica_frac[i], self.garnet_frac[i],
					self.sulphide_frac[i], self.graphite_frac[i], self.ol_frac[i], self.mixture_frac[i], self.other_frac[i]]
					m_list = [SEEL.quartz_m[i], SEEL.plag_m[i], SEEL.amp_m[i], SEEL.kfelds_m[i],
					SEEL.opx_m[i], SEEL.cpx_m[i], SEEL.mica_m[i], SEEL.garnet_m[i],
					SEEL.sulphide_m[i], SEEL.graphite_m[i], SEEL.ol_m[i], SEEL.mixture_m[i], SEEL.other_m[i]]
					
				frac_abundant = max(phase_list) #fraction of abundant mineral
				idx_max_ph = phase_list.index(frac_abundant) #index of the abundant mineral
				del phase_list[idx_max_ph] #deleting the abundant mineral form local list
				del m_list[idx_max_ph] #deleting the exponent of the abundant mineral from local list
				connectedness = np.asarray(phase_list)**np.asarray(m_list) #calculating the connectedness of the rest

				if sum(phase_list) != 0.0:
					m_abundant = np.log(1.0 - np.sum(connectedness)) / np.log(frac_abundant) #analytic solution to the problem
				else:
					m_abundant = 1

				if SEEL.solid_phase_method == 1:
					
					if idx_max_ph == 0:
						SEEL.granite_m[idx_node] = m_abundant
					elif idx_max_ph == 1:
						SEEL.granulite_m[idx_node] = m_abundant
					elif idx_max_ph == 2:
						SEEL.sandstone_m[idx_node] = m_abundant
					elif idx_max_ph == 3:
						SEEL.gneiss_m[idx_node] = m_abundant
					elif idx_max_ph == 4:
						SEEL.amphibolite_m[idx_node] = m_abundant
					elif idx_max_ph == 5:
						SEEL.basalt_m[idx_node] = m_abundant
					elif idx_max_ph == 6:
						SEEL.mud_m[idx_node] = m_abundant
					elif idx_max_ph == 7:
						SEEL.gabbro_m[idx_node] = m_abundant
					elif idx_max_ph == 8:
						SEEL.other_rock_m[idx_node] = m_abundant
					
					self.bulk_cond[idx_node] = (self.granite_cond[idx_node]*(self.granite_frac[idx_node]**SEEL.granite_m[idx_node])) +\
					(self.granulite_cond[idx_node]*(self.granulite_frac[idx_node]**SEEL.granulite_m[idx_node])) +\
					(self.sandstone_cond[idx_node]*(self.sandstone_frac[idx_node]**SEEL.sandstone_m[idx_node])) +\
					(self.gneiss_cond[idx_node]*(self.gneiss_frac[idx_node]**SEEL.gneiss_m[idx_node])) +\
					(self.amphibolite_cond[idx_node]*(self.amphibolite_frac[idx_node]**SEEL.amphibolite_m[idx_node])) +\
					(self.basalt_cond[idx_node]*(self.basalt_frac[idx_node]**SEEL.basalt_m[idx_node])) +\
					(self.mud_cond[idx_node]*(self.mud_frac[idx_node]**SEEL.mud_m[idx_node])) +\
					(self.gabbro_cond[idx_node]*(self.gabbro_frac[idx_node]**SEEL.gabbro_m[idx_node])) +\
					(self.other_rock_cond[idx_node]*(self.other_rock_frac[idx_node]**SEEL.other_rock_m[idx_node]))
				
				elif SEEL.solid_phase_method == 2:
					if idx_max_ph == 0:
						SEEL.quartz_m[idx_node] = m_abundant
					elif idx_max_ph == 1:
						SEEL.plag_m[idx_node] = m_abundant
					elif idx_max_ph == 2:
						SEEL.amp_m[idx_node] = m_abundant
					elif idx_max_ph == 3:
						SEEL.kfelds_m[idx_node] = m_abundant
					elif idx_max_ph == 4:
						SEEL.opx_m[idx_node] = m_abundant
					elif idx_max_ph == 5:
						SEEL.cpx_m[idx_node] = m_abundant
					elif idx_max_ph == 6:
						SEEL.mica_m[idx_node] = m_abundant
					elif idx_max_ph == 7:
						SEEL.garnet_m[idx_node] = m_abundant
					elif idx_max_ph == 8:
						SEEL.sulphide_m[idx_node] = m_abundant
					elif idx_max_ph == 9:
						SEEL.graphite_m[idx_node] = m_abundant
					elif idx_max_ph == 10:
						SEEL.ol_m[idx_node] = m_abundant
					elif idx_max_ph == 11:
						SEEL.mixture_m[idx_node] = m_abundant
					elif idx_max_ph == 12:
						SEEL.other_m[idx_node] = m_abundant
						
					self.bulk_cond[idx_node] = (self.quartz_cond[idx_node]*(self.quartz_frac[idx_node]**SEEL.quartz_m[idx_node])) +\
					(self.plag_cond[idx_node]*(self.plag_frac[idx_node]**SEEL.plag_m[idx_node])) +\
					(self.amp_cond[idx_node]*(self.amp_frac[idx_node]**SEEL.amp_m[idx_node])) +\
					(self.kfelds_cond[idx_node]*(self.kfelds_frac[idx_node]**SEEL.kfelds_m[idx_node])) +\
					(self.opx_cond[idx_node]*(self.opx_frac[idx_node]**SEEL.opx_m[idx_node])) +\
					(self.cpx_cond[idx_node]*(self.cpx_frac[idx_node]**SEEL.cpx_m[idx_node])) +\
					(self.mica_cond[idx_node]*(self.mica_frac[idx_node]**SEEL.mica_m[idx_node])) +\
					(self.garnet_cond[idx_node]*(self.garnet_frac[idx_node]**SEEL.garnet_m[idx_node])) +\
					(self.sulphide_cond[idx_node]*(self.sulphide_frac[idx_node]**SEEL.sulphide_m[idx_node])) +\
					(self.graphite_cond[idx_node]*(self.graphite_frac[idx_node]**SEEL.graphite_m[idx_node])) +\
					(self.ol_cond[idx_node]*(self.ol_frac[idx_node]**SEEL.ol_m[idx_node])) +\
					(self.mixture_cond[idx_node]*(self.mixture_frac[idx_node]**SEEL.mixture_m[idx_node])) +\
					(self.other_cond[idx_node]*(self.other_frac[idx_node]**SEEL.other_m[idx_node]))
					
		elif method == 1:
			
			if indexing_method == 'array':
				
				self.bulk_cond = np.zeros(len(self.ol_cond))
				start_idx = 0
				end_idx = len(self.T)
			elif indexing_method == 'index':
				start_idx = sol_idx
				end_idx = sol_idx + 1
				
			for i in range(start_idx,end_idx):
				
				if SEEL.solid_phase_method == 1:
					list_i = [self.granite_cond[i], self.granulite_cond[i], self.sandstone_cond[i],
					self.gneiss_cond[i],self.amphibolite_cond[i], self.basalt_cond[i], self.mud_cond[i],
					  self.gabbro_cond, self.other_rock_cond[i]]
				elif SEEL.solid_phase_method == 2:
					list_i = [self.quartz_cond[i], self.plag_cond[i], self.amp_cond[i],
					self.kfelds_cond[i],self.opx_cond[i],self.cpx_cond[i],self.mica_cond[i],
					self.garnet_cond[i],self.sulphide_cond[i],self.graphite_cond[i],self.ol_cond[i], self.mixture_cond[i], self.other_cond[i]]				
					
				while True:
				
					#while loop for deleting the zero arrays that could be encountered due to non-existence of the mineral.
					
					min_local = np.amin(np.asarray(list_i))
				
					if (min_local != 0.0):
						
						break
					
					else:
					
						list_i = np.delete(list_i, np.argwhere(list_i == 0))
						
				if SEEL.solid_phase_method == 1:
				
					self.bulk_cond[i] = (((self.granite_frac[i] / (self.granite_cond[i] + (2*min_local))) +\
					(self.granulite_frac[i] / (self.granulite_cond[i] + (2*min_local))) +\
					(self.sandstone_frac[i] / (self.sandstone_cond[i] + (2*min_local))) +\
					(self.gneiss_frac[i] / (self.gneiss_cond[i] + (2*min_local))) +\
					(self.amphibolite_frac[i] / (self.amphibolite_cond[i] + (2*min_local))) +\
					(self.basalt_frac[i] / (self.basalt_cond[i] + (2*min_local))) +\
					(self.mud_frac[i] / (self.mud_cond[i] + (2*min_local))) +\
					(self.gabbro_frac[i] / (self.gabbro_cond[i] + (2*min_local))) +\
					(self.other_rock_frac[i] / (self.other_rock_cond[i] + (2*min_local))))**(-1.0)) -\
					2.0*min_local
						
				elif SEEL.solid_phase_method == 2:
				
					self.bulk_cond[i] = (((self.quartz_frac[i] / (self.quartz_cond[i] + (2*min_local))) +\
					(self.plag_frac[i] / (self.plag_cond[i] + (2*min_local))) +\
					(self.amp_frac[i] / (self.amp_cond[i] + (2*min_local))) +\
					(self.kfelds_frac[i] / (self.kfelds_cond[i] + (2*min_local))) +\
					(self.opx_frac[i] / (self.opx_cond[i] + (2*min_local))) +\
					(self.cpx_frac[i] / (self.cpx_cond[i] + (2*min_local))) +\
					(self.mica_frac[i] / (self.mica_cond[i] + (2*min_local))) +\
					(self.garnet_frac[i] / (self.garnet_cond[i] + (2*min_local))) +\
					(self.sulphide_frac[i] / (self.sulphide_cond[i] + (2*min_local))) +\
					(self.graphite_frac[i] / (self.graphite_cond[i] + (2*min_local))) +\
					(self.ol_frac[i] / (self.ol_cond[i] + (2*min_local))) +\
					(self.mixture_frac[i] / (self.mixture_cond[i] + (2*min_local))) +\
					(self.other_frac[i] / (self.other_cond[i] + (2*min_local)))
					)**(-1.0)) -\
					2.0*min_local
					
		elif method == 2:
		
			if indexing_method == 'array':
				
				self.bulk_cond = np.zeros(len(self.ol_cond))
				start_idx = 0
				end_idx = len(self.T)
			elif indexing_method == 'index':
				start_idx = sol_idx
				end_idx = sol_idx + 1
				
			for i in range(start_idx,end_idx):
				
				if SEEL.solid_phase_method == 1:
					list_i = [self.granite_cond[i], self.granulite_cond[i], self.sandstone_cond[i],
					self.gneiss_cond[i],self.amphibolite_cond[i], self.basalt_cond[i], self.mud_cond[i],
					  self.gabbro_cond, self.other_rock_cond[i]]
				elif SEEL.solid_phase_method == 2:
					list_i = [self.quartz_cond[i], self.plag_cond[i], self.amp_cond[i],
					self.kfelds_cond[i],self.opx_cond[i],self.cpx_cond[i],self.mica_cond[i],
					self.garnet_cond[i],self.sulphide_cond[i],
					self.graphite_cond[i],self.ol_cond[i], self.mixture_cond[i], self.other_cond[i]]				
					
				while True:
				
					#while loop for deleting the zero arrays that could be encountered due to non-existence of the mineral.
					
					max_local = np.amax(np.asarray(list_i))
				
					if (max_local != 0.0):
						
						break
					
					else:
					
						list_i = np.delete(list_i, np.argwhere(list_i == 0))
						
				if SEEL.solid_phase_method == 1:
				
					self.bulk_cond[i] = (((self.granite_frac[i] / (self.granite_cond[i] + (2*max_local))) +\
					(self.granulite_frac[i] / (self.granulite_cond[i] + (2*max_local))) +\
					(self.sandstone_frac[i] / (self.sandstone_cond[i] + (2*max_local))) +\
					(self.gneiss_frac[i] / (self.gneiss_cond[i] + (2*max_local))) +\
					(self.amphibolite_frac[i] / (self.amphibolite_cond[i] + (2*max_local))) +\
					(self.basalt_frac[i] / (self.basalt_cond[i] + (2*max_local))) +\
					(self.mud_frac[i] / (self.mud_cond[i] + (2*max_local))) +\
					(self.gabbro_frac[i] / (self.gabbro_cond[i] + (2*max_local))) +\
					(self.other_rock_frac[i] / (self.other_rock_cond[i] + (2*max_local))))**(-1.0)) -\
					2.0*max_local
						
				elif SEEL.solid_phase_method == 2:
				
					self.bulk_cond[i] = (((self.quartz_frac[i] / (self.quartz_cond[i] + (2*max_local))) +\
					(self.plag_frac[i] / (self.plag_cond[i] + (2*max_local))) +\
					(self.amp_frac[i] / (self.amp_cond[i] + (2*max_local))) +\
					(self.kfelds_frac[i] / (self.kfelds_cond[i] + (2*max_local))) +\
					(self.opx_frac[i] / (self.opx_cond[i] + (2*max_local))) +\
					(self.cpx_frac[i] / (self.cpx_cond[i] + (2*max_local))) +\
					(self.mica_frac[i] / (self.mica_cond[i] + (2*max_local))) +\
					(self.garnet_frac[i] / (self.garnet_cond[i] + (2*min_local))) +\
					(self.sulphide_frac[i] / (self.sulphide_cond[i] + (2*min_local))) +\
					(self.graphite_frac[i] / (self.graphite_cond[i] + (2*min_local))) +\
					(self.ol_frac[i] / (self.ol_cond[i] + (2*min_local))) +\
					(self.mixture_frac[i] / (self.mixture_cond[i] + (2*min_local))) +\
					(self.other_frac[i] / (self.other_cond[i] + (2*min_local)))					
					)**(-1.0)) -\
					2.0*max_local
					
		elif method == 3:
		
			#Parallel model for maximum, minimum bounds and neutral w/o errors
			
			if SEEL.solid_phase_method == 1:
				self.bulk_cond[idx_node] = (self.granite_frac[idx_node]*self.granite_cond[idx_node]) +\
				(self.granulite_frac[idx_node]*self.granulite_cond[idx_node]) +\
				(self.sandstone_frac[idx_node]*self.sandstone_cond[idx_node]) +\
				(self.gneiss_frac[idx_node]*self.gneiss_cond[idx_node]) +\
				(self.amphibolite_frac[idx_node]*self.amphibolite_cond[idx_node]) +\
				(self.basalt_frac[idx_node]*self.basalt_cond[idx_node]) +\
				(self.mud_frac[idx_node]*self.mud_cond[idx_node]) +\
				(self.gabbro_frac[idx_node]*self.gabbro_cond[idx_node]) +\
				(self.other_rock_frac[idx_node]*self.other_rock_cond[idx_node])
				
			elif SEEL.solid_phase_method == 2:
			
				self.bulk_cond[idx_node] = (self.quartz_frac[idx_node]*self.quartz_cond[idx_node]) +\
				(self.plag_frac[idx_node]*self.plag_cond[idx_node]) +\
				(self.amp_frac[idx_node]*self.amp_cond[idx_node]) +\
				(self.kfelds_frac[idx_node]*self.kfelds_cond[idx_node]) +\
				(self.opx_frac[idx_node]*self.opx_cond[idx_node]) +\
				(self.cpx_frac[idx_node]*self.cpx_cond[idx_node]) +\
				(self.mica_frac[idx_node]*self.mica_cond[idx_node]) +\
				(self.garnet_frac[idx_node]*self.garnet_cond[idx_node]) +\
				(self.sulphide_frac[idx_node]*self.sulphide_cond[idx_node]) +\
				(self.graphite_frac[idx_node]*self.graphite_cond[idx_node]) +\
				(self.ol_frac[idx_node]*self.ol_cond[idx_node]) +\
				(self.mixture_frac[idx_node]*self.mixture_cond[idx_node]) +\
				(self.other_frac[idx_node]*self.other_cond[idx_node])
				
				
		elif method == 4:
		
			if indexing_method == 'array':
				self.bulk_cond = np.zeros(len(self.ol_cond))
				start_idx = 0
				end_idx = len(self.T)
			elif indexing_method == 'index':
				start_idx = sol_idx
				end_idx = sol_idx + 1

			#Perpendicular model for maximum, minimum bounds and neutral w/o errors				
			if SEEL.solid_phase_method == 1:
				for i in range(start_idx,end_idx):
					if self.granite_frac[i] == 0.0:
						self.granite_cond[i] = -999
					if self.granulite_frac[i] == 0.0:
						self.granulite_cond[i] = -999
					if self.sandstone_frac[i] == 0.0:
						self.sandstone_cond[i] = -999
					if self.gneiss_frac[i] == 0.0:
						self.gneiss_cond[i] = -999
					if self.amphibolite_frac[i] == 0.0:
						self.amphibolite_cond[i] = -999
					if self.basalt_frac[i] == 0.0:
						self.basalt_cond[i] = -999
					if self.mud_frac[i] == 0.0:
						self.mud_cond[i] = -999
					if self.gabbro_frac[i] == 0.0:
						self.gabbro_cond[i] = -999
					if self.other_rock_frac[i] == 0.0:
						self.other_rock_cond[i] = -999
	
				self.bulk_cond[idx_node] = 1.0 / ((self.granite_frac[idx_node] / self.granite_cond[idx_node]) +\
				(self.granulite_frac[idx_node] / self.granulite_cond[idx_node]) +\
				(self.sandstone_frac[idx_node] / self.sandstone_cond[idx_node]) +\
				(self.gneiss_frac[idx_node] / self.gneiss_cond[idx_node]) +\
				(self.amphibolite_frac[idx_node] / self.amphibolite_cond[idx_node]) +\
				(self.basalt_frac[idx_node] / self.basalt_cond[idx_node]) +\
				(self.mud_frac[idx_node] / self.mud_cond[idx_node]) +\
				(self.gabbro_frac[idx_node] / self.gabbro_cond[idx_node]) +\
				(self.other_rock_frac[idx_node] / self.other_rock_cond[idx_node]))
				
			elif SEEL.solid_phase_method == 2:
			
				for i in range(start_idx,end_idx):
					if self.quartz_frac[i] == 0.0:
						self.quartz_cond[i] = -999
					if self.plag_frac[i] == 0.0:
						self.plag_cond[i] = -999
					if self.amp_frac[i] == 0.0:
						self.amp_cond[i] = -999
					if self.kfelds_frac[i] == 0.0:
						self.kfelds_cond[i] = -999
					if self.opx_frac[i] == 0.0:
						self.opx_cond[i] = -999
					if self.cpx_frac[i] == 0.0:
						self.cpx_cond[i] = -999
					if self.mica_frac[i] == 0.0:
						self.mica_cond[i] = -999
					if self.garnet_frac[i] == 0.0:
						self.garnet_cond[i] = -999
					if self.sulphide_frac[i] == 0.0:
						self.sulphide_cond[i] = -999
					if self.graphite_frac[i] == 0.0:
						self.graphite_cond[i] = -999
					if self.ol_frac[i] == 0.0:
						self.ol_cond[i] = -999
					if self.mixture_frac[i] == 0.0:
						self.mixture_cond[i] = -999
					if self.other_frac[i] == 0.0:
						self.other_cond[i] = -999
	
				self.bulk_cond[idx_node] = 1.0 / ((self.quartz_frac[idx_node] / self.quartz_cond[idx_node]) +\
				(self.plag_frac[idx_node] / self.plag_cond[idx_node]) +\
				(self.amp_frac[idx_node] / self.amp_cond[idx_node]) +\
				(self.kfelds_frac[idx_node] / self.kfelds_cond[idx_node]) +\
				(self.opx_frac[idx_node] / self.opx_cond[idx_node]) +\
				(self.cpx_frac[idx_node] / self.cpx_cond[idx_node]) +\
				(self.mica_frac[idx_node] / self.mica_cond[idx_node]) +\
				(self.garnet_frac[idx_node] / self.garnet_cond[idx_node]) +\
				(self.sulphide_frac[idx_node] / self.sulphide_cond[idx_node]) +\
				(self.graphite_frac[idx_node] / self.graphite_cond[idx_node]) +\
				(self.ol_frac[idx_node] / self.ol_cond[idx_node]) +\
				(self.mixture_frac[idx_node] / self.mixture_cond[idx_node]) +\
				(self.other_frac[idx_node] / self.other_cond[idx_node]))
				
		elif method == 5:
		
			#Random model for maximum, minimum bounds and neutral w/o errors
			
			if SEEL.solid_phase_method == 1:
				
				self.bulk_cond[idx_node] = (self.granite_cond[idx_node]**self.granite_frac[idx_node]) *\
				(self.granulite_cond[idx_node]**self.granulite_frac[idx_node]) *\
				(self.sandstone_cond[idx_node]**self.sandstone_frac[idx_node]) *\
				(self.gneiss_cond[idx_node]**self.gneiss_frac[idx_node]) *\
				(self.amphibolite_cond[idx_node]**self.amphibolite_frac[idx_node]) *\
				(self.basalt_cond[idx_node]**self.basalt_frac[idx_node]) *\
				(self.mud_cond[idx_node]**self.mud_frac[idx_node]) *\
				(self.gabbro_cond[idx_node]**self.gabbro_frac[idx_node]) *\
				(self.other_rock_cond[idx_node]**self.other_rock_frac[idx_node]) 
				
			elif SEEL.solid_phase_method == 2:

				self.bulk_cond[idx_node] = (self.quartz_cond[idx_node]**self.quartz_frac[idx_node]) *\
				(self.plag_cond[idx_node]**self.plag_frac[idx_node]) *\
				(self.amp_cond[idx_node]**self.amp_frac[idx_node]) *\
				(self.kfelds_cond[idx_node]**self.kfelds_frac[idx_node]) *\
				(self.opx_cond[idx_node]**self.opx_frac[idx_node]) *\
				(self.cpx_cond[idx_node]**self.cpx_frac[idx_node]) *\
				(self.mica_cond[idx_node]**self.mica_frac[idx_node]) *\
				(self.garnet_cond[idx_node]**self.garnet_frac[idx_node]) *\
				(self.sulphide_cond[idx_node]**self.sulphide_frac[idx_node]) *\
				(self.graphite_cond[idx_node]**self.graphite_frac[idx_node]) *\
				(self.ol_cond[idx_node]**self.ol_frac[idx_node]) *\
				(self.mixture_cond[idx_node]**self.mixture_frac[idx_node]) *\
				(self.other_cond[idx_node]**self.other_frac[idx_node])
				
		elif method == -1:
			
			#In case the bulk conductivity is determined by a solid phase conductivity entry...
			
			self.bulk_cond = self.bckgr_res

		self.solid_phase_cond = np.array(self.bulk_cond)
			
		#Calculations regarding solid phases and fluid phases mixing take place after this.
		#checking if there's any melt/fluid on the list at all.
		if np.mean(self.melt_fluid_mass_frac) != 0:
			
			if SEEL.fluid_or_melt_method == 0:
				
				dens = iapws.iapws08.SeaWater(T = self.T[idx_node], P = self.p[idx_node], S = 0)
				self.dens_melt_fluid[idx_node] = dens.rho / 1e3
				
			elif SEEL.fluid_or_melt_method == 1:
				
				self.dens_melt_dry = float(self.dens_mat[1][SEEL.melt_cond_selection]) / 1e3 #index 1 is equate to melt
				#Determining xvol, first have to calculate the density of the melt from Sifre et al. (2014)
				
				self.dens_melt_fluid[idx_node] = (((self.h2o_melt[idx_node] * 1e-4) / 1e2) * 1.4) +\
				(((self.co2_melt[idx_node] * 1e-4) / 1e2) * 2.4) + (1 - (((self.h2o_melt[idx_node] * 1e-4) +\
				(self.co2_melt[idx_node] * 1e-4)) / 1e2)) * self.dens_melt_dry #calculating how much volatiles changed its density
				
			if indexing_method == 'array':
				self.melt_fuid_frac = np.zeros(len(self.melt_fluid_mass_frac))
				start_idx = 0
				end_idx = len(self.T)
			elif indexing_method == 'index':
				start_idx = sol_idx
				end_idx = sol_idx + 1

			self.melt_fluid_frac = np.zeros(len(self.melt_fluid_mass_frac))

			for i in range(start_idx,end_idx):
				if self.melt_fluid_mass_frac[i] != 0.0:
					self.melt_fluid_frac[i] = 1.0 / (1 + (((1.0/self.melt_fluid_mass_frac[i]) - 1) * (self.dens_melt_fluid[i] / (self.density_solids[i] / 1e3))))
	
	
			if melt_method == 0:

				#Modified Archie's Law taken from Glover et al. (2000) from eq. 8

				for i in range(start_idx,end_idx):

					if self.melt_fluid_mass_frac[i] != 0.0:

						p = np.log10(1.0 - self.melt_fluid_frac[i]**SEEL.melt_fluid_m[i]) / np.log10(1.0 - self.melt_fluid_frac[i])

						self.bulk_cond[i] = (self.bulk_cond[i] * (1.0 - self.melt_fluid_frac[i])**p) + (self.melt_fluid_cond[i] * (self.melt_fluid_frac[i]**SEEL.melt_fluid_m[i]))
							
			elif melt_method == 1:

				#Tubes model for melt and solid mixture from ten Grotenhuis et al. (2005) eq.5

				self.bulk_cond[idx_node] = ((1.0/3.0) * self.melt_fluid_frac[idx_node] * self.melt_fluid_cond[idx_node]) + ((1.0 - self.melt_fluid_frac[idx_node]) * self.bulk_cond[idx_node])
				
			elif melt_method == 2:

				#Spheres model for melt ans solid mixture got from ten Grotenhuis et al. (2005), eq.3

				self.bulk_cond[idx_node] = self.melt_fluid_cond[idx_node] + ((1.0 - self.melt_fluid_frac[idx_node]) / ((1.0 / (self.bulk_cond[idx_node] - self.melt_fluid_cond[idx_node])) +\
				 	(self.melt_fluid_frac[idx_node] / (3.0 * self.melt_fluid_cond[idx_node]))))
			
			elif melt_method == 3:
			
				#Modified brick-layer model from Schilling et al. (1997)

				ones = (1.0 - self.melt_fluid_frac[idx_node])
				two_thirds = (1.0 - self.melt_fluid_frac[idx_node])**(2.0/3.0)

				self.bulk_cond[idx_node] = self.melt_fluid_cond[idx_node] * (((self.melt_fluid_cond[idx_node] * (two_thirds - 1.0)) - (self.bulk_cond[idx_node] * two_thirds)) /\
				((self.bulk_cond[idx_node] * (ones - two_thirds)) + (self.melt_fluid_cond[idx_node] * (two_thirds - ones - 1.0))))
				
			elif melt_method == 4:
			
				#Hashin-shtrikman upper bound from Glover et al. (2000)
				vol_matrix = 1.0 - self.melt_fluid_frac[idx_node]

				self.bulk_cond[idx_node] = self.melt_fluid_cond[idx_node] * (1 -\
				((3 * vol_matrix * (self.melt_fluid_cond[idx_node] - self.bulk_cond[idx_node])) /\
				(3 * self.melt_fluid_cond[idx_node] - (self.melt_fluid_frac[idx_node] * (self.melt_fluid_cond[idx_node] - self.bulk_cond[idx_node])))))
				
			elif melt_method == 5:
			
				#Hashin-shtrikman lower bound from Glover et al. (2000)
				vol_matrix = 1.0 - self.melt_fluid_frac[idx_node]

				self.bulk_cond[idx_node] = self.bulk_cond[idx_node] * (1 +\
				((3 * self.melt_fluid_frac[idx_node] * (self.melt_fluid_cond[idx_node] - self.bulk_cond[idx_node])) /\
				(3 * self.bulk_cond[idx_node] + (vol_matrix * (self.melt_fluid_cond[idx_node] - self.bulk_cond[idx_node])))))
								
	def calculate_conductivity(self, method = None):

		if method == 'index':

			index = 0

		elif method == 'array':

			index = None
		
		if np.mean(self.melt_fluid_mass_frac) != 0.0:
			if SEEL.fluid_or_melt_method == 0:
				self.melt_fluid_cond = self.calculate_fluids_conductivity(method= method, sol_idx = index)
			elif SEEL.fluid_or_melt_method == 1:
				self.melt_fluid_cond = self.calculate_melt_conductivity(method = method, sol_idx = index)
	
		
		if SEEL.solid_phase_method == 1:
		
			if np.mean(self.granite_frac) != 0:
				self.granite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 2, sol_idx = index)
			else:
				self.granite_cond = np.zeros(len(self.T))
				
			if np.mean(self.granulite_frac) != 0:
				self.granulite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 3, sol_idx = index)
			else:
				self.granulite_cond = np.zeros(len(self.T))
				
			if np.mean(self.sandstone_frac) != 0:
				self.sandstone_cond = self.calculate_rock_conductivity(method = method, rock_idx= 4, sol_idx = index)
			else:
				self.sandstone_cond = np.zeros(len(self.T))
				
			if np.mean(self.gneiss_frac) != 0:
				self.gneiss_cond = self.calculate_rock_conductivity(method = method, rock_idx= 5, sol_idx = index)
			else:
				self.gneiss_cond = np.zeros(len(self.T))
				
			if np.mean(self.amphibolite_frac) != 0:
				self.amphibolite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 6, sol_idx = index)
			else:
				self.amphibolite_cond = np.zeros(len(self.T))

			if np.mean(self.basalt_frac) != 0:
				self.basalt_cond = self.calculate_rock_conductivity(method = method, rock_idx= 7, sol_idx = index)
			else:
				self.basalt_cond = np.zeros(len(self.T))

			if np.mean(self.mud_frac) != 0:
				self.mud_cond = self.calculate_rock_conductivity(method = method, rock_idx= 8, sol_idx = index)
			else:
				self.mud_cond = np.zeros(len(self.T))

			if np.mean(self.gabbro_frac) != 0:
				self.gabbro_cond = self.calculate_rock_conductivity(method = method, rock_idx= 9, sol_idx = index)
			else:
				self.gabbro_cond = np.zeros(len(self.T))
				
			if np.mean(self.other_rock_frac) != 0:
				self.other_rock_cond = self.calculate_rock_conductivity(method = method, rock_idx= 10, sol_idx = index)
			else:
				self.other_rock_cond = np.zeros(len(self.T))
				
		
			self.phase_mixing_function(method == SEEL.phs_mix_method, melt_method = SEEL.phs_melt_mix_method, indexing_method= method, sol_idx = index)
			
		elif SEEL.solid_phase_method == 2:
		
			if np.mean(self.quartz_frac) != 0:
				self.quartz_cond = self.calculate_mineral_conductivity(method = method, min_idx= 11, sol_idx = index)
			else:
				self.quartz_cond = np.zeros(len(self.T))
				
			if np.mean(self.plag_frac) != 0:
				self.plag_cond = self.calculate_mineral_conductivity(method = method, min_idx= 12, sol_idx = index)
			else:
				self.plag_cond = np.zeros(len(self.T))
				
			if np.mean(self.amp_frac) != 0:
				self.amp_cond = self.calculate_mineral_conductivity(method = method, min_idx= 13, sol_idx = index)
			else:
				self.amp_cond = np.zeros(len(self.T))
				
			if np.mean(self.kfelds_frac) != 0:
				self.kfelds_cond = self.calculate_mineral_conductivity(method = method, min_idx= 14, sol_idx = index)
			else:
				self.kfelds_cond = np.zeros(len(self.T))
				
			if np.mean(self.opx_frac) != 0:
				self.opx_cond = self.calculate_mineral_conductivity(method = method, min_idx= 15, sol_idx = index)
			else:
				self.opx_cond = np.zeros(len(self.T))
				
			if np.mean(self.cpx_frac) != 0:
				self.cpx_cond = self.calculate_mineral_conductivity(method = method, min_idx= 16, sol_idx = index)
			else:
				self.cpx_cond = np.zeros(len(self.T))
				
			if np.mean(self.mica_frac) != 0:
				self.mica_cond = self.calculate_mineral_conductivity(method = method, min_idx= 17, sol_idx = index)
			else:
				self.mica_cond = np.zeros(len(self.T))
				
			if np.mean(self.garnet_frac) != 0:
				self.garnet_cond = self.calculate_mineral_conductivity(method = method, min_idx= 18, sol_idx = index)
			else:
				self.garnet_cond = np.zeros(len(self.T))

			if np.mean(self.sulphide_frac) != 0:
				self.sulphide_cond = self.calculate_mineral_conductivity(method = method, min_idx= 19, sol_idx = index)
			else:
				self.sulphide_cond = np.zeros(len(self.T))

			if np.mean(self.graphite_frac) != 0:
				self.graphite_cond = self.calculate_mineral_conductivity(method = method, min_idx= 20, sol_idx = index)
			else:
				self.graphite_cond = np.zeros(len(self.T))

			if np.mean(self.ol_frac) != 0:
				self.ol_cond = self.calculate_mineral_conductivity(method = method, min_idx= 21, sol_idx = index)
			else:
				self.ol_cond = np.zeros(len(self.T))
				
			if np.mean(self.mixture_frac) != 0:
				self.mixture_cond = self.calculate_mineral_conductivity(method = method, min_idx= 22, sol_idx = index)
			else:
				self.mixture_cond = np.zeros(len(self.T))
	
			if np.mean(self.other_frac) != 0:
				self.other_cond = self.calculate_mineral_conductivity(method = method, min_idx= 23, sol_idx = index)
			else:
				self.other_cond = np.zeros(len(self.T))
	
				
			self.phase_mixing_function(method == SEEL.phs_mix_method, melt_method = SEEL.phs_melt_mix_method, indexing_method= method, sol_idx = index)
		
		self.cond_calculated = True
		
	def calculate_fugacity(self,mode):

		#Function that calculates oxygen fugacity buffers from selection.

		self.A_list = [-999,-27489.0,-999,-24930.0,-30650.0]
		self.B_list = [-999,6.702,-999,9.36,8.92]
		self.C_list = [-999,0.055,-999,0.046,0.054]

		#OXYGEN FUGACITY BUFFER CONSTANTS in the lists above(self.A_list ...)
		#Index 0: FMQ:
		#Index 1: IW: Hirsch (1991)
		#Index 2: QIF:
		#Index 3: NNO: Li et al. (1998)
		#Index 4 MMO: Xu et al. (2000)

		self.A_FMQ_low = -26455.3
		self.A_FMQ_high = -25096.3
		self.B_FMQ_low = 10.344
		self.B_FMQ_high = 8.735
		self.C_FMQ_low = 0.092
		self.C_FMQ_high = 0.11
		self.T_crit = 846.0

		self.A_QIF_low = -29435.7
		self.A_QIF_high = -29520.8
		self.B_QIF_low = 7.391
		self.B_QIF_high = 7.492
		self.C_QIF_low = 0.044
		self.C_QIF_high = 0.05


		if (mode == 0):

			self.fo2 = np.zeros(len(self.T))

			for i in range(0,len(self.T)):

				if self.T[i] < self.T_crit:

					self.fo2[i] = 10**((self.A_FMQ_low / self.T[i]) + self.B_FMQ_low + ((self.C_FMQ_low * ((self.p[i]*1e4) - 1)) / self.T[i]))

				else:

					self.fo2[i] = 10**((self.A_FMQ_high / self.T[i]) + self.B_FMQ_high + ((self.C_FMQ_high * ((self.p[i]*1e4) - 1)) / self.T[i]))


		elif (mode == 2):

			self.fo2 = np.zeros(len(self.T))

			for i in range(0,len(self.T)):

				if self.T[i] < self.T_crit:

					self.fo2[i] = 10**((self.A_QIF_low / self.T[i]) + self.B_QIF_low + ((self.C_QIF_low * ((self.p[i]*1e4) - 1)) / self.T[i]))

				else:

					self.fo2[i] = 10**((self.A_QIF_high / self.T[i]) + self.B_QIF_high + ((self.C_QIF_high * ((self.p[i]*1e4) - 1)) / self.T[i]))

		else:

			self.fo2 = 10**((self.A_list[mode] / self.T) + self.B_list[mode] + ((self.C_list[mode] * ((self.p*1e4) - 1)) / self.T))

		#self.fo2 is in bars multiply by 1e5 for Pa and 1e-4 for GPa

		return self.fo2

	def savetextfile(self):

		if self.write_file_save_name != '':

			lines = ['Parameter,Value,Objects,Type,Description,Unit\n']

			for i in range(1,len(self.init_params)):

				if self.init_params[i][2] == 'SEEL':
					val = getattr(SEEL,self.init_params[i][0])
				elif self.init_params[i][2] == 'self':
					val = getattr(self,self.init_params[i][0])
				try:
					if ('frac' in self.init_params[i][0]) == True:
						lines.append(','.join((self.init_params[i][0],str(val[0] * 1e2),self.init_params[i][2],self.init_params[i][3],self.init_params[i][4], self.init_params[i][5] + '\n')))
					else:
						lines.append(','.join((self.init_params[i][0],str(val[0]),self.init_params[i][2],self.init_params[i][3],self.init_params[i][4], self.init_params[i][5] + '\n')))
				except TypeError:
					lines.append(','.join((self.init_params[i][0],str(int(val)),self.init_params[i][2],self.init_params[i][3],self.init_params[i][4], self.init_params[i][5] + '\n')))

			filesave_composition = open(self.write_file_save_name ,'w')
			filesave_composition.writelines(lines)
			filesave_composition.close()

			print("Files are saved at the chosen location ")

	def read_pt_file(self):

		data_pt = self.read_csv(filename = self.pt_file, delim = ',')

		self.p = np.zeros(len(data_pt) - 1)
		self.T = np.zeros(len(data_pt) - 1)
		self.depth = np.zeros(len(data_pt) - 1)

		for i in range(1,len(data_pt)):

			self.p[i-1] = data_pt[i][0]
			self.depth[i-1] = data_pt[i][1]
			self.T[i-1] = data_pt[i][2]

	def duplicate_composition_for_pt_file(self):

		for i in range(0,len(self.T)):

			if SEEL.solid_phase_method == 1:

				self.granite_frac = np.ones(len(self.T)) * self.granite_frac[0]
				self.granulite_frac = np.ones(len(self.T)) * self.granulite_frac[0]
				self.sandstone_frac = np.ones(len(self.T)) * self.sandstone_frac[0]
				self.gneiss_frac = np.ones(len(self.T)) * self.gneiss_frac[0]
				self.amphibolite_frac = np.ones(len(self.T)) * self.amphibolite_frac[0]
				self.basalt_frac = np.ones(len(self.T)) * self.basalt_frac[0]
				self.mud_frac = np.ones(len(self.T)) * self.mud_frac[0]
				self.gabbro_frac = np.ones(len(self.T)) * self.gabbro_frac[0]
				self.other_rock_frac = np.ones(len(self.T)) * self.other_rock_frac[0]

				SEEL.granite_m = np.ones(len(self.T)) * SEEL.granite_m[0]
				SEEL.granulite_m = np.ones(len(self.T)) * SEEL.granulite_m[0]
				SEEL.sandstone_m = np.ones(len(self.T)) * SEEL.sandstone_m[0]
				SEEL.gneiss_m = np.ones(len(self.T)) * SEEL.gneiss_m[0]
				SEEL.amphibolite_m = np.ones(len(self.T)) * SEEL.amphibolite_m[0]
				SEEL.basalt_m = np.ones(len(self.T)) * SEEL.basalt_m[0]
				SEEL.mud_m = np.ones(len(self.T)) * SEEL.mud_m[0]
				SEEL.gabbro_m = np.ones(len(self.T)) * SEEL.gabbro_m[0]
				SEEL.other_rock_m = np.ones(len(self.T)) * SEEL.other_rock_m[0]

			elif SEEL.solid_phase_method == 2:

				self.quartz_frac = np.ones(len(self.T)) * self.quartz_frac[0]
				self.plag_frac = np.ones(len(self.T)) * self.plag_frac[0]
				self.amp_frac = np.ones(len(self.T)) * self.amp_frac[0]
				self.kfelds_frac = np.ones(len(self.T)) * self.kfelds_frac[0]
				self.garnet_frac = np.ones(len(self.T)) * self.garnet_frac[0]
				self.pyx_frac = np.ones(len(self.T)) * self.pyx_frac[0]
				self.mica_frac = np.ones(len(self.T)) * self.mica_frac[0]
				self.clay_frac = np.ones(len(self.T)) * self.clay_frac[0]
				self.carbonate_frac = np.ones(len(self.T)) * self.carbonate_frac[0]
				self.graphite_frac = np.ones(len(self.T)) * self.graphite_frac[0]
				self.sulphide_frac = np.ones(len(self.T)) * self.sulphide_frac[0]
				self.other_frac = np.ones(len(self.T)) * self.other_frac[0]

				SEEL.quartz_m = np.ones(len(self.T)) * SEEL.quartz_m[0]
				SEEL.plag_m = np.ones(len(self.T)) * SEEL.plag_m[0]
				SEEL.amp_m = np.ones(len(self.T)) * SEEL.amp_m[0]
				SEEL.kfelds_m = np.ones(len(self.T)) * SEEL.kfelds_m[0]
				SEEL.garnet_m = np.ones(len(self.T)) * SEEL.garnet_m[0]
				SEEL.pyx_m = np.ones(len(self.T)) * SEEL.pyx_m[0]
				SEEL.mica_m = np.ones(len(self.T)) * SEEL.mica_m[0]
				SEEL.clay_m = np.ones(len(self.T)) * SEEL.clay_m[0]
				SEEL.carbonate_m = np.ones(len(self.T)) * SEEL.carbonate_m[0]
				SEEL.graphite_m = np.ones(len(self.T)) * SEEL.graphite_m[0]
				SEEL.sulphide_m = np.ones(len(self.T)) * SEEL.sulphide_m[0]
				SEEL.other_m = np.ones(len(self.T)) * SEEL.other_m[0]

			SEEL.melt_fluid_m = np.ones(len(self.T)) * SEEL.melt_fluid_m[0]
			self.melt_fluid_mass_frac = np.ones(len(self.T)) * self.melt_fluid_mass_frac[0]
			self.salinity_fluid = np.ones(len(self.T)) * self.salinity_fluid[0]
			self.co2_melt = np.ones(len(self.T)) * self.co2_melt[0]
			self.h2o_melt = np.ones(len(self.T)) * self.h2o_melt[0]
			self.na2o_melt = np.ones(len(self.T)) * self.na2o_melt[0]
			self.k2o_melt = np.ones(len(self.T)) * self.k2o_melt[0]

		
	def write_results(self):

		lines = ['Bulk_Cond[S/m],Melt_Fluid_Cond[S/m],Solid_Cond[S/m],T[K],Depth[km],P[GPa]\n']

		for i in range(0,len(self.T)):
			lines.append(','.join((str(self.bulk_cond[i]),str(self.melt_fluid_cond[i]),str(self.solid_phase_cond[i]), str(self.T[i]), str(self.depth[i]), str(self.p[i]) + '\n')))
		filesave_results = open(self.write_file_save_name ,'w')
		filesave_results.writelines(lines)
		filesave_results.close()
		print("Files are saved at the chosen location...")
