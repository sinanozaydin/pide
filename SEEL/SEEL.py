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

	def check_composition(self):

		continue_adjusting = True

		if SEEL.solid_phase_method == 0:

			pass

		elif SEEL.solid_phase_method == 1:

			tot = self.granite_frac[0] + self.granulite_frac[0] + self.sandstone_frac[0] +\
			self.gneiss_frac[0] + self.amphibolite_frac[0] + self.basalt_frac[0] + self.mud_frac[0] +\
				 self.gabbro_frac[0] + self.other_rock_frac[0]

			if (tot <= 0.99) and (tot >= 1.01):
				QMessageBox.about(self, "Warning!", "The total number of does not add up to 100%. Currently it is:  " + str(tot*1e2))
				continue_adjusting = False
		elif SEEL.solid_phase_method == 2:

			tot = self.quartz_frac[0] + self.plag_frac[0] + self.amp_frac[0] + self.kfelds_frac[0] +\
			self.opx_frac[0] + self.cpx_frac[0] + self.mica_frac[0] + self.garnet_frac[0] + self.sulphide_frac[0] + self.graphite_frac[0] +\
			self.ol_frac[0] + self.mixture_frac[0] + self.other_frac[0]

			if (tot <= 0.99) and (tot > 1.01):
				QMessageBox.about(self, "Warning!", "The total number of does not add up to 100%. Currently it is:  " + str(tot*1e2))
				continue_adjusting = False

		return continue_adjusting
	
	def calculate_arrhenian_single(self, T, sigma, E, r=0, water=1):

		cond = (10.0**sigma) * (water**r) * np.exp((-E) / (self.R * T))

		return cond
	
	def calculate_fluids_conductivity(self, method, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond_fluids = np.zeros(len(self.T))

		if SEEL.type[0][SEEL.fluid_cond_selection] == '0':

			self.melt_fluid_cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[0][SEEL.fluid_cond_selection],
								   E = self.h_pol[0][SEEL.fluid_cond_selection])
			
		elif SEEL.type[0][SEEL.fluid_cond_selection] == '1':

			self.melt_fluid_cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[0][SEEL.fluid_cond_selection],
								   E = self.h_pol[0][SEEL.fluid_cond_selection]) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[0][SEEL.fluid_cond_selection],
								   E = self.h_p[0][SEEL.fluid_cond_selection])
			
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

		if SEEL.type[1][SEEL.melt_cond_selection] == '0':

			self.melt_fluid_cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[1][SEEL.melt_cond_selection],
								   E = self.h_pol[1][SEEL.melt_cond_selection])
			
		elif SEEL.type[1][SEEL.melt_cond_selection] == '1':

			self.melt_fluid_cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[1][SEEL.melt_cond_selection],
								   E = self.h_pol[1][SEEL.melt_cond_selection]) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[1][SEEL.melt_cond_selection],
								   E = self.h_p[1][SEEL.melt_cond_selection])
			
		elif SEEL.type[1][SEEL.melt_cond_selection] == '3':

			if ('*' in SEEL.name[1][SEEL.melt_cond_selection]) == True:

				melt_odd_function = SEEL.name[1][SEEL.melt_cond_selection].replace('*','')

			else:

				melt_odd_function = SEEL.name[1][SEEL.melt_cond_selection]

			cond_melt[idx_node] = eval(melt_odd_function + '(T = self.T[idx_node], P = self.p[idx_node], Melt_H2O = self.h2o_melt[idx_node],' +
			'Melt_CO2 = self.co2_melt, Melt_Na2O = self.na2o_melt[idx_node], Melt_K2O = self.k2o_melt[idx_node], method = method)')
		
		return cond_melt

	def calculate_rock_conductivity(self, method, rock_idx = None, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond = np.zeros(len(self.T))

		rock_sub_idx = rock_idx - self.fluid_num

		if SEEL.type[rock_idx][self.rock_cond_selections[rock_sub_idx]] == '0':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[rock_idx][self.rock_cond_selections[rock_sub_idx]],
								   E = self.h_pol[rock_idx][self.rock_cond_selections[rock_sub_idx]])
			
		elif SEEL.type[rock_idx][self.rock_cond_selections[rock_sub_idx]] == '1':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[rock_idx][self.rock_cond_selections[rock_sub_idx]],
								   E = self.h_pol[rock_idx][self.rock_cond_selections[rock_sub_idx]]) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[rock_idx][self.rock_cond_selections[rock_sub_idx]],
								   E = self.h_p[rock_idx][self.rock_cond_selections[rock_sub_idx]])
			
		elif SEEL.type[rock_idx][self.rock_cond_selections[rock_sub_idx]] == '3':

			if ('*' in SEEL.name[rock_idx][self.rock_cond_selections[rock_sub_idx]]) == True:

				odd_function = SEEL.name[rock_idx][self.rock_cond_selections[rock_sub_idx]].replace('*','')

			else:

				odd_function = SEEL.name[rock_idx][self.rock_cond_selections[rock_sub_idx]]

			cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], method = method)')

		return cond
	
	def calculate_mineral_conductivity(self, method, min_idx = None, sol_idx = None):

		if method == 'array':
			idx_node = None
		elif method == 'index':
			idx_node = sol_idx

		cond = np.zeros(len(self.T))

		min_sub_idx = min_idx - self.fluid_num - self.rock_num

		if SEEL.type[min_idx][self.minerals_cond_selections[min_sub_idx]] == '0':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[min_idx][self.minerals_cond_selections[min_sub_idx]],
								   E = self.h_pol[min_idx][self.minerals_cond_selections[min_sub_idx]])
			
		elif SEEL.type[min_idx][self.minerals_cond_selections[min_sub_idx]] == '1':

			cond[idx_node] = self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_pol[min_idx][self.minerals_cond_selections[min_sub_idx]],
								   E = self.h_pol[min_idx][self.minerals_cond_selections[min_sub_idx]]) + self.calculate_arrhenian_single(T = self.T[idx_node],
								   sigma = self.sigma_p[min_idx][self.minerals_cond_selections[min_sub_idx]],
								   E = self.h_p[min_idx][self.minerals_cond_selections[min_sub_idx]])
			
		elif SEEL.type[min_idx][self.minerals_cond_selections[min_sub_idx]] == '3':

			if ('*' in SEEL.name[min_idx][self.minerals_cond_selections[min_sub_idx]]) == True:

				odd_function = SEEL.name[min_idx][self.minerals_cond_selections[min_sub_idx]].replace('*','')

			else:

				odd_function = SEEL.name[min_idx][self.minerals_cond_selections[min_sub_idx]]

			cond[idx_node] = eval(odd_function + '(T = self.T[idx_node], P = self.p[idx_node], method = method)')

		return cond
	
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

		if SEEL.fluid_or_melt_method == 0:
			self.melt_fluid_cond = self.calculate_fluids_conductivity(method= method, sol_idx = index)
		elif SEEL.fluid_or_melt_method == 1:
			self.melt_fluid_cond = self.calculate_melt_conductivity(method = method, sol_idx = index)
	
		if SEEL.solid_phase_method == 0:

			self.phase_mixing_function(method == -1, melt_method = SEEL.phs_melt_mix_method, indexing_method= method, sol_idx = index)
			
		elif SEEL.solid_phase_method == 1:
		
			if self.granite_frac[0] != 0:
				self.granite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 2, sol_idx = index)
			else:
				self.granite_cond = np.array([0])
				
			if self.granulite_frac[0] != 0:
				self.granulite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 3, sol_idx = index)
			else:
				self.granulite_cond = np.array([0])
				
			if self.sandstone_frac[0] != 0:
				self.sandstone_cond = self.calculate_rock_conductivity(method = method, rock_idx= 4, sol_idx = index)
			else:
				self.sandstone_cond = np.array([0])
				
			if self.gneiss_frac[0] != 0:
				self.gneiss_cond = self.calculate_rock_conductivity(method = method, rock_idx= 5, sol_idx = index)
			else:
				self.gneiss_cond = np.array([0])
				
			if self.amphibolite_frac[0] != 0:
				self.amphibolite_cond = self.calculate_rock_conductivity(method = method, rock_idx= 6, sol_idx = index)
			else:
				self.amphibolite_cond = np.array([0])

			if self.basalt_frac[0] != 0:
				self.basalt_cond = self.calculate_rock_conductivity(method = method, rock_idx= 7, sol_idx = index)
			else:
				self.basalt_cond = np.array([0])

			if self.mud_frac[0] != 0:
				self.mud_cond = self.calculate_rock_conductivity(method = method, rock_idx= 8, sol_idx = index)
			else:
				self.mud_cond = np.array([0])

			if self.gabbro_frac[0] != 0:
				self.gabbro_cond = self.calculate_rock_conductivity(method = method, rock_idx= 9, sol_idx = index)
			else:
				self.gabbro_cond = np.array([0])
				
			if self.other_rock_frac[0] != 0:
				self.other_rock_cond = self.calculate_rock_conductivity(method = method, rock_idx= 10, sol_idx = index)
			else:
				self.other_rock_cond = np.array([0])
				
		
			self.phase_mixing_function(method == SEEL.phs_mix_method, melt_method = SEEL.phs_melt_mix_method, indexing_method= method, sol_idx = index)
			
		elif SEEL.solid_phase_method == 2:
		
			if self.quartz_frac[0] != 0:
				self.quartz_cond = self.calculate_mineral_conductivity(method = method, min_idx= 11, sol_idx = index)
			else:
				self.quartz_cond = np.array([0])
				
			if self.plag_frac[0] != 0:
				self.plag_cond = self.calculate_mineral_conductivity(method = method, min_idx= 12, sol_idx = index)
			else:
				self.plag_cond = np.array([0])
				
			if self.amp_frac[0] != 0:
				self.amp_cond = self.calculate_mineral_conductivity(method = method, min_idx= 13, sol_idx = index)
			else:
				self.amp_cond = np.array([0])
				
			if self.kfelds_frac[0] != 0:
				self.kfelds_cond = self.calculate_mineral_conductivity(method = method, min_idx= 14, sol_idx = index)
			else:
				self.kfelds_cond = np.array([0])
				
			if self.opx_frac[0] != 0:
				self.opx_cond = self.calculate_mineral_conductivity(method = method, min_idx= 15, sol_idx = index)
			else:
				self.opx_cond = np.array([0])
				
			if self.cpx_frac[0] != 0:
				self.cpx_cond = self.calculate_mineral_conductivity(method = method, min_idx= 16, sol_idx = index)
			else:
				self.cpx_cond = np.array([0])
				
			if self.mica_frac[0] != 0:
				self.mica_cond = self.calculate_mineral_conductivity(method = method, min_idx= 17, sol_idx = index)
			else:
				self.mica_cond = np.array([0])
				
			if self.garnet_frac[0] != 0:
				self.garnet_cond = self.calculate_mineral_conductivity(method = method, min_idx= 18, sol_idx = index)
			else:
				self.garnet_cond = np.array([0])

			if self.sulphide_frac[0] != 0:
				self.sulphide_cond = self.calculate_mineral_conductivity(method = method, min_idx= 19, sol_idx = index)
			else:
				self.sulphide_cond = np.array([0])

			if self.graphite_frac[0] != 0:
				self.graphite_cond = self.calculate_mineral_conductivity(method = method, min_idx= 20, sol_idx = index)
			else:
				self.graphite_cond = np.array([0])

			if self.ol_frac[0] != 0:
				self.ol_cond = self.calculate_mineral_conductivity(method = method, min_idx= 21, sol_idx = index)
			else:
				self.ol_cond = np.array([0])
				
			if self.mixture_frac[0] != 0:
				self.mixture_cond = self.calculate_mineral_conductivity(method = method, min_idx= 22, sol_idx = index)
			else:
				self.mixture_cond = np.array([0])
	
			if self.other_frac[0] != 0:
				self.other_cond = self.calculate_mineral_conductivity(method = method, min_idx= 23, sol_idx = index)
			else:
				self.other_cond = np.array([0])
	
				
			self.phase_mixing_function(method == SEEL.phs_mix_method, melt_method = SEEL.phs_melt_mix_method, indexing_method= method, sol_idx = index)
		
		self.cond_calculated = True

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
