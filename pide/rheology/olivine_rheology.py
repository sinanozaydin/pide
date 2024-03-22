#!/usr/bin/env python3

import numpy


class olivine_rheology(object):

	def __init__(self, difffusion_model = None, dislocation_model = None, GBS_model = None):
	
		self.R_const = 8.3144621
		self.Water_cutoff_rate = 5.0 #ppm
		self.Water_cutoff_rate_hsi = self.Water_cutoff_rate * (4.39e4) / 2695.0 
