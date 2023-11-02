#!/usr/bin/env python3

import os, csv
import numpy as np

def convert_model2MARE2DEM(filename, conductivity_array, mesh, **kwargs):

	cond_type = kwargs.pop('cond_type', 'conductivity')
	
	if cond_type == 'conductivity':
		conductivity_array = 1.0 / conductivity_array
	elif cond_type == 'resistivity':
		pass
	else:
		raise ValueError('Please enter a valid response for cond_type: "conductivity" or "resistivity".')
	
	
	