#!/usr/bin/env python3
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def plot_2D_underworld_Field_scatter(x_array = None, y_array = None, Field = None,cblimit_up = None, cblimit_down = None, log_bool = False,cb_name = 'coolwarm',**kwargs):
	
	"""A function to plot 2D underworld fields
	
	Input:
	array: x_array - X array || np.meshgrid
	array: y_array - Y arrray || np.meshgrid
	array: Field - value array.
	float: cblimit_up - colorbar upper limit
	float: cblimit_down - colobar lower limit
	bool: log_bool - boolean to determine whether the colobar will be logarithmic.
	str: cb_name - colorbar to use
	bool: plot_save - boolean to save the figure or not.
	str: label - name of the figure when saving it.
	
	"""
	
	plot_save = kwargs.pop('plot_save', False)
	label = kwargs.pop('label', 'Interpolated_UW_Figure.png')

	fig = plt.figure(figsize = (12,7))
	ax = plt.subplot(111)
	ax.set_ylim(np.amax(y_array),np.amin(y_array))
	ax.set_xlim(np.amin(x_array),np.amax(x_array))
	ax.set_ylabel('Depth [km]')
	ax.set_xlabel('Distance [km]')
	if log_bool == True:
		cax = ax.scatter(x_array, y_array ,c = Field, cmap = cb_name, norm=colors.LogNorm(), marker = 's', linewidth = 0.005, edgecolor = 'k')
	elif log_bool == False:
		cax = ax.scatter(x_array, y_array ,c = Field, cmap = cb_name, marker = 's', linewidth = 0.005, edgecolor = 'k')
	cax.set_clim(cblimit_down,cblimit_up)

	if log_bool == True:
		bondary = np.logspace(np.log10(cblimit_down),np.log10(cblimit_up))
		tick_array = np.arange(np.log10(cblimit_down),np.log10(cblimit_up)+1, 1)
		tick_array_list = 10.0**tick_array
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05,
		ticks = tick_array_list, ax = ax)
		
	elif log_bool == False:
		bondary = np.linspace(cblimit_down, cblimit_up)
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05, ax = ax)
		
	if plot_save == False:
		plt.show()
	elif plot_save == True:
		plt.savefig(label, dpi = 300)
		print('The file is saved as: ' + label + ' at location: ' + os.getcwd())
		

def plot_2D_underworld_Field(xmesh = None, ymesh = None, Field = None,cblimit_up = None, cblimit_down = None, log_bool = False,cb_name = 'coolwarm',**kwargs):

	"""
	plots the field of 2D underworld field
	
	Parameters
	----------
	
	xmesh: mesh in x direction in np.array
	ymesh: mesh in y direction in np.array
	Field: values in UW format 2D tuple of np.arrays
	cblimit_up: colorbar upper limit in float
	cblimit_down: colorbar lower limit in float
	log_bool: True-Logarithmic colorbar, False-normal colorbar in boolean
	cb_name: colorbar style, default:'coolwarm' in str
	
	**kwargs:
	
	plot_save: False-showing the plot, True saving image
	label: label of the saved image. try to put something.png in str
	cbar_label: colorbar label in str
	
	"""

	
	
	plot_save = kwargs.pop('plot_save', False)
	label = kwargs.pop('label', 'Interpolated_UW_Figure.png')
	cbar_label = kwargs.pop('cbar_label',None)
	
	fnew = Field
	
	contains_nan = np.isnan(fnew).any()
	
	if contains_nan == True:
	
		from scipy.interpolate import griddata
		#interpolating 
		
		xi = xmesh
		yi = ymesh
		x_i,  y_i = np.meshgrid(xmesh,ymesh)
		
		mask = np.isnan(fnew)
		points = np.column_stack((x_i[~mask], y_i[~mask]))
		values = fnew[~mask]
				
		fnew = griddata(points, values, (x_i, y_i),method = 'linear')
		
	else:
	
		xi = xmesh
		yi = ymesh

	fig = plt.figure(figsize = (12,7))
	ax = plt.subplot(111)
	if log_bool == True:
		cax = ax.pcolormesh(xi, yi, fnew, cmap = cb_name, norm=colors.LogNorm())
	elif log_bool == False:
		cax = ax.pcolor(xi, yi, fnew, cmap = cb_name)
		
	cax.set_clim(cblimit_down,cblimit_up)
	ax.set_ylim(np.amax(ymesh),np.amin(ymesh))
	ax.set_xlim(np.amin(xmesh),np.amax(xmesh))
	ax.set_ylabel('Depth [km]')
	ax.set_xlabel('Distance [km]')
	
	if log_bool == True:
		bondary = np.logspace(np.log10(cblimit_down),np.log10(cblimit_up))
		tick_array = np.arange(np.log10(cblimit_down),np.log10(cblimit_up)+1, 1)
		tick_array_list = 10.0**tick_array
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05,
		ticks = tick_array_list, ax = ax,label = cbar_label)
	elif log_bool == False:
		bondary = np.linspace(cblimit_down, cblimit_up)
		cbar_cax = fig.colorbar(cax,boundaries=bondary ,orientation="vertical", pad=0.05, ax = ax,label = cbar_label)
		
	if plot_save == False:
		plt.show()
	elif plot_save == True:
		plt.savefig(label, dpi = 300)
		print('The file is saved as: ' + label + ' at location: ' + os.getcwd())
	