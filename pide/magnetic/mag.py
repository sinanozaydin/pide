#!/usr/bin/env python3

import os, psutil
import numpy as np
import harmonica as hm

def calculate_3DModel_magnetic(sus_array, mesh, height, **kwargs):

	"""
	This is a function that converts the constructed 3D geodynamic model into prisms
	that can be used by Harmonica to calculate the magnetic response from the 3D susceptibility field
	
	Input Paramaters
	-----------------
	sus_array: calculated susceptibility field in np.ndarray()
	mesh: mesh of the associated susceptibility field.
	height: observation height to calculate magnetic field    
	
	mesh_unit: 'kilometres' or 'metres'
	
	Out
	-----------------
	Vertical magnetic field in nT
	"""
	
	#optional kwarg calls
	mesh_unit = kwargs.pop('mesh_unit', 'kilometres')

	#converting mesh in kilometers into meters.
	if mesh_unit == 'kilometres':
		mesh[0] = mesh[0] * 1e3
		mesh[1] = mesh[1] * 1e3
		mesh[2] = mesh[2] * 1e3
	elif mesh_unit == 'metres':
		pass
	else:
		raise ValueError('Please enter a valid response for mesh_unit: "kilometres" or "metres".')
		
	#get z vertical negtive downward
	mesh[2] = -mesh[2]
	
	#Match sus array with mesh
	sus_mesh = sus_array.reshape(mesh[0].shape[2],mesh[0].shape[0],mesh[0].shape[1])
	sus_mesh_t = np.transpose(sus_mesh,[1,2,0])
	
	#Find the mean value of sus
	sus_mean=sus_mesh_t.mean(axis=(0,1))
	sus_mean=np.broadcast_to(sus_mean[np.newaxis, np.newaxis,:], (sus_mesh_t.shape[0],sus_mesh_t.shape[1],sus_mean.size))
	sus_3d=sus_mesh_t-sus_mean
	#To do change sus to magnetization
	
	#Get grid size
	half_res_x=np.abs((np.unique(mesh[0].flatten())[0]-np.unique(mesh[0].flatten())[1])/2)
	half_res_y=np.abs((np.unique(mesh[1].flatten())[0]-np.unique(mesh[1].flatten())[1])/2)
	half_res_z=np.abs((np.unique(mesh[2].flatten())[0]-np.unique(mesh[2].flatten())[1])/2)

	#Build prism layer for Harmonica input
	prisms_rock = np.array([
		mesh[0].flatten() - half_res_x,
		mesh[0].flatten() + half_res_x,
		mesh[1].flatten() - half_res_y,
		mesh[1].flatten() + half_res_y,
		mesh[2].flatten() - half_res_z,
		mesh[2].flatten() + half_res_z]
	)
	
	#Calcuate gravity using harmonica
	pred_coords = (mesh[0][:,:,0], mesh[1][:,:,0], mesh[2][:,:,0]-mesh[2][:,:,0]+height)
	mag = hm.prism_magnetic(pred_coords, prisms_rock.T,np.array([sus_3d.flatten(),sus_3d.flatten(),sus_3d.flatten()]).T)
	return mag[2]