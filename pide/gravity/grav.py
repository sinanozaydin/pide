#!/usr/bin/env python3


import os, psutil
import numpy as np
import harmonica as hm

def calculate_3DModel_gravity(density_array, mesh, height, num_cpu, gravity_type='free_air', **kwargs):


	#This is a function that converts the constructed 3D geodynamic model into prisms
	#that can be used by Harmonica to calculate the gravity response from the 3D density field
	"""
	Input Paramaters
	-----------------
	density_array: calculated density field in np.ndarray()
	mesh: mesh of the associated density field.
	height: observation height to calculate gravity	
	
	gravity_type: 'free_air' or 'bouguer'
	mesh_unit: 'kilometres' or 'metres'
	"""
	
	#optional kwarg calls
	mesh_unit = kwargs.pop('mesh_unit', 'kilometres')

	#converting mesh in kilometers into meters.
	if mesh_unit == 'kilometres':
		mesh[0] = mesh[0] * 1e3
		mesh[1] = mesh[1] * 1e3
		mesh[2] = mesh[2] * 1e3
		
		height = height * 1e3
	elif mesh_unit == 'metres':
		pass
	else:
		raise ValueError('Please enter a valid response for mesh_unit: "kilometres" or "metres".')
	#get z vertical negtive downward
	mesh[2] = -mesh[2]
	
	#Match density array with mesh
	density_mesh = density_array.reshape(mesh[0].shape[2],mesh[0].shape[0],mesh[0].shape[1])
	density_mesh_t = np.transpose(density_mesh,[1,2,0])
	
	#Find the mean value of density
	density_mean=density_mesh_t.mean(axis=(0,1))
	density_mean=np.broadcast_to(density_mean[np.newaxis, np.newaxis,:], (density_mesh_t.shape[0],density_mesh_t.shape[1],density_mean.size))
	density_mean_copy=np.copy(density_mean)
	#Chose to calculate Bouguer or Free-air gravity
	if gravity_type == 'bouguer':
		density_mean_copy[mesh[2] >= 0]=2670
		density_3d=density_mesh_t-density_mean_copy
	elif gravity_type == 'free_air':
		density_mean_copy[mesh[2] >= 0]=0
		density_3d=density_mesh_t-density_mean_copy
	else:
		raise ValueError("Invalid datatype. Accepted values are 'free_air' or 'bouguer'.")
	
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
	
	print('Started calculation of gravity using harmonica...')
	
	if num_cpu == 1:
		g_z = hm.prism_gravity(pred_coords, prisms_rock.T, density_3d, field='g_z', parallel=False)
		
	else:
		#Get the process ID of the current process
		pid = os.getpid()
		
		p = psutil.Process(pid)
	
		total_cpu_count = psutil.cpu_count()
		
		cpu_affinity = list(range(num_cpu))
		
		cpu_affinity = cpu_affinity[:total_cpu_count]
		p.cpu_affinity(cpu_affinity)
		
		g_z = hm.prism_gravity(pred_coords, prisms_rock.T, density_3d, field='g_z')
		
	return g_z