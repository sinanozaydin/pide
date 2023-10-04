#!/usr/bin/env python3

import os,sys
import h5py

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL/sel_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

import SEL
from sel_src.geodyn.read_uw_model import *

#setting up source folder of the files
source_folder = os.path.join('.','example_data','uwconversion')

#setting up filename folders for the h5 files.
temp_fnm = os.path.join(source_folder,'temperature-501.h5')
pstrain_fnm = os.path.join(source_folder,'projPlasticStrain-501.h5')
material_fnm = os.path.join(source_folder,'projMaterialField-501.h5')
stress_fnm = os.path.join(source_folder,'projStressField-501.h5')
melt_fnm = os.path.join(source_folder,'projMeltField-501.h5')
strain_rate_fnm = os.path.join(source_folder,'strainRateField-501.h5')
pressure_fnm = os.path.join(source_folder,'pressureField-501.h5')
mesh_fnm = os.path.join(source_folder,'mesh.h5')

py_start_fnm = os.path.join(source_folder,'Col262.py')

#reading h5 files related to 2D UW model
temp_data, pressure_data, mesh_data, material_data, pstrain_data, stress_data, melt_data, strain_rate_data = read_h5_files(temp_h5 = temp_fnm, pressure_h5 = pressure_fnm, mesh_h5 = mesh_fnm,
material_h5= material_fnm, strain_h5 = pstrain_fnm, stress_h5 = stress_fnm, melt_h5 = melt_fnm, strain_rate_h5= strain_rate_fnm)

#reading startup py file to get material order and properties.

material_names = read_uw_material_names_from_py_input(py_start_fnm)

#reading 2d mesh params mesh itself, mesh_centers mesh, array in x direction, array in y direction, borders of the mesh[max_x, min_x, max_y, min_y]
mesh, mesh_center, x_mesh, y_mesh, borders_mesh = setup_2d_mesh(mesh_data)

#getting material_array
material_array, air_material_idx = setup_material(material_data, material_names)

#getting strain array
pstrain_array = setup_uw_data_array_PROJ(pstrain_data)
stress_array = setup_uw_data_array_PROJ(stress_data)
melt_array = setup_uw_data_array_PROJ(melt_data)
strain_rate_array = setup_uw_data_array_PROJ(strain_rate_data)
temp_array = setup_uw_data_array_PROJ(temp_data) #mesh
pressure_array = setup_uw_data_array_PROJ(pressure_data) / 1e9 #converting to gigapascal
stress_array = setup_uw_data_array_PROJ(stress_data) 


#Plotting fields if want to.

# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = material_array,cblimit_up = len(material_names), cblimit_down = 0, log_bool=False, cb_name = 'tab20b')
# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = pstrain_array, cblimit_up = 1e2, cblimit_down = 1e-2, log_bool=True, cb_name = 'binary')
# plot_2D_underworld_Field(x_array = mesh[0], y_array = mesh[1], Field = temp_array, cblimit_up = 1500, cblimit_down = 600, log_bool=False, cb_name = 'coolwarm')
# plot_2D_underworld_Field(x_array = mesh_center[0], y_array = mesh_center[1], Field = pressure_array, cblimit_up = 10, cblimit_down = 0, log_bool=False, cb_name = 'coolwarm')


print(len(temp_array))
