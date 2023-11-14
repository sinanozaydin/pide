#!/usr/bin/env python3

import os,sys
import h5py
import numpy as np

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL/sel_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

import SEL
from material import Material
from model import Model
from geodyn.read_uw_model import *
from geodyn.interpolate_fields import interpolate_2d_fields
from geodyn.plot_models import *
from geodyn.material_process import *
from mt.mt_model_conversion import convert_2DModel_2_MARE2DEM

source_folder = "/media/sinan/Geodyn_HDD/Geodyn/Test000q/Test000q"

#setting up filename folders for the h5 files.
temp_fnm = os.path.join(source_folder,'temperature-0.h5')
pstrain_fnm = os.path.join(source_folder,'projPlasticStrain-0.h5')
material_fnm = os.path.join(source_folder,'projMaterialField-0.h5')
stress_fnm = os.path.join(source_folder,'projStressField-0.h5')
melt_fnm = os.path.join(source_folder,'projMeltField-0.h5')
strain_rate_fnm = os.path.join(source_folder,'strainRateField-0.h5')
pressure_fnm = os.path.join(source_folder,'pressureField-0.h5')
mesh_fnm = os.path.join(source_folder,'mesh.h5')

py_start_fnm = os.path.join(source_folder,'..','Test000q.py')

temp_data = read_h5_file(temp_fnm)
pressure_data = read_h5_file(pressure_fnm)
mesh_data = read_h5_file(mesh_fnm)
material_data = read_h5_file(material_fnm)
pstrain_data = read_h5_file(pstrain_fnm)
stress_data = read_h5_file(stress_fnm)
melt_data = read_h5_file(melt_fnm)
strain_rate_data = read_h5_file(strain_rate_fnm)

#reading startup py file to get material order and properties.

material_names = read_uw_material_names_from_py_input(py_start_fnm)
print(material_names)

#reading 2d mesh params mesh itself, mesh_centers mesh, array in x direction, array in y direction, borders of the mesh[max_x, min_x, max_y, min_y]
mesh, mesh_center, x_mesh, y_mesh, x_mesh_centers, y_mesh_centers, borders_mesh = setup_2d_mesh(mesh_data)

#getting material_array
material_array, air_material_idx = setup_material(material_data, material_names)
sys.exit()
#getting strain array
pstrain_array = setup_uw_data_array_PROJ(pstrain_data)
stress_array = setup_uw_data_array_PROJ(stress_data)
melt_array = setup_uw_data_array_PROJ(melt_data)
strain_rate_array = setup_uw_data_array_PROJ(strain_rate_data)

temp_array = setup_uw_data_array_PROJ(temp_data) #mesh
pressure_array = setup_uw_data_array_PROJ(pressure_data) / 1e9 #converting to gigapascal
