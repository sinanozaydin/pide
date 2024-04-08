#!/usr/bin/env python3

import os,sys
import h5py
import numpy as np

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')
core_path_ext_2 = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide/pide_src')
sys.path.append(core_path_ext)
sys.path.append(core_path_ext_2)

from geodyn.read_uw_model import *
from geodyn.interpolate_fields import interpolate_3d_fields
from gravity.calculate_3DModel_gravity import calculate_3DModel_gravity

import pide

source_folder = "/home/sinan/Desktop/Research/SEL/3d_work_dir/ProtPA300r"

mesh_fnm = os.path.join(source_folder,'mesh.h5')
density_fnm = os.path.join(source_folder,'projDensityField-158.h5')

mesh_data = read_h5_file(mesh_fnm)
density_data = read_h5_file(density_fnm)

mesh, mesh_center, x_mesh, y_mesh, z_mesh, x_mesh_centers, y_mesh_centers, z_mesh_centers, borders_mesh = setup_3d_mesh(mesh_data)
density_array = setup_uw_data_array_PROJ_3D(density_data)

g = calculate_3DModel_gravity(density_array = density_array, mesh = mesh, height = 30, num_cpu = 4)




