import os,sys
import time
import numpy as np
import pide
from pide.inversion import conductivity_metropolis_hastings_two_param, conductivity_solver_single_param
from pide.imaging.plot_distribution import plot_posterior_distribution_two_params, plot_posterior_distribution_heatmap_two_params
from pide.utils.utils import save_h5_files
import matplotlib.pyplot as plt
import seaborn as sns

# Observed value (target output)
T = np.array([1200.0])
P = np.array([2])

p_obj = pide.pide()
p_obj.set_temperature(T)
p_obj.set_pressure(P)
p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.4)
p_obj.set_melt_fluid_frac(0.15)

p_obj.set_bulk_water(25000.0)
p_obj.mantle_water_distribute()
cond = p_obj.calculate_conductivity()
print(cond)
import ipdb
ipdb.set_trace()
