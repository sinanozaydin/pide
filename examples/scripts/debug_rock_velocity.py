import pide
import numpy as np
import matplotlib.pyplot as plt

#Setting the temperature and pressure arrays
temp = np.arange(600,1500,10) #setting up temperature array
pressure = np.arange(1,6,0.1)
#Creating a meshgrid
T,P = np.meshgrid(temp,pressure)
#flattening the arrays to load into pide
T_array = T.ravel()
P_array = P.ravel()

p_obj = pide.pide() #creating the initial object

#Setting the temperature and pressure arrays.
p_obj.set_temperature(T_array) #setting temperature array in K
p_obj.set_pressure(P_array) 

p_obj.set_solid_phase_method('rock')
p_obj.set_composition_solid_rock(granite = 1)
v_bulk_granite,v_p_granite,v_s_granite = p_obj.calculate_seismic_velocities()

fig = plt.figure(figsize = (20,10))
ax = plt.subplot(121)
cax = ax.tricontourf(T_array,P_array, v_p_granite, cmap = 'coolwarm_r',levels = 100)
cax.set_clim(np.amin(v_p_granite),np.amax(v_p_granite))
ax.set_xlabel(r'Temperature [$K^{\circ}$]')
ax.set_ylabel(r'Pressure [GPa]')
cax_cb = fig.colorbar(cax, orientation="horizontal", pad=0.1,
			 ax = ax, label = r'$V_P$ [km/s]')

ax = plt.subplot(122)
cax = ax.tricontourf(T_array,P_array, v_s_granite, cmap = 'coolwarm_r',levels = 100)
cax.set_clim(np.amin(v_s_granite),np.amax(v_s_granite))
ax.set_xlabel(r'Temperature [$K^{\circ}$]')
ax.set_ylabel(r'Pressure [GPa]')
cax_cb = fig.colorbar(cax, orientation="horizontal", pad=0.1,
			ax = ax, label = r'$V_S$ [km/s]')
fig.suptitle('Lherzolite with 8% Phlogopite - Composition')

plt.show()