import pide
import numpy as np
import unittest

class test_mixing_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		
		#tolerance of acceptance value
		self.atol = 1e-10
	
		temp = np.array([1400])
		pres = np.array(2.0)
		
		self.p_obj = pide.pide()
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.cond_check_solid = [0.012005345480507728, 0.019389690076153224, 0.03185792193441639, 0.0378625380743522, 0.015849050568884276, 0.020373418961111167]
		
		self.cond_check_melt_solid = [0.01938969007624769, 0.028016271558536296, 0.03752488975921886, 0.037592677815007, 0.03752488975921913, 0.020698902992675225]
		
	def test_1_mixing_function(self):
	
		self.p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.2, cpx = 0.1,garnet = 0.1)
		self.p_obj.set_phase_interconnectivities(ol = 1.0,opx = 2.5,cpx = 3,garnet = 3)
		self.p_obj.set_bulk_water(100)
		self.p_obj.mantle_water_distribute()
		
		mixing_methods = self.p_obj.list_phs_mix_methods()

		cond_list = []
		for i in range(0,len(mixing_methods)):
			
			self.p_obj.set_solid_phs_mix_method(i)
			cond = self.p_obj.calculate_conductivity()
			cond_list.append(cond[0])
			
		assert np.allclose(np.array(cond_list),np.array(self.cond_check_solid),atol = self.atol)
	
	def test_2_mixing_function_melt(self):
	
		self.p_obj.set_composition_solid_mineral(ol = 0.6,opx = 0.2, cpx = 0.1,garnet = 0.1)
		self.p_obj.set_phase_interconnectivities(ol = 1.0,opx = 2.5,cpx = 3,garnet = 3)
		self.p_obj.set_solid_phs_mix_method(1)
		self.p_obj.set_melt_fluid_frac(0.02)
		self.p_obj.set_melt_properties(water = 10000,co2 = 100)
		melt_mixing_methods = self.p_obj.list_phs_melt_fluid_mix_methods()
		
		cond_melt_list = []
		for i in range(0,len(melt_mixing_methods)):
		
			self.p_obj.set_solid_melt_fluid_mix_method(i)
			cond = self.p_obj.calculate_conductivity()
			cond_melt_list.append(cond[0])
			
		assert np.allclose(np.array(cond_melt_list), np.array(self.cond_check_melt_solid),atol = self.atol)
	
if __name__ == "__main__":
	
	unittest.main()