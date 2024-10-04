import pide
import numpy as np
import unittest

class test_mixing_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		
		#tolerance of acceptance value
		self.atol = 1e-10
	
		temp = np.array([1000])
		pres = np.array(2.0)
		
		self.p_obj = pide.pide()
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.cond_check_solid = [0.0003442257078239972, 0.0006463744196063008, 0.0021093477280029055, 0.0026937719092410975, 0.0004587858438986189, 0.0007677865611706103]
		
		self.cond_check_melt_solid = [0.0006463744196074318, 0.000726540376246293, 0.0008333378894745669, 0.0008340140844370986, 0.0008333378894745661, 0.0006864496718682038]
		
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