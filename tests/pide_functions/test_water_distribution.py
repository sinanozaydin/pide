import pide
import numpy as np
import unittest


class test_water_distribution_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
	
		#Test written for water distribution functions at October 2 2024. 
		
		#tolerance of acceptance value
		self.atol = 1e-5
		
		temp = np.array([1000])
		pres = np.array(2.0)
		
		self.p_obj = pide.pide()
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.lithosphere_water_results = [148.9998866855862, 290.10277937683634, 352.18477416347935, 119.19990934846896]
		
		self.lithosphere_water_melt_results = [251.76250201369436, 490.1815914206629, 595.0804519846848, 201.4100016109555, 83441.11143560433]
		
		self.tz_water_results = [1357.6349998302956, 703.2549299120931, 357.2616502053423, 90.4184909886977]
		
	def test_mantle_water_distribute(self):
	
		self.p_obj.set_composition_solid_mineral(ol = 0.6, opx = 0.3, cpx = 0.05,garnet = 0.05)
		self.p_obj.set_bulk_water(200)
		
		self.p_obj.mantle_water_distribute()
		
		water_list = [self.p_obj.ol_water[0],self.p_obj.opx_water[0],self.p_obj.cpx_water[0],self.p_obj.garnet_water[0]]
		
		assert np.allclose(np.array(water_list),np.array(self.lithosphere_water_results),atol = self.atol)
		
	def test_mantle_water_distribute_with_melt(self):
	
		self.p_obj.set_melt_fluid_frac(0.02)
		self.p_obj.set_bulk_water(2000.0)
		self.p_obj.mantle_water_distribute()
		
		water_list = [self.p_obj.ol_water[0],self.p_obj.opx_water[0],self.p_obj.cpx_water[0],self.p_obj.garnet_water[0], self.p_obj.melt_water[0]]
		
		assert np.allclose(np.array(water_list),np.array(self.lithosphere_water_melt_results),atol = self.atol)
		
	def test_mantle_transition_zone_water_distribute(self):
	
		self.p_obj.reset()
		
		temp = np.array([2000])
		pres = np.array(13.0)
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.p_obj.set_composition_solid_mineral(rwd_wds = 0.6, cpx = 0.2, garnet = 0.1, perov = 0.1)
		self.p_obj.set_bulk_water(1000.0)
		self.p_obj.transition_zone_water_distribute()
		
		water_list = [self.p_obj.rwd_wds_water[0],self.p_obj.cpx_water[0],self.p_obj.garnet_water[0],self.p_obj.perov_water[0]]
		
		assert np.allclose(np.array(water_list),np.array(self.tz_water_results),atol = self.atol)
				
if __name__ == "__main__":
	
	unittest.main()