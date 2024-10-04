import pide
import numpy as np
import unittest


class test_water_solubility_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		
		#tolerance of acceptance value
		self.atol = 1e-3
	
		temp = np.array([1000])
		pres = np.array(2.0)
		
		self.p_obj = pide.pide()
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.mineral_names_sol = ['ol','opx','cpx','garnet','rwd_wds']
		self.p_obj.set_parameter('ti_ol',0.01)
		
		self.sol_check = [[45.38739975258276, 50.29536192118538, 51.08186026610852, 12.200000000000003, 3.726218152020939,
		168.7764517659536, 68.82584537321634, 142.7, 61.0, 858.6846951266635], [1671.859101411614, 47671.4231294835, 1665.4231253174394,
		1671.859101411614], [1042.4432198837694, 2029.6369491136993, 14906.486], [686.9477561013309, 684.303287238804, 14906.486],
		[28777.0, 28777.0, 28777.0]]
		
		self.bulk_sol_check = [370.74528244]
		
	def test_1_mineral_water_solubility(self):
	
		all_sols = []
		for mineral_name in self.mineral_names_sol:
			
			list_sols = eval(f"self.p_obj.list_mantle_water_solubilities(mineral_name='{mineral_name}')")
			sols = []
			for j in range(0,len(list_sols)):
				eval(f"self.p_obj.set_mantle_water_solubility({mineral_name} = {j})")
				max_water = eval(f"self.p_obj.calculate_mineral_water_solubility('{mineral_name}')")
				sols.append(max_water[0])
			all_sols.append(sols)
		
		for mineral in range(0,len(self.sol_check)):
			assert np.allclose(np.array(all_sols[mineral]),self.sol_check[mineral], atol = self.atol)
			
	
	def test_2_bulk_water_solubility(self):
	
		self.p_obj.set_composition_solid_mineral(ol = 0.6, opx = 0.2, garnet = 0.1, cpx = 0.1)
		max_water_lithosphere = self.p_obj.calculate_bulk_mantle_water_solubility()
				
		assert np.allclose(np.array(max_water_lithosphere),np.array(self.bulk_sol_check),atol = self.atol)

if __name__ == "__main__":
	
	unittest.main()