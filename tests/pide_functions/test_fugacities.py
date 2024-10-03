import pide
import numpy as np
import unittest


class test_fugacity_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
	
		#tolerance of acceptance value
		self.atol = 1e-10
	
		temp = np.array([1000])
		pres = np.array(2.0)
		
		self.p_obj = pide.pide()
		
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		
		self.water_fugacity_calculated = [13.05021266]
		
		self.o2_fugacity_calculated = [6.895884836192688e-15, 2.0556302500027243e-20, 9.357288076216961e-22, 2.238484028213431e-15, 2.23844279419359e-21]
		
	def test_1_water_fugacity(self):
	
		self.p_obj.calculate_water_fugacity()
		
		assert np.allclose(np.array(self.p_obj.water_fugacity), np.array(self.water_fugacity_calculated),atol = self.atol)
		
	def test_2_o2_fugacity(self):
	
		mode_list = [0,1,2,3,4]

		fugs = []
		for i in mode_list:
		
			o2_fug = self.p_obj.calculate_o2_fugacity(mode = i)
			fugs.append(o2_fug[0])
			
		assert np.allclose(np.array(fugs),np.array(self.o2_fugacity_calculated),atol = self.atol)
		
if __name__ == "__main__":
	
	unittest.main()