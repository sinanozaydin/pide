from pide.rheology.olivine_rheology import olivine_rheology
import numpy as np
import unittest


class test_viscosity(unittest.TestCase):

	@classmethod
	def setUpClass(self):
	
		#Test written for viscosity and strain functions at October 31 2024. With newly added functions, this test file has to change
		#since there are new studies added to the arrays.
		
		#tolerance of acceptance value
		self.atol = 1e-5
		
		self.temp = np.arange(800,1200,100)
		self.pres = 3 * np.ones(len(self.temp))
		self.ol_water = 100 * np.ones(len(self.temp))
		
		self.results_visc = [3.199506412566774e+28, 7.789609000791768e+25, 5.414339311434777e+23, 6.88292256505474e+21]
		
	def test_viscosity_functions(self):
		
		rheol_obj = olivine_rheology(T = self.temp,P = self.pres,water = self.ol_water, xFe = 0.1)
		stress_input = rheol_obj.Stress_from_grainSize_vanderWAL1977(grain_size = 1) 
		
		diff_strain = rheol_obj.Hirth_Kohlstedt_2003_diff_fugacity(gr_sz = 1,stress = stress_input, melt = 0.0, fugacity_model= 'Zhao2004', calibration_model="Withers2012")
		disl_strain = rheol_obj.Hirth_Kohlstedt_2003_dislocation_fugacity(stress = stress_input, melt = 0.0, fugacity_model= 'Zhao2004', calibration_model="Withers2012")
		gbs_strain = rheol_obj.Ohuchi_et_al_2014_GBS(gr_sz=1, stress = stress_input, fugacity_model= 'Zhao2004', calibration_model="Withers2012")
		
		eff_visc = rheol_obj.calculate_effective_viscosity(stress = stress_input, strain_diff=diff_strain,strain_disl=disl_strain,strain_GBS=gbs_strain)
		
		assert np.allclose(np.array(eff_visc), np.array(self.results_visc),atol = self.atol)
		
		
if __name__ == "__main__":
	
	unittest.main()