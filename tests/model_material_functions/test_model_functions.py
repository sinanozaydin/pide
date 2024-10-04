import pide
import unittest
import numpy as np

class test_model_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
	
		#Test written for electrical conductivity models at October 2 2024. With newly added functions, this test file has to change
		#since there are new studies added to the arrays.
		
		#tolerance of acceptance value
		self.atol = 1e-5
		
		self.check_deformcond = 0.007827479821995546
		
		self.check_model_cond = [8.219948950770384e-05, 0.00012332868124250433, 0.0001804442443035664, 0.000258126207567983,
		0.00036184266230103684, 0.0004980549434863003, 0.0006743349570385356, 0.0008994988983728286, 0.0011837628186641453, 
		0.001538926883429566, 0.0019785966390192923, 0.0025184509933538182, 0.003176567704138736, 0.00397381768407388,
		0.0049343390964659934, 0.006086100767774966, 0.007461561699491269, 0.009098429362685826, 0.011040514136305383,
		0.013338671080386638]
		
	def test_1_deformcond(self):
	
		from pide.geodyn.deform_cond import plastic_strain_2_conductivity
	
		deformcond = plastic_strain_2_conductivity(strain = 0.5, low_cond = 1e-6,high_cond = 10, low_strain = 0.01, high_strain=1, function_method = "exponential",
		conductivity_decay_factor = 0.7, conductivity_decay_factor_2 = 0.2, strain_decay_factor = 0.3)
		
		assert np.isclose(deformcond,self.check_deformcond, atol = self.atol)
		
	def test_2_material_and_model(self):
	
		from pide.model import Model
		from pide.material import Material
		
		material_array = np.ones(20)
		temp_array = np.linspace(900,1500,20)
		pressure_array = np.ones(20)
		melt_array = np.zeros(20)
		
		mantleLithosphere = Material(name = 'mantleLithosphere', material_index = 1,
		calculation_type = 'mineral', composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
		interconnectivities = {'ol':1,'opx':2,'garnet':5, 'cpx':5}, 
		el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0}, water_distr = True, water = {'bulk':100},
		xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1}, solid_phase_mixing_idx = 1)
		
		material_object_list = [mantleLithosphere]
		material_skip_list = [None]
		
		model_object = Model(material_list = material_object_list, material_array = material_array, T = temp_array, P = pressure_array, model_type = 'underworld_3d', melt = melt_array,
		material_node_skip_rate_list = material_skip_list)
		cond_model = model_object.calculate_model(type = "conductivity",num_cpu=2)
		
		assert np.allclose(cond_model,np.array(self.check_model_cond),atol = self.atol)
		
if __name__ == "__main__":
	
	unittest.main()