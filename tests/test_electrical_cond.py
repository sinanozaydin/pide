import pide
import unittest
import numpy as np

class test_conductivity_functions(unittest.TestCase):

	@classmethod
	def setUpClass(self):
	
		#Test written for electrical conductivity models at October 2 2024. With newly added functions, this test file has to change
		#since there are new studies added to the arrays.
		
		#tolerance of acceptance value
		self.atol = 1e-5
		
		#Setting up temperature and pressure
		temp = np.array([1000])
		pres = np.array(2.0)
		

		self.p_obj = pide.pide()
		
		#setting up temperature pressure and equivalent parameters used in calculations
		self.p_obj.set_temperature(temp)
		self.p_obj.set_pressure(pres)
		self.p_obj.set_param1_mineral(mica = 0.2, plag = 0.1)
		self.p_obj.set_param1_rock(sandstone = 10)
		self.p_obj.set_fluid_properties(salinity = 0.1)
		self.p_obj.set_melt_properties(co2 = 100, water = 1000, na2o = 0.1, k2o = 0.1)
		
		#calculated values for all minerals and rocks
		self.mineral_cond_list_check = [[0.0002147666215751601, 0.0019628863936549724, 0.0019628863936549724, 1.0190583930014192e-05,
						8.208132334583077e-05, 0.00020788083347927274, 0.0001440326145468413, 6.855909754207005e-07, 0.0010714138960102133,
						7.377942552234648e-06, 1.7927554236795322e-06, 3.297803999685759e-06, 1.7484794023483996e-06,
						0.00018056600809558101, 1.6910366224498637e-06, 1.325877798709181e-05, 5.684421480001739e-06, 3.826141179766439e-06,
						5.065540593063711e-06, 4.670308887813639e-05, 1.0328167482905911e-07], [0.000781880478677682,
						0.0041784297728858005, 7.933807799716656e-05, 6.5886975752919734e-06, 2.079446696555967e-06, 8.528539400246339e-07,
						2.465168638439833e-06], [0.0004962607047489137, 0.002462870542595579, 0.007108929846228112, 2.77195158939179e-06,
						0.0009367836943396204, 0.00020649961515201444, 0.0003657934479451687, 5.342668461994264e-07, 3.036074531806225e-07,
						0.0006794702021050355, 0.002166081366067568, 0.0020843655959144006], [0.017628516537043188, 3.419744010057346e-05,
						1.2886396201172216e-10, 1.9792124574723175e-11, 3.8305090578826584e-06, 0.00036438748674276397, 0.20570613002486426,
						2.623717614437895, 6.958467616527965, 0.0027019883214688574, 0.0005859226310618902, 0.0034794253202184602,
						5.2657573700746e-05, 7.373852865647446e-05, 0.00012507246648229607, 0.0016101702596428337,
						2.2867640880756637e-07, 0.0011355912598567696, 7.758182189447044e-05], [0.008479833615942118,
						0.1640862682924139], [0.0070481722066036735,
						0.02080240802198666, 13667.57423820947, 0.03391151823805323, 0.14736974502790565,
						0.04356001857421698], [0.0003188940114855732, 0.013619445169725835, 0.00012394697219022993,
						2.506389402139377e-05, 3.689650732899158e-05, 3.2249459603471138e-06, 6.596016408306563,
						0.025493674803551688, 0.004856571995839248], [0.0004128716087684925, 2.197819422218023e-05,
						0.0018054801889687924, 4.18379313046424e-06, 0.30574923240698204,
						0.5276860781809379, 0.02360645107435132], [0.0910880250289498,
						0.06145860932937158, 0.05067569206031278, 0.0989669940613521], 
						[3332.728995761379, 344.3499307633384, 153.81546403030342],
						[100000.0], [15.668221458353507, 232.0239164777426], 
						[0.00038867129908321056, 6.617775585469027e-05, 7.855161191415685e-05,
						5.96958471942372e-05, 0.000725938263780609, 0.0001791219710105094, 0.00029143475739094724, 1.7437241829326952e-05,
						1.4652025419612293e-05, 1.3456087012655397e-05, 7.869745619559822e-05, 0.0003563383662745001],
						[9.136938593750724e-05, 0.00031829235011104895, 35.07677858258606, 0.016354842407704934, 0.007612290451624446,
						46.83657179784756], [0.00019107660645723495, 0.0009147662433918977, 0.0025793801343671945, 0.19739108420395143, 2.170179463231433,
						0.09033385965305728, 0.023804816701909138, 0.031559219353188686, 0.026456434639873783, 0.15243823410196836,
						0.6470873295769766, 8.420299371544626, 8.63964005700373, 122.06699678239542, 0.0003513802736255081, 0.0003931961527606333,
						0.0005317099241385271, 0.0025914554862907784, 0.04564406736998614, 0.14886276382806088,
						0.7179351351044376, 2.4126573631693167], [244.21164927494843, 0.0027744528159345428, 0.003604655509112567,
						0.0003936708988150554]]
						
		self.rock_cond_list_check = [[6.850646712231316e-06, 4.289840479223029e-05, 0.01770295328980004, 0.0002448778106938003, 0.00021047417616575027, 0.00012485445378909127,
		4.200470588638835e-06, 0.00033903506939415547, 0.0003958561083780382, 8.315216946910493e-05, 1.2486261946646642e-05, 0.00024114209299575372],
		[0.024589951680813342, 3.63368464619912, 0.0036273508031219375, 0.0069424914042069885, 0.04346239059278374, 0.006664214605720846], 
		[0.031188895840939354], [0.002073930843385905, 0.00745376397205312, 0.01218077065864558, 0.0775612069832346, 0.46737810548631803],
		[0.14736974502790565, 3.9464210219680207, 58.54288802397429, 13667.540326691233, 0.008290972861156798], [0.0012790228002886473, 
		0.06626895617758016, 0.012726584640203181, 0.004621835463107331, 0.0013870339125108635, 0.003985402950653599, 0.15773933062288584], 
		[0.0068899523886569655, 0.04174000094758388, 0.11470690410306332], [0.021357794074482208, 0.0], [0.00015341929653709275, 
		0.0003606185841661195, 0.0013154565463910474, 100000.0]]
		
		self.fluid_cond_list_check = [25.081286700000014, 20.65512240933996, 1.847006803747041, 2.347808955788135, 415.3305805225229, 0.05541742774428419]
		
		self.melt_cond_list_check = [0.0018385859285729245, 0.0399511738130767, 39396802.36751015, 0.016063600072673018, 0.0022231159686968128, 0.0023407173410140434,
		0.00014588209386490485, 1.843408259149654e-14, 0.0001037151459511539, 0.0015140360467432313, 0.008814750009255202, 0.00021485904252158006, 595703291.8333696]


	def test_mineral_conds(self):
	
		#testing mineral conductivities
		
		mineral_list = self.p_obj.list_available_minerals()
		all_cond_list = []
		for i in range(0,len(mineral_list)):
		
			list_models = self.p_obj.list_mineral_econd_models(mineral_list[i])
			mineral_cond_list = []
			for j in range(0,len(list_models)):
				
				eval('self.p_obj.set_mineral_conductivity_choice(' + mineral_list[i] + '=' + str(j) + ')')
				eval(f"self.p_obj.set_mineral_water({mineral_list[i]}= 50)")
				cond = eval(f"self.p_obj.calculate_mineral_conductivity(min_idx = '{mineral_list[i]}')")
				mineral_cond_list.append(cond[0])
			all_cond_list.append(mineral_cond_list)
		
		for mineral in range(0,len(all_cond_list)):
			
			assert np.allclose(np.array(all_cond_list[mineral]), np.array(self.mineral_cond_list_check[mineral]), atol=1e-5)
			
	def test_rock_conds(self):
		
		#testing rock conductivities
		
		rock_list = self.p_obj.list_available_rocks()
		all_cond_rock_list = []
		for i in range(0,len(rock_list)):
		
			list_models = self.p_obj.list_rock_econd_models(rock_list[i])
			rock_cond_list = []
			for j in range(0,len(list_models)):
				
				eval('self.p_obj.set_rock_conductivity_choice(' + rock_list[i] + '=' + str(j) + ')')
				eval(f"self.p_obj.set_rock_water({rock_list[i]}= 50)")
				cond = eval(f"self.p_obj.calculate_rock_conductivity(rock_idx = '{rock_list[i]}')")
				rock_cond_list.append(cond[0])
			all_cond_rock_list.append(rock_cond_list)
		
		for rock in range(0,len(all_cond_rock_list)):
		
			assert np.allclose(np.array(all_cond_rock_list[rock]), np.array(self.rock_cond_list_check[rock]), atol=self.atol)
			
	def test_fluid_conds(self):
	
		#testing fluid conductivities
	
		list_models = self.p_obj.list_fluid_econd_models()
		
		
		fluid_cond_list = []
		for j in range(0,len(list_models)):
			
			eval(f"self.p_obj.set_melt_fluid_conductivity_choice(fluid={str(j)})")
			
			cond = self.p_obj.calculate_fluids_conductivity()
			fluid_cond_list.append(cond[0])
			
		assert np.allclose(np.array(fluid_cond_list), np.array(self.fluid_cond_list_check), atol = self.atol)
		
	def test_melt_conds(self):
		
		#testing melt conductivities
	
		list_models = self.p_obj.list_melt_econd_models()
		
		melt_cond_list = []
		for j in range(0,len(list_models)):
			
			eval(f"self.p_obj.set_melt_fluid_conductivity_choice(melt={str(j)})")
			
			cond = self.p_obj.calculate_melt_conductivity()
			melt_cond_list.append(cond[0])
			
		assert np.allclose(np.array(melt_cond_list), np.array(self.melt_cond_list_check), atol = self.atol)
				
if __name__ == "__main__":
	
	unittest.main()