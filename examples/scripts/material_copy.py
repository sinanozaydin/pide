from pide.material import Material
from pide.model import Model
from pide.utils.utils import sort_through_external_list
#For building a geotherm.
from pide.geodyn.geotherm import calculate_hasterok2011_geotherm
import numpy as np
import sys, ipdb,time
start = time.time()

moho = 38 #km

T, depth, p, idx_LAB = calculate_hasterok2011_geotherm(SHF = 40, T_0 =25.0,max_depth = 250,moho = 38.0)
#Here, other than te compositional and geometrical parameters, the user should also indicate the top and the bottom of the layers in km.
Layer_1_Upper_Crust = Material(name = 'Layer_1_Upper_Crust', calculation_type = 'rock', composition = {'granite':1.0},
							   el_cond_selections = {'granite': 10},solid_phase_mixing_idx = 1,top = 0, bottom = 15)

#Instead of using 15, the user here can also use another objects bottom variable.
Layer_2_Lower_Crust = Material(name = 'Layer_2_Lower_Crust', calculation_type = 'rock', composition = {'granulite':1.0},
							   el_cond_selections = {'granulite': 2},solid_phase_mixing_idx = 1,top = Layer_1_Upper_Crust.bottom, bottom = 38)

Layer_3_Upper_Mantle = Material(name = 'Layer_3_Upper_Mantle', calculation_type = 'mineral',
								composition = {'ol':0.65,'opx':0.25,'garnet':0.05,'cpx':0.05},
								solid_phase_mixing_idx = 1,
								el_cond_selections = {'ol':4, 'opx':0, 'garnet':0,'cpx':0},
								water_distr = True, water = {'bulk':0}, xfe = {'ol':0.1,'opx':0.1,'garnet':0.1, 'cpx':0.1},
								top = Layer_2_Lower_Crust.bottom, bottom = 75)

#A layer's attributes can also be copied by the copy_attributes method in a material object:
Layer_4_Upper_Mantle = Material()
Layer_3_Upper_Mantle.copy_attributes(Layer_4_Upper_Mantle)
#Now change just only the water content and layer position from Layer_3
Layer_4_Upper_Mantle.water = {'bulk':200}
Layer_4_Upper_Mantle.top = Layer_3_Upper_Mantle.bottom
Layer_4_Upper_Mantle.bottom = 150

#
Layer_5_Upper_Mantle = Material()
Layer_4_Upper_Mantle.copy_attributes(Layer_5_Upper_Mantle)
Layer_5_Upper_Mantle.composition = {'ol':0.60,'opx':0.20,'garnet':0.05,'mica':0.1,'cpx':0.05}
Layer_5_Upper_Mantle.top = Layer_4_Upper_Mantle.bottom
Layer_5_Upper_Mantle.bottom = 250

#Now creating a model object to put these layers in alongside with our thermodynamic conditions (T,P,Depth)
Layers = Model(material_list = [Layer_1_Upper_Crust,Layer_2_Lower_Crust,Layer_3_Upper_Mantle,Layer_4_Upper_Mantle,Layer_5_Upper_Mantle],
			   T = T, P = p, depth = depth)

cond = Layers.calculate_geothermal_block(type = 'conductivity')