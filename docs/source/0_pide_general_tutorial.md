![Example Image](../../docs/figures/pide_logo.png)
## <span style="color:green"> Notebook - Starting with pide </span>
<span style="color:purple">Sinan Ozaydin, School of Geosciences, The University of Sydney, NSW 2006, Australia <br/> </span>
<span style="color:purple">sinan.ozaydin@protonmail.com</span>

pide is a python3 library to calculate geophysical observableslike electrical conductivity and seismic velocity for the given compositional and thermodynamic environment. To derive this information it uses library of experimental models on electrical conductivity and thermoelastic constants. pide is constituted by three main classes:

<span style="color:red;">-pide:</span> This is the class where calculations related to electrical conductivity and seismic velocities are performed. pide object is central to the library and all material and model classes utilizes pide object to calculate electrical conductivities and seismic velocities. 

<span style="color:red;">-material:</span> A material holder class where the user can define a compositional environment and choices relevant to how the geophysical observables are going to be calculated. This can be done in the pide class by directly setting these parameters on a pide class, however, material class should create a tidy and interoperable environment for more complex calculations that can be streamlined through other processes.

<span style="color:red;">-model:</span> The model class represents a collection of materials that are geospatially indexed and thermodynamically defined. The users can create model objects from informations derived from thermomechanical numerical (like underworld or aspect outputs), or geological models.

**Other tools included in the toolkit of pide are:**

<span style="color:green;">-geodyn:</span> A collection of scripts and methods to deal with input/output/processing of thermomechanical models.

<span style="color:green;">-gravity and mag:</span> A collection of methods to calculate gravity and magnetic anomalies from given 3D environment set by a model object.

<span style="color:green;">-mt:</span> A collection of tools to convert the outputs of pide into MT input scripts for forward and inverse modelling of these synthetic models. The codes currently implemented are ModEM and MARE2DEM.

<span style="color:green;">-rheology:</span> A collection of methods to perform rheological calculations from the outputs of pide (e.g., mantle viscosity from MT driven water contents).

<span style="color:green;">-inversion:</span> A collection of methods to invert for compositional parameterisations used in pide (e.g., mineral/melt contents, interconnectivity of phases, water contents) from real data.

To understand the functions used in pide, one can utilise the following workflow:


```python
import pide

#creating the pide object
p_object = pide.pide()

#Now, we can list all the methods used in pide class by:
p_object.list_methods()


```

    - calculate_arrhenian_single
    - calculate_bulk_mantle_water_solubility
    - calculate_conductivity
    - calculate_density_solid
    - calculate_fluids_conductivity
    - calculate_melt_conductivity
    - calculate_mineral_conductivity
    - calculate_mineral_water_solubility
    - calculate_o2_fugacity
    - calculate_rock_conductivity
    - calculate_seismic_velocities
    - calculate_transition_zone_water_solubility
    - calculate_water_fugacity
    - get_method_manual
    - get_mineral_index
    - get_rock_index
    - list_available_minerals
    - list_available_rocks
    - list_fluid_econd_models
    - list_mantle_water_partitions_melt
    - list_mantle_water_partitions_solid
    - list_mantle_water_solubilities
    - list_melt_econd_models
    - list_methods
    - list_mineral_econd_models
    - list_phs_melt_fluid_mix_methods
    - list_phs_mix_methods
    - list_rock_econd_models
    - list_transition_zone_water_partitions_solid
    - mantle_water_distribute
    - reset
    - revalue_arrays
    - set_alopx
    - set_bulk_water
    - set_composition_solid_mineral
    - set_composition_solid_rock
    - set_depth
    - set_fluid_properties
    - set_grain_boundary_H_Diffusion
    - set_grain_boundary_water_partitioning
    - set_grain_size
    - set_mantle_transition_zone_water_partitions
    - set_mantle_water_partitions
    - set_mantle_water_solubility
    - set_melt_fluid_conductivity_choice
    - set_melt_fluid_frac
    - set_melt_fluid_interconnectivity
    - set_melt_or_fluid_mode
    - set_melt_properties
    - set_mineral_conductivity_choice
    - set_mineral_water
    - set_o2_buffer
    - set_param1_mineral
    - set_param1_rock
    - set_parameter
    - set_phase_interconnectivities
    - set_pressure
    - set_rock_conductivity_choice
    - set_rock_water
    - set_seismic_velocity_properties
    - set_solid_melt_fluid_mix_method
    - set_solid_phase_method
    - set_solid_phs_mix_method
    - set_temperature
    - set_watercalib
    - set_xfe_mineral
    - transition_zone_water_distribute
    - write_data



```python
#In order to get manual of a certain method one can do:
p_object.get_method_manual('calculate_mineral_conductivity')
```

    [91mThe Manual for the [93mcalculate_mineral_conductivity:[0m
    A method to calculate mineral conductivity with the environment set up.
    		
    		Input:
    		int: min_idx - Name of the mineral or index of the mineral chosen.
    		
    		str: method - 'array' or 'index'|| Default - 'array
    		int: sol_idx: index parameter if method|index is chosen.
    		
    		Output:
    		float/array: conductivity || in (S/m)
    		
    		Example:
    		calculate_mineral_conductivity(min_idx = 'ol')


A lot of environment variables to set in pide can be access through the **list methods**. An example can be:


```python
#To list all available electrical electrical conductivity models and get their indexes.
p_object.list_mineral_econd_models('ol')
```

    [91mElectrical conductivity models for the given mineral: ol[0m
    0.   Dai2014_DryandWetOlivine_fo2  -----  polaron+proton
    1.   Dai2020_WetOlivine_200ppmTi_fo2  -----  polaron+proton
    2.   Dai2020_WetOlivine_683ppmTi_fo2  -----  polaron+proton
    3.   Fei2020_WetOlivineIonic_Isotropic  -----  ionicWet
    4.   Gardes2014_DryandWetOlivine  -----  ionic+polaron+proton
    5.   Jones2012_WetOlivine  -----  proton
    6.   Liu2021_DryOlivine_NNOBuffer  -----  polaron
    7.   Poe2010_DryandWetOlivine  -----  polaron+proton
    8.   Wang2006_DryandWetOlivine  -----  polaron+proton
    9.   Yoshino2009_DryandWetOlivine  -----  ionic+polaron+proton
    10.   Constable2006_dryOlivine_fo2  -----  polaron
    11.   Dai2014_DryOlivine_xFe  -----  polaron
    12.   Fullea2011_DryOlivine_xFe  -----  polaron
    13.   Pommier2018_ShearedDryOlivine  -----  polaron
    14.   Xu1998_DryOlivine  -----  polaron
    15.   Yoshino2012_DryOlivine_xFe  -----  polaron
    16.   Novella2017_HDiffusion  -----  protondiffusion
    17.   Sun2019_HDiffusion  -----  protondiffusion
    18.   DuFrane2012_HDiffusion  -----  protondiffusion
    19.   Kohlstedt1998_HDiffusion  -----  protondiffusion
    20.   Demouchy2006_HDiffusion  -----  protondiffusion





    ['Dai2014_DryandWetOlivine_fo2',
     'Dai2020_WetOlivine_200ppmTi_fo2',
     'Dai2020_WetOlivine_683ppmTi_fo2',
     'Fei2020_WetOlivineIonic_Isotropic',
     'Gardes2014_DryandWetOlivine',
     'Jones2012_WetOlivine',
     'Liu2021_DryOlivine_NNOBuffer',
     'Poe2010_DryandWetOlivine',
     'Wang2006_DryandWetOlivine',
     'Yoshino2009_DryandWetOlivine',
     'Constable2006_dryOlivine_fo2',
     'Dai2014_DryOlivine_xFe',
     'Fullea2011_DryOlivine_xFe',
     'Pommier2018_ShearedDryOlivine',
     'Xu1998_DryOlivine',
     'Yoshino2012_DryOlivine_xFe',
     'Novella2017_HDiffusion',
     'Sun2019_HDiffusion',
     'DuFrane2012_HDiffusion',
     'Kohlstedt1998_HDiffusion',
     'Demouchy2006_HDiffusion']



Then we can use this information to set our olivine conduction model as Gardes2014_DryandWetOlivine:


```python
p_object.set_mineral_conductivity_choice(ol = 4)
```

The pide object can be reset to default values can be done as follows:


```python
p_object.reset()
```
