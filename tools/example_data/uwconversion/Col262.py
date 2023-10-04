#!/usr/bin/env python3
# coding: utf-8

##from underworld.UWGeodynamics.surfaceProcesses import SedimentationThreshold
##from UWGeodynamics import visualisation as vis


import underworld as uw
from underworld import UWGeodynamics as GEO
from underworld import function as fn
import numpy as np
##from UWGeodynamics import visualisation as vis

#import scipy
#import os.path
#from mpi4py import MPI
#import argparse
#import math

GEO.__version__

u = GEO.UnitRegistry

plotting_bool = False

resolution = (768, 384)
#resolution_setting_up = (768,384)
half_rate = 1.25 * u.centimeter / u.year
model_length = 1536e3 * u.meter
surfaceTemp = 293.15 * u.kelvin
baseModelTemp = 1673.15 * u.kelvin
bodyforce = 3370 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT

Model = GEO.Model(elementRes = resolution, 
				  minCoord = (0. * u.kilometer,  -720. * u.kilometer), 
				  maxCoord = (1536. * u.kilometer, 48. * u.kilometer),
				  periodic = [False, False],
				  gravity = (0.0, -9.81 * u.meter / u.second**2))


Model.outputDir = "Col262"

Model.diffusivity = 1e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)

Model.minViscosity = 2e18 * u.pascal * u.second
Model.maxViscosity = 3e22 * u.pascal * u.second

#Calling registry objects.
rh = GEO.ViscousCreepRegistry()
pl = GEO.PlasticityRegistry()
solidii = GEO.SolidusRegistry()
liquidii = GEO.LiquidusRegistry()

#Defining Air Layer and its properties
Air_Object = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom = (24.0 * u.kilometer)))
Air_Object.density = 1. * u.kilogram / u.metre**3
Air_Object.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Air_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Air_Object.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
Air_Object.viscosity = 1e18 * u.pascal * u.second
Air_Object.compressibility = 1e4  # 1.e4 seems to be a good value, nb: Compressibility should be zero when using Lecode isostasy


#Defining Air Layer and its properties
Sticky_Air_Object = Model.add_material(name="Sticky Air", shape=GEO.shapes.Layer(top=Air_Object.bottom, bottom = (-5.0 * u.kilometer)))
Sticky_Air_Object.density = 1. * u.kilogram / u.metre**3
Sticky_Air_Object.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Sticky_Air_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Sticky_Air_Object.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
Sticky_Air_Object.viscosity = 3e18 * u.pascal * u.second
Sticky_Air_Object.compressibility = 1e4  # 1.e4 seems to be a good value, nb: Compressibility should be zero when using Lecode isostasy

# +
# #Eclogitized oceanic material phase change object
Eclogite_Object = Model.add_material(name = "Eclogitized Oceanic Material")
Eclogite_Object.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Eclogite_Object.density = 3500 * u.kilogram / u.metre**3
Eclogite_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Eclogite_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Eclogite_Object.viscosity = 5e19 * u.pascal * u.second #0.1 * rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
# -

#Lithospheric Mantle Layer
Lithospheric_Mantle_Object = Model.add_material(name="Lithospheric Mantle",
shape=GEO.shapes.Layer(top=Sticky_Air_Object.bottom, bottom = (-155.0 * u.kilometer)))
#Lithospheric_Mantle_Object.density = 3395 * u.kilogram / u.metre**3
Lithospheric_Mantle_Object.density = GEO.LinearDensity(reference_density=3395. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Lithospheric_Mantle_Object.radiogenicHeatProd = 0.02e-6 * u.watt / u.meter**3
Lithospheric_Mantle_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Lithospheric_Mantle_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
#Lithospheric_Mantle_Object.viscosity = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
#Lithospheric_Mantle_Object.viscosity = rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003

Lithospheric_Mantle_Object.viscosity = GEO.CompositeViscosity([rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003, rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003])
Lithospheric_Mantle_Object.add_melt_modifier(solidii.Mantle_Solidus, liquidii.Mantle_Liquidus,
                                                                                        latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                                                        meltFraction=0.,
                                                                                        meltFractionLimit=0.02,
                                                                                        meltExpansion=0.13,
                                                                                        viscosityChangeX1 = 0.00,
                                                                                        viscosityChangeX2 = 0.02,
                                                                                        viscosityChange = 1e-2)
#Lithospheric_Mantle_Object.plasticity = pl.Rey_et_al_2014_LithosphericMantle
Lithospheric_Mantle_Object.plasticity = GEO.DruckerPrager(name="Lithospheric mantle",
                                           cohesion=10. * u.megapascal, 
                                           cohesionAfterSoftening=2 * u.megapascal,
                                           frictionCoefficient=0.577,
                                           frictionAfterSoftening=0.11,
                                           epsilon1=0.0, epsilon2=0.15)

Lithospheric_Mantle_Object.elasticity = GEO.Elasticity(shear_modulus=40e9 * u.pascal, observation_time=20000 * u.year)

Lithospheric_Mantle_Object.stressLimiter = 200. * u.megapascal #does this work same here? Not sure what to enter on stress= ?

#Asthenospheric Mantle
Asthenospheric_Mantle_Object = Model.add_material(name="Asthenospheric Mantle", shape=GEO.shapes.Polygon([
								(0. * u.kilometer,-140. * u.kilometer),
								(42. * u.kilometer,-139.5 * u.kilometer),
								(90. * u.kilometer, -139 * u.kilometer),
								(112. * u.kilometer,-138.5 * u.kilometer),
								(132. * u.kilometer,-137.7 * u.kilometer),
								(200. * u.kilometer,-136.5 * u.kilometer),
								(268. * u.kilometer,-134.1 * u.kilometer),
								(284. * u.kilometer,-129.9 * u.kilometer),
								(298. * u.kilometer,-123.5 * u.kilometer),
								(400. * u.kilometer,-42.0 * u.kilometer),
                                                                (544. * u.kilometer,-42.0 * u.kilometer),
								(535. * u.kilometer,-52.0 * u.kilometer),
                                                                (525. * u.kilometer,-62.0 * u.kilometer),
								(500. * u.kilometer,-88.0 * u.kilometer),
								(503.3 * u.kilometer,-95. * u.kilometer),
								(510.1 * u.kilometer,-107. * u.kilometer),
								(520. * u.kilometer,-117.1 * u.kilometer),
								(534. * u.kilometer,-124 * u.kilometer),
								(550. * u.kilometer,-127.4 * u.kilometer),
								(570. * u.kilometer,-127.3 * u.kilometer),
								(586. * u.kilometer,-123.6 * u.kilometer),
								(600. * u.kilometer,-118.1 * u.kilometer),
								(614. * u.kilometer,-111.4 * u.kilometer),
								(628. * u.kilometer,-107. * u.kilometer),
								(648. * u.kilometer,-105.2 * u.kilometer),
								(678. * u.kilometer,-104.7 * u.kilometer),
								(1522. * u.kilometer,-104.7 * u.kilometer),
								(1536. * u.kilometer,-104.7 * u.kilometer),
								(1536. * u.kilometer,-720. * u.kilometer),
								(0. * u.kilometer,-720. * u.kilometer),
								(0. * u.kilometer,-140. * u.kilometer)]))
#Asthenospheric_Mantle_Object.density = 3395 * u.kilogram / u.metre**3
Asthenospheric_Mantle_Object.density = GEO.LinearDensity(reference_density=3395. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Asthenospheric_Mantle_Object.radiogenicHeatProd = 0.02e-6 * u.watt / u.meter**3
Asthenospheric_Mantle_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Asthenospheric_Mantle_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Asthenospheric_Mantle_Object.viscosity = 3e20 * u.pascal * u.second #Newtonian?
Asthenospheric_Mantle_Object.add_melt_modifier(solidii.Mantle_Solidus, liquidii.Mantle_Liquidus, 
											latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
											meltFraction=0.,
											meltFractionLimit=0.02,
											meltExpansion=0.13, 
											viscosityChangeX1 = 0.01,
											viscosityChangeX2 = 0.03,
											viscosityChange = 1e-2)
#Asthenospheric_Mantle_Object.plasticity = pl.Rey_et_al_2014_LithosphericMantle
Asthenospheric_Mantle_Object.plasticity = GEO.DruckerPrager(name="Asthenosphere",
                                           cohesion=10. * u.megapascal,
                                           cohesionAfterSoftening=2 * u.megapascal,
                                           frictionCoefficient=0.577,
                                           frictionAfterSoftening=0.11,
                                           epsilon1=0.0, epsilon2=0.15)

Asthenospheric_Mantle_Object.stressLimiter = 200. * u.megapascal

Mantle_Wedge_Object = Model.add_material(name="Mantle Wedge", shape=GEO.shapes.Polygon([
                                                                (400. * u.kilometer,-42.0 * u.kilometer),
                                                                (544. * u.kilometer,-42.0 * u.kilometer),
                                                                (535. * u.kilometer,-52.0 * u.kilometer),
                                                                (525. * u.kilometer,-62.0 * u.kilometer),
								(584 * u.kilometer,-5 * u.kilometer),
								(414 * u.kilometer,-5 * u.kilometer)]))
#Mantle_Wedge_Object.density = 3370 * u.kilogram / u.metre**3
Mantle_Wedge_Object.density = GEO.LinearDensity(reference_density=3350. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Mantle_Wedge_Object.radiogenicHeatProd = 0.02e-6 * u.watt / u.meter**3
Mantle_Wedge_Object.diffusivity = 7.5e-7 * u.metre**2 / u.second 
Mantle_Wedge_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Mantle_Wedge_Object.viscosity = GEO.CompositeViscosity([0.01 * rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003, 0.01 *rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003])
#Mantle_Wedge_Object.viscosity = 6e20 * u.pascal * u.second
Mantle_Wedge_Object.add_melt_modifier(solidii.Mantle_Solidus, liquidii.Mantle_Liquidus, 
											latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
											meltFraction=0.,
											meltFractionLimit=0.08,
											meltExpansion=0.13, 
											viscosityChangeX1 = 0.02,
											viscosityChangeX2 = 0.08,
											viscosityChange = 1e-2)
#Mantle_Wedge_Object.plasticity = pl.Rey_et_al_2014_LithosphericMantle
Mantle_Wedge_Object.plasticity = GEO.DruckerPrager(name="Lithospheric mantle",
                                           cohesion=10. * u.megapascal,
                                           cohesionAfterSoftening=2 * u.megapascal,
                                           frictionCoefficient=0.577,
                                           frictionAfterSoftening=0.11,
                                           epsilon1=0.0, epsilon2=0.15)

Mantle_Wedge_Object.elasticity = GEO.Elasticity(shear_modulus=40e9 * u.pascal, observation_time=20000 * u.year)

Mantle_Wedge_Object.stressLimiter = 200. * u.megapascal


#Oceanic Sediment Layer
Oceanic_Sediment = Model.add_material(name="Oceanic Sediment", shape=GEO.shapes.Polygon([
                                                                (384 * u.kilometer,-5 * u.kilometer),
                                                                (852 * u.kilometer,-5 * u.kilometer),
                                                                (922 * u.kilometer,-7 * u.kilometer),
                                                                (380 * u.kilometer,-7 * u.kilometer)]))
Oceanic_Sediment.density = 2800 * u.kilogram / u.metre**3
Oceanic_Sediment.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Oceanic_Sediment.diffusivity = 1e-6 * u.metre**2 / u.second
Oceanic_Sediment.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Oceanic_Sediment.viscosity = 5e19 * u.pascal * u.second
Oceanic_Sediment.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                                                        latentHeatFusion= 250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                                                        meltFraction=0,
                                                                                        meltFractionLimit=1.,
                                                                                        meltExpansion=0.13,
                                                                                        viscosityChangeX1 = 0.15,
                                                                                        viscosityChangeX2 = 0.30,
                                                                                        viscosityChange = 1e-3)
Oceanic_Sediment.plasticity = GEO.DruckerPrager(name="Ocean sediments",
                                           cohesion=1.5 * u.megapascal,
                                           cohesionAfterSoftening=0.15 * u.megapascal,
                                           frictionCoefficient=0.577,
                                           frictionAfterSoftening=0.011,
                                           epsilon1=0.0, epsilon2=0.15)

Oceanic_Sediment.stressLimiter = 25. * u.megapascal

#Oceanic Upper-Crust
Oceanic_Upper_Crust = Model.add_material(name="Oceanic Upper-Crust", shape=GEO.shapes.Polygon([
                                                                (380 * u.kilometer,-7 * u.kilometer),
                                                                (922 * u.kilometer,-7 * u.kilometer),
                                                                (1022 * u.kilometer,-13 * u.kilometer),
                                                                (374 * u.kilometer,-13 * u.kilometer)]))
Oceanic_Upper_Crust.density = 2900 * u.kilogram / u.metre**3
Oceanic_Upper_Crust.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Oceanic_Upper_Crust.diffusivity = 1e-6 * u.metre**2 / u.second
Oceanic_Upper_Crust.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Oceanic_Upper_Crust.viscosity = 1 * rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998
Oceanic_Upper_Crust.plasticity = pl.Rey_et_al_2014_UpperCrust
Oceanic_Upper_Crust.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                                                        latentHeatFusion= 250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                                                        meltFraction=0,
                                                                                        meltFractionLimit=0.3,
                                                                                        meltExpansion=0.13,
                                                                                        viscosityChangeX1 = 0.15,
                                                                                        viscosityChangeX2 = 0.30,
                                                                                        viscosityChange = 1e-3)
Oceanic_Upper_Crust.stressLimiter = 300. * u.megapascal

# +
# #Oceanic Lower-Crust (This wasn't in the file?)
# Oceanic_Lower_Crust = Model.add_material(name="Oceanic Lower-Crust", shape=GEO.shapes.Polygon([
#                                                               (586 * u.kilometer,-13 * u.kilometer),
#                                                               (1536 * u.kilometer,-13 * u.kilometer),
#                                                               (1536 * u.kilometer,-20 * u.kilometer),
#                                                               (577 * u.kilometer,-20 * u.kilometer)]))
# Oceanic_Lower_Crust.density = 3200 * u.kilogram / u.metre**3
# Oceanic_Lower_Crust.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
# Oceanic_Lower_Crust.diffusivity = 1e-6 * u.metre**2 / u.second 
# Oceanic_Lower_Crust.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
# Oceanic_Lower_Crust.viscosity = 2e21 * u.pascal * u.second
# Oceanic_Lower_Crust.plasticity = pl.Rey_et_al_2014_LowerCrust
# Oceanic_Lower_Crust.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus, 
#                                                                                       latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
#                                                                                       meltFraction=0.2,
#                                                                                       meltFractionLimit=0.3,
#                                                                                       meltExpansion=0.13, 
#                                                                                       viscosityChangeX1 = 0.15,
#                                                                                       viscosityChangeX2 = 0.30,
#                                                                                       viscosityChange = 1e-3)
# Oceanic_Lower_Crust.stressLimiter = 150. * u.megapascal
# -

#Continental Sediments UP
Continental_Sediments_UP_Object = Model.add_material(name="Continental Upper Crust Upper Plate", shape=GEO.shapes.Polygon([
                                                                (0 * u.kilometer, 1 * u.kilometer),
                                                                (348 * u.kilometer, 1 * u.kilometer),
                                                                (384 * u.kilometer, -5 * u.kilometer),
                                                                (0 * u.kilometer,-5 * u.kilometer)]))

#Continental_Sediments_Object.density = 2580 * u.kilogram / u.metre**3
Continental_Sediments_UP_Object.density = GEO.LinearDensity(reference_density=2580. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Sediments_UP_Object.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Continental_Sediments_UP_Object.diffusivity = 1e-6 * u.metre**2 / u.second
Continental_Sediments_UP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Sediments_UP_Object.viscosity = 8e19 * u.pascal * u.second #1.0 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Continental_Sediments_UP_Object.plasticity = GEO.DruckerPrager(name="Continental Sediments UP",
                                                cohesion=1.5 * u.megapascal,
                                                cohesionAfterSoftening=0.15 * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Continental_Sediments_UP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                meltFraction=0.,
                                                meltFractionLimit=0.3,
                                                meltExpansion=0.13,
                                                viscosityChangeX1 = 0.15,
                                                viscosityChangeX2 = 0.30,
                                                viscosityChange = 1e-3)
Continental_Sediments_UP_Object.stressLimiter = 25. * u.megapascal

#Continental Upper-Crust_Layer
Continental_Upper_Crust_UP_Object = Model.add_material(name="Continental Upper Crust Upper Plate", shape=GEO.shapes.Polygon([
								(0 * u.kilometer, -5 * u.kilometer),
								(384 * u.kilometer, -5 * u.kilometer),
								(372 * u.kilometer, -17 * u.kilometer),
								(0 * u.kilometer, -17 * u.kilometer)]))

#Continental_Upper_Crust_Object.density = 2580 * u.kilogram / u.metre**3
Continental_Upper_Crust_UP_Object.density = GEO.LinearDensity(reference_density=2580. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Upper_Crust_UP_Object.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Continental_Upper_Crust_UP_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Continental_Upper_Crust_UP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Upper_Crust_UP_Object.viscosity = 1.0 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
#Continental_Upper_Crust_UP_Object.plasticity = pl.Rey_et_al_2014_UpperCrust
Continental_Upper_Crust_UP_Object.plasticity = GEO.DruckerPrager(name="Continental Upper Crust",
                                                cohesion=15. * u.megapascal,
                                                cohesionAfterSoftening=1.5 * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Continental_Upper_Crust_UP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus, 
											latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
											meltFraction=0.,
											meltFractionLimit=0.3,
											meltExpansion=0.13, 
											viscosityChangeX1 = 0.15,
											viscosityChangeX2 = 0.30,
											viscosityChange = 1e-3)
Continental_Upper_Crust_UP_Object.stressLimiter = 100. * u.megapascal


#Continental Lower Crust upper plate
Continental_Lower_Crust_UP_Object = Model.add_material(name="Continental Lower Crust Upper Plate", shape=GEO.shapes.Polygon([
								(0 * u.kilometer,-17 * u.kilometer),
								(372 * u.kilometer,-17 * u.kilometer),
								(350 * u.kilometer,-39 * u.kilometer),
								(0 * u.kilometer,-39 * u.kilometer)]))

#Continental_Lower_Crust_UP_Object.density = 2750 * u.kilogram / u.metre**3
Continental_Lower_Crust_UP_Object.density = GEO.LinearDensity(reference_density=2750. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Lower_Crust_UP_Object.radiogenicHeatProd = 0.5e-6 * u.watt / u.meter**3
Continental_Lower_Crust_UP_Object.diffusivity = 1e-6 * u.metre**2 / u.second 
Continental_Lower_Crust_UP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Lower_Crust_UP_Object.viscosity = 1 * rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012
Continental_Lower_Crust_UP_Object.plasticity = pl.Rey_et_al_2014_LowerCrust
Continental_Lower_Crust_UP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus, 
											latentHeatFusion= 250.0 * u.kilojoules / u.kilogram / u.kelvin,
											meltFraction=0.0,
											meltFractionLimit=0.3,
											meltExpansion=0.13, 
											viscosityChangeX1 = 0.15,
											viscosityChangeX2 = 0.30,
											viscosityChange = 1e-3)
Continental_Lower_Crust_UP_Object.stressLimiter = 150. * u.megapascal

# Decollement upper plate
Decollement_UP = Model.add_material(name="Decollement Upper Plate", shape=GEO.shapes.Polygon([
                                                                (0 * u.kilometer, -16 * u.kilometer),
                                                                (371 * u.kilometer, -16 * u.kilometer),
                                                                (373 * u.kilometer, -18 * u.kilometer),
                                                                (0 * u.kilometer, -18 * u.kilometer)]))
Decollement_UP.density = 2650 * u.kilogram / u.metre**3
Decollement_UP.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Decollement_UP.diffusivity = 1e-6 * u.metre**2 / u.second
Decollement_UP.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Decollement_UP.viscosity = 0.5 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Decollement_UP.plasticity = GEO.DruckerPrager(name="Decollement UP",
                                                cohesion=1. * u.megapascal,
                                                cohesionAfterSoftening=0.1 * u.megapascal,
                                                frictionCoefficient=0.2,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Decollement_UP.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                meltFraction=0.,
                                                meltFractionLimit=0.3,
                                                meltExpansion=0.13,
                                                viscosityChangeX1 = 0.15,
                                                viscosityChangeX2 = 0.30,
                                                viscosityChange = 1e-3)
Decollement_UP.stressLimiter = 30. * u.megapascal
# -

#Continental Sediments LP
Continental_Sediments_LP_Object = Model.add_material(name="Continental Sediments Lower Plate", shape=GEO.shapes.Polygon([
                                                                (852 * u.kilometer, -5 * u.kilometer),
                                                                (1062 * u.kilometer, 1 * u.kilometer),
                                                                (1536 * u.kilometer, 1 * u.kilometer),
                                                                (1536 * u.kilometer, -5 * u.kilometer)]))

#Continental_Sediments_Object.density = 2580 * u.kilogram / u.metre**3
Continental_Sediments_LP_Object.density = GEO.LinearDensity(reference_density=2580. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Sediments_LP_Object.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Continental_Sediments_LP_Object.diffusivity = 1e-6 * u.metre**2 / u.second
Continental_Sediments_LP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Sediments_LP_Object.viscosity = 8e19 * u.pascal * u.second #1.0 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Continental_Sediments_LP_Object.plasticity = GEO.DruckerPrager(name="Continental Sediments LP",
                                                cohesion=1.5 * u.megapascal,
                                                cohesionAfterSoftening=0.15 * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Continental_Sediments_LP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                meltFraction=0.,
                                                meltFractionLimit=0.3,
                                                meltExpansion=0.13,
                                                viscosityChangeX1 = 0.15,
                                                viscosityChangeX2 = 0.30,
                                                viscosityChange = 1e-3)
Continental_Sediments_LP_Object.stressLimiter = 25. * u.megapascal


# Upper crust lower plate

Continental_Upper_Crust_LP_Object = Model.add_material(name="Continental Upper Crust Lower Plate", shape=GEO.shapes.Polygon([
                                                                (852 * u.kilometer, -5 * u.kilometer),
                                                                (1536 * u.kilometer, -5 * u.kilometer),
                                                                (1536 * u.kilometer, -17 * u.kilometer),
                                                                (957 * u.kilometer,-17 * u.kilometer)]))

#Continental_Upper_Crust_LP_Object.density = 2580 * u.kilogram / u.metre**3
Continental_Upper_Crust_LP_Object.density = GEO.LinearDensity(reference_density=2580. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Upper_Crust_LP_Object.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Continental_Upper_Crust_LP_Object.diffusivity = 1e-6 * u.metre**2 / u.second
Continental_Upper_Crust_LP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Upper_Crust_LP_Object.viscosity = 1.0 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
#Continental_Upper_Crust_LP_Object.plasticity = pl.Rey_et_al_2014_UpperCrust
Continental_Upper_Crust_LP_Object.plasticity = GEO.DruckerPrager(name="Continental Upper Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Continental_Upper_Crust_LP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                                                        latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                                                        meltFraction=0.,
                                                                                        meltFractionLimit=0.3,
                                                                                        meltExpansion=0.13,
                                                                                        viscosityChangeX1 = 0.15,
                                                                                        viscosityChangeX2 = 0.30,
                                                                                        viscosityChange = 1e-3)
Continental_Upper_Crust_LP_Object.stressLimiter = 100. * u.megapascal

# Lower crust lower plate

Continental_Lower_Crust_LP_Object = Model.add_material(name="Continental Lower Crust Lower Plate", shape=GEO.shapes.Polygon([
                                                                (957 * u.kilometer, -17 * u.kilometer),
                                                                (1536 * u.kilometer,-17 * u.kilometer),
                                                                (1536 * u.kilometer,-39 * u.kilometer),
                                                                (1062 * u.kilometer, -39 * u.kilometer)]))
#Continental_Lower_Crust_LP_Object.density = 2750 * u.kilogram / u.metre**3
Continental_Lower_Crust_LP_Object.density = GEO.LinearDensity(reference_density=2750. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Continental_Lower_Crust_LP_Object.radiogenicHeatProd = 0.5e-6 * u.watt / u.meter**3
Continental_Lower_Crust_LP_Object.diffusivity = 1e-6 * u.metre**2 / u.second
Continental_Lower_Crust_LP_Object.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Continental_Lower_Crust_LP_Object.viscosity = 1 * rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998
Continental_Lower_Crust_LP_Object.plasticity = pl.Rey_et_al_2014_LowerCrust
Continental_Lower_Crust_LP_Object.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                                                        latentHeatFusion= 250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                                                        meltFraction=0.0,
                                                                                        meltFractionLimit=0.3,
                                                                                        meltExpansion=0.13,
                                                                                        viscosityChangeX1 = 0.15,
                                                                                        viscosityChangeX2 = 0.30,
                                                                                        viscosityChange = 1e-3)
Continental_Lower_Crust_LP_Object.stressLimiter = 150. * u.megapascal

# Decollement lower plate
Decollement_LP = Model.add_material(name="Decollement lower plate", shape=GEO.shapes.Polygon([
                                                                (922 * u.kilometer,-16 * u.kilometer),
                                                                (1536 * u.kilometer,-16 * u.kilometer),
                                                                (1536 * u.kilometer,-18 * u.kilometer),
                                                                (992 * u.kilometer,-18 * u.kilometer)]))
Decollement_LP.density = 2650 * u.kilogram / u.metre**3
Decollement_LP.radiogenicHeatProd = 1.5e-6 * u.watt / u.meter**3
Decollement_LP.diffusivity = 1e-6 * u.metre**2 / u.second
Decollement_LP.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Decollement_LP.viscosity = 0.5 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Decollement_LP.plasticity = GEO.DruckerPrager(name="Decollement LP",
                                                cohesion=1. * u.megapascal,
                                                cohesionAfterSoftening=0.1 * u.megapascal,
                                                frictionCoefficient=0.2,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)
Decollement_LP.add_melt_modifier(solidii.Crustal_Solidus, liquidii.Crustal_Liquidus,
                                                latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                meltFraction=0.,
                                                meltFractionLimit=0.3,
                                                meltExpansion=0.13,
                                                viscosityChangeX1 = 0.15,
                                                viscosityChangeX2 = 0.30,
                                                viscosityChange = 1e-3)
Decollement_LP.stressLimiter = 30. * u.megapascal

#Benioff faultt
Fault = Model.add_material(name="Faulted Crust", shape=GEO.shapes.Polygon([
						(579 * u.kilometer,-11 * u.kilometer),
					        (589 * u.kilometer,-11 * u.kilometer),
                                                (504 * u.kilometer,-96 * u.kilometer),
                                                (499 * u.kilometer,-91 * u.kilometer)]))
#Faulted_Crust.density = 3395 * u.kilogram / u.metre**3

Fault.density = GEO.LinearDensity(reference_density=3395. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 3.0e-5 * u.kelvin**-1)
Fault.radiogenicHeatProd = 0.0 * u.watt / u.meter**3
Fault.diffusivity = 1e-6 * u.metre**2 / u.second
Fault.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
Fault.viscosity = 5e19 * u.pascal * u.second
#Fault.viscosity = 0.02 * rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003

#Fault.plasticity = pl.Rey_et_al_2014_LithosphericMantle
Fault.plasticity = GEO.DruckerPrager(name="Faulted Mantle",
                                           cohesion=1.5 * u.megapascal,
                                           cohesionAfterSoftening=0.15 * u.megapascal,
                                           frictionCoefficient=0.577,
                                           frictionAfterSoftening=0.011,
                                           epsilon1=0.0, epsilon2=0.15)

Fault.stressLimiter = 25. * u.megapascal

Fault.add_melt_modifier(solidii.Mantle_Solidus, liquidii.Mantle_Liquidus,
                                                latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                                                meltFraction=0.,
                                                meltFractionLimit=0.02,
                                                meltExpansion=0.13,
                                                viscosityChangeX1 = 0.00,
                                                viscosityChangeX2 = 0.02,
                                                viscosityChange = 1e-3)

#===========================
# Damage oceanic lithosphere
#===========================

#def gaussian(zz, centre, widthx):
#    return ( np.exp( -(zz - centre)**2 / widthx ))


#widthz = GEO.nd(48.0 * u.kilometer) # width of the damaged region along z, for a normal distribution twice as large
#widthx = GEO.nd(480. * u.kilometer) # width of the damaged region along x, for a normal distribution twice as large

#maxDamage = 0.5
#angle = 0.0 # tan(degree*2*pi/360): (10>0.176327, 20>0.3639702, 30>0.5773503, 40>0.8390996, 45>1, 50>1.191754, 60>1.732051, 70>2.747477, 80>5.671282)

#centre = angle * Model.swarm.particleCoordinates.data[:,1] + (((GEO.nd(640. * u.kilometer)-angle*(widthz/2))))

#Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
#Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centre , widthx)
#Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], GEO.nd(-3. * u.kilometer), widthz*100)

#Air_Object_mask = Model.swarm.particleCoordinates.data[:,1] > GEO.nd(-5 * u.kilometer)
#Model.plasticStrain.data[Air_Object_mask] = 0.0

#noncrust_mask2 = Model.swarm.particleCoordinates.data[:,1] <= GEO.nd(-70 * u.kilometer)
#Model.plasticStrain.data[noncrust_mask2] = 0.0


#===========================

#====================================
# Damage oceanic lithosphere
#====================================

def gaussian(zz, centre, width):
    return ( np.exp( -(zz - centre)**2 / width ))

maxDamage = 0.25
centre = (GEO.nd(550. * u.kilometer), GEO.nd(-24. * u.kilometer))
width = GEO.nd(50 * u.kilometer)  # this gives a normal distribution
                                  # from about -100 km to 100 km

Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centre[0], width)
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], centre[1], width/50)

# The following lines make the random damage only apply to the crust
noncrust_mask = Model.swarm.particleCoordinates.data[:,1] <= GEO.nd(-70 * u.kilometer)
Air_Object_mask = Model.swarm.particleCoordinates.data[:,1] > GEO.nd(-5 * u.kilometer)
LateralExclusion0_mask = Model.swarm.particleCoordinates.data[:,0] < GEO.nd(150 * u.kilometer)
LateralExclusion1_mask = Model.swarm.particleCoordinates.data[:,0] > GEO.nd(1000 * u.kilometer)

Model.plasticStrain.data[noncrust_mask] = 0.0
Model.plasticStrain.data[Air_Object_mask] = 0.0

Model.plasticStrain.data[LateralExclusion0_mask] = 0.0
Model.plasticStrain.data[LateralExclusion1_mask] = 0.0

#==========================
# Damage Benioff plane
#==========================

#def gaussian(xx, centre, widthz):
#    return ( np.exp( -(xx - centre)**2 / widthz ))


#widthx = GEO.nd(3.0 * u.kilometer) # width of the damaged region along x, for a normal distribution twice as large
#widthz = GEO.nd(48. * u.kilometer) # width of the damaged region along z, for a normal distribution twice as large

#maxDamage = 0.5
#angle = 1.0 # tan(degree*2*pi/360): (10>0.176327, 20>0.3639702, 30>0.5773503, 40>0.8390996, 45>1, 50>1.191754, 60>1.732051, 70>2.747477, 80>5.671282)

#centre = angle * Model.swarm.particleCoordinates.data[:,1] + (((GEO.nd(618. * u.kilometer)-angle*(widthx/2))))

#Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
#Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centre , widthz)
#Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], GEO.nd(-3. * u.kilometer), widthx*100)

#Air_Object_mask = Model.swarm.particleCoordinates.data[:,1] > GEO.nd(-5 * u.kilometer)
#Model.plasticStrain.data[Air_Object_mask] = 0.0

#noncrust_mask2 = Model.swarm.particleCoordinates.data[:,1] <= GEO.nd(-100 * u.kilometer)
#Model.plasticStrain.data[noncrust_mask2] = 0.0


#"""
#  Benioff random damage
#====================================
#"""


#===========================


pts_lithosphere = GEO.circles_grid(radius=2.5*u.kilometer,
					minCoord=[0.0 * u.kilometer, Lithospheric_Mantle_Object.bottom], 
					maxCoord=[1636. * u.kilometer, Sticky_Air_Object.bottom])

Model.add_passive_tracers(name="Lithosphere", vertices=pts_lithosphere)

# +
# #Adding the surface tracers
# x_surface = np.array([GEO.nd(0*u.kilometer),GEO.nd(566*u.kilometer),GEO.nd(584*u.kilometer),
#               GEO.nd(1436*u.kilometer),GEO.nd(1536*u.kilometer)])
# y_surface = np.array([GEO.nd(1*u.kilometer),GEO.nd(1*u.kilometer),
#               GEO.nd(-5*u.kilometer),GEO.nd(-5*u.kilometer),GEO.nd(-3*u.kilometer)])

# +
# coords_surface = np.ndarray((len(x_surface),2))
# coords_surface[:,0] = x_surface
# coords_surface[:,1] = y_surface
# Model.add_passive_tracers(name = 'Surface', vertices = coords_surface)

# +
# #Adding the moho tracers
# x_moho = np.array([GEO.nd(0*u.kilometer),GEO.nd(550*u.kilometer), GEO.nd(584*u.kilometer)])
# y_moho = np.array([GEO.nd(-39*u.kilometer),GEO.nd(-39*u.kilometer), GEO.nd(-5*u.kilometer)])

# +
# coords_moho = np.ndarray((len(x_moho),2))
# coords_moho[:,0] = x_moho
# coords_moho[:,1] = y_moho
# Model.add_passive_tracers(name = 'Moho', vertices = coords_moho)
# -

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400))
		Fig.Points(Model.swarm, Model.materialField, colours='#d3fefb #F0F8FF #19497c #96c556 #b5dca8 #1a6615 #a58436 #e3e137 #fffb00 #dfa0b7 #6c6c6c #663547 #171717 #ac0000', fn_size=2.0, discrete=True)
	#     Fig.Points(Model.Lithosphere_tracers, pointSize=1.0, colour="black")
		Fig.show()


# Setting up Temperature Boundary conditions.

Model.set_temperatureBCs(top = 293.15 * u.kelvin,
                        bottom = 1673.15 * u.kelvin,
                        materials=[(Air_Object, 293.15 * u.kelvin),(Sticky_Air_Object, 293.15 * u.kelvin),(Asthenospheric_Mantle_Object, 1603.15 * u.kelvin)])

Model.swarm.allow_parallel_nn = True

Model.init_model()

# +
Oceanic_Sediment.phase_changes = GEO.PhaseChange((Model.temperature > GEO.nd(523 * u.kelvin)), Eclogite_Object.index)
Oceanic_Upper_Crust.phase_changes = GEO.PhaseChange((Model.temperature > GEO.nd(523 * u.kelvin)), Eclogite_Object.index)
#Oceanic_Lower_Crust.phase_changes = GEO.PhaseChange((Model.temperature > GEO.nd(523 * u.kelvin)), Eclogite_Object.index)


vel_left_in = 0.30 * u.centimeter / u.year
vel_left_out = -0.066 * u.centimeter / u.year
depth_left_z1 = 20. * u.kilometre
depth_left_z2 = 10. * u.kilometre
depth_left_z3 = -80. * u.kilometre
depth_left_z4 = -140. * u.kilometre
depth_left_z5 = -160. * u.kilometre
grad_left_z1z2 = (GEO.nd(0 * u.centimeter / u.year) - GEO.nd(vel_left_in)) / (GEO.nd(depth_left_z1) - GEO.nd(depth_left_z2))
grad_left_z3z4 = (GEO.nd(vel_left_in) - GEO.nd(0 * u.centimeter / u.year)) / (GEO.nd(depth_left_z3) - GEO.nd(depth_left_z4))
grad_left_z4z5 = (GEO.nd(0 * u.centimeter / u.year) - GEO.nd(vel_left_out)) / (GEO.nd(depth_left_z4) - GEO.nd(depth_left_z5))

# +
# conditions_left_wall = [(Model.y > GEO.nd(depth_left_top), GEO.nd(vel_left_top)),
#                                               (Model.y > GEO.nd(depth_left_down),
#                                               ((m_left * Model.y) + b_left)),
#                                               (True, GEO.nd(vel_left_down))]
# -

conditions_left_wall = [(Model.y > GEO.nd(depth_left_z1), GEO.nd(0 * u.centimeter / u.year)),
                        (Model.y > GEO.nd(depth_left_z2), grad_left_z1z2 * (Model.y - GEO.nd(depth_left_z1)) + GEO.nd(0 * u.centimeter / u.year)),
                                                (Model.y > GEO.nd(depth_left_z3), GEO.nd(vel_left_in)),
                                                (Model.y > GEO.nd(depth_left_z4), grad_left_z3z4 * (Model.y - GEO.nd(depth_left_z3)) + GEO.nd(vel_left_in)),
                                                (Model.y > GEO.nd(depth_left_z5), grad_left_z4z5 * (Model.y - GEO.nd(depth_left_z4))),
                                                (True, GEO.nd(vel_left_out))]

function_left_wall = fn.branching.conditional(conditions_left_wall)

# +
# vel_right_top = -1.25 * u.centimeter / u.year
# vel_right_down = 0.219541 * u.centimeter / u.year
# depth_right_top = -110. * u.kilometre
# depth_right_down = -120. * u.kilometre
# m_right = (GEO.nd(vel_right_top) - GEO.nd(vel_right_down)) / (GEO.nd(depth_right_top) - GEO.nd(depth_right_down))
# b_right = GEO.nd(vel_right_top) - (GEO.nd(m_right)*(GEO.nd(depth_right_top)))
# -

vel_right_in = -0.30 * u.centimeter / u.year
vel_right_out = 0.066 * u.centimeter / u.year #0.6818181818
depth_right_z1 = 20. * u.kilometre
depth_right_z2 = 10. * u.kilometre
depth_right_z3 = -80. * u.kilometre
depth_right_z4 = -140. * u.kilometre
depth_right_z5 = -160. * u.kilometre
grad_right_z1z2 = (GEO.nd(0 * u.centimeter / u.year) - GEO.nd(vel_right_in)) / (GEO.nd(depth_right_z1) - GEO.nd(depth_right_z2))
grad_right_z3z4 = (GEO.nd(vel_right_in) - GEO.nd(0 * u.centimeter / u.year)) / (GEO.nd(depth_right_z3) - GEO.nd(depth_right_z4))
grad_right_z4z5 = (GEO.nd(0 * u.centimeter / u.year) - GEO.nd(vel_right_out)) / (GEO.nd(depth_right_z4) - GEO.nd(depth_right_z5))

# +
# conditions_right_wall = [(Model.y > GEO.nd(depth_right_top), GEO.nd(vel_right_top)),
#                                               (Model.y > GEO.nd(depth_right_down),
#                                               ((m_right * Model.y) + b_right)),
#                                               (True, GEO.nd(vel_right_down))]
# -

conditions_right_wall = [(Model.y > GEO.nd(depth_right_z1), GEO.nd(0 * u.centimeter / u.year)),
                        (Model.y > GEO.nd(depth_right_z2), grad_right_z1z2 * (Model.y - GEO.nd(depth_right_z1)) + GEO.nd(0 * u.centimeter / u.year)),
                                                (Model.y > GEO.nd(depth_right_z3), GEO.nd(vel_right_in)),
                                                (Model.y > GEO.nd(depth_right_z4), grad_right_z3z4 * (Model.y - GEO.nd(depth_right_z3)) + GEO.nd(vel_right_in)),
                                                (Model.y > GEO.nd(depth_right_z5), grad_right_z4z5 * (Model.y - GEO.nd(depth_right_z4))),
                                                (True, GEO.nd(vel_right_out))]

function_right_wall = fn.branching.conditional(conditions_right_wall)

# +
# Model.set_velocityBCs(left=[function_left_wall, None],
#                       right=[function_right_wall, None],
#                       bottom=GEO.LecodeIsostasy(reference_mat=Asthenospheric_Mantle_Object, average=True))
# -

Model.set_velocityBCs(left=[function_left_wall, None],
                      right=[function_right_wall, None],
                      bottom=[None,0.],
                      top=[None,0.])

# +
#Model.plasticStrain.data[air_mask] = 0.0

# +
# #Surface processes
# Model.surfaceProcesses = GEO.surfaceProcesses.SedimentationThreshold(air=[Air_Object], sediment=[Continental_Sediment_1, Oceanic_Sediment], threshold = -5. * u.kilometre)
# -

#Model.surfaceProcesses = GEO.surfaceProcesses.ErosionThreshold(sediment=[Oceanic_Sediment, Continental_Sediments_UP_Object, Continental_Sediments_LP_Object], air=[Sticky_Air_Object],threshold = 10. * u.kilometre)


Model.init_model(temperature="steady-state", pressure="lithostatic")

Model.set_temperatureBCs(top = 293.15 * u.kelvin,
                        bottom = 1673.15 * u.kelvin,
                        materials=[(Air_Object, 293.15 * u.kelvin), (Sticky_Air_Object, 293.15 * u.kelvin)])

# +
# npoints = 1536 # This is the number of points used to define the surface
# coords = np.ndarray((npoints, 2))

# coords[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints)
# coords[:, 1] = 0.

# surface_tracers = Model.add_passive_tracers(name="Surface", vertices=coords)
# coords[:, 1] -= -GEO.nd(uppermantle.top)
# moho_tracers = Model.add_passive_tracers(name="Moho", vertices=coords)
# -

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400), title="Temperature Field (Kelvin)", quality=3)
		Fig.Surface(Model.mesh, GEO.dimensionalise(Model.temperature, u.kelvin), colours="coolwarm", drawSides = "xyzXYZ")
	#     Fig.Points(Model.Surface_tracers, pointSize=10, colour="black")
	#     Fig.Points(Model.Moho_tracers, pointSize=10, colour="red")
		Fig.show()

if plotting_bool == True:
	distances, temperature = GEO.extract_profile(Model.temperature, line = [(100.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])
	distances, pressure = GEO.extract_profile(Model.pressureField, line = [(100.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])
	distances, viscosity = GEO.extract_profile(Model.viscosityField, line = [(100.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])


# +
import matplotlib.pyplot as plt

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15,7))
		ax1.plot(GEO.dimensionalise(temperature , u.celsius), GEO.dimensionalise(distances, u.kilometer))
		ax1.set_xlabel("Temperature [$C^{\circ}$]")
		ax1.set_ylabel("Depth [km]")
		ax1.set_ylim(250, 40)
		ax1.set_xlim(0,1700)
		ax1.set_title("Temperature profile")


		ax2.plot(GEO.dimensionalise(pressure, u.gigapascal), GEO.dimensionalise(distances, u.kilometer))
		ax2.set_xlabel("Pressure [GPa]")
		ax2.set_title("Pressure profile")
		ax2.set_ylim(250, 40)
		ax2.set_xlim(0,10)

		ax3.plot(GEO.dimensionalise(viscosity, u.pascal * u.second), GEO.dimensionalise(distances, u.kilometer))
		ax3.set_xlabel("Viscosity [Pa.s]")
		ax3.set_title("Viscosity profile")
		ax3.set_ylim(250, 40)
		ax3.set_xlim(1e18,1.5e23)
# -


if plotting_bool == True:
	distances2, temperature2 = GEO.extract_profile(Model.temperature, line = [(2.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])
	distances2, pressure2 = GEO.extract_profile(Model.pressureField, line = [(2.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])
	distances2, viscosity2 = GEO.extract_profile(Model.viscosityField, line = [(2.* u.kilometer, Model.top), (100.* u.kilometer, Model.bottom)])

# +
import matplotlib.pyplot as plt

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15,7))
		ax1.plot(GEO.dimensionalise(1.01*temperature, u.celsius), GEO.dimensionalise(distances, u.kilometer))
		ax1.plot(GEO.dimensionalise(temperature2, u.celsius), GEO.dimensionalise(distances2, u.kilometer))
		ax1.set_xlabel("Temperature [$C^{\circ}$]")
		ax1.set_ylabel("Depth [km]")
		ax1.set_ylim(250, 40)
		ax1.set_xlim(0,1700)
		ax1.set_title("Temperature profile")


		ax2.plot(GEO.dimensionalise(pressure2, u.gigapascal), GEO.dimensionalise(distances2, u.kilometer))
		ax2.set_xlabel("Pressure [GPa]")
		ax2.set_title("Pressure profile")
		ax2.set_ylim(250, 40)
		ax2.set_xlim(0,10)

		ax3.plot(GEO.dimensionalise(viscosity2, u.pascal * u.second), GEO.dimensionalise(distances2, u.kilometer))
		ax3.set_xlabel("Viscosity [Pa.s]")
		ax3.set_title("Viscosity profile")
		ax3.set_ylim(250, 40)
		ax3.set_xlim(1e18,1.5e23)
# -

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400))
		Fig.Surface(Model.mesh, GEO.dimensionalise(Model.pressureField, u.gigapascal))
		Fig.show()


if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400), title="Viscosity Field [Pa.s]", quality=3)
		Fig.Points(Model.swarm, 
				   GEO.dimensionalise(Model.viscosityField, u.pascal * u.second),
				   logScale=True,
				   fn_size=3.0)
		Fig.show()


if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400), title="Strain Field", quality=3)
		Fig.Surface(Model.mesh, GEO.dimensionalise(Model.strainRateField, u.second**-1))
		Fig.show()

if plotting_bool == True:
	if GEO.nProcs == 1:
		Fig = vis.Figure(figsize=(1200,400), title="Velocity (cm/yr)", quality=2)
		Fig.Surface(Model.mesh, Model.velocityField[0])
		Fig.show()


def post_hook():  
#    Stop any brittle yielding near the edges of the model
    coords = fn.input()
    zz = (coords[0] - GEO.nd(Model.minCoord[0])) / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*20.0) + fn.math.tanh((1.0-zz)*20.0) - fn.math.tanh(20.0), 4)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm)


# +
# To check default parameters run: GEO.rcParamsDefault

# +
GEO.rcParams["initial.nonlinear.tolerance"] = 1e-3
GEO.rcParams["nonlinear.tolerance"] =  1e-3  # 5e-4
GEO.rcParams["nonlinear.min.iterations"] = 1
GEO.rcParams["nonlinear.max.iterations"] = 100
GEO.rcParams["CFL"] = 0.1 # Courant factor

GEO.rcParams["advection.diffusion.method"] = "SLCN" #SUPG or SLCN
GEO.rcParams["shear.heating"] = True
GEO.rcParams["surface.pressure.normalization"] = True  # Make sure the top of the model is approximately 0 Pa

GEO.rcParams["swarm.particles.per.cell.2D"] = 60
GEO.rcParams["popcontrol.split.threshold"] = 0.15 #0.95
GEO.rcParams["popcontrol.max.splits"] = 100
GEO.rcParams["popcontrol.particles.per.cell.2D"] = 60
# -

GEO.rcParams["default.outputs"].append("projMeltField")
GEO.rcParams["default.outputs"].append("projStressTensor")
GEO.rcParams["default.outputs"].append("projSurfaceProcess")

# +
# population_control = uw.swarm.PopulationControl(swarm, 
#                                                 aggressive=True,splitThreshold=0.15, 
#                                                 maxDeletions=2,maxSplits=10,
#                                                 particlesPerCell=20)
# -

# This is a bunch of solver options. You can try playing with them, but these should be good enough.
# For more details: https://underworld2.readthedocs.io/en/v2.14.0b/build/user_guide/08_StokesSolver.html
solver = Model.solver

# Decide whether to use mumps or multigrid
if resolution[0] * resolution[1] < 1e6:
    print("Using mumps")
    solver.set_inner_method("mumps") # "mumps" is a good alternative for "lu" 
    #solver.set_inner_method("lu") # use "lu" direct solve and large penalty (if running in serial)
    #solver.set_inner_method("superludist")
else:
    print("Using multigrid with coarse mumps")
    #solver.options.mg.levels = 4
    solver.options.A11.mg_coarse_pc_factor_mat_solver_package = "mumps"
    solver.options.A11.mg_coarse_pc_type = "lu"
    solver.options.A11.mg_coarse_ksp_type = "preonly"
    solver.options.A11.mg_coarse_ksp_view = ""

solver.set_penalty(1e5) #1e6
#GEO.rcParams["initial.nonlinear.tolerance"] = 1e-3

solver.options.A11.ksp_rtol=1e-8 #Try 1e-4
solver.options.A11.ksp_set_min_it_converge = 10
solver.options.A11.use_previous_guess = True
solver.options.scr.ksp_rtol=1e-6 #Try 1e-3
solver.options.scr.use_previous_guess = True
solver.options.scr.ksp_set_min_it_converge = 10
solver.options.scr.ksp_type = "cg"

# +
#solver.options.main.help = ""
solver.options.main.remove_constant_pressure_null_space=True
#solver.options.main.Q22_pc_type = "uwscale"
#solver.set_penalty(1e6)
#GEO.rcParams["initial.nonlinear.tolerance"] = 1e-3

# Penalty method with mult grid (mg)
#solver.set_penalty(1.0)
#solver.solve(nonLinearIterate=True, nonLinearTolerance=0.01) 
# -

Model.run_for(50000000. * u.years, checkpoint_interval = 5e4 * u.year, restartStep=488, restartDir="Col262")





