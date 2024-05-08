# %%
import matplotlib.pyplot as plt
import underworld as uw
from underworld import UWGeodynamics as GEO
import underworld.function as fn
from underworld.UWGeodynamics.surfaceProcesses import SedimentationThreshold
import numpy as np
import scipy
import os.path
from mpi4py import MPI
import argparse
import math

# %%
u = GEO.UnitRegistry

GEO.rcParams["CFL"] = 0.2
GEO.rcParams["shear.heating"] = True
GEO.rcParams["popcontrol.split.threshold"] = 0.95
GEO.rcParams["popcontrol.max.splits"] = 100
GEO.rcParams["swarm.particles.per.cell.2D"] = 60
GEO.rcParams["popcontrol.particles.per.cell.2D"] = 60
GEO.rcParams["advection.diffusion.method"] = "SLCN"


resolution = (1360,440) #500 m resolution

half_rate = 1 * u.centimeter / u.year
model_length = 680e3 * u.meter
model_height = 210e3 * u.meter
surfaceTemp = 293.15 * u.degK
baseModelTemp = 1603.15 * u.degK
bodyforce = 3150 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT

# %%
Model = GEO.Model(elementRes=resolution, 
                  minCoord=(0. * u.kilometer, -190. * u.kilometer), 
                  maxCoord=(680. * u.kilometer, 30. * u.kilometer), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))


Model.outputDir="NRiftExtInv"

Model.diffusivity = 9e-7 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)

# %%
air = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0 * u.kilometer))
air.density = 1. * u.kilogram / u.metre**3
air.diffusivity = 1e-5 * u.metre**2 / u.second
air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
air.compressibility = None

# %%
sediment = Model.add_material(name="Sediment")
sediment.density           = GEO.LinearDensity(reference_density=2300. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
sediment.radiogenicHeatProd   = 1.15e-6 * u.microwatt / u.meter**3

Sediment1 = Model.add_material(name="Sediment Layer 1", shape=GEO.shapes.Layer(top=0. * u.kilometer, bottom=-1. * u.kilometer))
Sediment1.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment1.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment2 = Model.add_material(name="Sediment Layer 2", shape=GEO.shapes.Layer(top=-1. * u.kilometer, bottom=-2. * u.kilometer))
Sediment2.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment2.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment3 = Model.add_material(name="Sediment Layer 3", shape=GEO.shapes.Layer(top=-2. * u.kilometer, bottom=-3. * u.kilometer))
Sediment3.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment3.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment4 = Model.add_material(name="Sediment Layer 4", shape=GEO.shapes.Layer(top=-3. * u.kilometer, bottom=-4. * u.kilometer))
Sediment4.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment4.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment5 = Model.add_material(name="Sediment Layer 5", shape=GEO.shapes.Layer(top=-4. * u.kilometer, bottom=-5. * u.kilometer))
Sediment5.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment5.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment6 = Model.add_material(name="Sediment Layer 6", shape=GEO.shapes.Layer(top=-5. * u.kilometer, bottom=-6. * u.kilometer))
Sediment6.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
Sediment6.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

# %%
continentalcrustL3 = Model.add_material(name="Continental Crust Layer3", shape=GEO.shapes.Layer(top=-6. * u.kilometer, bottom=-9. * u.kilometer))
continentalcrustL3.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL3.density  = GEO.LinearDensity(reference_density  = 2675. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL4 = Model.add_material(name="Continental Crust Layer4", shape=GEO.shapes.Layer(top=-9. * u.kilometer, bottom=-12. * u.kilometer))
continentalcrustL4.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL4.density  = GEO.LinearDensity(reference_density  = 2675. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL5 = Model.add_material(name="Continental Crust Layer5", shape=GEO.shapes.Layer(top=-12. * u.kilometer, bottom=-15. * u.kilometer))
continentalcrustL5.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL5.density  = GEO.LinearDensity(reference_density  = 2700. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL6 = Model.add_material(name="Continental Crust Layer6", shape=GEO.shapes.Layer(top=-15. * u.kilometer, bottom=-18. * u.kilometer))
continentalcrustL6.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL6.density  = GEO.LinearDensity(reference_density  = 2725. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL7 = Model.add_material(name="Continental Crust Layer7", shape=GEO.shapes.Layer(top=-18. * u.kilometer, bottom=-21. * u.kilometer))
continentalcrustL7.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL7.density  = GEO.LinearDensity(reference_density  = 2750. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL8 = Model.add_material(name="Continental Crust Layer8", shape=GEO.shapes.Layer(top=-21. * u.kilometer, bottom=-24. * u.kilometer))
continentalcrustL8.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL8.density  = GEO.LinearDensity(reference_density  = 2775. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL9 = Model.add_material(name="Continental Crust Layer9", shape=GEO.shapes.Layer(top=-24. * u.kilometer, bottom=-27. * u.kilometer))
continentalcrustL9.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL9.density  = GEO.LinearDensity(reference_density  = 2800. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

continentalcrustL10 = Model.add_material(name="Continental Crust Layer10", shape=GEO.shapes.Layer(top=-27. * u.kilometer, bottom=-36. * u.kilometer))
continentalcrustL10.radiogenicHeatProd = 1.15e-6 * u.watt / u.meter**3
continentalcrustL10.density  = GEO.LinearDensity(reference_density  = 2825. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

# %%
uppermantle = Model.add_material(name="Upper Mantle", shape=GEO.shapes.Layer(top=-36. * u.kilometer, bottom=-154. * u.kilometer))
uppermantle.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
uppermantle.temperatureLimiter = 1603. * u.kelvin

asthenosphere = Model.add_material(name="Asthenosphere", shape=GEO.shapes.Layer(top=-154. * u.kilometer, bottom=Model.bottom))
asthenosphere.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3,
                                         thermalExpansivity= 2.8e-5 * u.kelvin**-1)
asthenosphere.temperatureLimiter = 1603. * u.kelvin

# %%
rh = GEO.ViscousCreepRegistry()

Model.minViscosity = 1e18 * u.pascal * u.second
Model.maxViscosity = 5e23 * u.pascal * u.second

air.viscosity = 1e18 * u.pascal * u.second

sediment.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990

Sediment1.viscosity = 1e20 * u.pascal *u.second
Sediment2.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Sediment3.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Sediment4.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Sediment5.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Sediment6.viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990

continentalcrustL3.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
continentalcrustL4.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
continentalcrustL5.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
continentalcrustL6.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
continentalcrustL7.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
continentalcrustL8.viscosity = rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012
continentalcrustL9.viscosity = rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012
continentalcrustL10.viscosity = rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012

uppermantle.viscosity = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
asthenosphere.viscosity = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003

# %%
sediment.plasticity = GEO.DruckerPrager(cohesion=1.0 * u.megapascal,
                               cohesionAfterSoftening=0.1 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
sediment.stressLimiter = 100 * u.megapascal

Sediment1.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment1.stressLimiter = 100 * u.megapascal

Sediment2.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment2.stressLimiter = 100 * u.megapascal

Sediment3.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment3.stressLimiter = 100 * u.megapascal

Sediment4.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment4.stressLimiter = 100 * u.megapascal

Sediment5.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment5.stressLimiter = 100 * u.megapascal

Sediment6.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
Sediment6.stressLimiter = 100 * u.megapascal

continentalcrustL3.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.25)
continentalcrustL3.stressLimiter = 125 * u.megapascal

continentalcrustL4.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.25)
continentalcrustL4.stressLimiter = 150 * u.megapascal

continentalcrustL5.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.25)
continentalcrustL5.stressLimiter = 150 * u.megapascal

continentalcrustL6.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.25)
continentalcrustL6.stressLimiter = 150 * u.megapascal

continentalcrustL7.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.25)
continentalcrustL7.stressLimiter = 150 * u.megapascal

continentalcrustL8.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
continentalcrustL8.stressLimiter = 150 * u.megapascal

continentalcrustL9.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
continentalcrustL9.stressLimiter = 150 * u.megapascal

continentalcrustL10.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.1154,
                                                epsilon1=0.0, epsilon2=0.15)
continentalcrustL10.stressLimiter = 150 * u.megapascal


uppermantle.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.0577,
                                                epsilon1=0.0, epsilon2=0.15)
uppermantle.stressLimiter = 250 * u.megapascal

asthenosphere.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=20. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.577,
                                                frictionAfterSoftening=0.0577,
                                                epsilon1=0.0, epsilon2=0.15)
asthenosphere.stressLimiter = 250 * u.megapascal



# ## Melt

# %%


solidii = GEO.SolidusRegistry()
my_crust_solidus = GEO.Solidus(A1=923 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.2e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mid_crust_solidus = GEO.Solidus(A1=1263 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.2e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mantle_solidus = solidii.Mantle_Solidus

liquidii = GEO.LiquidusRegistry()
my_crust_liquidus = GEO.Liquidus(A1=1423 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.6e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mid_crust_liquidus = GEO.Liquidus(A1=1763 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.6e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mantle_liquidus = liquidii.Mantle_Liquidus


continentalcrustL3.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

continentalcrustL4.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

continentalcrustL5.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

continentalcrustL6.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )


continentalcrustL7.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )


continentalcrustL8.add_melt_modifier(my_crust_solidus, my_crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

continentalcrustL9.add_melt_modifier(my_crust_solidus, my_crust_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

continentalcrustL10.add_melt_modifier(my_crust_solidus, my_crust_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-3
                        )

uppermantle.add_melt_modifier(mantle_solidus, mantle_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.02,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.001,
                         viscosityChangeX2 = 0.03,
                         viscosityChange = 1e-2
                        )

asthenosphere.add_melt_modifier(mantle_solidus, mantle_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.02,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.001,
                         viscosityChangeX2 = 0.03,
                         viscosityChange = 1e-2
                        )



# %%

# %%

Model.set_temperatureBCs(top=293.15 * u.degK, nodeSets = [(air.shape, 293.15 * u.degK)])

Model.set_heatFlowBCs(bottom=(0.020 * u.watt / u.metre**2, asthenosphere))

Model.init_model()

Model.set_temperatureBCs(nodeSets = [(air.shape, 293.15 * u.degK), (asthenosphere.shape, 1603.15 * u.degK)])

Model.init_model()

#Model.set_temperatureBCs(top=293.15 * u.degK, bottom=1603.15 * u.degK)

# %%
import underworld.function as fn

velocity = 2*1.073 * u.centimeter / u.year

Model.set_velocityBCs(left=[0.0 * u.centimeter / u.year, None],
                      right=[-velocity, None],
                      top=[None, None],
                     bottom=GEO.LecodeIsostasy(reference_mat=asthenosphere, average=False))

# Model.set_velocityBCs(left=[0., None],
#                       right=[0., None],
#                       top=[None, None],
#                      bottom=GEO.LecodeIsostasy(reference_mat=asthenosphere, average=False))

# Model.set_stressBCs(bottom=[0., 6.48174e9 * u.pascal])


#velocity = 0.375 * u.centimeter / u.year

#transition = 10. * u.kilometer

#conditions_left = [(Model.y > GEO.nd(uppermantle.bottom), GEO.nd(velocity)),
#                   (Model.y > GEO.nd(uppermantle.bottom) -  GEO.nd(transition),
#                   GEO.nd(velocity)*(GEO.nd(transition)-(Model.y-GEO.nd(uppermantle.bottom)))/GEO.nd(transition)),
#                   (True, GEO.nd(0.6793355 * u.centimeter / u.year))]

#fn_condition_left = fn.branching.conditional(conditions_left)
#
#conditions_right = [(Model.y > GEO.nd(uppermantle.bottom), GEO.nd(-velocity)),
#                   (Model.y > GEO.nd(uppermantle.bottom) -  GEO.nd(transition),
#                   GEO.nd(-velocity)*(GEO.nd(transition)-(Model.y-GEO.nd(uppermantle.bottom)))/GEO.nd(transition)),
#                   (True, GEO.nd(-0.6793355 * u.centimeter / u.year))]
#
#fn_condition_right = fn.branching.conditional(conditions_right)

#Model.set_velocityBCs(left=[fn_condition_left, 0.0 * u.centimeter / u.year],
#                      right=[fn_condition_right, 0.0 * u.centimeter / u.year],
#                      bottom=GEO.LecodeIsostasy(reference_mat=uppermantle, average=False))





# %%
def gaussian(xx, centre, width):
    return ( np.exp( -(xx - centre)**2/width))

maxDamage = 0.25
centre = (GEO.nd(340. * u.kilometer), GEO.nd(-40. * u.kilometer))
width = GEO.nd(75. * u.kilometer)  # this gives a normal distribution

Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centre[0], width)
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], centre[1], width*100)

air_mask = Model.swarm.particleCoordinates.data[:,1] > GEO.nd(0 * u.kilometer)

Model.plasticStrain.data[air_mask] = 0.0

npoints = int(Model.length.magnitude)
coords = np.ndarray((npoints, 2))
coords[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints)
coords[:, 1] = GEO.nd(Sediment1.top)

#coords_surface2 = np.ndarray((npoints_surface, 2))
#coords_surface2[:, 0] = np.linspace(Model.minCoord[0], Model.maxCoord[0], npoints_surface)
#coords_surface2[:, 1] = 0. * u.kilometre

#Model.surfaceProcesses = GEO.surfaceProcesses.velocitySurface_2D(airIndex = air.index,sedimentIndex = sediment.index, sedimentationRate = 0.6 * u.millimeter / u.year, erosionRate =  0.1 * u.millimeter / u.year, surfaceElevation = 0.0*u.kilometer, surfaceArray = coords)
npoints_surface = 1000
coords_surface = np.ndarray((npoints_surface, 2))
coords_surface[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints_surface)
coords_surface[:, 1] = GEO.nd(0. * u.kilometre)
surface_tracers = Model.add_passive_tracers(name="Topography", vertices=coords_surface)


npoints_moho = 1000
coords_moho = np.ndarray((npoints_moho, 2))
coords_moho[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints_moho)
coords_moho[:, 1] = GEO.nd(0. * u.kilometre) - GEO.nd(continentalcrustL10.bottom)
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=coords_moho)


coords_FSE_Crust = GEO.circles_grid(radius = 1.5 * u.kilometer,
                           minCoord=[Model.minCoord[0], continentalcrustL10.bottom],
                           maxCoord=[Model.maxCoord[0], 0.*u.kilometer])

FSE_Crust = Model.add_passive_tracers(name="FSE_Crust", vertices=coords_FSE_Crust)


coords_FSE_Mantle = GEO.circles_grid(radius=1.5 * u.kilometer, 
                    minCoord=[Model.minCoord[0], -154.*u.kilometer], 
                    maxCoord=[Model.maxCoord[0], continentalcrustL10.bottom])

FSE_Mantle = Model.add_passive_tracers(name="FSE_Mantle", vertices=coords_FSE_Mantle)



uppermantle.phase_changes = GEO.PhaseChange((Model.temperature > GEO.nd(1523.*u.kelvin)), asthenosphere.index)



#Model.set_heatFlowBCs(bottom=(0.036 * u.watt / u.metre**2, asthenosphere))

# %%
Model.surfaceProcesses = GEO.surfaceProcesses.velocitySurface_2D( 
    airIndex=air.index, sedimentIndex=sediment.index,
    sedimentationRate= 1. *u.millimeter / u.year, erosionRate= 1. * u.millimeter / u.year,
    surfaceElevation=0.*u.kilometer,
    surfaceArray = coords_surface, updateSurfaceRB = 5. * u.kilometer, updateSurfaceLB = 5. * u.kilometer)


Model.set_temperatureBCs(top=293.15 * u.degK, nodeSets = [(air.shape, 293.15 * u.degK)])


solver = Model.solver

# Decide whether to use mumps or multigrid
if resolution[0] * resolution[1] < 1e6:
    print("Using mumps")
    solver.set_inner_method("mumps")
else:
    print("Using multigrid with coarse mumps")
    #solver.options.mg.levels = 4
    solver.options.A11.mg_coarse_pc_factor_mat_solver_package = "mumps"
    solver.options.A11.mg_coarse_pc_type = "lu"
    solver.options.A11.mg_coarse_ksp_type = "preonly"
    #solver.options.A11.mg_coarse_ksp_view = ""

solver.options.A11.ksp_rtol=1e-8  #1e-4
solver.options.A11.ksp_set_min_it_converge = 10
solver.options.A11.use_previous_guess = True
solver.options.scr.ksp_rtol=1e-6  #1e-3
solver.options.scr.use_previous_guess = True
solver.options.scr.ksp_set_min_it_converge = 10
solver.options.scr.ksp_type = "cg"

#solver.options.main.help = ""
solver.options.main.remove_constant_pressure_null_space=True
#solver.options.main.Q22_pc_type = "uwscale"
solver.set_penalty(1.0e2)
#Model.solver = solver

def post_hook():  
    
    coords = fn.input()
    zz = (coords[0] - GEO.nd(Model.minCoord[0])) / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*20.0) + fn.math.tanh((1.0-zz)*20.0) - fn.math.tanh(20.0), 4)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm)

Model.run_for(9000000 * u.years, checkpoint_interval=100000. * u.year, restartStep=152, restartDir="NRiftExtInv")
