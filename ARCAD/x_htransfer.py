import math
import numpy as np
from ARC_settings import *

from x_upperliquidreservoir import UpperReservoirLiquidSectionLength, UpperReservoirLiquidVolume
from x_definitions import Coolant_CoolantOutlet, Coolant_CoolantOutlet, Coolant_CoolantOutlet, Coolant_CoolantOutlet, Expander_CoolantOutlet, Steel_CoolantOutlet
from x_flowconversion import RelativeTubeArea

# Intitial ARC temperature
T0ARC   = CoolantOutletTemperature
TARC    = T0ARC
T       = T0ARC
T0ARC_C = T0ARC - 273

# Outlet tube outer diameter (m)
Dh = CoolantOutletTubeDiameter/100

# Outlet tube inner diameter (m)
Dhi = Dh - 2*TubeThickness/1000

# Outlet tube length (m)
Length = UpperReservoirLiquidSectionLength/100

# ARC reservoir volume (m^3)
ExpanderVolume = UpperReservoirLiquidVolume/100 ** 3

# Coolant outlet tube volume (m^3)
SteelVolume = math.pi * ((Dh/2)**2 - (Dhi/2)**2) * Length

# Heat transfer area between sodium and steel
Area1 = Length * math.pi * (Dhi/2) * 2

# Heat transfer area between steel and potassium
Area2 = Length * math.pi * (Dh/2) * 2

# Sodium
CoolantDensity          = Coolant_CoolantOutlet.density
CoolantConductivity     = Coolant_CoolantOutlet.conductivity
CoolantHeatCapacity     = Coolant_CoolantOutlet.heatcapacity
CoolantDynamicViscosity = Coolant_CoolantOutlet.viscosity

# Potassium
ExpanderDensity      = 1000 * Expander_CoolantOutlet.density
ExpanderHeatCapacity = Expander_CoolantOutlet.heatcapacity

# Steel
SteelDensity      = 1000 * Steel_CoolantOutlet.density
SteelHeatCapacity = Steel_CoolantOutlet.heatcapacity
SteelConductivity = Steel_CoolantOutlet.conductivity

# Derived
Prandtl  = CoolantHeatCapacity * CoolantDynamicViscosity / CoolantConductivity

vlist = [FullFlowCoolantVelocity*RelativeTubeArea, 0.5*FullFlowCoolantVelocity*RelativeTubeArea, 0.25*FullFlowCoolantVelocity*RelativeTubeArea, NatFlowCoolantVelocity*RelativeTubeArea]

# For the loop
tau1inversecollect = []
h1collect = []
h2collect = []
tau2inversecollect = []
tau3inversecollect = []
reynoldscollect    = []

for v in vlist:

	Reynolds = CoolantDensity * v * Dh / CoolantDynamicViscosity
	Peclet   = Reynolds * Prandtl

	############ NUSSELT OF SODIUM #################################

	if Reynolds < 1e6:

		Nusselt = 4.8 + 0.0156 * (Reynolds ** 0.85) * (Prandtl ** 0.93)

	else:

		Nusselt = 4.8 + 0.0156 * (Reynolds ** 0.85) * (Prandtl ** 0.86)	

	############ SODIUM TO STEEL ###################################

	# Heat transfer coefficient between sodium and steel per square meter [UNIT: W/m^2/K]
	h1 = Nusselt * CoolantConductivity / Dh
	
	# Heat effect (watts) transfer per kelvin of temperature difference between sodium and steel [UNIT: W/K]
	HeatTransfer_1 = h1 * Area1
	
	# Heat capacity steel [J/K]
	SteelHeatCap   = SteelHeatCapacity * SteelDensity * SteelVolume
	
	# Time delay between sodium and steel [UNIT: s^-1]
	tau1 = HeatTransfer_1/SteelHeatCap

	############ STEEL TO POTASSIUM ################################

	#  Heat effect (watts) transfer per kelvin of temperature difference between steel and potassium [UNIT: W/K]
	HeatTransfer_2 = 2 * math.pi * SteelConductivity * Length / (math.log(Dh/Dhi))

	# Heat capacity potassium [UNIT: J/K]
	ExpanderHeatCap   = ExpanderHeatCapacity * ExpanderDensity * ExpanderVolume

	# Time delay between steel and potassium  [UNIT: s^-1]
	tau2 = HeatTransfer_2/ExpanderHeatCap

	############ COLLECT VALUES ####################################

	# Approx combined
	tau3 = 1/tau1 + 1/tau2

	# Collect values
	tau1inversecollect.append(1/tau1)
	tau2inversecollect.append(1/tau2)
	tau3inversecollect.append(tau3)
	reynoldscollect.append(Reynolds)
	h1collect.append(h1)
	h2collect.append(HeatTransfer_2/Area2)
