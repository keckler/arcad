from x_definitions import * 
from ARC_settings import *
from x_react import CoreExpansionVolume
import math

## UPPER LIQUID RESERVOIR CALCULATIONS

# Liquid expansion parameter (A)
A = 1 / ( (Expander_AtCore.density/Expander_FullActuation.density) -1)

# Temperature reduction when going from upper to lower reservoir compensation parameter (B)
B = Expander_CoolantInlet.density / Expander_MidActuation.density

# The liquid volume required to provide the needed expansion (assuming no liquid is lost)
FirstTerm = A * B * CoreExpansionVolume

# The make-up term for lost liquid
SecondTerm = B * CoreExpansionVolume / 2

# The required total volume of the liquid part of the upper reservoir in order to provide "CoreExpansion" over the determined temperature range [cm^3]
UpperReservoirLiquidVolume = FirstTerm + SecondTerm 

# The total volume for expansion from operation to at core temperature in order to reach the bottom of the active core at "TemperatureRiseToCore" [cm^3]  -- EQ. 14
ExpansionVolumeBelowCore = UpperReservoirLiquidVolume  * (Expander_CoolantOutlet.density/Expander_AtCore.density - 1)/(Expander_CoolantInlet.density/Expander_HalfWayToCore.density)

# The total expansion from operation to full actuation [cm^3]
TotalExpansion = CoreExpansionVolume + ExpansionVolumeBelowCore

# The area of the liquid-filled section of the upper reservoir [cm^2]
UpperReservoirLiquidSectionArea = ((math.sqrt(3)/2) * InnerFlatToFlat ** 2) - math.pi*((CoolantOutletTubeDiameter/2)**2)

# The length of the liquid-filled section of the upper reservoir [cm]
UpperReservoirLiquidSectionLength = UpperReservoirLiquidVolume / UpperReservoirLiquidSectionArea

