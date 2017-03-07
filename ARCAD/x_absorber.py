from x_definitions import * 
from ARC_settings import *
from x_upperliquidreservoir import *
from x_belowcoreregion import *
import math

# The volume of absorber liquid at full actuation (eq. 21)
AbsorberVolume_FullActuation = CoreExpansionVolume + UpperReservoirLiquidVolume * (Expander_FullActuation.density/Expander_Maximum.density - 1) + AbsorberVolumeMargin

# The volumetric average temperature of the absorber liquid at full actuation (eq. 22)
AbsorberAverageTemperature_FullActuation = (CoolantAverageTemperature * CoreExpansionVolume + CoolantInletTemperature * (UpperReservoirLiquidVolume * (Expander_FullActuation.density/Expander_Maximum.density- 1)) + AbsorberVolumeMargin) / AbsorberVolume_FullActuation

# Get the data for the absorber at the volumetric average actuation temperature
AbsorberData_VolumetricAverageActuation = ARCliquids(material=Absorber, pressure=1, temperature=AbsorberAverageTemperature_FullActuation)

# The total mass of absorber per assembly (eq. 23)
AbsorberMass = AbsorberData_VolumetricAverageActuation.density * AbsorberVolume_FullActuation