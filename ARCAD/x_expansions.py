from x_definitions import * 
from ARC_settings import *
from x_upperliquidreservoir import *
from x_belowcoreregion import *
from x_absorber import *
from x_react import InnerTubeInnerRadius
import math

# EXPANDER LIQUID ALWAYS IN LOWER RESERVOIR

# The area available for liquid in the lower reservoir (lower part) Eq. 25
AL0 = (math.sqrt(3)/2 * InnerFlatToFlat ** 2 - math.pi * (CoolantInletTubeDiameter/2)**2)

# The volume of expander liquid always present at the bottom of the lower reservoir (eq. 24)
VolumeL0_RoomTemperature = ExpansionLiquidLowerLevel * AL0

# Volume of expansion of the expansion liquid always present at bottom of lower reservoir (VolumeL0_RoomTemperature) between room temperature and operating temperature (eq. 26)
dVolumeL0_RoomTemperature_to_OperatingTemperature = VolumeL0_RoomTemperature * (Expander_Roomtemperature.density/Expander_CoolantInlet.density -1)


# ABSORBER LIQUID

# Volume of expansion of the absorber liquid going from room temperature to operating temperature
dVolumeAbsorber_RoomTemperature_to_OperatingTemperature = AbsorberMass * (1/Absorber_CoolantInlet.density - 1/Absorber_Roomtemperature.density)


# EXPANDER LIQUID IN INNER ARC-TUBE

# Total length of the inner ARC-tube
InnerARCTubeLength = (InitialGuessLowerReservoirLength - ExpansionLiquidLowerLevel) + BelowCoreLength + CoreLength + AboveCoreLength

# Total volume of liquid inside the inner ARC-tube (eq. 26)
InnerARCTubeVolume = math.pi * InnerARCTubeLength * InnerTubeInnerRadius ** 2

# Total mass of liquid in inner ARC tube at room temperature (eq. 27)
InnerARCTubeLiquidMass_RoomTemperature = InnerARCTubeVolume * Expander_Roomtemperature.density

# Volumetrically averaged temperature of liquid in inner ARC tube at operational temperature (eq. 28)
InnerARCTubeLiquidTemperature_Operation = (CoolantInletTemperature * (BelowCoreLength + InitialGuessLowerReservoirLength - ExpansionLiquidLowerLevel) + CoolantAverageTemperature * CoreLength + CoolantOutletTemperature * AboveCoreLength) / InnerARCTubeLength

# Collect density at this temperature
Expander_InnerARCOperation = ARCliquids(material=Expander, pressure=1, temperature=InnerARCTubeLiquidTemperature_Operation)

# Volume expansion of inner ARC tube liquid from room temperature to operating temperature (eq. 29)
dVolumeInnerARCTube_RoomTemperature_To_OperatingTemperature = InnerARCTubeLiquidMass_RoomTemperature * (1/Expander_InnerARCOperation.density - 1/Expander_Roomtemperature.density)


# EXPANDER LIQUID IN UPPER RESERVOIR

# Mass of expander liquid in upper reservoir at room temperature (eq. 30)
LiquidMassUpperReservoir_RoomTemperature = UpperReservoirLiquidVolume * Expander_Roomtemperature.density

# Volume expansion of expander liquid in upper reservoir to lower reservoir when going from room temperature to operating temperature (eq. 31)
dVolumeUpperReservoirLiquid_RoomTemperature_to_OperatingTemperature = LiquidMassUpperReservoir_RoomTemperature * (1/Expander_CoolantOutlet.density - 1/Expander_Roomtemperature.density)



# TOTAL VOLUME EXPANSION, ROOM TEMPERATURE TO OPERATING TEMPERATURE (eq. 32)
TotalExpansion_RoomTemperature_To_OperatingTemperature = dVolumeL0_RoomTemperature_to_OperatingTemperature + dVolumeAbsorber_RoomTemperature_to_OperatingTemperature + dVolumeInnerARCTube_RoomTemperature_To_OperatingTemperature + dVolumeUpperReservoirLiquid_RoomTemperature_to_OperatingTemperature

# Lower reservoir Section 1 available volume (eq. 33)
LowerReservoirFreeVolume = VolumeL0_RoomTemperature + TotalExpansion_RoomTemperature_To_OperatingTemperature

# Axial length of the lower section of the lower reservoir (eq. 34)
LowerReservoirLowerSectionLength = LowerReservoirFreeVolume / AL0
