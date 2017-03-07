import math
import os

print("-- Reading ARC material properties database")

from ARC_properties import *

print("-- Reading ARC user settings")

from ARC_settings import *

print("-- Performing pre-calculations")

# Make output folder
respath = "Output/" + Name
mlab = respath + "/htrans"
if not os.path.exists("Output"): os.makedirs("Output")
if not os.path.exists(respath):  os.makedirs(respath)
if not os.path.exists(mlab):  os.makedirs(mlab)

# Core geometry definition
BelowCoreLength = BelowCoreLength + LowerReservoirUpperConnectorLength # [cm] The length from the bottom of the active core to the mid-section of the lower ARC reservoir

# ARC-rod geometry definition (1/2)
OuterTubeOuterRadius = FuelRodOuterRadius                                               # [cm] Should be identical to the fuel rod outer radius
OuterTubeInnerRadius = OuterTubeOuterRadius - OuterTubeThickness                        # [cm] The inner radius of the outer ARC-tube

# Inlet and outlet tube definitions
CoolantOutletTubeInnerDiameter = CoolantOutletTubeDiameter - 2*TubeThickness/10
CoolantInletTubeInnerDiameter  = CoolantInletTubeDiameter - 2*TubeThickness/10

# Temperature definitions
ActuationTemperatureSpan     = CoolantTemperatureLimit - (TemperatureRiseToCore + CoolantOutletTemperature) # [deg. C] The temperature rise in the upper ARC reservoir for the absorber to travels across the core axially (change in temperature, deg. C)
RoomTemperature              = 20  # [deg. C] The definition of room-temperature (absolute temperature in Celcius)
CoolantAverageTemperature    = (CoolantInletTemperature   + CoolantOutletTemperature)/2   # [deg. C] The average coolant temperature in the core (absolute value, deg. C)
AbsorberAtCoreTemperature    = CoolantOutletTemperature   + TemperatureRiseToCore         # [deg. C] The temperature in the upper reservoir when the absorber reaches the core (absolute value, deg. C)
HalfWayToCoreTemperature     = (AbsorberAtCoreTemperature + CoolantOutletTemperature)/2   # [deg. C] The temperature in the upper reservoir when the absorber is half-way to the active  core (absolute value, deg. C)
AverageAtCoreTemperature     = (AbsorberAtCoreTemperature + CoolantInletTemperature)/2    # [deg. C] The average coolant temperature when absorber is at core (absolute value, dec. C)
FullActuationTemperature     = AbsorberAtCoreTemperature  + ActuationTemperatureSpan      # [deg. C] The temperature in the upper reservoir at full system actuation (absolute value, deg. C)
MidActuationTemperature      = (FullActuationTemperature  + CoolantOutletTemperature)/2   # [deg. C] The temperature at the middle of the actuation (absolute value, deg. C)
AverageActuationTemperature  = (FullActuationTemperature  + CoolantInletTemperature)/2    # [deg. C] The average coolant temperature at full actuation (absolute value, dec. C)

# Assembly area
AssemblyArea = math.sqrt(3)/2 * AssemblyHexagonFlatToFlat ** 2

# Liquid densitites
Expander_Roomtemperature  = ARCliquids(material=Expander, pressure=1, temperature=RoomTemperature)
Expander_CoolantInlet     = ARCliquids(material=Expander, pressure=1, temperature=CoolantInletTemperature)
Expander_CoolantAverage   = ARCliquids(material=Expander, pressure=1, temperature=CoolantAverageTemperature)
Expander_CoolantOutlet    = ARCliquids(material=Expander, pressure=1, temperature=CoolantOutletTemperature)
Expander_HalfWayToCore    = ARCliquids(material=Expander, pressure=1, temperature=HalfWayToCoreTemperature)
Expander_AtCore           = ARCliquids(material=Expander, pressure=1, temperature=AbsorberAtCoreTemperature)
Expander_AverageAtCore    = ARCliquids(material=Expander, pressure=1, temperature=AverageAtCoreTemperature)
Expander_MidActuation     = ARCliquids(material=Expander, pressure=1, temperature=MidActuationTemperature)
Expander_FullActuation    = ARCliquids(material=Expander, pressure=1, temperature=FullActuationTemperature)
Expander_AverageActuation = ARCliquids(material=Expander, pressure=1, temperature=AverageActuationTemperature)
Expander_Maximum          = ARCliquids(material=Expander, pressure=1, temperature=SystemMaximumTemperature)
Expander_MeltingDensity   = 0.828 # [g/cm^3]
 
Absorber_Roomtemperature  = ARCliquids(material=Absorber, pressure=1, temperature=RoomTemperature)
Absorber_CoolantInlet     = ARCliquids(material=Absorber, pressure=1, temperature=CoolantInletTemperature)
Absorber_CoolantAverage   = ARCliquids(material=Absorber, pressure=1, temperature=CoolantAverageTemperature)
Absorber_CoolantOutlet    = ARCliquids(material=Absorber, pressure=1, temperature=CoolantOutletTemperature)
Absorber_AtCore           = ARCliquids(material=Absorber, pressure=1, temperature=AbsorberAtCoreTemperature)
Absorber_AverageAtCore    = ARCliquids(material=Absorber, pressure=1, temperature=AverageAtCoreTemperature)
Absorber_FullActuation    = ARCliquids(material=Absorber, pressure=1, temperature=FullActuationTemperature)
Absorber_AverageActuation = ARCliquids(material=Absorber, pressure=1, temperature=AverageActuationTemperature)
Absorber_Maximum          = ARCliquids(material=Absorber, pressure=1, temperature=SystemMaximumTemperature)

Coolant_CoolantInlet     = ARCliquids(material=Coolant, pressure=1, temperature=CoolantInletTemperature)
Coolant_CoolantAverage   = ARCliquids(material=Coolant, pressure=1, temperature=CoolantAverageTemperature)
Coolant_CoolantOutlet    = ARCliquids(material=Coolant, pressure=1, temperature=CoolantOutletTemperature)

Steel_CoolantInlet     = ARCliquids(material=Tubes, pressure=1, temperature=CoolantInletTemperature)
Steel_CoolantAverage   = ARCliquids(material=Tubes, pressure=1, temperature=CoolantAverageTemperature)
Steel_CoolantOutlet    = ARCliquids(material=Tubes, pressure=1, temperature=CoolantOutletTemperature)