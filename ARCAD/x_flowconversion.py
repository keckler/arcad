from ARC_settings import *
from x_definitions import * 
import math

# Calculate the flow area in the assembly
BarePinArea   = PinsPerAssembly * math.pi * FuelRodOuterRadius ** 2

if WirePitch > 0:

	# Calculate the area of wire in the assembly
	wirewidth = FuelRodPitch - FuelRodOuterRadius*2
	costheta = WirePitch / math.sqrt( WirePitch ** 2 + (math.pi*(FuelRodOuterRadius*2 + wirewidth))**2 )
	WireArea = PinsPerAssembly * (math.pi * wirewidth ** 2) / 8 * costheta

else:

	WireArea = 0

# The pin + wire area
TotalPinArea = BarePinArea + WireArea

TotalInnerArea = math.sqrt(3)/2 * InnerFlatToFlat ** 2
BundleFlowArea = TotalInnerArea - TotalPinArea

# Flow area in the outlet tube
TubeFlowArea = math.pi * (CoolantOutletTubeInnerDiameter/2) ** 2

# Flow area in the outlet tube
InTubeFlowArea = math.pi * (CoolantInletTubeInnerDiameter/2) ** 2

# Contraction part


# Relative area
RelativeTubeArea   =  BundleFlowArea / TubeFlowArea
InRelativeTubeArea =  BundleFlowArea / InTubeFlowArea
ContRelativeArea   =  BundleFlowArea / TotalInnerArea
#print(BundleFlowArea/PinsPerAssembly)

