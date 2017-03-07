from x_definitions import * 
from ARC_settings import *
from x_upperliquidreservoir import *
from x_react import InnerTubeOuterRadius
import math

### CALCULATIONS OF BELOW CORE REGION

# Clean volume below core at same diameter
CleanBelowCoreVolume = BelowCoreLength * math.pi * (OuterTubeInnerRadius ** 2 - InnerTubeOuterRadius ** 2)

#print(CleanBelowCoreVolume)

# Excess volume in below core region
BelowCoreExcessVolume = CleanBelowCoreVolume - ExpansionVolumeBelowCore

# Set the "excess" volume below the core (which must be removed in some way)
CalcExcessVolume        = BelowCoreExcessVolume

# Define the outer diameter of the inner ARC tube in the below core region as equal to that in the core region (as an initial value)
InnerTubeOuterRadius_BC = InnerTubeOuterRadius

# While there is too much volume in the below core region, adjust (increase) the inner ARC tube outer diameter
while CalcExcessVolume > 0:

	InnerTubeOuterRadius_BC += 1e-5
	CalcExcessVolume         = BelowCoreLength * math.pi * (OuterTubeInnerRadius ** 2 - InnerTubeOuterRadius_BC ** 2) - ExpansionVolumeBelowCore

	#print(CalcExcessVolume)

# Assume that we don't need wire (until we know better)
Wire = 0

# Check if the calculated required diameter of the inner ARC tube is allowed
MaxInnerTubeOuterRadius_BC = OuterTubeInnerRadius - MinimumTubeSeperation

if InnerTubeOuterRadius_BC > MaxInnerTubeOuterRadius_BC:

	# If not allowed, set it to the maximum value allowed
	InnerTubeOuterRadius_BC = MaxInnerTubeOuterRadius_BC

	# Since its not allowed, we need wire
	Wire = 1

# Since the tube diameter is set now, define the inner tube inner diameter
InnerTubeInnerRadius_BC = InnerTubeOuterRadius_BC - InnerTubeThickness

# Define the new below-core volume
BelowCoreFreeVolume = BelowCoreLength * math.pi * (OuterTubeInnerRadius ** 2 - InnerTubeOuterRadius_BC ** 2)

NWires = NumberOfWires
HD = 0

# If wires are needed, do wire-calcalations
if Wire == 1:

	# Define the required volume to be removed by the wires
	WireVolume_Required = BelowCoreFreeVolume - ExpansionVolumeBelowCore 

	# Define diameter of the wire around the inner ARC-tube
	WireDiameter = OuterTubeInnerRadius - InnerTubeOuterRadius_BC
	
	# Set a very high initial value for the pitch of the wires
	HD = 100
	
	# Set initial value for the calculated volume taken up by wires
	WireVolume_Calc = 0
	
	# Loop until the appropriate wire configuration is found
	while WireVolume_Calc < WireVolume_Required:
	
		# Calculate the volume taken up by wires
		WireVolume_Calc = ((NWires * math.pi * WireDiameter ** 2)/4) * math.sqrt((((math.pi ** 2) * (BelowCoreLength ** 2)) / (HD**2)) + (BelowCoreLength ** 2))
		
		# Reduce the relative wire pitch until the wires take up enough space
		HD -= 0.01

		# If too much pitch reductin is required, instead increase the number of wires
		if HD < MinimumRelativeWirePitch:
	
			# Step up number of wires by 1
			NWires += 1

			# Re-set the pitch
			HD = 100
