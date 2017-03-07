from x_definitions import * 
from ARC_settings import *

UARC = -100*1e-5

if RunQSRB == "on":

	from x_ABC import UARC

import math

ReactivityWorth        = UARC # in raw nominal unit
ReactivityWorth        = 100 * 1e5 * ReactivityWorth/Beta # In cents
CalcWorth              = 0.0 # Initial value of calculated reactivity worth
RelativeExpansionSpace = 0.0 # Initial value of the relative expansion space
ARCTubesPerAssembly    = 1.0 # Initial value for number of ARC tubes per assembly

# Assembly area
AssemblyArea = math.sqrt(3)/2 * AssemblyHexagonFlatToFlat ** 2

while abs(CalcWorth) < abs(ReactivityWorth):

	# The outer radius of the inner ARC-tube in the core region and above
	InnerTubeOuterRadius = math.sqrt(OuterTubeInnerRadius** 2 *(1-RelativeExpansionSpace))  

	# The inner radius of the inner ARC-tube in the core region and above
	InnerTubeInnerRadius = InnerTubeOuterRadius - InnerTubeThickness

	# The area of the ARC annulus
	CoreExpansionArea = ARCTubesPerAssembly * math.pi * (OuterTubeInnerRadius ** 2  - InnerTubeOuterRadius ** 2)

	# ARC-expansion volume fraction
	VolumeFractionARCExpansion = CoreExpansionArea / AssemblyArea

	# Calculated reactivity worth
	CalcWorth = VolumeFractionARCExpansion*100*dkARCA*100

	# Increase the relative space for expansion
	RelativeExpansionSpace += 0.00005

	if InnerTubeInnerRadius < MinimumInnerTubeRadius:

		InnerTubeInnerRadius   = OuterTubeInnerRadius - InnerTubeThickness
		RelativeExpansionSpace = 0.0
		ARCTubesPerAssembly += 1

# ARC-expansion area
CoreExpansionArea = ARCTubesPerAssembly * math.pi * (OuterTubeInnerRadius ** 2  - InnerTubeOuterRadius ** 2)

# ARC-expansion volume
CoreExpansionVolume = CoreExpansionArea * CoreLength

# ARC-expansion volume fraction
VolumeFractionARCExpansion = CoreExpansionArea / AssemblyArea