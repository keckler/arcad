from ARC_settings import *
from x_flowconversion import *
from x_definitions import * 
from x_expansions import LowerReservoirLowerSectionLength
from x_upperliquidreservoir import UpperReservoirLiquidSectionLength
from x_gases import GasReservoirLength
from x_htransfer import reynoldscollect, vlist
import math

# Inlet tube straight part
vlistin = [FullFlowCoolantVelocity*InRelativeTubeArea, 0.5*FullFlowCoolantVelocity*InRelativeTubeArea, 0.25*FullFlowCoolantVelocity*InRelativeTubeArea, NatFlowCoolantVelocity*InRelativeTubeArea]

inreynoldscollect = []

for v in vlistin:

	InReynolds = Coolant_CoolantInlet.density * v * (CoolantOutletTubeInnerDiameter/100) / Coolant_CoolantInlet.viscosity
	inreynoldscollect.append(InReynolds)

# Inlet tube straight part
i = 0
inlettube_pdlist = []

for reynolds in inreynoldscollect:

	lambdapd = 1 / (1.8 * math.log10(reynolds) - 1.62) ** 2
	inlettube_pdlist.append(lambdapd * (LowerReservoirLowerSectionLength/CoolantInletTubeInnerDiameter)/2 * (Coolant_CoolantInlet.density) * vlistin[i]**2)

	i += 1

# Outlet tube straight part
i = 0
outlettube_pdlist = []

for reynolds in reynoldscollect:

	lambdapd = 1 / (1.8 * math.log10(reynolds) - 1.62) ** 2

	outlettube_pdlist.append(lambdapd * (UpperReservoirLiquidSectionLength/CoolantOutletTubeInnerDiameter)/2 * (Coolant_CoolantOutlet.density) * vlist[i]**2)

	i += 1

# Area ratio in and out
n0 = 1/RelativeTubeArea

# Angle frustrum
a0 = math.atan(0.5 * (InnerFlatToFlat-CoolantOutletTubeInnerDiameter)/(GasReservoirLength))
ar = 0.01745 * a0
# Reynolds at start of contraction
# Inlet tube straight part
vlistcont = [FullFlowCoolantVelocity*ContRelativeArea, 0.5*FullFlowCoolantVelocity*ContRelativeArea, 0.25*FullFlowCoolantVelocity*ContRelativeArea, NatFlowCoolantVelocity*ContRelativeArea]

pdcont = []
pdtotal = []

i = 0
for v in vlistcont:

	ContReynolds = ContReynolds = Coolant_CoolantOutlet.density * v * (CoolantOutletTubeInnerDiameter/100) / Coolant_CoolantOutlet.viscosity
	SmoothContLambda = 1 / (1.8 * math.log10(ContReynolds) - 1.62) ** 2
	xifr = (SmoothContLambda / (8 * math.sin(a0/2))) * (1 - 1/n0)
	xi = (-0.0125 * n0 ** 4 + 0.0224 * n0 ** 3 - 0.00723 * n0 ** 2 + 0.00444 * n0 - 0.00745) * (ar ** 3 - 2 * math.pi * ar ** 2 - 10 * ar) + xifr	

	pdcont.append(xi * Coolant_CoolantOutlet.density/2 * vlist[i] ** 2)

	pdtotal.append(pdcont[i] + outlettube_pdlist[i] + inlettube_pdlist[i])

	i += 1
