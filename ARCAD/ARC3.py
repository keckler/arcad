import math

print("-- Starting ARC-reservoir calculations")

from x_definitions import *

if RunQSRB == "on":

	print("-- Calculating QSRB")
	from x_ABC import *

print("-- Setting required absorber volume")

from x_react import *

print("-- Calculating the upper ARC reservoir geometry")

from x_upperliquidreservoir import *

print("-- Performing calculations in the below-core region")

from x_belowcoreregion import *

print("-- Performing absorber liquid calculations")

from x_absorber import *

print("-- Calculating liquid expansions")

from x_expansions import *

print("-- Calculating gas production and pressure")

from x_gases import *

print("-- Performing flow velocity conversion")

from x_flowconversion import *

print("-- Calculating pressure drop")

from x_pressuredrop import *

print("-- Validating SASSYS")

from x_ssvalid import *

print("-- Estimating time-delay parameters")

from x_statespace import *

print("-- Calculations DONE")

print("")
print("---- ARC ACTUATION -------------------------")
print("Vol. frac. in core: {0:3.3f}".format(VolumeFractionARCExpansion*100) + "%")
print("Actuation worth:   {0:3.1f}".format(VolumeFractionARCExpansion*100*dkARCA*100) + " cents ({0:3.1f}".format(VolumeFractionARCExpansion*100*dkARCA*Beta) + " pcm) ")
print("Actuation start:    {0:3.1f}".format(TemperatureRiseToCore+CoolantOutletTemperature) + " deg. C ({0:3.1f}".format(TemperatureRiseToCore+CoolantOutletTemperature+273) + "K) ")
print("Actuation end:      {0:3.1f}".format(TemperatureRiseToCore+CoolantOutletTemperature+ActuationTemperatureSpan) + " deg. C ({0:3.1f}".format(TemperatureRiseToCore+CoolantOutletTemperature+ActuationTemperatureSpan+273) + "K) ")
print("--------------------------------------------")
print("")
print("---- ARC TUBES -----------------------------")
print("Tubes per assembly:        {0:3.1f}".format(ARCTubesPerAssembly) + " ")
print("Outer tube outer diameter: {0:3.3f}".format(OuterTubeOuterRadius*2) + " cm")
print("Outer tube inner diameter: {0:3.3f}".format(OuterTubeInnerRadius*2) + " cm")
print("Inner tube outer diameter: {0:3.3f}".format(InnerTubeOuterRadius*2) + " cm (in-core)")
print("Inner tube inner diameter: {0:3.3f}".format(InnerTubeInnerRadius*2) + " cm (in-core)")
print("Tube annulus width:        {0:3.3f}".format(OuterTubeInnerRadius-InnerTubeOuterRadius) + " cm (in-core)")
print("Inner tube outer diameter: {0:3.3f}".format(InnerTubeOuterRadius_BC*2) + " cm (below-core)")
print("Inner tube inner diameter: {0:3.3f}".format(InnerTubeInnerRadius_BC*2) + " cm (below-core)")
print("Tube annulus width:        {0:3.3f}".format(OuterTubeInnerRadius-InnerTubeOuterRadius_BC) + " cm (below-core)")

if NWires == 0:

	print("Nr of wires below core:     No wire")

else:

	print("Nr of wires below core:    {0:3.1f}".format(NWires))

if NWires == 0:

	print("Relative wire-pitch:       No wire")

else:

	print("Wire-pitch:                {0:3.1f}".format(HD*InnerTubeOuterRadius_BC*2) + " cm")
print("--------------------------------------------")
print("")
print("---- UPPER RESERVOIR -----------------------")
print("Total length:               {0:3.2f}".format(UpperReservoirLiquidSectionLength+GasReservoirLength) + " cm")
print("For expansion liquid:       {0:3.2f}".format(UpperReservoirLiquidSectionLength) + " cm")
print("For gas:                    {0:3.2f}".format(GasReservoirLength) + " cm")
print("Volume of liquid:           {0:3.1f}".format(UpperReservoirLiquidVolume) + " cm^3")
print("Volume of gas reservoir:    {0:3.1f}".format(GasVolumeReservoir) + " cm^3")
print("Outlet tube outer diameter: {0:3.2f}".format(CoolantOutletTubeDiameter) + " cm")
print("Outlet tube inner diameter: {0:3.2f}".format(CoolantOutletTubeInnerDiameter*100) + " cm")
print("Inner ass. flat-to-flat:    {0:3.2f}".format(InnerFlatToFlat) + " cm")
print("EOL gas pressure at full actuation: {0:3.2f}".format(EOLGasPressure_Actuation/1e6) + " MPa")
print("Full-flow coolant v. through res:   {0:3.2f}".format(FullFlowCoolantVelocity*RelativeTubeArea) + " m/s")
print("--------------------------------------------")
print("")
print("---- LOWER RESERVOIR -----------------------")
print("Total length:              {0:3.2f}".format(LowerReservoirLowerSectionLength + LowerReservoirUpperConnectorLength) + " cm")
print("For liquids:               {0:3.2f}".format(LowerReservoirLowerSectionLength) + " cm")
print("For shielding:             {0:3.2f}".format(LowerReservoirUpperConnectorLength) + " cm")
print("Volume of liquid:          {0:3.1f}".format(LowerReservoirFreeVolume) + " cm^3")
print("Inlet tube outer diameter: {0:3.2f}".format(CoolantInletTubeDiameter) + " cm")
print("Inlet tube inner diameter: {0:3.2f}".format(CoolantInletTubeInnerDiameter) + " cm")
print("Mass of absorber:          {0:3.2f}".format(AbsorberMass) + " g")
print("--------------------------------------------")
print("")

if FlowTauEquation != "on":

	print("---- PRESSURE DROP -------------------------")
	print("Full flow: {0:3.0f}".format(pdtotal[0]) + " Pa ({0:3.2f}".format(vlist[0]/RelativeTubeArea) + " m/s in rod bundle)")
	print("50% flow:  {0:3.0f}".format(pdtotal[1]) + " Pa ({0:3.2f}".format(vlist[1]/RelativeTubeArea) + " m/s in rod bundle)")
	print("25% flow:  {0:3.0f}".format(pdtotal[2]) + " Pa ({0:3.2f}".format(vlist[2]/RelativeTubeArea) + " m/s in rod bundle)")
	print("Nat. circ: {0:3.0f}".format(pdtotal[3]) + " Pa ({0:3.2f}".format(vlist[3]/RelativeTubeArea) + " m/s in rod bundle)")
	print("--------------------------------------------")
	print("")
print("---- TIME DELAY (tau (s)) ------------------")

if FlowTauEquation == "on":

	print("")
	print("The lag be described as:")
	print("")
	print("            {0:3.3f}".format(flowtau[0]) + " * f + {0:3.3f}".format(flowtau[1]) + "")
	print("tau(f) =   -------------------- [s]")
	print("                {0:3.3f}".format(flowtau[2]) + " + f")
	print("")
	print("Where f is the normalized flow rate")
	print("(f=1.0 is full flow)")

else:

	for x in np.arange(len(vlistv)):
	
		print("{:6.1f}".format(vlistv[x]*100) + "% flow: {0:3.3f}".format(taulist[x]) + " s ")

if RunQSRB == "on":

	print("")
	print("---- QSRB-RESULTS --------------------------")
	print("Required ARC worth: ")
	print("- ULOF:  {:6.1f}".format(UARC_ULOF*1e5)  + " pcm ")
	print("- ULOHS: {:6.1f}".format(UARC_ULOHS*1e5) + " pcm ")
	print("- UTOP:  {:6.1f}".format(UARC_UTOPc*1e5)  + " pcm (coolant)")
	print("- UTOP:  {:6.1f}".format(UARC_UTOPf*1e5)  + " pcm (fuel)")
	print("")
	print(">>" + WhatLimits + "<< determines ARC-worth")
	print("")
	print("Coolant outlet temperatures (QS): ")

	X1 = ""
	X2 = ""
	X3 = ""

	if ULOF_Tout_NARC > 883:

		X1 = ">Boiling!"

	if ULOHS_Tout_NARC > 883:

		X2 = ">Boiling!"

	if UTOP_Tout_NARC > 883:

		X3 = ">Boiling!"

	print("__________w._ARC_______|___w/o._ARC___")
	print("- ULOF:  {:6.1f}".format(ULOF_Tout)   + " deg. C |  {:6.1f}".format(ULOF_Tout_NARC)    + " deg.c " + X1 + "")
	print("- ULOHS: {:6.1f}".format(ULOHS_Tout)  + " deg. C |  {:6.1f}".format(ULOHS_Tout_NARC)   + " deg.c " + X2 + "")
	print("- UTOP:  {:6.1f}".format(UTOP_Tout)   + " deg. C |  {:6.1f}".format(UTOP_Tout_NARC)    + " deg.c " + X3 + "")
	print("  fuel:  {:6.1f}".format(UTOP_TFuel)  + " deg. C |  {:6.1f}".format(UTOP_TFuel_NARC)   + " deg.c")
	print("")
	print("Quasi-static ARC reactivity: ")
	print("- ULOF:  {:6.1f}".format(ULOF_Reactivity*1e5)  + " pcm ({:3.0f}".format(ULOF_Relworth*100)  + "% worth / {:3.0f}".format(ULOF_Relact*100)  + "% dist.)")
	print("- ULOHS: {:6.1f}".format(ULOHS_Reactivity*1e5) + " pcm ({:3.0f}".format(ULOHS_Relworth*100) + "% worth / {:3.0f}".format(ULOHS_Relact*100) + "% dist.)")
	print("- UTOP:  {:6.1f}".format(UTOP_Reactivity*1e5)  + " pcm ({:3.0f}".format(UTOP_Relworth*100)  + "% worth / {:3.0f}".format(UTOP_Relact*100)  + "% dist.)")
	print("--------------------------------------------")

#print("--------------------------------------------------------")
#print("")
#print("---- HEAT TRANSFER -------------------------------------")
#print("Full flow: {0:3.1f}".format(h1collect[0]/1e3) + " kW/(m^2*K) [COT] + {0:3.1f}".format(h2collect[0]/1e3) + " kW/(m^2*K) [UR2]")
#print("50% flow:  {0:3.1f}".format(h1collect[1]/1e3) + " kW/(m^2*K) [COT] + {0:3.1f}".format(h2collect[1]/1e3) + " kW/(m^2*K) [UR2]")
#print("25% flow:  {0:3.2f}".format(h1collect[2]/1e3) + " kW/(m^2*K) [COT] + {0:3.1f}".format(h2collect[2]/1e3) + " kW/(m^2*K) [UR2]")
#print("Nat. circ: {0:3.2f}".format(h1collect[3]/1e3) + " kW/(m^2*K) [COT] + {0:3.1f}".format(h2collect[3]/1e3) + " kW/(m^2*K) [UR2]")
#print("--------------------------------------------------------")
print("")

