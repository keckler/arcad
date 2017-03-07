import math
import numpy as np
from ARC_settings import *

################################ INPUT ################

Cent   = 100.0
Tin0   = CoolantInletTemperature
Tout0  = CoolantOutletTemperature
dTc    = Tout0 - Tin0
dTac   = TemperatureRiseToCore
pcm    = 1e-5
dTf    = DeltaTFuel
TLC    = CoolantTemperatureLimit
TLF    = FuelTemperatureLimit
T0ARC  = Tout0 + dTac
T1ARC  = CoolantTemperatureLimit
Pd     = ULOHS_Power
Fn     = ULOF_Flow
rhoExt = UTOP_External * Beta * pcm / Cent # The external reactivity 
pr     = FuelTempPeaking
A      = -0.008919768
B      = -0.002043634
C      = -0.000021120

a0 =      0.4959
a1 =     -0.5025
b1 =      0.0482
w =       3.004

#######################################################

UARC = -0.005311

## Calculate the ARC reactivity in ULOHS

Uarc0 = 0.0
Uarc  = Uarc0
Udiff = 1.0
Tval  = 0.0
relworth = 0.00
z = 0.0

print("   - ULOHS")

while Udiff > pcm:

	Reactivity = 100
	dTin = -0.1

	while Reactivity > pcm:
	
		dTin += 0.1
		Reactivity = (Pd-1) * (A+B) + C * dTin + Uarc
	
	Tval = Tin0 + dTin

	if Tval < T0ARC:
	
		UarcC = 0
	
	if Tval > T1ARC:
	
		UarcC = UARC
	
	if(Tval > T0ARC and Tval < T1ARC):
	
		UpFrom = Tval - T0ARC
		Tot    = T1ARC - T0ARC
	
		z = UpFrom/Tot
	
		relworth = a0 + a1*math.cos(z*w) + b1*math.sin(z*w)
		UarcC = relworth * UARC
	
	Udiff = abs(Uarc-UarcC)
	Uarc -= 1e-6

ULOHS_Tout       = Tval
ULOHS_Reactivity = UarcC
ULOHS_Relworth   = ULOHS_Reactivity/UARC
ULOHS_Relact     = z

print(ULOHS_Tout)

## Calculate the ARC reactivity in ULOF

print("   - ULOF")

Uarc0 = 0.0
Uarc  = Uarc0
Udiff = 1.0
Tval  = Tout0
relworth = 0.0
z = 0.0

while Udiff > pcm:

	Reactivity = 100
	dTout = -0.1

	while Reactivity > pcm:
	
		dTout += 0.1
		Reactivity = (Fn * dTout / dTc + Fn - 1) * A + (dTout/dTc) * B + Uarc

	Tval = Tout0 + dTout

	if Tval < T0ARC:
	
		UarcC = 0
	
	if Tval > T1ARC:
	
		UarcC = UARC
	
	if(Tval > T0ARC and Tval < T1ARC):
	
		UpFrom = Tval - T0ARC
		Tot    = T1ARC - T0ARC
	
		z = UpFrom/Tot
	
		relworth = a0 + a1*math.cos(z*w) + b1*math.sin(z*w)
		UarcC = relworth * UARC

	Udiff = abs(Uarc-UarcC)
	Uarc -= 1e-6

ULOF_Tout       = Tval
ULOF_Reactivity = UarcC
ULOF_Relworth   = ULOF_Reactivity/UARC
ULOF_Relact     = z

print(ULOF_Tout)

## Calculate the ARC reactivity in UTOP

print("   - UTOP")

Uarc0 = 0.0
Uarc  = Uarc0
Udiff = 1.0
Tval  = Tout0
relworth = 0.0
z = 0.0

while Udiff > pcm:

	Reactivity = 100
	dTout = -0.1

	while Reactivity > pcm:
	
		dTout += 0.1
		Reactivity = (dTout/dTc) * A + (dTout/dTc) * B + rhoExt + Uarc

	Tval = Tout0 + dTout

	if Tval < T0ARC:
	
		UarcC = 0
	
	if Tval > T1ARC:
	
		UarcC = UARC
	
	if(Tval > T0ARC and Tval < T1ARC):
	
		UpFrom = Tval - T0ARC
		Tot    = T1ARC - T0ARC
	
		z = UpFrom/Tot
	
		relworth = a0 + a1*math.cos(z*w) + b1*math.sin(z*w)
		UarcC = relworth * UARC

	Udiff = abs(Uarc-UarcC)
	Uarc -= 1e-6

UTOP_Tout       = Tval
UTOP_TFuel      = Tin0 + (UTOP_Tout - Tin0)/2 + dTf
UTOP_Reactivity = UarcC
UTOP_Relworth   = UTOP_Reactivity/UARC
UTOP_Relact     = z

## Without ARC system results
ULOF_Tout_NARC  = Tout0 + ((A+B)/(A*Fn + B) - 1)*dTc
ULOHS_Tout_NARC = Tin0 + (1-Pd)*(A+B)/C
UTOP_Tout_NARC  = Tin0 + (1 - rhoExt/(A+B))*dTc
UTOP_TFuel_NARC = Tin0 + (UTOP_Tout_NARC - Tin0)/2 + dTf
