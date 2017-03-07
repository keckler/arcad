from x_definitions import * 
from ARC_settings import *
from x_upperliquidreservoir import *
from x_belowcoreregion import *
from x_absorber import *
from x_expansions import *
from x_react import CoreExpansionArea
import math

EOLGasPressure_Actuation = 100e6
GasReservoirLength = 1e-3

# Avogrados number
Avogadro = 6.022141 * 1e23

# Gas constant
GasConstant = 8.3145

# Gas volume at room temperature

while EOLGasPressure_Actuation > GasPressureLimit:

	# Gas volume in lower reservoir at room temperature (equal to the expansion of liquids in to the lower reservoir going to operating)
	GasVolumeLowerReservoir_RoomTemperature = TotalExpansion_RoomTemperature_To_OperatingTemperature
	
	# Gas volume in ARC tube room temperature
	GasVolumeARCTubes_RoomTemperature = CoreExpansionArea * (AboveCoreLength + CoreLength) + ExpansionVolumeBelowCore
	
	# Gas volume in ARC tube full-actuation temperature
	GasVolumeARCTubes_ActuationTemperature = CoreExpansionArea * AboveCoreLength
	
	# Gas volume in gas reservoir
	# Without tube
	GasVolumeExTube = GasReservoirLength * (math.sqrt(3)/2) * InnerFlatToFlat ** 2
	
	# The tube removal
	TubeRemoval = (math.pi * GasReservoirLength / 3) * ( ((CoolantOutletTubeDiameter/2)**2) + (InnerFlatToFlat/2)*(CoolantOutletTubeDiameter/2) + ((InnerFlatToFlat/2)**2))
	
	# Net volume for gas
	GasVolumeReservoir = GasVolumeExTube - TubeRemoval
	
	# Total volume of gas room temperature
	GasVolumeRoomTemperature = GasVolumeLowerReservoir_RoomTemperature + GasVolumeARCTubes_RoomTemperature + GasVolumeReservoir
	
	# Total volume of gas actuation temperature
	GasVolumeActuationTemperature = GasVolumeARCTubes_ActuationTemperature + GasVolumeReservoir
	
	
	## FILL GAS PRESSURE
	
	# Hydrostatic head at room-temperature conditions [cm]
	HydrostaticHead = (LowerReservoirLowerSectionLength - ExpansionLiquidLowerLevel) + BelowCoreLength + CoreLength + AboveCoreLength + GasReservoirLength + UpperReservoirLiquidSectionLength
	
	# Gas fill pressure [Pa]
	GasFillPressure = GasPressureMargin * (HydrostaticHead/100 * 9.82 * Expander_MeltingDensity*1000)
	
	# Initial moles of gas
	InitialMoles = GasFillPressure * (GasVolumeRoomTemperature/100**3) / (GasConstant * (RoomTemperature+273)) 
	
	
	## GAS PRODUCTION
	
	# Absorber liquid atomic density
	AbsorberAtomicDensity = Absorber_CoolantInlet.density * Avogadro / (Enrichment * 6.015 + (1-Enrichment) * 7.016)
	
	# Gas atoms produced during operation
	GasAtomsProduced = AbsorberAtomicDensity * Flux * (Enrichment * Li6XS + (1-Enrichment) * Li7XS) * (ResidenceTime * 24 * 60 * 60) * AbsorberMass/Absorber_CoolantInlet.density
	
	# Moles of gas produced
	GasMolesProduced = GasAtomsProduced / Avogadro
	
	# Moles of helium produced
	HeliumMolesProduced = GasMolesProduced / 2
	
	# Moles of T2 gas produced
	DiTritiumMolesProduced = HeliumMolesProduced / 2
	
	# Moles of gas that adds pressure
	PressureAddingGasMoles = HeliumMolesProduced + DiTritiumMolesProduced
	
	# Add the produced gas
	AddedMoles   = PressureAddingGasMoles
	GasMoles     = InitialMoles + AddedMoles
	
	# Baseline gas pressure at end of life at room temperature
	EOLGasPressure = GasMoles * GasConstant * (RoomTemperature + 273) / (GasVolumeRoomTemperature/100**3)
	
	# Gas pressure at EOL at full actuation
	EOLGasPressure_Actuation = GasMoles * GasConstant * (FullActuationTemperature + 273) / (GasVolumeActuationTemperature/100**3)

	GasReservoirLength += 0.001
	
	#print(EOLGasPressure_Actuation/1e6)
	