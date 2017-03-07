################## Import from core design #############################################################
CoreLength                = 137.16 # [cm] The length of the active core (fuelled region)
FuelRodOuterRadius        = 0.404  # |cm] The outer fuel rod cladding radius
FuelRodPitch              = 1.099 * (FuelRodOuterRadius*2) # [cm] The pitch between fuel rods
InnerFlatToFlat           = 15.300 # [cm] The flat-to-flat distance between two inner duct walls
AssemblyHexagonFlatToFlat = 16.142 # [cm] The flat-to-flat distance of the fuel assembly unit cell (including half the inter-assembly space)
CoolantInletTemperature   = 355.0  # [deg. C] The inlet temperature of the coolant
CoolantOutletTemperature  = 510.0  # [deg. C] The mixed-mean coolant outlet temperature
AboveCoreLength           = 170.82 # [cm] The length from the top of the active core to the bottom of the upper ARC reservoir, 
                                   #      Plenum is 170.82 cm, total pin is 422.28 cm, leaving 114.3 cm for upper or lower structures
BelowCoreLength           = 114.3  # [cm] The length from the bottom of the active core to the bottom of the rods
CladdingThickness         = 0.0635 # [cm] The thickness of the fuel rod cladding
PinsPerAssembly           = 271.0  # [#] The total number of rods (including ARC rods) in the assembly
Beta                      = 330.0  # [pcm] Value of delayed neutrons in pcm
FullFlowCoolantVelocity   = 9.065   # [m/s] The coolant velocity in the fuel bundle at full flow
NatFlowCoolantVelocity    = 0.05 * FullFlowCoolantVelocity  # [m/s] The coolant velocity in the fuel bundle at natural circulation
WirePitch                 = 20.0 * (FuelRodOuterRadius*2)   # [cm] The axial height of one wire rotation in the bundle
DuctThickness             = 0.22  # [cm]
##########################################################################################################

# Name of run
Name = "ABR_MOX_075"
 
# Material selection for ARC
Expander = "K"
Absorber = "Li"
Gas      = "Ar"
Coolant  = "Na"
Tubes    = "HT9"

# Absorber enrichment (fraction 6-Li)
Enrichment = 0.9

# dkARCA
dkARCA = -16.5 # Dollars per % of absorber in core

# Gas pressure limit in system [Pa]
GasPressureLimit = 1.0e6

# Temperature input to ARC
TemperatureRiseToCore    = 50.0  # [deg. C] The temperature rise in the upper ARC reservoir for the absorber to reach the bottom of the active core (change in temperature, deg. C)
SystemMaximumTemperature = 750.0 # [deg. C] The highest absolute temperature for which the system is designed (regarding absorber liquid coverage and gas pressures etc.)

# ARC-tubes input
OuterTubeThickness         = CladdingThickness  # [cm] Recommended to be identical to the fuel rod cladding thickness
InnerTubeThickness         = OuterTubeThickness # [cm] Recommended to be identical to the outer tube thickness
MinimumTubeSeperation      = 0.1                # [cm] Determines the distance between tubes in the below-core region (recommended 1 mm = 0.1 cm)
MinimumInnerTubeRadius     = 0.1                # [cm] Determines the minimum acceptable inner ARC tube inner radius
NumberOfWires              = 0.0			    # [-] The number of spacer wires or fins used in the below-core region (if not needed at all, this number will be adjusted to 0, if more are needed, the value will be increased)
MinimumRelativeWirePitch   = 10                 # [-] The smallest axial distance of one wire rotation (measured in the unit of inner ARC-tube diameters)

# ARC-reservoir input
CoolantOutletTubeDiameter          = InnerFlatToFlat-1.5 # [cm] The _outer_ diameter of the outlet tube! Must be smaller than InnerFlatToFlat
CoolantInletTubeDiameter           = 10.0                # [cm] Must, and should, be (significantly) smaller than InnerFlatToFlat
LowerReservoirUpperConnectorLength = 10.0                # [cm] The length of the upper section of the lower reservoir (not critical, only used for shielding)
TubeThickness                      = 2.0                 # [mm] The thickness of the inlet and outlet tubes

# State-space heat transfer analysis settings
MaxTime              = 50.0                                   # [s] Time total "simulation"
TimeStep             = 0.1                                    # [s] Time between calculations
TubeRegions          = 10                                     # [#] Number of regions inside the coolant outlet tube (for heat transfer)
ReservoirRegions     = 50                                     # [#] Number of regions inside the upper reservoir (for heat transfer)
DuctRegions          = 10                                      # [#] Number of regions inside the duct (for heat transfer)
T_Na                 = CoolantOutletTemperature + 50 + 273.0  # [K] Initial temperature of coolant (put to a value above or below steady-state to see a difference)
T_HT9                = 792.071 #CoolantOutletTemperature + 273.0       # [K] Initial temperature of tube       7.921E+02
T_K                  = CoolantOutletTemperature + 273.0       # [K] Initial temperature of reservoir 
FlowTauEquation      = "off"                                  # [on/off] Make a correlation for flow dependent tau-values
FlowList             = [1.0, 0.5, 0.25, 0.05]                 # [#] Fractional flow at which time-lag is calculated (if FlowTauEquation is off)
NusseltCorrelation   = "CC"                                   # [CC/NS] The type of Nusselt correlation used for outlet tube flow
Step1Time            = 30                                     # [s] Times at which coolant temperature
Step2Time            = 30                                     #     is set at values given below
Step3Time            = 30                                     # ___
TemperatureStep1     = T_Na                                   # [K] Coolant temperature definitions 
TemperatureStep2     = T_Na + 50                              #     at times corresponding to above
TemperatureStep3     = T_Na                                   # 
TemperatureStep4     = T_Na                                   # ___
CalculateTau         = "on"                                   # Estimate the volume-averaged upper reservoir tau
MatlabPlotFile       = "on"                                   # [on/off] Create a files for Matlab plotting of results
PlotAll              = "on"                                   # [on/off] Plot all regions of temperature calculations

# Quasi-static evaluation
RunQSRB                 = "on"         # [on/off]
Beta                    = 330.0        # [pcm] Effective delayed neutron fraction
QSRB_A                  = -0.008919768 # [reactivity (raw)]
QSRB_B                  = -0.002043634 # [reactivity (raw)]
QSRB_C                  = -0.000021120 # [reactivity (raw)]/K

DeltaTFuel              = 519.8       # [K] Radial temperature increase across the fuel rod
CoolantTemperatureLimit = 750.0       # [deg. C] The maximum allowable LONG-TERM coolant outlet temperature
FuelTemperatureLimit    = 2000.0      # [deg. C] The maximum allowable LONG-TERM fuel centerline temperature
ULOHS_Power             = 0.005       # [#] Fraction of nominal power that remains in critical state of ULOHS (= passive heat removal systems)
ULOF_Flow               = 0.03        # [#] Fraction of nominal flow following ULOF (= Natural circulation flow) 
UTOP_External           = 28          # [cents] The reactivity added in UTOP * Beta * pcm / 100 # The external reactivity in cents
FuelTempPeaking         = 1.4         # [#] Peak-to-average fuel temperature in the core

# Special settings

# Absorber liquid microscopic capture cross-section
Li6XS = 3 * 1e-24 # [barns]
Li7XS = 0         # [barns]

# Volumetric average flux in lower ARC reservoir absorber liquid
Flux = 1e13 # [n/cm^2*s]

# Residence time in core [days]
ResidenceTime = 2118

# [-] An additional multiplier (or absolute value) for the total volume of absorber liquid to account for uncertainties and the expansion of liquids in the lower reservoir,
#     to ensure that the system never loses negative reactivity worth all the way up to temperature = SystemMaximumTemperature
AbsorberVolumeMarginMultiplier = 1.00
AbsorberVolumeMargin           = 0

# [cm] The lowermost level of the expansion liquid from the bottom of the lower reservoir (at room temperature, when all is solid)
ExpansionLiquidLowerLevel = 4.0 

# Initial guess for the length of the lower reservoir [cm]
InitialGuessLowerReservoirLength = 15 

# Load pressure multiplier [-]
GasPressureMargin = 1.0
