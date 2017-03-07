import numpy as np
from ARC_settings import *
from ARC_properties import *
from scipy import signal as ss
from scipy.optimize import curve_fit
from x_flowconversion import RelativeTubeArea, BundleFlowArea
from x_upperliquidreservoir import UpperReservoirLiquidSectionLength, UpperReservoirLiquidVolume
from x_definitions import CoolantOutletTubeInnerDiameter, Coolant_CoolantOutlet, mlab

############# INPUT #############################################################

ReservoirVolume                = UpperReservoirLiquidVolume/(100**3)
TubeLength                     = UpperReservoirLiquidSectionLength/100
CoolantOutletTubeInnerDiameter = CoolantOutletTubeInnerDiameter/100
TubeThickness                  = TubeThickness/1000
DuctThickness                  = DuctThickness/100
vlistv                         = FlowList

if FlowTauEquation == "on":

	a = np.arange(101)
	a = a/100
	vlistv = a[101:0:-1]
	vlistv=np.append(vlistv, 0)

Printing = "off"

################################################################################

TotalRegions             = TubeRegions + ReservoirRegions + DuctRegions
ReservoirRegionNR        = TubeRegions + ReservoirRegions
TubeRegionThickness      = TubeThickness/TubeRegions

TubeOuterDiameter        = CoolantOutletTubeInnerDiameter + 2 * TubeThickness
TubeOuterRadius          = TubeOuterDiameter/2

ReservoirOuterRadius     = math.sqrt(TubeLength * math.pi * (ReservoirVolume + TubeLength * math.pi * TubeOuterRadius ** 2)) / TubeLength / math.pi     # meters (5 cm)

UnDuctedArea    = math.sqrt(3)/2 * ((InnerFlatToFlat/100) ** 2)
DuctedArea      = math.sqrt(3)/2 * (((InnerFlatToFlat/100)+2*DuctThickness) ** 2)
DuctArea        = DuctedArea - UnDuctedArea
DuctVolume1     = DuctArea * TubeLength
DuctOuterRadius = ReservoirOuterRadius + DuctThickness

TubeVolume               = TubeLength * math.pi * (TubeOuterRadius ** 2 - (CoolantOutletTubeInnerDiameter/2) ** 2)
DuctVolume               = TubeLength * math.pi * (DuctOuterRadius ** 2 - ReservoirOuterRadius ** 2)

ReservoirThickness       = ReservoirOuterRadius - TubeOuterRadius
ReservoirRegionThickness = ReservoirThickness/ReservoirRegions
DuctRegionThickness      = DuctThickness/DuctRegions

Na_Conductivity  = ARCliquids(material="Na", pressure=1, temperature=T_Na).conductivity
HT9_Conductivity = ARCliquids(material="HT9", pressure=1, temperature=T_HT9).conductivity
K_Conductivity   = ARCliquids(material="K", pressure=1, temperature=T_K).conductivity

HT9_Density = 1000 * ARCliquids(material="HT9", pressure=1, temperature=T_HT9).density
K_Density   = 1000 * ARCliquids(material="K", pressure=1, temperature=T_K).density

CP_HT9 = ARCliquids(material="HT9", pressure=1, temperature=T_HT9).heatcapacity
HT9_Conductivity_Duct = HT9_Conductivity
CP_K   = ARCliquids(material="K", pressure=1, temperature=T_K).heatcapacity

################################################################################

vx = [1.0]
vlist = np.zeros(len(vx))

filenameslisti = []
filenameslistm = []
filenameslistx = []
filenameslisty = []

for x in np.arange(len(vx)):

	vlist[x]          = vx[x]*FullFlowCoolantVelocity*RelativeTubeArea
	filenameslistm.append(mlab + "/m_" + str(int(vlistv[x]*100)) + "_flowSAS.m")

################################################################################

iii = 0

taulist = np.zeros(len(vlist))

#reyfile = open("reynolds.m", 'w')
#n1file = open("nsch.m", 'w')
#n2file = open("nchi.m", 'w')
#
#reyfile.write("Rey=[")
#n1file.write("N1=[")
#n2file.write("N2=[")

for Velocity in vx:

	if FlowTauEquation == "on":

		if iii == 0:

			print("   - 0%  done")

		elif iii == 24:

			print("   - 25% done")

		elif iii == 49:

			print("   - 50% done")

		elif iii == 74:

			print("   - 75% done")

		elif iii == 99:

			print("   - 100% done")

	else:

		print("   - At " + str(vlistv[iii]*100) + "% flow")

	H  = np.zeros([TotalRegions])
	A  = np.zeros([TotalRegions])
	V  = np.zeros([TotalRegions])
	M  = np.zeros([TotalRegions])
	D  = np.zeros([TotalRegions+1])
	
	AS  = np.zeros([TotalRegions,TotalRegions])
	BS  = np.zeros([TotalRegions,TotalRegions])

	
	##### DIAMETER #################################################################
	
	# First diameter (by definition)
	
	D[0] = CoolantOutletTubeInnerDiameter 
	
	# Diameters within tube
	
	for x in np.arange(1,TubeRegions+1,1):
	
		D[x] = D[x-1] + TubeRegionThickness*2
	
	# Diameters within reservoir
	
	for x in np.arange(TubeRegions+1,ReservoirRegionNR+1,1):
	
		D[x] = D[x-1] + ReservoirRegionThickness*2

	# Diameters within duct
	
	for x in np.arange(ReservoirRegionNR+1,TotalRegions+1,1):
	
		D[x] = D[x-1] + DuctRegionThickness*2

	Dmid = np.zeros(TotalRegions)

	for x in np.arange(TotalRegions):

		Dmid[x] = (D[x] + D[x+1])/2

	OuterDiameter = Dmid[TotalRegions-1]
	DMultiPlier = MaxTime/OuterDiameter

	if Printing == "on":

		print("D")
		print(D)
		
	################################################################################
	
	for x in np.arange(0,TotalRegions,1):
	
		A[x] = math.pi * D[x] * TubeLength


	if Printing == "on":

		print("A")
		print(A)
	
	A[0] = 5*A[0]
	
	################################################################################
	
	for x in np.arange(0,TotalRegions,1):
	
		V[x] = TubeLength * math.pi * ( (D[x+1]/2) ** 2 - (D[x]/2) ** 2)

	if Printing == "on":

		print("V")
		print(V)

	################################################################################
	
	TubeMass = 0
	ReservoirMass =0
	DuctMass = 0

	for x in np.arange(0,TubeRegions,1):
	
		M[x] = HT9_Density * V[x]
		TubeMass +=  HT9_Density * V[x]

	for x in np.arange(TubeRegions,ReservoirRegionNR,1):
	
		M[x] = K_Density * V[x]
		ReservoirMass +=  K_Density * V[x]
	
	for x in np.arange(ReservoirRegionNR,TotalRegions,1):
	
		M[x] = HT9_Density * V[x]	
		DuctMass +=  HT9_Density * V[x]

	if Printing == "on":

		print("masses")
		print(TubeMass)
		print(ReservoirMass)
		print(DuctMass)
	
		print("M")	
		print(M)

	###### COOLANT TO TUBE STEEL ###############################################################
	
	CoolantDensity          = ARCliquids(material="Na", pressure=1, temperature=T_Na).density
	CoolantConductivity     = ARCliquids(material="Na", pressure=1, temperature=T_Na).conductivity
	CoolantHeatCapacity     = ARCliquids(material="Na", pressure=1, temperature=T_Na).heatcapacity
	CoolantDynamicViscosity = ARCliquids(material="Na", pressure=1, temperature=T_Na).viscosity
	Prandtl                 = CoolantHeatCapacity * CoolantDynamicViscosity / CoolantConductivity
	Reynolds                = CoolantDensity * Velocity * CoolantOutletTubeInnerDiameter / CoolantDynamicViscosity
	Peclet                  = Reynolds * Prandtl

	############ NUSSELT OF COOLANT #################################

	if NusseltCorrelation == "NS":
	
		Nusselt = 4.8 + 0.0156 * (Reynolds ** 0.85) * (Prandtl ** 0.93)
	
	else:
	
		Nusselt = 4.8 + 0.0156 * (Reynolds ** 0.85) * (Prandtl ** 0.86)	


	#reyfile.write("\n")
	#reyfile.write(str(Reynolds))

	#n1file.write("\n")
	#n1file.write(str(Nusselt))

	#n2file.write("\n")
	#n2file.write(str(NusseltChi))


	h_Na     = Nusselt * Na_Conductivity / CoolantOutletTubeInnerDiameter # W/m**2/K
	H[0]    = A[0] / (1/h_Na + (TubeRegionThickness/HT9_Conductivity))
	
	#############################################################################################
	
	for x in np.arange(1,TubeRegions,1):
	
		H[x] = HT9_Conductivity * A[x] / TubeRegionThickness
	
	# Heat transfer from tube to reservoir
	H[TubeRegions] = A[TubeRegions] /( (TubeRegionThickness/2)/HT9_Conductivity + (ReservoirRegionThickness/2)/K_Conductivity)
	
	for x in np.arange(TubeRegions,ReservoirRegionNR,1):
	
		H[x] = K_Conductivity * A[x] / ReservoirRegionThickness

	# Heat transfer from reservoir to duct
	H[ReservoirRegionNR] = A[ReservoirRegionNR] /( (DuctThickness/2)/HT9_Conductivity_Duct + (ReservoirRegionThickness/2)/K_Conductivity)

	for x in np.arange(ReservoirRegionNR,TotalRegions,1):
	
		H[x] = HT9_Conductivity_Duct * A[x] / DuctRegionThickness

	if Printing == "on":

		print("H")
		print(H)	
	
	# First entry in AS state-space matrix (set)
	
	AS[0,0] = -(H[0]+H[1])/CP_HT9/M[0]
	AS[0,1] = H[1]/CP_HT9/M[0]

	for x in np.arange(1,TotalRegions-1,1):
	
		CP = CP_K

		if x < TubeRegions:
	
			CP = CP_HT9
	
		if x > ReservoirRegionNR-1:
	
			CP = CP_HT9

		AS[x,x-1] = H[x]/CP/M[x]
		AS[x,x] = -(H[x]+H[x+1])/M[x]/CP
		AS[x,x+1] = H[x+1]/M[x]/CP
	
	# Last
	AS[TotalRegions-1,TotalRegions-2] = H[TotalRegions-1]/CP_HT9/M[TotalRegions-1]
	AS[TotalRegions-1,TotalRegions-1] = -H[TotalRegions-1]/CP_HT9/M[TotalRegions-1]
	
	if Printing == "on":	
	
		print(AS)
	
	BS[0,0] = H[0]/CP_HT9/M[0]

	from x_convtoss import *

	Step     = 0.5
	LastTime = TimeX[-1]
	T = np.arange(0,LastTime,Step)
	
	TimeTemp = {k:v for k,v in zip(TimeX, TempX)}
	u1 = []
	
	for t in T:

		for k,v in TimeTemp.items():
	
			if k == t:
	
				u1.append(v)
		
	x=1
	
	U1 = np.column_stack((u1,u1))

	while x < TotalRegions-1:
	
		U1 = np.column_stack((U1,u1))
				
		x += 1
	
	U = U1

	C = np.identity(TotalRegions)
	X0 = TempX[0] * np.ones(TotalRegions)
	DS = np.zeros([TotalRegions,TotalRegions])
	
	sys1 = ss.lti(AS,BS,C,DS)	
	t, yout, xout = ss.lsim(sys1,U,T,X0)

	## VOLUME AVERAGED TEMPERATURE
	
	ReservoirVolumeAveragedTemperature = np.arange(0,LastTime,Step)
	TubeVolumeAveragedTemperature = np.arange(0,LastTime,Step)
	DuctVolumeAveragedTemperature = np.arange(0,LastTime,Step)
	
	for T in np.arange(0,LastTime,Step):
	
		T = int(T)	
		TempVol = 0.0
	
		#print("----------------------------------")
	
		for x in np.arange(0,TubeRegions,1):
	
			# Temperature multiplied by volume
			TempVol += yout[T,x] * V[x]
	
	
		TubeVolumeAveragedTemperature[T] = (TempVol/TubeVolume)
	
		TempVol = 0.0
	
		#print("----------------------------------")
	
		for x in np.arange(TubeRegions,ReservoirRegionNR,1):
	
			# Temperature multiplied by volume
			TempVol += yout[T,x] * V[x]

		ReservoirVolumeAveragedTemperature[T] = (TempVol/ReservoirVolume)

		TempVol = 0.0
	
		#print("----------------------------------")
	
		for x in np.arange(ReservoirRegionNR,TotalRegions,1):
	
			# Temperature multiplied by volume
	
			TempVol += yout[T,x] * V[x]
	
		DuctVolumeAveragedTemperature[T] = (TempVol/DuctVolume)
	
	if MatlabPlotFile == "on":
	
		mfile = open(filenameslistm[iii], 'w')
		
		mfile.write("TimeX=[")
		
		for Tt in TimeX:
		
			mfile.write(str(Tt))
			mfile.write("\n")
		
		mfile.write("];\n")

		mfile.write("TempSAS=[")
		
		for Tt in TempSAS:
		
			mfile.write(str(Tt))
			mfile.write("\n")
		
		mfile.write("];\n")

		mfile.write("Time=[")
		
		for Tt in t:
		
			mfile.write(str(Tt))
			mfile.write("\n")
		
		mfile.write("];\n")
		
		mfile.write("Tout=[")
		
		for Tt in U[:,0]:
		
			mfile.write(str(Tt-273.0))
			mfile.write("\n")
		
		mfile.write("];\n")
		
		mfile.write("T_Reservoir=[")
		
		for Tt in ReservoirVolumeAveragedTemperature:
		
			mfile.write(str(Tt-273.0))
			mfile.write("\n")
		
		mfile.write("];\n")
	
		mfile.write("T_Tube=[")
		
		for Tt in TubeVolumeAveragedTemperature:
		
			mfile.write(str(Tt-273.0))
			mfile.write("\n")
		
		mfile.write("];\n")
	
		mfile.write("T_Duct=[")
		
		for Tt in DuctVolumeAveragedTemperature:
		
			mfile.write(str(Tt-273.0))
			mfile.write("\n")
		
		mfile.write("];\n")

		mfile.write("RadLocation=[\n")

		for Dia in Dmid:
		
			mfile.write(str(Dia*50))
			mfile.write("\n")
		
		mfile.write("];\n")

		mfile.write("RadYLocation=[0\n")
		mfile.write(str(CoolantOutletTubeInnerDiameter*100) + "\n")

		for Dia in Dmid:
		
			mfile.write(str(Dia*100))
			mfile.write("\n")
		
		mfile.write("];\n")
		mfile.write("figure1 = figure;\n")
		mfile.write("axes1 = axes('Parent',figure1);\n")
		mfile.write("hold(axes1,'on');\n")
		mfile.write("xlim(axes1,[-1 " + str(MaxTime+1) + "]);\n")
		mfile.write("ylim(axes1,[" + str(T_HT9-273-10) + " " + str(T_Na-273+10) + "]);\n")
		mfile.write("plot1 = plot(Time,Tout,'b',Time,T_Tube,'k',Time,T_Reservoir,'r',Time,T_Duct,'m','LineWidth',3,'Parent',axes1);\n")
		mfile.write("box(axes1,'on');\n")
		mfile.write("set(axes1,'FontSize',24,'GridLineStyle','--','XGrid','on','YGrid','on');\n")
		mfile.write("xlabel('Elapsed time (s)'); \n")
		mfile.write("ylabel('Temperature (deg. C)'); \n")
		mfile.write("legend1 = legend('Coolant','Tube','Reservoir','Duct');\n")
		mfile.write("set(legend1,'Location','southeast');\n")
		mfile.close()

	iii += 1
