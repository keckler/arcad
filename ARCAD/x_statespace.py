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

vlist = np.zeros(len(vlistv))

filenameslisti = []
filenameslistm = []
filenameslistx = []
filenameslisty = []

for x in np.arange(len(vlistv)):

	vlist[x]          = vlistv[x]*FullFlowCoolantVelocity*RelativeTubeArea
	filenameslisti.append(mlab + "/i_" + str(int(vlistv[x]*100)) + "_flow.txt")
	filenameslistm.append(mlab + "/m_" + str(int(vlistv[x]*100)) + "_flow.m")
	filenameslistx.append(mlab + "/x_" + str(int(vlistv[x]*100)) + "_flow.m")
	filenameslisty.append(mlab + "/y_" + str(int(vlistv[x]*100)) + "_flow.m")

################################################################################

iii = 0
taulist = np.zeros(len(vlist))

for Velocity in vlist:

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

	ofile = open(filenameslisti[iii], 'w')
	
	ofile.write("Information by region \n")
	ofile.write("\n")
	
	for x in np.arange(0,TotalRegions,1):
	
		Material = "Reservoir (K)"

		if x < TubeRegions:
	
			Material = "Tube (HT9)"
	
		if x > ReservoirRegionNR-1:
	
			Material = "Duct (HT9)"
	
		if x == TubeRegions:
	
			ofile.write("\n")
			ofile.write("############################ \n")
			ofile.write("Tube --> Reservoir \n")
			ofile.write("############################ \n")
			ofile.write("\n")
			ofile.write("\n")

		if x == ReservoirRegionNR:
	
			ofile.write("\n")
			ofile.write("############################ \n")
			ofile.write("Reservoir --> Duct \n")
			ofile.write("############################ \n")
			ofile.write("\n")
			ofile.write("\n")
	
		ofile.write("Region ID: " + str(int(x)) + "\n")
		ofile.write("Material: " + Material + "\n")
		ofile.write("Inner diameter:     {0:.3f}".format(D[x]*100) + " cm \n")
		ofile.write("Outer diameter:     {0:.3f}".format(D[x+1]*100) + " cm\n")
		ofile.write("Radial thickness:   {0:.3f}".format((D[x+1]-D[x])*1000/2) + " mm \n")
		ofile.write("Volume:             {0:.2f}".format(V[x] * 100**3) + " cm^3 \n")
		ofile.write("Mass:               {0:.2f}".format(M[x]*1000) + " g \n")
		ofile.write("Outer area:         {0:.2f}".format(A[x] * 100**2) + " cm^2 \n")
		ofile.write("H-trans rate (in):  {0:.2f}".format(H[x]) + " W/K \n")
	
		if x < TotalRegions-1:
	
			ofile.write("H-trans rate (out): {0:.2f}".format(H[x+1]) + " W/K \n")
	
		else:
	
			ofile.write("H-trans rate (out): 0.0000 W/K \n")
	
		ofile.write("\n")
	
	ofile.close()
	
	
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

	if (Step1Time+Step2Time+Step3Time) > MaxTime:
	
		Step1Time = int(MaxTime/Step1Time)
		Step2Time = int(MaxTime/Step2Time)
		Step3Time = int(MaxTime/Step3Time)
	
	Step4Time = MaxTime-(Step1Time+Step2Time+Step3Time) # Seconds
	
	Step1 = int(Step1Time/TimeStep)
	Step2 = int(Step2Time/TimeStep)
	Step3 = int(Step3Time/TimeStep)
	Step4 = int(Step4Time/TimeStep)
	
	u1 = int(TemperatureStep1) * np.ones(Step1)
	u2 = int(TemperatureStep2) * np.ones(Step2)
	u3 = int(TemperatureStep3) * np.ones(Step3)
	u4 = int(TemperatureStep4) * np.ones(Step4)
	
	U1 = np.column_stack((u1,u1))
	U2 = np.column_stack((u2,u2))
	U3 = np.column_stack((u3,u3))
	U4 = np.column_stack((u4,u4))
	
	x=1
	
	while x < TotalRegions-1:
	
		U1 = np.column_stack((U1,u1))
		U2 = np.column_stack((U2,u2))
		U3 = np.column_stack((U3,u3))
		U4 = np.column_stack((U4,u4))
				
		x += 1
	
	U = np.append(U1,U2, axis=0)
	U = np.append(U,U3,  axis=0)
	U = np.append(U,U4,  axis=0)

	T = np.arange(0,MaxTime,TimeStep)
	C = np.identity(TotalRegions)

	X0 = T_HT9 * np.ones(TotalRegions)
	DS = np.zeros([TotalRegions,TotalRegions])
	
	sys1          = ss.lti(AS,BS,C,DS)
	t, yout, xout = ss.lsim(sys1,U,T,X0)
	
	## VOLUME AVERAGED TEMPERATURE
	
	ReservoirVolumeAveragedTemperature = np.arange(0,MaxTime,TimeStep)
	TubeVolumeAveragedTemperature = np.arange(0,MaxTime,TimeStep)
	DuctVolumeAveragedTemperature = np.arange(0,MaxTime,TimeStep)
	
	for T in np.arange(0,MaxTime/TimeStep,1):
	
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
	
	
	if CalculateTau == "on":
	
		def func(T, tau_K):
		  return  T_Na+(T_K-T_Na)*np.exp(-t/tau_K)
		
		tau_K, pcov = curve_fit(func, t, ReservoirVolumeAveragedTemperature, p0=2.0)

	taulist[iii] = tau_K

	if MatlabPlotFile == "on":
	
		mfile = open(filenameslistm[iii], 'w')
		xfile = open(filenameslistx[iii], 'w')
		yfile = open(filenameslisty[iii], 'w')
		
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

		if PlotAll == "on":
	
			for x in np.arange(0,TotalRegions,1):
			
				mfile.write("\n")
				mfile.write("T" + str(int(x)) + "=[")
				
				for Tt in yout[:,x]:
				
					mfile.write(str(Tt-273.0))
					mfile.write("\n")
				
				mfile.write("];\n")
	
			# X-File for 3D mesh plot

			xfile.write("TempMatrix=[\n")

			for x in np.arange(TotalRegions):

				xfile.write("T" + str(int(x)))

				if x == TotalRegions:
				
					xfile.write("")

				else:

					xfile.write(",")						

			xfile.write("];\n")
			xfile.write("figure1 = figure;\n")
			xfile.write("axes1 = axes('Parent',figure1,'Position',[0.0994449106852214 0.11 0.678471755981446 0.840889444694245]);\n")
			xfile.write("hold(axes1,'on');\n")
			xfile.write("xlabel('Radial position (mm)','FontWeight','bold');\n")
			xfile.write("ylabel('Time (s)','FontWeight','bold');\n")
			xfile.write("zlabel('Temperature (deg. C)','FontWeight','bold');\n")
			xfile.write("mesh(RadLocation*10,Time,TempMatrix,'Parent',axes1,'Tag','meshz');\n")
			xfile.write("view(axes1,[-70.8 40.4]);\n")
			xfile.write("grid(axes1,'on');\n")
			xfile.write("colorbar('peer',axes1,'Position',[0.881840277777778 0.11 0.0138888888888888 0.815]);\n")

			# Y-File for 2D circle plot

			yfile.write("CircMatrix=[")

			for TimeY in np.arange(1,MaxTime/TimeStep,1):

				TimeY = int(TimeY)

				yfile.write(str(U[TimeY,0]-273.0) + ", " + str(U[TimeY,0]-273.0) + ",") 
	
				for region in np.arange(TotalRegions):
	
					yfile.write("T" + str(int(region)) + "(" + str(int(TimeY)) + ")")
	
					if region == TotalRegions-1:
					
						yfile.write(";\n")

					else:
	
						yfile.write(",")
	
			yfile.write("];\n")
			yfile.write("CircMatrix=CircMatrix';\n")

		mfile.write("figure1 = figure;\n")
		mfile.write("axes1 = axes('Parent',figure1);\n")
		mfile.write("hold(axes1,'on');\n")
		mfile.write("xlim(axes1,[-1 " + str(MaxTime+1) + "]);\n")
		mfile.write("ylim(axes1,[" + str(T_HT9-273-10) + " " + str(T_Na-273+10) + "]);\n")
		
		mfile.write("plot1 = plot(Time,Tout,'b',Time,T_Tube,'k',Time,T_Reservoir,'r',Time,T_Duct,'m','LineWidth',3,'Parent',axes1);\n")
	
		if PlotAll == "on":
		
			for x in np.arange(0,TubeRegions,1):
			
				mfile.write("hold on; \n")
				mfile.write("plot(Time,T"+str(int(x)) + ",'--k','LineWidth',0.5,'Parent',axes1); \n")
		
			for x in np.arange(TubeRegions,ReservoirRegionNR,1):
			
				mfile.write("hold on; \n")
				mfile.write("plot(Time,T"+str(int(x)) + ",'--r','LineWidth',0.5,'Parent',axes1); \n")

			for x in np.arange(ReservoirRegionNR,TotalRegions,1):
			
				mfile.write("hold on; \n")
				mfile.write("plot(Time,T"+str(int(x)) + ",'--m','LineWidth',0.5,'Parent',axes1); \n")
		
		mfile.write("box(axes1,'on');\n")
		mfile.write("set(axes1,'FontSize',24,'GridLineStyle','--','XGrid','on','YGrid','on');\n")
		mfile.write("xlabel('Elapsed time (s)'); \n")
		mfile.write("ylabel('Temperature (deg. C)'); \n")
		mfile.write("legend1 = legend('Coolant','Tube','Reservoir','Duct');\n")
		mfile.write("set(legend1,'Location','southeast');\n")
		
		mfile.close()
		xfile.close()
		yfile.close()

	iii += 1

if FlowTauEquation == "on":

	def func(v, f1, f2, f3):
	  return  (f1 * v + f2)/(v + f3)
	
	flowtau, pcov = curve_fit(func, vlistv, taulist, p0=[4.0, 90.0, 11.0])

