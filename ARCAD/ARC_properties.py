import math

class ARCliquids(object):

	def __init__(self, material, temperature, pressure):

		self.material        = material
		self.temperature     = temperature
		self.pressure        = pressure

		if self.material == "K":

			if self.temperature < 63.5:

				self.density = 0.864 - 2.4162 * 1e-4 * self.temperature

			elif self.temperature >= 63.5:
			
				self.density = 0.8415 - 2.172 * 1e-4 * self.temperature - 2.70 * 1e-8 * self.temperature ** 2 + 4.77 * 1e-12 * self.temperature ** 3

			self.heatcapacity = 838.47 - 0.3672 * self.temperature + 4.5899 * 1e-4 * self.temperature ** 2
			self.conductivity = 56.16 * math.exp(-7.958*1e-4*(self.temperature-273.0))

		if self.material == "Li":

			if self.temperature < 180.7:

				self.density = 0.534

			elif self.temperature >= 180.7:
			
				self.density = (562 - 0.1 * (self.temperature + 273))/1000

			self.conductivity = 35.0 - 0.019 * (self.temperature + 273)
			self.heatcapacity = 4640 - 0.885 + 3.91 * 1e-4 * (self.temperature + 273) ** 2
			self.viscosity    = 0.0312 * math.exp(253/(self.temperature + 273))/((self.temperature + 273) ** 0.7368)

		if self.material == "He":

			self.density      = 4.8 * 1e-4 * self.pressure / (self.temperature + 273) # [g/cm^3]
			self.conductivity = 0.0640 + 3.23 * 1e-4 * (self.temperature + 273) - 3.13 * 1e-8 * (self.temperature + 273) ** 2 # [W/m*K]
			self.heatcapacity = 5190 # [J/kg*K]
			self.viscosity    = 8.33 * 1e-6 + 4.16 * 1e-8 * (self.temperature + 273) - 5.30 * 1e-12 * (self.temperature + 273) ** 2

		if self.material == "Na":

			self.density      = 1014 - 0.235 * (self.temperature + 273)
			self.conductivity = 104 - 0.047 * (self.temperature + 273)
			self.heatcapacity = -3.001*1E-6 * (self.temperature + 273) ** -2 + 1658 - 0.8479 * (self.temperature + 273) + 4.454 * 1E-4 * (self.temperature + 273) ** 2
			self.viscosity    = math.exp(556.835/(self.temperature + 273) - 0.3958 * math.log((self.temperature + 273)) - 6.4406)

		if self.material == "HT9":

			self.density      = 7.874 - 3.23 * 1e-4 * (self.temperature + 273)
			self.conductivity = 17.622 + 2.42 * 1e-2 * (self.temperature+273) - 1.606 * 1e-5 * (self.temperature+273) ** 2

			if self.temperature+273 < 800:
	
				self.heatcapacity = ((self.temperature+273)-500)/6 + 500
	
			else:
	
				self.heatcapacity = 3*((self.temperature+273)-800)/5 + 550
	
			
