#import libraries
from django.shortcuts import render
from numpy import arange, mean
import fiona
from shapely.geometry import mapping, shape, Point, Polygon, LineString
import json, requests
import math
from math import radians
from datetime import datetime
from pyproj import Proj, transform


#offset from GMT (hrs)
offsetGMT = 1

## CANAL PARAMETERS (constants)
# Canal Depth - this will be kept constant
CanalD = 2.1
# Canal section length (m) - kept constant
CanalL = 1
# Bank Height - kept constant
BankH = 0.5

# Latitude and Longitude (Degrees)
lat = 53.478
lng = -2.24


# Total canal length (km)
total_CanalL = 2

#section depth - split depth into 9 sections
section = CanalD/9
sectionD = arange(0, 2, 0.2)



# vertical extinction coefficient
VEC = 3

# relative light intensity from the shade
shade = 0.2

# Gravity (N/kg)
gravity = 9.807

# Latent Heat of Vapourisation (J/g)
L = 2454.9

# Gas Constant for Water Vapour
Rv = 461.5
	
# L/Rv
LRv = (L*1000)/Rv

# Refractive index of water
riW = 1.33

# refractive index of air
riA = 1

# overhead sun intensity (W/m^2)
sun = 476

# Emissivity of water
e = 0.95

# Stefan-Boltzmann Constant
SB = 5.67*10**(-8)

# Triple point of water (K)
To = 273.0 

# Specific heat capacity of water
capW_kJ = 4.186     # kJ/kg K
capW = 4186         # J/kg K

# heat condicutivity of water (W/m K)
condW = 0.58

# Thermal conductivity (kc)
conductivity = 0.58

#Convection coefficent (hc)
convection = 200

# Baseline Pressure at Triple Point (kPa)
Eo = 0.611

# constants to calculate dew point (in evaporation)
b = 17.67
c = 243.5

# Specific weight of water (g/m^3)
p = 999286

        
'''
These functions are used to calculate the suns positioning at different times of the day
i.e. the sun altitude and azimuth
in both radians and degrees
'''
def HRA(DOY, hr, mins):

	#latitude (rad)
	latR = lat*math.pi/180

	#longitude (rad)
	lngR = lng*math.pi/180

	# local standard time meridian (deg)
	LSTM = 15*offsetGMT

	# B (rad)
	B = (360.0/365.0)*(DOY-81)*math.pi/180
	

	# Equation of time (minutes)
	# Equation of Time Corrects for Eccentricity of Earth's Orbit and Axial Tilt
	EoT = 9.87*math.sin(2*B)-7.53*math.cos(B)-1.5*math.sin(B)

	# Declination angle (rad)
	# Declination Angle is that between the Equator and the line from Earth to the Sun created due to the Earth's Tilt
	dA = 23.45*math.sin((360.0/365.0)*(DOY-81)*math.pi/180)*math.pi/180

	# Time correction factor (min)
	# Accounts for Variations of Local Solar Time within the time zone
	TC = 4*(lng-LSTM)+EoT

	# local time (hours)
	LT = hr+(mins/60.0)
  

	# local solar time (hrs)
	LST = LT + TC/60.0

	# hour angle (rad)
	HRA = 15*(LST-12)*math.pi/180

	return HRA

def SunAltitude_R(DOY, hr, mins):
	#latitude (rad)
	latR = lat*math.pi/180

	#longitude (rad)
	lngR = lng*math.pi/180

	# Declination angle (rad)
	# Declination Angle is that between the Equator and the line from Earth to the Sun created due to the Earth's Tilt
	dA = 23.45*math.sin((360.0/365.0)*(DOY-81)*math.pi/180)*math.pi/180            

	# Altitude
	Alt_r = math.asin(math.sin(dA)*math.sin(latR)+math.cos(dA)*math.cos(latR)*math.cos(HRA(DOY, hr, mins)))  #radians

	return Alt_r


def SunAltitude_D(DOY, hr, mins):
	Alt_d = SunAltitude_R(DOY, hr, mins)*180/math.pi                                        #degrees
	return Alt_d


def SunAzimuth_R(DOY, hr, mins):
	#latitude (rad)
	latR = lat*math.pi/180

	#longitude (rad)
	lngR = lng*math.pi/180

	# Declination angle (rad)
	# Declination Angle is that between the Equator and the line from Earth to the Sun created due to the Earth's Tilt
	dA = 23.45*math.sin((360.0/365.0)*(DOY-81)*math.pi/180)*math.pi/180

	#Azimuth (rad)
	Az_r = math.acos((math.sin(dA)*math.cos(latR)-math.cos(dA)*math.sin(latR)*math.cos(HRA(DOY, hr, mins)))/math.cos(SunAltitude_R(DOY, hr, mins)))
	
	#Azimuth correction (rad)
	#Accounts for artifact of calculation whereby original only correct for Solar morning
	if HRA(DOY, hr, mins)>0:
		Azimuth_r = (2*math.pi)-Az_r
	elif HRA(DOY, hr, mins)<=0:
		Azimuth_r = Az_r

	return Azimuth_r

def SunAzimuth_D(DOY, hr, mins):
	#azimuth (deg)
	Az_d = SunAzimuth_R(DOY, hr, mins)*180/math.pi

	#azimuth correction (deg)
	if HRA(DOY, hr, mins)>0:
		Azimuth_d = 360-Az_d
	elif HRA(DOY, hr, mins)<=0:
		Azimuth_d = Az_d

	return Azimuth_d

'''
The following functions are used to calculate the solar radiation
outputs are absorbed and relfected radiation due to sunlight
'''
def Absorptivity(DOY, hr, mins):

	## Reflectivity Angle Dependence
	## Calculates the Absoptivity of the Water given the direct sunlight angle utilising the Fresnel Equations
	# zenith angle (rad)
	zen = 90 - SunAltitude_D(DOY, hr, mins)
	zen_r = zen*math.pi/180
	
	x = math.sqrt(1-(((riA/riW)*math.sin(zen_r))**2))
	y1 = riA*math.cos(zen_r)-riW*x
	y2 = riA*math.cos(zen_r)+riW*x
	Rs = (y1/y2)**2
	z1 = riA*x-riW*math.cos(zen_r)
	z2 = riA*x+riW*math.cos(zen_r)
	Rp = (z1/z2)**2
	reflectivity = (Rs+Rp)/2
	absorptivity = 1-reflectivity

	return absorptivity

def SunIntensity(DOY, hr, mins):

	## Intensity of the Sun based on Angle (W/m^2)
	#intensity
	intensity = sun*math.cos(math.radians(90 - SunAltitude_D(DOY, hr, mins)))

	return intensity

def AbsorbedLight(DOY, hr, mins):
	#absorped light (W/m^2)
	if SunIntensity(DOY, hr, mins) > 0:
		AbsorbLight = SunIntensity(DOY, hr, mins)*Absorptivity(DOY, hr, mins)
	elif SunIntensity(DOY, hr, mins) <= 0:
		AbsorbLight = 0

	return AbsorbLight

def ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
	
	## shading - Accounts for Shading due to Banking and Buildings
	# Length of Shadow Cast by Building 1
	B1_shadow = (B1_h*math.sin(SunAzimuth_R(DOY, hr, mins) - AngleR)/(math.tan(SunAltitude_R(DOY, hr, mins)))) - B1_d

	# The shaded water width from building1 (m)
	if B1_shadow < 0:
		B1_shadow = 0

	# Length of Shadow Cast by Building 2
	B2_shadow = (-B2_h*math.sin(SunAzimuth_R(DOY, hr, mins) - AngleR)/(math.tan(SunAltitude_R(DOY, hr, mins)))) - B2_d
	
	# The shaded water width from building 2 (m)
	if B2_shadow < 0:
		B2_shadow = 0

	# total building shadow
	if B1_shadow > B2_shadow:
		building_shadow = B1_shadow
	elif B1_shadow <= B2_shadow:
		building_shadow = B2_shadow
		
	# shading from bank
	bankShade = abs(BankH*math.sin(SunAzimuth_R(DOY, hr, mins) - AngleR)/(math.tan(SunAltitude_R(DOY, hr, mins))))
	
	#Total shading width (m)
	if building_shadow > bankShade:
		totalShade = building_shadow
	elif building_shadow < bankShade:
		totalShade = bankShade
		
	# Limiting Shade width to that of canal
	if totalShade > CanalW:
		totalShade = CanalW


	return totalShade

def UnshadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
	
	# Unshaded Canal Width (m)
	return(CanalW-ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR))


def AbsorbedRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
	
	#Absorbed radiation due to sunlight
	absorbedRad = (UnshadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)+
				   ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)*
				   shade)*AbsorbedLight(DOY, hr, mins)

	return absorbedRad

def ReflectedRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
	reflectRad = (UnshadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR) +
				  ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR) *
				  shade)*SunIntensity(DOY, hr, mins)*(1-Absorptivity(DOY, hr, mins))

	return reflectRad

'''
The following are used to calculate change in water temperature at different depths
the canal is split into 0.2m sections
energy gained due to sunlight is calculated
temperature rise and new temperature are calculated from this
'''

#energy gained from sunlight, J
def SunlightEnergy(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
        sun = []
        for d in range(9):
                sun.append((AbsorbedRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)*math.exp
                               (-VEC*(sectionD[d]))-AbsorbedRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)
                               *math.exp(-VEC*(sectionD[d+1])))*(3600/4))
        return sun



def NetSunlightEnergy(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy):
        net = []
        for i in range(9):
                net.append(SunlightEnergy(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)[i] - InitialEnergy[i] + InitialEnergy[i+1])
        return net 


#temp rise
def TempRise(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy):
        temp = []
        for s in NetSunlightEnergy(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy):
                temp.append(s/(capW*(section*CanalL*CanalW)*1000))
        return temp     
        
def NewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy, Tw_C):
        temp = []
        for t, w in zip(TempRise(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy), Tw_C):
                temp.append(w+t)
                
        return temp

def NewTempK(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy, Tw_C):
        temp = []
        for t in NewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, InitialEnergy, Tw_C):
                temp.append(t+273.15)
        return temp

        
'''
Functions to calculate evaporation rate
'''
def Ew(Tw_C):
        
        # Saturated Vapour Pressure at Surface Water Temperature
        Ew_K = math.exp(77.345+0.0057*(Tw_C[0]+273.15)-7235/(Tw_C[0]+273.15))/math.pow((Tw_C[0]+273.15), 8.2)    

        # Saturated Vapour Pressure at Surface Water Temperature (mbar)
        Ew = Ew_K/100

        return Ew

def Ea(Ta, RH):

	# Dew point temperature (oC)
	y = math.log(RH*0.01)+ b*Ta/(c+Ta)
	Tdp_C = c*y/(b-y)

	# Dew point temperature (K)
	Tdp = Tdp_C + 273.15

	# Vapour Pressure at Air Temperature in kelvin

	Ea_K = Eo*(math.exp(LRv*((1.0/To)-(1/Tdp))))
	
	# Vapour Pressure at Air Temperature (mbar)
	Ea = Ea_K*10

	return Ea

def Ev(Ta, U, RH, Tw_C):

	# Evaporation rate (mm/day)
	Ev_mm = 0.165*(0.8+0.864*U)*(Ew(Tw_C)-Ea(Ta, RH))

	# Evaporation rate (m/day)
	Ev = Ev_mm/1000

	return Ev

def EvaporationEnergy(Ta, U, RH, Tw_C):

	# Energy loss to evaporation (J/m^2/day)
	Qe = Ev(Ta, U, RH, Tw_C)*L*p

	return Qe

def EvaporationPower(Ta, U, RH, CanalW, Tw_C):

	# Power Loss to Evaporation per metre length (W/m)
	Pe = EvaporationEnergy(Ta, U, RH, Tw_C)*CanalW*(1.0/86400.0)   

	return Pe

def EvapoativeCooling(Ta, U, RH, CanalW, Tw_C):
	air = EvaporationPower(Ta, U, RH, CanalW, Tw_C)/2
	return air    

# Radiative Emissions From Water
def WaterEmmissions(CanalW, Tw_C):
	# W/m^2
	p = e*SB*((Tw_C[0]+273.15)**4)
	# W/m
	Pwat = p*CanalW
	return Pwat

# Thermal Absorption of Water from Sky
def ThermalAbsorption(CanalW, Ta, RH, CloudH, CloudC):
	#Pth (W/m^2)
	thermalA = (1+ CloudH*CloudC**2)*8.78*(10**-13)*(Ta+273.15)**5.852*RH**0.07195
	#W/m
	pth = CanalW*thermalA
	return pth

def NetThermalEmissions(CanalW, Ta, RH, CloudH, CloudC, Tw_C):
	net = WaterEmmissions(CanalW, Tw_C)-ThermalAbsorption(CanalW, Ta, RH, CloudH, CloudC)
	return net


'''
sensible heat transfer
'''
def SensibleHeat(Ta, U, RH, AirP, CanalW, Tw_C):

        
        # Bowen Ratio
        BR = (0.61*AirP*(Tw_C[0]-Ta))/(Ew(Tw_C)-Ea(Ta, RH))/1000


        # Sensible heat (W/m)
        Qh = EvaporationPower(Ta, U, RH, CanalW, Tw_C)*BR

        return Qh




'''
heat transfer between depths
sections slipt into 0.2m intervals
'''

# conductive/convective effects
def SurfaceHeatTransfer(CanalW, Ta, U, RH, AirP, CloudC, CloudH, Tw_C):
        heat = (EvapoativeCooling(Ta, U, RH, CanalW, Tw_C) + NetThermalEmissions(CanalW, Ta, RH, CloudH, CloudC, Tw_C)
                + SensibleHeat(Ta, U, RH, AirP, CanalW, Tw_C))*(3600/4)
        return heat


def HeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, U, RH, AirP, CloudC, CloudH, InitialEnergy, Tw_C):

        heatTransfer = []
        
        for t in range(8):

                conduct = conductivity*(CanalW*CanalL)*(Tw_C[t+1]-Tw_C[t])/section

                if Tw_C[t+1] > Tw_C[t]:
                        heat = (conduct + convection*(CanalW*CanalL)*(Tw_C[t+1]-Tw_C[t]))*(3600/4)
                        
                elif Tw_C[t+1] <= Tw_C[t]:
                        heat = conduct*(3600/4)
                        
                heatTransfer.append(heat)

        #insert calculated surface heat transfer to beginning of list
        heatTransfer.insert(0, SurfaceHeatTransfer(CanalW, Ta, U, RH, AirP, CloudC, CloudH, Tw_C))
        heatTransfer.append(0)
        
        return heatTransfer

#calculate energy gained/lost form the canal
def EnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C):
        energylost = -EvapoativeCooling(Ta, U, RH, CanalW, Tw_C)-SensibleHeat(Ta, U, RH, AirP, CanalW, Tw_C)-NetThermalEmissions(CanalW, Ta, RH, CloudH, CloudC, Tw_C)
        return(-energylost-AbsorbedRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR))

#from energy gain, temperature change can be calculated from various assumptions
def TempChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C):
        tempchange = EnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C)/852.04

        return tempchange


'''
concrete equivalent
for comparison
'''



##CONCRETE EQUIVELENT
concreteDepth = 0.05 #m, section depth
concreteEm = 0.85 #emissivity
concreteAlbedo = 0.15
concreteAbsorb = 0.85 #absorbivity
concreteConduct = 1.21 #W/(m*K) thermal conductivity
concreteCapacity = 921 #J/(kg*K) specific heat capacity
concreteDensity = 2238 #kg/m3
#ground
groundConduct = 1
groundCapacity = 1900
groundDensity = 1500


#solar radiation W/m
def ConcreteSolarRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR):
	SR = (UnshadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)+
              ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)
              *shade)*SunIntensity(DOY, hr, mins)*concreteAbsorb
	if SR < 0:
		SR = 0

	return SR

	
##thermal radiation
#emissions W/m
def ConcreteEmissions(CanalW, concreteTemp):         
	return(CanalW*concreteEm*SB*(concreteTemp[0]**4))

#Absorbed W/m
def ConcreteAbsorb(CanalW, Ta, RH, CloudH, CloudC):
	return(ThermalAbsorption(CanalW, Ta, RH, CloudH, CloudC)*concreteAbsorb)

#Net radiated W/m
def ConcreteNetRadiated(CanalW, Ta, RH, CloudH, CloudC, concreteTemp):
	return(ConcreteEmissions(CanalW, concreteTemp)-ConcreteAbsorb(CanalW, Ta, RH, CloudH, CloudC))


#temp after sunlight
def ConcreteAfterSunlight(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, concreteTemp):
        temp = []
        for s in range(1, 10):
                temp.append(concreteTemp[s])
        sun = (ConcreteSolarRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)*
               (3600/4))/(concreteCapacity*(CanalW*CanalL*concreteDepth*concreteDensity))+concreteTemp[0]
        temp.insert(0, sun)
        return temp            


#Convection to air  W/m
def ConcreteEnergytoAir(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, U, concreteTemp):
	h = 0.664*0.025*(0.713**0.3)*((15*(10**(-6)))**(-0.5))*(U**0.5)
	energytoAir = h*CanalW*(ConcreteAfterSunlight(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, concreteTemp)[0]-(Ta+273.15))
	if energytoAir < 0:
		energytoAir = 0
	
	return energytoAir

##change in concrete temps
def ConcHeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp):
	return(ConcreteNetRadiated(CanalW, Ta, RH, CloudH, CloudC, concreteTemp) + ConcreteEnergytoAir(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, U, concreteTemp))

def ConcreteHeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp):
        heat = [] 
        for t in range(9):
                temp = CanalW*(concreteTemp[t+1]-concreteTemp[t])/concreteDepth
                if t < 4:
                        transfer = (concreteConduct*temp)
                elif t >= 4:
                        transfer = (groundConduct*temp)
                heat.append(transfer)
                
        #insert surface temperature at the beginning of this list
        heat.insert(0, ConcHeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp))
        
        return heat

#calculate change in energy       
def ConcreteEnergyChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp):
        energy = []
        for e in range(9):
                energy.append((ConcreteHeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp)[e+1]
                              -ConcreteHeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp)[e])*(3600/4))
        energy.append(0)
        return energy

#temperature rise of concrete
def ConcreteTempRise(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp):
        temp = []
        for t in range(10):
                if t <= 4:
                        rise=ConcreteEnergyChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp)[t]/(concreteCapacity*(CanalW*CanalL*concreteDepth*concreteDensity))
                if t > 4:
                        rise=ConcreteEnergyChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp)[t]/(groundCapacity*(CanalW*CanalL*concreteDepth*groundDensity))
                temp.append(rise)
        return temp

#the new temperature of concrete 
def ConcreteNewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, RH, CloudH, CloudC, U, concreteTemp):
        temp = []
        for s, t in zip(ConcreteAfterSunlight(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, concreteTemp), ConcreteTempRise(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, RH, CloudH, CloudC, U, concreteTemp)):
                temp.append(s+t)
        return temp

#calculate energy lost/gained from the concrete equivilent 
def ConcreteEnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, RH, CloudH, CloudC, U, concreteTemp):
        energy = (ConcreteNetRadiated(CanalW, Ta, RH, CloudH, CloudC, concreteTemp) + ConcreteEnergytoAir(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, AngleR, CanalW, Ta, U, concreteTemp)) - ConcreteSolarRadiation(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR)
        return energy

#calculate concrete temperature change fromthe energy gain
def ConcreteTempChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, RH, CloudH, CloudC, U, concreteTemp):
        tempchange = ConcreteEnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, RH, CloudH, CloudC, U, concreteTemp)/852.04
        return tempchange

# Final net energy gain i.e. difference in canal and concrete energy gain
def NetEnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C, concreteTemp):
    netEnergy = EnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C)-ConcreteEnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, CanalW, AngleR, Ta, RH, CloudH, CloudC, U, concreteTemp)
    return netEnergy

# Final net temperaure change - THIS IS THE FINAL FIGURE TO USE IN MODEL
def NetTempChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, concreteTemp, InitialEnergy, Tw_C):
        netTemp = ((NetEnergyGain(B1_h, B1_d, B2_h, B2_d, DOY, hr, mins, Ta, U, RH, AirP, CanalW, CloudH, CloudC, AngleR, InitialEnergy, Tw_C, concreteTemp)*3600/1000)/852.04)/6
        return netTemp

def offsetPoint(point, distance, direction):
        """
        * Offset a point by a given distance and direction along a plane
        """
        return Point(point.x + math.sin(radians(direction)) * distance, point.y + math.cos(math.radians(direction)) * distance)

def wrap(bearing):
        """
        * Force a number to wrap into the range 0-360
        """
        return (bearing + 360) % 360

def createLine(point, length, angle):
        """
        * Create a line of a given length and direction that passes through a given point at its centre
        """
        l = length / 2
        return LineString([offsetPoint(point, l, angle), offsetPoint(point, l, wrap(angle - 180))])

def getWidthAndDirection(point, geom):
        """
        * Calculate the width and direction of a canal at a given point
        """
        minWidth = float('inf')
        minBearing = None
        
        # test all directions at 3 degree intervals
        for i in range(0, 360, 3):
                
                # construct a line
                line = createLine(point, 1000, i).intersection(geom)
                
                # record direction and width if smaller than previous smallest
                if line.length < minWidth:
                        minWidth = line.length
                        minBearing = i

        # calculate direction of canal (perpendicular to width direction)               
        # and force number to wrap into bearing 180-360
        bearing = wrap(minBearing+90)
        if bearing < 180:
                bearing = bearing+180
                        
        return minWidth, bearing
        
def getDirection(currentPoint, nextPoint):
        """
        * Calculate the direction between two points
        """
        return wrap(270 - math.degrees(math.atan2(currentPoint.y - nextPoint.y, currentPoint.x - nextPoint.x)))


def getCanalPoints(path, interval):
        """
        * Convert a lineString into a load of dots at a specified interval
        """
        with fiona.open(path) as line2:

                # get geometry of line
                line = line2[0]['geometry']

                # blank list to append points to
                points = []
                
                # for every coorindinate in the line
                for l in range(len(line['coordinates'])-1):

                        #get coorindinates of current and next points 
                        currentPoint = Point(line['coordinates'][l][0], line['coordinates'][l][1])
                        nextPoint = Point(line['coordinates'][l+1][0], line['coordinates'][l+1][1])

                        # get distance and direction between current and next points                        
                        distance = currentPoint.distance(nextPoint)
                        direction = getDirection(currentPoint, nextPoint)

                        # divide distance by the interval of canal points
                        n = int(distance/interval)

                        # for each point
                        for i in range(n):

                                # offset it by the specified interval and in the direction of the line
                                # append to list
                                points.append(offsetPoint(currentPoint, interval*i, direction))
                
                return points
        
def getDistanceAndHeight(path, point, canalWidth):
    '''
    gets the distance and height of the nearest building
    to a specified point
    '''
    minDist = float('inf')

    # open buildings shapefile
    with fiona.open(path) as buildings:

        # for each building
        for b in buildings:

            # get the geometry
            building = Polygon(shape(b['geometry']))

            # calculate distance from canal point to exterior of each building
            dist = building.exterior.distance(point)

            # if distance measured is smallest
            if (dist < minDist):

                # save min distance
                # minus half the canal width to get distance from edge of the canal
                minDist = dist - canalWidth/2

                # save height of building
                height = b['properties']['relhmax']
                
        return minDist, height

def chunks(l, n):
    #splits list into list of lists
    for i in range(0, len(l), n):
        yield l[i:i+n]
        
# this is the view for the homepage (where the user inputs the data)
def index(request):
    return render(request, 'cooling/index.html', {})

# this is the view for the results page (where the results are calculated)
def results(request):

    '''
    This is where you want to do all of your data processing,
     and then put your results in the context object to be passed
     to the results page.
    '''

    # the number a user puts in on the website will appear here as the interval between points along the canal
    interval = int(request.POST['interval'])

    # get a set of points representing the canal centreline         
    canalPoints = getCanalPoints('cooling/data/canal_line.shp', interval)

    with fiona.open('cooling/data/canal.shp') as canal:                
            # open the canal feature
            geom = shape(canal[0]['geometry'])

    # open file to write canal points to
    with fiona.open('cooling/data/canal_points.shp', 'w', driver='ESRI Shapefile', crs='+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs ', schema={'geometry': 'Point', 'properties': {}} ) as o:

        #empty lists to append results to
        canalWidth = []
        canalDirection = []
        canalDirectionDegrees = []
        B1_dist = []
        B1_height = []
        B2_dist = []
        B2_height = []
        xs = []
        ys = []

        #for each point along the canal
        for point in canalPoints:
                
            # if the point intersects the canal (this is needed as the canal polygon contains gaps):
            if point.intersection(geom):

                # write a feature to it (convert geometry to fiona style in the process)
                o.write({'geometry': mapping(point),'properties': {}})

                # get width and direction of the canal at each point
                canalW, canalD = getWidthAndDirection(point, geom)

                # get northing and easting of each point and append them to their appropriate list (created above)
                xs.append(mapping(point)['coordinates'][0])
                ys.append(mapping(point)['coordinates'][1])
                
                # append the results to lists
                canalDirectionDegrees.append(canalD)
                canalDirection.append(radians(canalD))
                canalWidth.append(canalW)

                #in canal width
                for width in canalWidth:
                        
                        # get distance and height of nearest building on each side of the canal
                        distNorth, heightNorth = getDistanceAndHeight('cooling/data/buildings_north.shp', point, width)
                        distSouth, heightSouth = getDistanceAndHeight('cooling/data/buildings_south.shp', point, width)

                # append results to lists
                B1_dist.append(distNorth)
                B1_height.append(heightNorth)
                B2_dist.append(distSouth)
                B2_height.append(heightSouth)

    # load weather data from openweathermap
    response = requests.get('http://history.openweathermap.org//storage/e0b7d3e03f39a638bd829614a59ca103.json')
    weather = json.loads(response.text)

    # define blank lists for weather data
    airTemp = []
    humidity = []
    wind = []
    cloudCover = []
    airPressure = []

    # extract weather data from jan 2018 0:00
    for w in range(8725, len(weather)-1):
            
        #get temperature data and convert to celcius
        airTemp.append(weather[w]['main']['temp']-273.15)

        #get humidity and append
        humidity.append(weather[w]['main']['humidity'])

        #get air pressure and append
        airPressure.append(weather[w]['main']['pressure'])

        #get wind speed and append
        wind.append(weather[w]['wind']['speed'])

        #get cloud cover, convert to between 0-1, and append
        cloudCover.append((weather[w]['clouds']['all'])/100)

    # fill in missing weather data
    wind.insert(344, 2)
    airTemp.insert(344, 6.5)
    cloudCover.insert(344, 0.8)
    airPressure.insert(344,1000)
    humidity.insert(344, 90)

    wind.insert(392, 2)
    airTemp.insert(392, 5.5)
    cloudCover.insert(392, 0.78)
    airPressure.insert(392, 1000)
    humidity.insert(392, 85)

    wind.insert(729, 7)
    airTemp.insert(729, 5)
    cloudCover.insert(729, 0.36)
    airPressure.insert(729, 1001)
    humidity.insert(729, 80)

    wind.insert(1736, 4)
    airTemp.insert(1736, 5)
    cloudCover.insert(1736, 0.36)
    airPressure.insert(1736, 1002)
    humidity.insert(1736, 75)

    wind.insert(1760, 3)
    airTemp.insert(1760, 7)
    cloudCover.insert(1760, 0.85)
    airPressure.insert(1760, 990)
    humidity.insert(1760, 85)

    wind.insert(1952, 4)
    airTemp.insert(1952, 8)
    cloudCover.insert(1952, 0.2)
    airPressure.insert(1952, 995)
    humidity.insert(1952, 75)

    wind[1969:1969] = 2, 2, 2, 3, 3, 3
    airTemp[1969:1969] = 5, 6, 7, 7, 8, 8
    cloudCover[1969:1969] = 0.75, 0.8, 0.8, 0.8, 0.8, 0.9
    airPressure[1969:1969] = 1000, 1000, 1000, 1000, 1000, 1000
    humidity[1969:1969] = 86, 83, 80, 77, 74, 71

    wind.insert(2623, 3)
    airTemp.insert(2623, 11.5)
    cloudCover.insert(2623, 0.85)
    airPressure.insert(2623, 1025)
    humidity.insert(2623, 80)

    wind.insert(2719, 4)
    airTemp.insert(2719, 11)
    cloudCover.insert(2719, 0.9)
    airPressure.insert(2719, 1013)
    humidity.insert(2719, 76)

    #split each weather parameter into 24 hour periods
    wind = list(chunks(wind, 24))
    airTemp = list(chunks(airTemp, 24))
    humidity = list(chunks(humidity, 24))
    cloudCover = list(chunks(cloudCover, 24))
    airPressure = list(chunks(airPressure, 24))
    
    # hour to loop through
    hour = range(24)

    # average yearly initial temperatures and energy values
    Tw_C = [9.5, 9.4, 9.1, 8.8, 8.4, 8.2, 8, 7.9, 7.8]
    InitialEnergy = [519808, 334464, 177913, 69102, -6227, -4970, -3635, -2464, -1272, 0]
    concreteTemp = [266.3171482, 267.0560374, 267.508365, 267.7275081, 267.7840425, 267.7377899, 267.6432967, 267.5625629, 267.5196339, 267.5196339]


    # calculate the average building distances and heights and canal width and direction.
    #these will be used to estimate initial temperature and energy values by iterating them
    B1_d = mean(B1_dist)
    B2_d = mean(B2_dist)
    B1_h = mean(B1_height)
    B2_h = mean(B2_height)
    canalW = mean(canalWidth)
    angle = mean(canalDirection)
    

    # get the entered date
    selectedDate = datetime.strptime(request.POST['date'], "%Y-%m-%d").date()

    # change format to day of year
    DOY = selectedDate.toordinal() - datetime.strptime(("2017-12-31"), "%Y-%m-%d").date().toordinal()

    # blank lists for weather data over specified day
    airT = []
    Wind = []
    cCover = []
    rHumidity = []
    airPress = []

    # loop through temperature for specified day
    for t in airTemp[DOY-1]:
            
            # add user-defined value onto actual air temperature and append to list
            airT.append(t+int(request.POST['airTemp']))

    # loop through wins speed for specified day
    for w in wind[DOY-1]:

            # add user-defined value onto actual wins speed
            win = w+int(request.POST['wind'])
            
            # if wind is smaller than 0 then =0 as can't have a negative wind speed
            if win < 0:
                    win = 0
                    
            # append to list
            Wind.append(win)
            
    # loop through cloud cover for specified day
    for c in cloudCover[DOY-1]:

            # add user-defined value onto actual cloud cover
            cloud = c+int(request.POST['cloudC'])
            
            #if cloud cover is greater than 1 then =1, if smaller than 0 =0, cloud cover must be between 0-100%
            if cloud >1:
                          cloud = 1
            if cloud <0:
                          cloud = 0

            # append to list
            cCover.append(cloud)
            
    # loop through humidity for specified day
    for h in humidity[DOY-1]:

            # add user-defined value onto actual humidity
            hum = h+int(request.POST['humidity'])
            
            #humidity must be between 0 and 100%
            if hum >100:
                    hum = 100
            if hum < 0:
                    hum = 0

            # append to list
            rHumidity.append(hum)
            
    # loop through air pressure for specified day                 
    for p in airPressure[DOY-1]:

            # add user-defined value onto actual air pressure and append to list
            airPress.append(p+int(request.POST['airP']))

    # get user inputted cloud height
    cloudH = float(request.POST['cloudH'])
    

    '''
    to get water temp, concrete temp and energy for particular DOY
    '''
    count = 0
    # iterate 5 times
    while count < 5:
        count = count+1
        # for every hour in the day and hourly weather inputs
        for hr, RH, Ta, U, airP, cloudC in zip(hour, rHumidity, airT, Wind, airPress, cCover):
                
            #calculate new concrete temperature
            concreteTemp = ConcreteNewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, Ta, RH, cloudH, cloudC, U, concreteTemp)
            
            #calculate new water temperature
            Tw_C = NewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, InitialEnergy, Tw_C)
            
            # calculate new energy
            InitialEnergy = HeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, Ta, U, RH, airP, cloudC, cloudH, InitialEnergy, Tw_C)


    # to calculate results to use for the map:
    mapTempChange = []
    shade = []

    # for different sections of the canal
    for B1_h, B2_h, B1_d, B2_d, canalW, angle in zip(B1_height, B2_height, B1_dist, B2_dist, canalWidth, canalDirection):

        #for every hour in a day and hourly weather
        for hr, RH, Ta, U, airP, cloudC in zip(hour, rHumidity, airT, Wind, airPress, cCover):

            # net cooling result every hour
            tempChange = NetTempChange(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, Ta, U, RH, airP, canalW, cloudH, cloudC, angle, concreteTemp, InitialEnergy, Tw_C)
        
            #append result to the list
            mapTempChange.append(tempChange)
            
            #calculate new concrete temperature
            concreteTemp = ConcreteNewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, Ta, RH, cloudH, cloudC, U, concreteTemp)
            
            #calculate new water temperature
            Tw_C = NewTemp(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, InitialEnergy, Tw_C)
            
            # calculate new energy
            InitialEnergy = HeatTransfer(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle, Ta, U, RH, airP, cloudC, cloudH, InitialEnergy, Tw_C)

            # calculate shaded width of canal and append to a list
            shade.append(ShadedWidth(B1_h, B1_d, B2_h, B2_d, DOY, hr, 0, canalW, angle))



    # this is for use in the cloud cover graph - change back to % cloud cover - easier to read for user
    cloudC = []
    for c in cCover:
            cloudC.append(c*100)

    
    '''
    output 3 - map
    '''
    # Split list into list of lists, each list at a different point along the canal
    mapTempChange = list(chunks(mapTempChange, 24))
    shade = list(chunks(shade, 24))

    # blank lists
    maxTempChange = []
    position = []

    # get the minimum temperature change in each list(i.e. the maximum urban cooling)
    for t in mapTempChange:

            # min temp change 
            m = min(t)

            #append to list(convert to positive number to represent urban cooling)
            maxTempChange.append(m*-1)

            #get position of this value in the list
            #use min() function to get one value (in case more than one position is selected)
            pos = min([i for i, j in enumerate(t) if j==m])
            
            #...and append it to a list
            position.append(pos)


    # get shade to display in extra info window on map
    maxShade = []
    
    # get the shaded width of canal at time of maximum cooling and append to a list
    for s, p, width in zip(shade, position, canalWidth):
            maxShade.append(s[p])
            
    
    ## Transform coordinates from WGS84 to BNG so markers can be placed on map
    # create a 'Proj' object for WGS84
    p1 = Proj(init='epsg:4326')

    # create a 'Proj' object for British National Grid
    p2 = Proj(init='epsg:27700')

    mapData = []
    for i in range(len(xs)):
            
        # transform the coordinates
        x, y = transform(p2,p1,xs[i],ys[i])

        # append lat and long coordinates and max temeprature change to list of lists   
        mapData.append([x, y, maxTempChange[i]])


    '''
    output 1 - figure of max cooling
    '''
    # get the maximum and minimum cooling through the day and round to 3 decimal places
    # also used to set the colour of the dots on on the map
    maxCool = round((max(maxTempChange)), 3)
    minCool = round((min(maxTempChange)), 3) 


    '''
    output 2 - line graph
    '''
    # get the average of each list i.e. average cooling at 0:00 at various points, average at 1:00, etc.
    graphTempChange = mean(mapTempChange, axis=0)

  

    #anything you want to send to the results page goes in the context
    context = {
        # output 1
        'maxCool': maxCool,

        #output 2
        'graphTempChange': graphTempChange,

        #output 3
        'mapTempChange': mapData,

        #date
        'date': selectedDate,

        #weather data
        'airTemp': airT,
        'wind': Wind,
        'humidity': rHumidity,
        'airP': airPress,
        'cloudC': cloudC,
        'cloudH': cloudH,

        # extra info window on map
        'canalWidth': canalWidth,
        'canalDirection': canalDirectionDegrees,
        'position': position,
        'shade': maxShade,

        # used to colour markers on map
        'minCool': minCool
    }

    return render(request, 'cooling/results.html', context)
