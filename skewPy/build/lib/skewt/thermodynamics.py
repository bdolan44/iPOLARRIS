from numpy import interp,pi,cos,sin,arctan2,sqrt,linspace,min,exp,log,where

#-----------------------------------------------------------------------
# Here we go. A set of functions that I use from time to time to calculate 
# the basic stuff that I'm sick of doing over and over! I'm going to 
# endeavour to include references and global constants to make it all nice 
# and legible.
#-----------------------------------------------------------------------

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218	      # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=R_s_da/R_s_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)

_devel="working"

def Theta(tempk,pres,pref=100000.):
    """Potential Temperature

    INPUTS: 
    tempk (K)
    pres (Pa)
    pref: Reference pressure (default 100000 Pa)

    OUTPUTS: Theta (K)

    Source: Wikipedia
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """

    try:
        minpres=min(pres)
    except TypeError:
        minpres=pres

    if minpres<2000:
        print "WARNING: P<2000 Pa; did you input a value in hPa?"

    return tempk*(pref/pres)**(Rs_da/Cp_da)

def TempK(theta,pres,pref=100000.):
    """Inverts Theta function."""

    try:
        minpres=min(pres)
    except TypeError:
        minpres=pres

    if minpres<2000:
        print "WARNING: P<2000 Pa; did you input a value in hPa?"

    return theta*(pres/pref)**(Rs_da/Cp_da)

def ThetaE():
    """Equivalent potential temperature"""
    raise NotImplementedError


def ThetaV(tempk,pres,e):
    """Virtual Potential Temperature

    INPUTS
    tempk (K)
    pres (Pa)
    e: Water vapour pressure (Pa) (Optional)
    """ 

    mixr=MixRatio(e,pres)
    theta=Theta(tempk,pres)

    return theta*(1+mixr/Epsilon)/(1+mixr)

def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """

    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=MixRatio(es,pres)

    if e is None:
    # assume saturated
        e=es

    w=MixRatio(e,pres)

    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma

def Density(tempk,pres,mixr):
    """Density of moist air"""

    virtualT=VirtualTempFromMixR(tempk,mixr)
    return pres/(Rs_da*virtualT)


def VirtualTemp(tempk,pres,e):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    e: vapour pressure (Pa)
    p: static pressure (Pa)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia)."""

    tempvk=tempk/(1-(e/pres)*(1-Epsilon))
    return tempvk


def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk*(1.0+0.6*mixr)

def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """

    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)

def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice


    RETURNS: e_sat  (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
    # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
    # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError

def RH(T, Td):
# calculate the RH given temperature and dew point

    T += 273.15
    Td += 273.15
    return 100.0*(np.exp((17.625*Td)/(243.04+Td))/np.exp((17.625*T)/(243.04+T)))

def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure

    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)

def MixR2VaporPress(qv,p):
    """Return Vapor Pressure given Mixing Ratio and Pressure
    INPUTS
    qv (kg kg^-1) Water vapor mixing ratio`
    p (Pa) Ambient pressure

    RETURNS
    e (Pa) Water vapor pressure
    """

    return qv*p/(Epsilon+qv)



def VaporPressure(dwpt):
    """Water vapor pressure
    INPUTS
    dwpt (C) Dew Point Temperature (for SATURATION vapor 
    pressure use tempc)

    RETURNS
    e (Pa) Water Vapor Pressure

    SOURCE:
    Bolton, Monthly Weather Review, 1980, p 1047, eq. (10)
    """

    return 611.2*exp(17.67*dwpt/(243.5+dwpt))

def DewPoint(e):
    """ Use Bolton's (1980, MWR, p1047) formulae to find tdew.
    INPUTS:
    e (Pa) Water Vapor Pressure
    OUTPUTS:
    Td (C) 
    """

    ln_ratio=log(e/611.2)
    Td=((17.67-ln_ratio)*degCtoK+243.5*ln_ratio)/(17.67-ln_ratio)
    return Td-degCtoK


