import numpy as np

# Multiply simulation units by the following to get real units 
class SI:
    def __init__(self):
        # Fundamental units
        self.pc   = 3.085677581e16           # [m]  m in a pc
        self.myr  = 365.25*24.*3600.*1e6     # [s]  s in a Myr
        self.msun = 1.9886e30                # [kg] Solar mass 
        self.m_h  = 1.6737236e-27            # [kg] Hydrogen mass
        self.k_b  = 1.380649e-23             # [J/K] Boltzmann's Constant
        self.c    = 2.99792458e8             # [m/s] Speed of light
        self.h    = 6.62607004e-34           # [J s] Planck Constant
        self.hbar = 1.0545718e-34            # [J s] reduced Planck constant
        # Derived units
            # Energy
        self.e  = (self.msun*self.pc**2
                        /self.myr**2)        # [J], Energy unit
        self.edens = self.e/self.pc**3       # [J/m^3], Energy dens unit
        self.esdens = self.edens*self.pc     # [J/m^2], Energy surface dens unit
            # Momentum
        self.mom  = (self.msun*self.pc
                        /self.myr)           # [kg m/s], momentum unit
        self.momdens = self.mom/self.pc**3   # [kg/(m^2 s)], mom dens unit
        self.momsdens = self.momdens*self.pc # [kg/(m s)], mom surface dens unit
            # Velocity 
        self.v   = self.pc/self.myr          # [m/s], velocity unit
            # Density   
        self.rho = self.msun/self.pc**3      # [kg/m^3], density unit
        self.rhos = self.msun/self.pc**2     # [kg/m^2], surface density unit
            # Boltzmann Constant
        self.k_b = self.e*1.38065e-23        # [J/K] 

class CGS:
    def __init__(self):
        # Fundamental units
        self.pc   = 3.085677581e18           # [cm]  cm in a pc
        self.myr  = 365.25*24.*3600.*1e6     # [s]  s in a Myr
        self.msun = 1.9886e33                # [g] Solar mass 
        self.m_h  = 1.6737236e-24            # [g] Hydrogen mass
        self.k_b  = 1.380649e-16             # [erg/K] Boltzmann's Constant
        self.c    = 2.99792458e10            # [cm/s] Speed of light
        self.h    = 6.62607004e-27           # [erg s] Planck Constant
        self.hbar = 1.0545718e-27            # [erg s] reduced Planck constant

        # Derived units
            # Energy
        self.e  = (self.msun*self.pc**2
                        /self.myr**2)        # [erg], Energy unit
        self.edens = self.e/self.pc**3       # [erg/cm^3], Energy dens unit
        self.esdens = self.edens*self.pc     # [erg/cm^2], Energy surface dens unit
            # Momentum
        self.mom  = (self.msun*self.pc
                        /self.myr)           # [g cm/s], momentum unit
        self.momdens = self.mom/self.pc**3   # [g/(cm^2 s)], mom dens unit
        self.momsdens = self.momdens*self.pc # [g/(cm s)], mom surface dens unit
            # Velocity 
        self.v   = self.pc/self.myr          # [cm/s], velocity unit
        # Density   
        self.rho = self.msun/self.pc**3      # [kg/m^3], density unit
        self.rhos = self.msun/self.pc**2     # [kg/m^2], surface density unit

class GAL:
    # "Galactic" units (pc, Myr, Msol)
    def __init__(self):
        # Fundamental units
        self.pc   = 1.0          # [pc]  
        self.myr  = 1.0          # [Myr]  
        self.msun = 1.0          # [Solar mass] 
        self.m_h  = 1.6737236e-24/1.9886e33         # Hydrogen mass in Solar mass unit
        self.k_b  = 1.3806e-23*(3.145e7*1e6)**2/ \
                        (1.9886e30*3.086e16**2)     # [M_sun*pc^2/(K*Myr^2)] Boltzmann's Constant
        self.G    = 4.4492e-3                       # Gravitational constant 
        self.c    = 2.99792458e10*(3.145e13/3.0866e18)  # [cm/s] Speed of light
        # Note this is incorrect -- do this later
        self.h    = 1.0                         # [erg s] Planck Constant
        self.hbar = 1.0                          # [erg s] reduced Planck constant

        # Derived units
            # Energy
        self.e  = (self.msun*self.pc**2
                        /self.myr**2)        # [erg], Energy unit
        self.edens = self.e/self.pc**3       # [erg/cm^3], Energy dens unit
        self.esdens = self.edens*self.pc     # [erg/cm^2], Energy surface dens unit
            # Momentum
        self.mom  = (self.msun*self.pc
                        /self.myr)           # [g cm/s], momentum unit
        self.momdens = self.mom/self.pc**3   # [g/(cm^2 s)], mom dens unit
        self.momsdens = self.momdens*self.pc # [g/(cm s)], mom surface dens unit
            # Velocity 
        self.v   = self.pc/self.myr          # [cm/s], velocity unit
        # Density   
        self.rho = self.msun/self.pc**3      # [kg/m^3], density unit
        self.rhos = self.msun/self.pc**2     # [kg/m^2], surface density unit

class ISM:
    # "ISM" units: n=1, T=1, G=1. Conversion to cgs. Astro units provided, but need to convert yourself.
    def __init__(self,mu=1):
        self.mu    = mu
        self.n0    = 1
        self.T0    = 1
        self.G0    = 1
        self.Myr   = 3.6e3*2.4e1*3.65e2*1e6
        self.Msun  = 1.9886e33
        self.Pc    = 3.0866e18
        self.mdns  = mu*1.6737236e-24*self.n0
        self.time  = 1/np.sqrt(6.67408e-8*self.mdns)
        self.velo  = np.sqrt(1.3806e-16*self.T0/(mu*1.6737236e-24))
        self.leng  = self.time*self.velo
        self.mass  = self.mdns*self.leng**3
        self.myrs  = self.time/self.Myr
        self.prsc  = self.leng/self.Pc
        self.kmsc  = self.velo*1e-5
        self.msol  = self.mass/self.Msun
        self.ener  = self.n0*self.T0*1.3806e-16*(self.leng)**3
        self.bfld  = np.sqrt(2.0*self.n0*self.T0*1.3806e-16)
        self.muga  = self.bfld*1e6

    def print(self):
        print("%25s%12s%20s"   % ("quantity"            , "comp unit","phys unit"))
        print("%25s%12s%20.8e" % ("density cm^(-3)"     , "1"        ,self.n0))
        print("%25s%12s%20.8e" % ("temperature K"       , "1"        ,self.T0))
        print("%25s%12s%20.8e" % ("Newton's constant"   , "1"        ,6.67408e-8))
        print("%25s%12s%20.8e" % ("time seconds"        , "1"        ,self.time))
        print("%25s%12s%20.8e" % ("time Myr"            , "1"        ,self.myrs))
        print("%25s%12s%20.8e" % ("velocity cm/s"       , "1"        ,self.velo))
        print("%25s%12s%20.8e" % ("velocity km/s"       , "1"        ,self.kmsc))
        print("%25s%12s%20.8e" % ("length cm"           , "1"        ,self.leng))
        print("%25s%12s%20.8e" % ("length pc"           , "1"        ,self.prsc))
        print("%25s%12s%20.8e" % ("mass g"              , "1"        ,self.mass))
        print("%25s%12s%20.8e" % ("mass Msol"           , "1"        ,self.msol))
        print("%25s%12s%20.8e" % ("energy erg"          , "1"        ,self.ener))
        print("%25s%12s%20.8e" % ("bfield muG"          , "1"        ,self.muga))

