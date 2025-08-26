# -*- coding: utf-8 -*-

# module constant

# Copyright (C) 2019 National Research Council Canada
# Author:  Harold Parks

# This file is part of MetroloPy.

# MetroloPy is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.

# MetroloPy is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
# details.

# You should have received a copy of the GNU General Public License along with 
# MetroloPy. If not, see <http://www.gnu.org/licenses/>.

"""
Constants used by both the siunits and the codata2018 modules.  Most constants
are from CODATA 2018
"""

import numpy as np
from warnings import warn
from .gummy import gummy
from .ummy import ummy,_getfinfo,_iinfo
from .unit import MFraction


def _rounding_u(x):
    if not ummy.rounding_u:
        try:
            fi,_ = _getfinfo(x)
            if fi is _iinfo:
                raise ValueError()
            return ummy(x,x*fi.rel_u)
        except:
            warn('numpy.finfo cannot get the floating point accuracy for float64 math constants')
            return ummy(x)
    else:
        return ummy(x)
    
pi = _rounding_u(np.pi)
euler = _rounding_u(np.e)
sqrt2 = _rounding_u(np.sqrt(2))


# constants from CODATA 2022:
    
(
alph, # fine structure constant
aral,  # alpha particle mass in u
ryd, # Rydberg constant (1/m)
are, # electron mass in u
mpsme, # proton-electron mass ratio
ae, # electron mag. mom. anomaly
rd, # deuteron rms charge radius (m)
ghn, # helion g factor
arh, # helion relative atomic mass
sigmah, # helion shielding shift
mmu, # muon mass (kg)
amu, # muon mag. mom. anomaly
arn, # neutron mass in u
gnn, # neutron g factor
ard, # deuteron relative atomic mass
gdn, # deuteron g factor
gp, # proton g factor
sigmapp, # proton sheilding factor
mtu, # triton relative atomic mass
gtn # triton g factor
) = gummy.create([7.2973525643e-3,       
           4.001506179129,        
           10973731.568157,      
           5.485799090441e-4,    
           1836.152673426,         
           1.15965218046e-3,    
           2.12778e-15,            
           -4.2552506995,          
           3.014932246932,        
           5.9967029e-5,           
           1.883531627e-28,       
           1.16592062e-3,         
           1.00866491606,         
           -3.82608552,            
           2.013553212544,        
           0.8574382335,           
           5.5856946893,           
           2.56715e-5,             
           3.01550071597,         
           5.957924930],
           u=[0.0000000011e-3,
              0.000000000062,
              0.000012,
              0.000000000097e-4,
              0.000000032,
              0.00000000018e-3,
              0.00027e-15 ,
              0.0000000034,
              0.000000000074,
              0.0000023e-5,
              0.000000042e-28,
              0.00000041e-3,
              0.00000000040,
              0.00000090,
              0.000000000015,
              0.000000002,
              0.0000000016,
              0.00041e-5,
              0.00000000010,
              0.000000012],
           correlation_matrix=[[       1, 0.00001, 0.00125,-0.03972, 0.03527, 0.97556, 0.00036, 0.00000,-0.00332, 0.00000,-0.00014, 0.00087, 0.00168, 0.00000,-0.01036, 0.00002, 0.00000,-0.00003,-0.00237, 0.00000],
                               [ 0.00001,       1, 0.00001,-0.00031, 0.00028, 0.00001,-0.00001, 0.00000,-0.00002, 0.00000, 0.00000,-0.00000, 0.00000,-0.00007, 0.00000, 0.00000, 0.00000, 0.00000,-0.00002, 0.00000],
                               [ 0.00125, 0.00001,       1,-0.03754,-0.04138, 0.00225, 0.60917, 0.00000, 0.00203, 0.00000, 0.00010, 0.00000, 0.00010, 0.00000, 0.00645, 0.00002, 0.00000,-0.00002, 0.00145, 0.00000],
                               [-0.03972,-0.00031,-0.03754,       1,-0.88879,-0.03875, 0.02306, 0.00009, 0.08307, 0.00000,-0.00001, 0.00000, 0.00397, 0.00001, 0.25967,-0.00035,-0.00008, 0.00052, 0.05946,-0.00001],
                               [ 0.03527, 0.00028,-0.04138,-0.88879,       1, 0.03441,-0.02598,-0.00010, 0.05140, 0.00000, 0.00001,-0.00001, 0.00130,-0.00001, 0.14561, 0.00040, 0.00009,-0.00059,-0.00232, 0.00000],
                               [ 0.97556, 0.00001, 0.00225,-0.03875, 0.03441,       1, 0.00035, 0.00000,-0.00323, 0.00000,-0.00013, 0.00085, 0.00164, 0.00000,-0.01011, 0.00002, 0.00000,-0.00003, 0.00116, 0.00000],
                               [ 0.00036,-0.00001, 0.60917, 0.02306,-0.02598, 0.00035,       1, 0.00000,-0.00134, 0.00000, 0.00006, 0.00000,-0.00003, 0.00000,-0.00375,-0.00001, 0.00000, 0.00002,-0.00095, 0.00000],
                               [ 0.00000, 0.00000, 0.00000, 0.00009,-0.00010, 0.00000, 0.00000,       1,-0.00001,-0.02844,-0.00004,-0.00199, 0.00000, 0.00296,-0.00001, 0.00000, 0.00000, 0.17032, 0.00000, 0.00000],
                               [-0.00332,-0.00002, 0.00203, 0.08307, 0.05140,-0.00323,-0.00134,-0.00001,       1, 0.00000, 0.00000,-0.00001, 0.00513, 0.00000, 0.29794, 0.00002, 0.00000,-0.00003, 0.71472, 0.00000],
                               [ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,-0.02844, 0.00000,       1, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
                               [-0.00014, 0.00000, 0.00010,-0.00001, 0.00001,-0.00013, 0.00006,-0.00004, 0.00000, 0.00000,       1, 0.08172, 0.00000, 0.00000, 0.00000,-0.00002,-0.00018,-0.00023, 0.00000,-0.00003],
                               [ 0.00087,-0.00000, 0.00000, 0.00000,-0.00001, 0.00085, 0.00000,-0.00199,-0.00001, 0.00000, 0.08172,       1, 0.00000,-0.00020,-0.00002, 0.00000,-0.00001,-0.01160, 0.00000, 0.00000],
                               [ 0.00168, 0.00000, 0.00010, 0.00397, 0.00130, 0.00164,-0.00003, 0.00000, 0.00513, 0.00000, 0.00000, 0.00000,       1, 0.00000, 0.01921, 0.00000, 0.00000, 0.00000, 0.00366, 0.00000],
                               [ 0.00000,-0.00007, 0.00000, 0.00001,-0.00001, 0.00000, 0.00000, 0.00296, 0.00000, 0.00000, 0.00000,-0.00020, 0.00000,       1, 0.00000, 0.00000, 0.00000, 0.01731, 0.00000, 0.00000],
                               [-0.01036, 0.00000, 0.00645, 0.25967, 0.14561,-0.01011,-0.00375,-0.00001, 0.29794, 0.00000, 0.00000,-0.00002, 0.01921, 0.00000,       1, 0.00006, 0.00001,-0.00009, 0.21297, 0.00000],
                               [ 0.00002, 0.00000, 0.00002,-0.00035, 0.00040, 0.00002,-0.00001, 0.00000, 0.00002, 0.00000,-0.00002, 0.00000, 0.00000, 0.00000, 0.00006,       1, 0.10706, 0.00770, 0.00001, 0.01581],
                               [ 0.00000, 0.00000, 0.00000,-0.00008, 0.00009, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,-0.00018,-0.00001, 0.00000, 0.00000, 0.00001, 0.10706,       1, 0.07194, 0.00000, 0.14767],
                               [-0.00003, 0.00000,-0.00002, 0.00052,-0.00059,-0.00003, 0.00002, 0.17032,-0.00003, 0.00000,-0.00023,-0.01160, 0.00000, 0.01731,-0.00009, 0.00770, 0.07194,       1,-0.00002, 0.01062],
                               [-0.00237,-0.00002, 0.00145, 0.05946,-0.00232, 0.00116,-0.00095, 0.00000, 0.71472, 0.00000, 0.00000, 0.00000, 0.00366, 0.00000, 0.21297, 0.00001, 0.00000,-0.00002,       1, 0.00000],
                               [ 0.00000, 0.00000, 0.00000,-0.00001, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,-0.00003, 0.00000, 0.00000, 0.00000, 0.00000, 0.01581, 0.14767, 0.01062, 0.00000,       1]]
)


c = 299792458 # speed of light in m/s
h = MFraction('6.62607015e-34') # planck constant in J s

    
e = MFraction('1.602176634e-19') # electron charge in C
k = MFraction('1.380649e-23') # Boltzmann constant in J/K

KJ = 2*e/h
RK = h/e**2

KJ90 = MFraction('483597.9e9') # Josephson constant Hz/V, 1990 conventional value
RK90 = MFraction('25812.807') # Von Klitzing constant in ohm, 1990 conventional value
    
G = ummy(6.67430e-11,0.00015e-11) # gravitational constant in m**3/kg s**2, CODATA 2022

dalton = 2*ryd*h/(are*c*alph**2)
me = 2*ryd*h/(c*alph**2) #electron mass in kg

a0 = alph/(4*pi*ryd) # bohr radius in m

mp = mpsme*me # in kg


# constants from IAU 2009:
earth_mass = ummy(3.986004418e14,8e5)/G # in kg
solar_mass = ummy(1.32712442099e20,1e10)/G # in kg
jupiter_mass = solar_mass/ummy(1.047348644e3,1.7e-5) # in kg

hbar = h/(2*pi)