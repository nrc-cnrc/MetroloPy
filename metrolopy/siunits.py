# -*- coding: utf-8 -*-

# module siunits

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
This module is loaded by the gummy.units module and is not intended be be
imported directly.  

The most of the units here are from the SI Brochure, draft 9th edition published 
by the BIPM.
"""

from .gummy import gummy, _lg10
from .ummy import MFraction
from .unit import _BuiltinLib,Conversion,Unit
from .prefixedunit import PrefixedUnit,BinaryPrefixedUnit
from .logunit import LogUnit,LogConversion
from .offsetunit import OffsetUnit,OffsetConversion
from .functions import sqrt
import numpy as np

pi = np.pi
e = np.e

# constants from draft of the 9 edition of the SI brochure dated 18 Dec. 2018
_const_c = 299792458 # speed of light in m/s
_const_h = MFraction('6.62607015e-34') # planck constant in J s
_const_hbar = _const_h/(2*pi)
_const_e = MFraction('1.602176634e-19') # electron charge in C
_const_k = MFraction('1.380649e-23') # Boltzmann constant in J/K
_const_dalton = gummy(1.660539040e-27,0.000000020e-27) # unified atomic mass unit in kg

# constants from CODATA 2014
_const_G = gummy(6.67408e-11,0.00031e-11) # gravitational constant in m**3/kg s**2
_const_me = gummy(9.10938356e-31,0.00000011e-31) # mass of the electron in kg
_const_alpha = gummy(7.2973525664e-3,0.0000000017e-3) # fine structure constant
_const_rydberg = gummy(10973731.568508,0.000065) # in 1/m
_const_a0 = _const_alpha/(4*pi*_const_rydberg) # bohr radius in m
_const_proton_mass = gummy(1.672621898e-27,0.000000021e-27) # in kg

# constants from IAU 2009
_const_earth_mass = gummy(3.986004418e14,8e5)/_const_G # in kg
_const_solar_mass = gummy(1.32712442099e20,1e10)/_const_G # in kg
_const_jupiter_mass = _const_solar_mass/gummy(1.047348644e3,1.7e-5) # in kg

with _BuiltinLib():
    #SI Base units
    _m = PrefixedUnit('metre','m',additional_names=('meter',),add_symbol=True,
                      order=1,description='SI unit of length',
                      base_description='SI base unit of length')
    Unit.alias('micron','um')
    _cm = Unit.unit('cm')
    _cm.description = 'SI unit of length and CGS base unit'
    
    _kg = PrefixedUnit('gram','g',add_symbol=True,order=1,base_prefix='kilo',
                      description='SI unit of mass',
                      base_description='SI base unit of mass')
    Unit.unit('g').description = 'SI unit of mass and CGS base unit'
    
    _s = PrefixedUnit('second','s',add_symbol=True,order=3,
                      description='SI unit of time',
                      base_description='SI and CGS base unit of time')
    
    _A = PrefixedUnit('ampere','A',add_symbol=True,order=4,
                      description='SI unit of electrical current',
                      base_description='SI base unit of electrical current')
    
    _K = PrefixedUnit('kelvin','K',add_symbol=True,order=5,
                      description='SI unit of temperature',
                      base_description='SI base unit of temperature')
    
    _mol = PrefixedUnit('mole','mol',add_symbol=True,order=6,
                      description='SI unit for amount of substance',
                      base_description='SI base unit for amount of substance')
    
    _cd = PrefixedUnit('candela','cd',add_symbol=True,order=7,
                      description='SI unit of luminosity',
                      base_description='SI base unit of luminosity')
    
    _ipk = Unit('international prototype kilogram','m(K)',
                Conversion(_kg,gummy(1,1.2e-8)),add_symbol=True,
                html_symbol='<i>m</i>(&#x1d4a6;)',latex_symbol='m({\mathcal {K}})',
                description='IPK, le grand k')
    Unit.alias('IPK',_ipk)
    
    #Coherent derived units in the SI with special names and symbols
    _rad = PrefixedUnit('radian','rad',Conversion(_m*_m**-1,1),add_symbol=True,order=0,
                        description='unit of angle, angular unit')
    _sr = PrefixedUnit('steradian','sr',Conversion(_m**2*_m**-2,1),add_symbol=True,order=0,
                       description='unit of solid angle, angular unit')
    _Hz = PrefixedUnit('hertz','Hz',Conversion(_s**-1,1),add_symbol=True,order=0,
                       description='frequency')
    _N = PrefixedUnit('newton','N',Conversion(_m*_kg*_s**-2,1),add_symbol=True,order=0,
                      description='SI derived unit of force')
    _Pa = PrefixedUnit('pascal','Pa',Conversion(_N*_m**-2,1),add_symbol=True,order=0,
                       description='SI derived unit of pressure')
    _J = PrefixedUnit('joule','J',Conversion(_N*_m,1),add_symbol=True,order=0,
                      description='SI derived unit of energy')
    _W = PrefixedUnit('watt','W',Conversion(_J*_s**-1,1),add_symbol=True,order=0,
                      description='SI derived unit of power')
    _C = PrefixedUnit('coulomb','C',Conversion(_s*_A,1),add_symbol=True,order=0,
                      description='SI derived unit of electrical charge')
    _V = PrefixedUnit('volt','V',Conversion(_W*_A**-1,1),add_symbol=True,order=0,
                      description='SI derived unit of electric potential, voltage')
    PrefixedUnit('farad','F',Conversion(_C*_V**-1,1),add_symbol=True,order=0,
                 description='SI derived unit of electrical capacitance')
    _ohm = PrefixedUnit('ohm','\u03A9',Conversion(_V*_A**-1,1),add_symbol=True,order=0,
                 html_symbol='&Omega;',latex_symbol=r'\Omega',ascii_symbol='ohm',
                 description='SI derived unit of electrical resistance')
    _siemens = PrefixedUnit('siemens','S',Conversion(_A*_V**-1,1),add_symbol=True,order=0,
                 description='SI derived unit of electrical conductance, susceptance and admittance')
    Unit.alias('mho',_siemens)
    _Wb = PrefixedUnit('weber','Wb',Conversion(_V*_s,1),add_symbol=True,order=0,
                       description='SI derived unit of magnetic flux')
    _T = PrefixedUnit('tesla','T',Conversion(_Wb*_m**-2,1),add_symbol=True,order=0,
                      description='SI derived unit of magnetic flux density')
    PrefixedUnit('henry','H',Conversion(_Wb*_A**-1,1),add_symbol=True,order=0,
                 description='SI derived unit of electrical inductance')
    
    _degC = OffsetUnit('degree Celsius','\u00B0C',OffsetConversion(_K,273.15),
               latex_symbol='^{\circ}C',ascii_symbol = 'degC',add_symbol=True,order=0,
               description='unit of temperature')
    Unit.alias('degree C',_degC)
    Unit.alias('deg C',_degC)
    
    _lm = PrefixedUnit('lumen','lm',Conversion(_cd*_sr,1),add_symbol=True,order=0,
                       description='SI derived unit for luminous flux')
    PrefixedUnit('lux','lx',Conversion(_lm*_m**-2,1),add_symbol=True,order=0,
                 description='SI derived unit for illuminance and luminous emittance')
    PrefixedUnit('becquerel','Bq',Conversion(_s**-1,1),add_symbol=True,order=0,
                 description='SI derived unit for radioactivity')
    PrefixedUnit('gray','Gy',Conversion(_J*_kg**-1,1),add_symbol=True,order=0,
                 description='SI derived unit for ionizing radiation')
    PrefixedUnit('sievert','Sv',Conversion(_J*_kg**-1,1),add_symbol=True,order=0,
                 description='SI derived unit for ionizing radiation dose')
    PrefixedUnit('katal','kat',Conversion(_mol*_s**-1,1),add_symbol=True,order=0,
                 description='SI unit for catalytic activity')
    
    # 1990 conventional unit definitions
    _KJ90 = MFraction('483597.9e9') # Josephson constant Hz/V, 1990 conventional value
    _KJ = 2*_const_e/_const_h
    _RK90 = MFraction('25812.807') # Von Klitzing constant in ohm, 1990 conventional value
    _RK = _const_h/_const_e**2
    _V90 = PrefixedUnit('volt 90','V(90)',Conversion(_V,_KJ90/_KJ),add_symbol=True,
                        html_symbol='V<sub>90</sub>',latex_symbol='V_{90}',
                        description='1990 conventional definition of voltage, electrical potential')
    _ohm90 = PrefixedUnit('ohm 90','\u03A9(90)',Conversion(_ohm,_RK/_RK90),add_symbol=True,
                          html_symbol='&Omega;<sub>90</sub>',latex_symbol=r'\Omega_{90}',ascii_symbol='ohm(90)',
                          description='1990 conventional definition of electrical resistance')
    _A90 = PrefixedUnit('ampere 90','A(90)',Conversion(_V90*_ohm90**-1,1),add_symbol=True,
                        html_symbol='A<sub>90</sub>',latex_symbol='A_{90}',
                        description='1990 conventional definition of electrical current')
    PrefixedUnit('coulomb 90','C(90)',Conversion(_s*_V90*_ohm90**-1,1),add_symbol=True,
                 html_symbol='C<sub>90</sub>',latex_symbol='C_{90}',
                 description='1990 conventional definition of electrical charge')
    PrefixedUnit('watt 90','W(90)',Conversion(_V90*_A90,1),add_symbol=True,
                 html_symbol='W<sub>90</sub>',latex_symbol='W_{90}',
                 description='1990 conventional definition of electrical power')
    PrefixedUnit('farad 90','F(90)',Conversion(_s*_ohm90**-1,1),add_symbol=True,
                 html_symbol='F<sub>90</sub>',latex_symbol='F_{90}',
                 description='1990 conventional definition of electrical capacitance')
    PrefixedUnit('henry 90','H(90)',Conversion(_s*_ohm90,1),add_symbol=True,
                 html_symbol='H<sub>90</sub>',latex_symbol='H_{90}',
                 description='1990 conventional definition of electrical inductance')
    
    # 1954 Kelvin definition
    PrefixedUnit('kelvin 54','K(54)',Conversion(_K,gummy(1,u=3.7e-7)),add_symbol=True,
                 html_symbol='K<sub>54</sub>',latex_symbol='K_{54}',
                 description = 'unit of temperature, 1954 SI definition based on the triple point of water')
    
    #Non-SI units accepted for use with the International System of Units
    _min = Unit('minute','min',Conversion(_s,60),add_symbol=True,order=0,
                description='unit of time')
    _h = Unit('hour','h',Conversion(_min,60),add_symbol=True,order=0,
              description='unit of time')
    _d = Unit('day','d',Conversion(_h,24),add_symbol=True,order=0,
              description='unit of time')
    Unit.alias('D',_d)
    _deg = Unit('degree','\t\u00B0',Conversion(_rad,pi/180),add_symbol=True,order=0,
                 latex_symbol='^{\circ}',ascii_symbol='deg',
                 description='unit of angle, angular unit')
    _arcmin = Unit('arcminute',"\t'",Conversion(_deg,MFraction(1,60)),add_symbol=True,order=0,
                   short_name='arcmin',description='unit of angle, angular unit')
    Unit('arcsecond','\t"',Conversion(_arcmin,MFraction(1,60)),add_symbol=True,order=0,
         short_name='arcsec',description='unit of angle, angular unit')
    Unit('hectare','ha',Conversion('hm**2',1),add_symbol=True,order=0,
         description='unit of area')
    PrefixedUnit('litre','L',Conversion('dm**3',1),add_symbol=True,order=0,
                 additional_names=('liter',),description='unit of volume')
    _tonne = PrefixedUnit('tonne','t',Conversion(_kg,1000),add_symbol=True,order=0,
                          prefixes=['centi','deci','deca','hecto','kilo','mega',
                                    'giga','tera','peta','exa','zetta','yotta'],
                          description='unit of mass')
    Unit.alias('metric ton',_tonne)
    
    _Da = PrefixedUnit('dalton','Da',Conversion(_kg,_const_dalton),
                 add_symbol=True,description='unified atomic mass unit')
                 
    #Other non-SI units
    PrefixedUnit('bar','bar',Conversion(_Pa,100000),add_symbol=True,order=0,
                 description='unit of pressure')
    Unit('millimetre of mercury','mmHg',Conversion(_Pa,133.322387415),add_symbol=True,
         order=0,description='unit of pressure')
    Unit('angstrom','\u212B',Conversion(_m,MFraction('1e-10')),add_symbol=True,order=0,
         html_symbol='&#8491;',ascii_symbol='angstrom',description='unit of length')
    
    _nm = Unit('nautical mile','M',Conversion(_m,1852),add_symbol=True,order=0,
               description='unit of length')
    Unit.alias('NM',_nm)
    Unit.alias('nmi',_nm)
    Unit.alias('Nm',_nm)
    
    Unit('barn','b',Conversion(_m**2,MFraction('1e-28')),add_symbol=True,order=0,
         description='unit of area')
    Unit('knot','kn',Conversion(_nm*_h**-1,1),add_symbol=True,order=0,
         description='unit of speed')
    _dbp = LogUnit('decibel','dB',LogConversion(1,10,10,_lg10),add_symbol=True)
    Unit.alias('dB power',_dbp)
    Unit.alias('dB-p',_dbp)
    
    _dbf = LogUnit('decibel field','dB',LogConversion(1,20,10,_lg10),short_name='dB-f')
    Unit.alias('dB field',_dbf)
    Unit.alias('dB root-power',_dbf)
    
    def _ln(x):
        try:
            return np.log(x)
        except:
            return np.log(float(x))
    LogUnit('neper','Np',LogConversion(1,1,e,_ln),add_symbol=True)
    
    #Non-SI units associated with the CGS and the CGS-Gaussian system of units
    PrefixedUnit('erg','erg',Conversion(_J,MFraction('1e-7')),add_symbol=True,order=0,
                 description='CGS unit of energy')
    _dyn = PrefixedUnit('dyne','dyn',Conversion(_N,MFraction('1e-5')),add_symbol=True,order=0,
                        description='CGS unit of force')
    PrefixedUnit('poise','P',Conversion('dyn s cm**-2',1),add_symbol=True,order=0,
                 description='CGS unit of dynamic viscosity')
    PrefixedUnit('stokes','St',Conversion('cm**2 s**-1',1),add_symbol=True,order=0,
                 description='CGS unit of kinematic viscosity')
    PrefixedUnit('stilb','sb',Conversion('cd cm**-2',1),add_symbol=True,order=0,
                 description='CGS unit of luminance')
    PrefixedUnit('phot','ph',Conversion('cd sr cm**-2',1),add_symbol=True,order=0,
                 description='CGS unit of illuminance or luminous flux')
    PrefixedUnit('galileo','Gal',Conversion('cm s**-2',1),add_symbol=True,order=0,
                 description='CGS unit of acceleration')
    PrefixedUnit('maxwell','Mx',Conversion(_Wb,MFraction('1e-8')),add_symbol=True,order=0,
                 description='CGS unit of magnetic flux')
    PrefixedUnit('gauss','G',Conversion(_T,MFraction('1e-4')),add_symbol=True,order=0,
                 description='CGS unit of magnetic flux density')
    PrefixedUnit('oersted','Oe',Conversion(_A*_m**-1,1000/(4*pi)),add_symbol=True,
                 order=0,description='CGS unit of auxiliary magnetic field')
    
    # decibel units
    _dBSPL = LogUnit('decibel sound pressure level','dB',
                 LogConversion(gummy(20,unit='uPa'),20,10,_lg10),
                 short_name='dB(SPL)',add_symbol=False,
                 description='sound pressure level in air')
    Unit.alias('dB(20uPa)',_dBSPL)
    
    _dBuPa = LogUnit('decibel \u03BCPa','dB',
                   LogConversion(gummy(1,unit='uPa'),20,10,_lg10),
                   short_name='dB(uPa)',add_symbol=False,
                   description='sound pressure level in water')
    
    _dBSIL = LogUnit('decibel sound intensity level','dB',
                 LogConversion(gummy(1e-12,unit=_W/_m**2),10,10,_lg10),
                 short_name='dB(SIL)',add_symbol=False)
    
    _dBSWL = LogUnit('decibel sound power level','dB',
                 LogConversion(gummy(1e-12,unit=_W),10,10,_lg10),
                 short_name='dB(SWL)',add_symbol=False)
    Unit.alias('dB(pW)',_dBSWL)
    
    _dBm = LogUnit('decibel milliwatt','dBm',LogConversion(gummy(1,unit='mW'),10,10,_lg10),
                   add_symbol=True,description='unit of power')
    Unit.alias('dB(mW)',_dBm)
    
    _dBV = LogUnit('decibel volt','dBV',LogConversion(gummy(1,unit=_V),20,10,_lg10),
                   add_symbol=True,description='unit of voltage')
    Unit.alias('dB(V)',_dBV)
    
    _dBv = LogUnit('decibel u','dBu',LogConversion(gummy(sqrt(0.6),unit=_V),20,10,_lg10),
                   add_symbol=True,
                   description='unit of voltage, volt, decibel relative to sqrt(0.6) V')
    Unit.alias('dBu',_dBv)
    
    _dBmV = LogUnit('decibel millivolt','dBmV',LogConversion(gummy(1,unit='mV'),20,10,_lg10),
                    add_symbol=True,description='unit of voltage, volt')
    Unit.alias('dB(mV)',_dBmV)
    
    _dBuV = LogUnit('decibel microvolt','dB\u03bcV',LogConversion(gummy(1,unit='uV'),20,10,_lg10),
                    add_symbol=True,ascii_symbol='dBuV',description='unit of voltage')
    Unit.alias('dB(uV)',_dBuV)
    
    _dBZ = LogUnit('decibel Z','dBZ',LogConversion(gummy(1,unit='mm**6 m**-3'),10,10,_lg10),
                    add_symbol=True,description='unit used with weather radar')
    Unit.alias('dB(Z)',_dBZ)
    
    _dBJ = LogUnit('decibel joule','dBJ',LogConversion(gummy(1,unit=_J),10,10,_lg10),
                    add_symbol=True,description='unit of energy')
    Unit.alias('dB(J)',_dBJ)
    
    _dBmu = LogUnit('decibel microvolt per metre','dB\u03bc',LogConversion(gummy(1,unit='uV/m'),20,10,_lg10),
                    add_symbol=True,description='unit of voltage')
    Unit.alias('dB(uV/m)',_dBmu)
    Unit.alias('decibel microvolt per meter',_dBmu)
    
    _dBf = LogUnit('decibel femtowatt','dBf',LogConversion(gummy(1,unit='fW'),10,10,_lg10),
                    add_symbol=True,description='unit of power')
    Unit.alias('dB(fW)',_dBf)
    
    _dBW = LogUnit('decibel watt','dBW',LogConversion(gummy(1,unit=_W),10,10,_lg10),
                    add_symbol=True,description='unit of power')
    Unit.alias('dB(W)',_dBW)
    
    _dBk = LogUnit('decibel kilowatt','dBk',LogConversion(gummy(1,unit='kW'),10,10,_lg10),
                    add_symbol=True,description='unit of power')
    Unit.alias('dB(kW)',_dBk)
    
    _dBsm = LogUnit('decibel square metre','dBsm',LogConversion(gummy(1,unit=_m**2),10,10,_lg10),
                    add_symbol=True,description='unit used to measure antenna effective area')
    Unit.alias('dB(m**2)',_dBsm)
    Unit.alias('decibel square meter',_dBsm)
    
    _dBmm1 = LogUnit('decibel reciprocal metre','dB(m\u207b\u00b9)',LogConversion(gummy(1,unit=_m**-1),10,10,_lg10),
                    add_symbol=True,html_symbol='dB(m<sup>-1</sup>)',
                    latex_symbol='\t\t\\mathrm{dB}(\\mathrm{m}^{-1})',
                    ascii_symbol='dB(m**-1)',
                    description='unit used to measure antenna factor')
    Unit.alias('decibel reciprocal meter',_dBmm1)
    
    _dBHz = LogUnit('decibel hertz','dB-Hz',LogConversion(gummy(1,unit=_Hz),10,10,_lg10),
                    add_symbol=True,description='frequency')
    Unit.alias('dB(Hz)',_dBHz)
    
    _dBK = LogUnit('decibel kelvin','dBK',LogConversion(gummy(1,unit=_K),10,10,_lg10),
                    add_symbol=True,description='unit used to measure noise temperature')
    Unit.alias('dB(K)',_dBK)
    
    _dBmm1 = LogUnit('decibel reciprocal kelvin','dB(K\u207b\u00b9)',LogConversion(gummy(1,unit=_K**-1),20,10,_lg10),
                    add_symbol=True,html_symbol='dB(K<sup>-1</sup>)',
                    latex_symbol='\t\t\\mathrm{dB}(\\mathrm{K}^{-1})',
                    ascii_symbol='dB(K**-1)',description='temperature')
    
    _mBm = LogUnit('millibel milliwatt','mBm',LogConversion(gummy(1,unit='mW'),1000,10,_lg10),
                    add_symbol=True,description='unit of power')
    Unit.alias('mB(mW)',_mBm)
    
    LogUnit('bel','B',LogConversion(1,1,10,_lg10),add_symbol=False)
    
    _bit = BinaryPrefixedUnit('bit','bit',order=0,add_symbol=True,
                              additional_names=('shannon',),
                              description='unit of information')
    Unit.alias('shannon',_bit)
    
    _nat = Unit('natural unit of information','nat',Conversion(_bit,1/np.log(2)),
         add_symbol=True,order=0,description='unit of information')
    Unit.alias('nit',_nat)
    Unit.alias('nepit',_nat)
    
    _Hart = Unit('hartley','Hart',Conversion(_bit,np.log2(10)),add_symbol=True,
                 order=0,description='unit of information')
    Unit.alias('ban',_Hart)
    Unit.alias('dit',_Hart)
    
    BinaryPrefixedUnit('byte','B',Conversion(_bit,8),add_symbol=True,order=0,
                       description='unit of information')
    
    _nibble = Unit('nibble','nibble',Conversion(_bit,4),add_symbol=True,order=0,
                   description='unit of information')
    Unit.alias('nybble',_nibble)
    Unit.alias('nyble',_nibble)
    
    PrefixedUnit('torr','Torr',Conversion(_Pa,MFraction(101325,760)),add_symbol=True,
                         order=0,prefixes=['milli'],description='unit of pressure')
    Unit('standard atmosphere','atm',Conversion(_Pa,101325),add_symbol=True,order=0,description='unit of pressure')
    PrefixedUnit('electronvolt','eV',Conversion(_J,_const_e),add_symbol=True,order=0,description='unit of energy, the energy gained by one electron moving across one volt')
                 
    # astronomical units
    _c = Unit('speed of light','c',Conversion(_m/_s,_const_c),add_symbol=True,order=0,
         html_symbol='<i>c</i>',latex_symbol='\t\tc',
         description='natural unit of velocity')
    _au = Unit('astronomical unit','au',Conversion(_m,149597870700),add_symbol=True,
               order=0,description='astronomical unit of length')
    Unit.alias('ua',_au)
    PrefixedUnit('parsec','pc',Conversion(_au,648000/pi),add_symbol=True,order=0,
                 prefixes=['kilo','mega','giga'],
                 description='astronomical unit of length')
    _a = PrefixedUnit('Julian year','a',Conversion(_d,MFraction('365.25')),add_symbol=True,order=0,
                       description='astronomical unit of time',prefixes=['kilo','mega',
                       'giga','tera'],additional_names=('annum','year'),
                        additional_short_names=('yr',))
    Unit('light year','ly',Conversion(_c*_a,1),add_symbol=True,order=0,
         description='astronomical unit of length')
    Unit('siriometer','Sm',Conversion(_au,1000000),add_unit=False,
         description='astronomical unit of length')
    Unit('light second','light-second',Conversion(_c*_s,1),add_symbol=True,
         order=0,description='unit of length',)
    Unit('light minute','light-minute',Conversion(_c*_min,1),add_symbol=True,
         order=0,description='astronomical unit of length')
    Unit('light hour','light-hour',Conversion(_c*_h,1),add_symbol=True,
         order=0,description='astronomical unit of length')
    _M_solar = Unit('solar mass','M\u2609',Conversion(_kg,_const_solar_mass),
                    add_symbol=True,order=0,html_symbol='<i>M</i><sub>&#x2609;</sub>',
                    latex_symbol='\t\tM_{\odot}',ascii_symbol='M(solar)',
                    description='astronomical unit of mass, mass of sun')
    Unit.alias('M(Sun)',_M_solar)
    _M_J = Unit('Jupiter mass','M(J)',Conversion(_kg,_const_jupiter_mass),add_symbol=True,
         order=0,html_symbol='<i>M</i><sub>J</sub>',latex_symbol='\t\tM_{\mathrm{J}}',
         description='astronomical unit of mass') #check
    Unit.alias('M(Jup)',_M_J)
    Unit('Earth mass','M\u2295',Conversion(_kg,_const_earth_mass),add_symbol=True,
         order=0,html_symbol='<i>M</i><sub>&#x2295;</sub>',latex_symbol='\t\tM_{\oplus}',
         ascii_symbol='M(E)',description='astronomical unit of mass')
    _Jy = Unit('jansky','Jy',Conversion(_W*_m**-2*_Hz**-1,MFraction('1e-26')),add_symbol=True,
         decription='astronomical unit, spectral flux density, spectral irradiance')
    LogUnit('monochromatic AB magnitude','m(AB)',
            LogConversion(gummy(3631,unit=_Jy),-2.5,10,_lg10),
            add_symbol=True,html_symbol='<i>m</i><sub>AB</sub>',
            latex_symbol='\t\tm_{\\mathrm{AB}}',
            description='astronomical unit, spectral flux density, spectral irradiance')
    
    # natural units
    _hbar = Unit('natural unit of action','\u210f',Conversion(_J*_s,_const_hbar),
         order=0,add_symbol=True,ascii_symbol='h-bar')
    _e = Unit('elementary charge','e',Conversion(_C,_const_e),
         order=0,add_symbol=True,html_symbol='<i>e</i>',latex_symbol='\t\te',
         description='natural unit of charge, natural unit, atomic unit of charge')
    _a0 = Unit('bohr','a0',Conversion(_m,_const_a0),
         order=0,add_symbol=True,html_symbol='<i>a</i><sub>0</sub>',latex_symbol='\t\ta_{0}',
         description='natural unit of length, atomic unit of length, Bohr radius')
    Unit.alias('a(0)',_a0)
    
    _me = Unit('electon mass','m(e)',Conversion(_kg,_const_me),
         order=0,add_symbol=True,html_symbol='<i>m</i><sub>e</sub>',latex_symbol='\t\tm_{\\mathrm{e}}',
         description='natural, atomic unit of mass')
    _Eh = Unit('hartree','E_h',Conversion(_hbar**2*_me**-1*_a0**-2,1),order=0,
               add_symbol=True,html_symbol='<i>E</i><sub>h</sub>',latex_symbol='\t\tE_{\\mathrm{h}}',
               description='natural, atomic unit of energy')
    Unit.alias('Ha',_Eh)
    Unit('reduced compton wavelength','\u019b(C)',
         Conversion(_hbar*_me**-1*_c**-1,1),order=0,add_symbol=True,
         html_symbol='\u019b<sub>C</sub>',latex_symbol='\t\t\u019b_{\\mathrm{C}}',
         ascii_symbol='lambda(C)',
         description='electron reduced compton wavelength, natural unit of length')
    Unit('proton mass','m(p)',Conversion(_kg,_const_proton_mass),
         order=0,add_symbol=True,html_symbol='<i>m</i><sub>p</sub>',latex_symbol='\t\tm_{\\mathrm{p}}')

    Unit('Planck length','l(P)',Conversion(_m,sqrt(_const_hbar*_const_G*_const_c**-3)),order=0,
         add_symbol=True,html_symbol='<i>l</i><sub>P</sub>',
         latex_symbol='\t\tl_{\\mathrm{P}}',description='natural unit of length')
    Unit('Planck mass','m(P)',Conversion(_kg,sqrt(_const_hbar*_const_G**-1*_const_c)),order=0,
         add_symbol=True,html_symbol='<i>m</i><sub>P</sub>',
         latex_symbol='\t\tm_{\\mathrm{P}}',description='natural unit of mass')
    Unit('Planck time','t(P)',Conversion(_s,sqrt(_const_hbar*_const_G*_const_c**-5)),order=0,
         add_symbol=True,html_symbol='<i>t</i><sub>P</sub>',
         latex_symbol='\t\tt_{\\mathrm{P}}',description='natural unit of time')
    Unit('Planck charge','q(P)',Conversion(_C,_const_e/sqrt(_const_alpha)),order=0,
         add_symbol=True,html_symbol='<i>q</i><sub>P</sub>',
         latex_symbol='\t\tq_{\\mathrm{P}}',description='natural unit of charge')
    Unit('Planck temperature','T(P)',Conversion(_K,sqrt(_const_hbar*_const_G**-1*_const_c**5*_const_k**-2)),order=0,
         add_symbol=True,html_symbol='<i>T</i><sub>P</sub>',
         latex_symbol='\t\tT_{\\mathrm{P}}',description='natural unit of temperature')
    
    PrefixedUnit('volt ampere reactive','var',Conversion(_V*_A,1),add_symbol=True,
                 order=0,description='volt-ampere reactive power')
    PrefixedUnit('volt ampere','VA',Conversion(_V*_A,1),add_symbol=True,
                 order=0,description='volt-ampere apparent power')
    
    _M_W = LogUnit('moment magnitude','M(W)',
                   LogConversion(gummy(1,unit=_dyn*_cm),MFraction(2,3),10,_lg10,offset=-10.7),
                   add_symbol=True,html_symbol='<i>M</i><sub>W</sub>',
                   latex_symbol='\\textit{M}_{W}',
                   desription='moment magnitude scale, earthquake intensity')
    Unit.alias('MMS',_M_W)
    
    Unit.alias('wavenumber',_cm**-1)
    
    Unit('root hertz','\u221aHz',Conversion(_Hz**0.5,1),add_symbol=True,order=0,
         html_symbol='&radic;<span style="text-decoration:overline;">Hz</span>',
         latex_symbol='\t\t\\sqrt{\\mathrm{Hz}}',ascii_symbol='sqrtHz')
    
