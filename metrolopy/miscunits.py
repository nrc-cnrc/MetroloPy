# -*- coding: utf-8 -*-

# module miscunits

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

Miscellaneous unit definitions: astronomical units, natural units, some 
non-SI units accepted for use with the SI not in the siunits module, some 
obsolete SI units, and a few other units.
"""

from numpy import sqrt,log,log2
from .gummy import gummy,_lg10
from .ummy import MFraction
from .unit import Unit,Conversion
from .prefixedunit import PrefixedUnit,BinaryPrefixedUnit
from .logunit import LogUnit,LogConversion
from .siunits import (_kg,_V,_ohm,_s,_K,_J,_m,_d,_min,_h,_W,_Hz,_C,_A,_N,_cm,
                      _Wb,_T,_Pa)
from .constcom import (dalton,KJ,RK,KJ90,RK90,e,c,solar_mass,jupiter_mass,
                       earth_mass,hbar,a0,me,mp,G,alph,k,pi,euler)

with Unit._builtin():
    
    # non-SI accepted for use with the SI, SI Brochure, 9th ed
    _Da = PrefixedUnit('dalton','Da',Conversion(_kg,dalton),
                       add_symbol=True,description='unified atomic mass unit')
    
    PrefixedUnit('electronvolt','eV',Conversion(_J,e),add_symbol=True,order=0,
                 description='unit of energy, the energy gained by one electron moving across one volt')


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
            return log(x)
        except:
            return log(float(x))
    LogUnit('neper','Np',LogConversion(1,1,euler,_ln),add_symbol=True)


    # astronomical units
    _c = Unit('speed of light','c',Conversion(_m/_s,c),add_symbol=True,order=0,
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
    _M_solar = Unit('solar mass','M\u2609',Conversion(_kg,solar_mass),
                    add_symbol=True,order=0,html_symbol='<i>M</i><sub>&#x2609;</sub>',
                    latex_symbol='\t\tM_{\\odot}',ascii_symbol='M(solar)',
                    description='astronomical unit of mass, mass of sun')
    Unit.alias('M(Sun)',_M_solar)
    _M_J = Unit('Jupiter mass','M(J)',Conversion(_kg,jupiter_mass),add_symbol=True,
         order=0,html_symbol='<i>M</i><sub>J</sub>',latex_symbol='\t\tM_{\\mathrm{J}}',
         description='astronomical unit of mass') #check
    Unit.alias('M(Jup)',_M_J)
    Unit('Earth mass','M\u2295',Conversion(_kg,earth_mass),add_symbol=True,
         order=0,html_symbol='<i>M</i><sub>&#x2295;</sub>',latex_symbol='\t\tM_{\\oplus}',
         ascii_symbol='M(E)',description='astronomical unit of mass')
    _Jy = Unit('jansky','Jy',Conversion(_W*_m**-2*_Hz**-1,MFraction('1e-26')),add_symbol=True,
         decription='astronomical unit, spectral flux density, spectral irradiance')
    LogUnit('monochromatic AB magnitude','m(AB)',
            LogConversion(gummy(3631,unit=_Jy),-2.5,10,_lg10),
            add_symbol=True,html_symbol='<i>m</i><sub>AB</sub>',
            latex_symbol='\t\tm_{\\mathrm{AB}}',
            description='astronomical unit, spectral flux density, spectral irradiance')
    
    # natural units
    _hbar = Unit('natural unit of action','\u210f',Conversion(_J*_s,hbar),
         order=0,add_symbol=True,ascii_symbol='h-bar')
    _e = Unit('elementary charge','e',Conversion(_C,e),
         order=0,add_symbol=True,html_symbol='<i>e</i>',latex_symbol='\t\te',
         description='natural unit of charge, natural unit, atomic unit of charge')
    _a0 = Unit('bohr','a0',Conversion(_m,a0),
         order=0,add_symbol=True,html_symbol='<i>a</i><sub>0</sub>',latex_symbol='\t\ta_{0}',
         description='natural unit of length, atomic unit of length, Bohr radius')
    Unit.alias('a(0)',_a0)
    
    _me = Unit('electon mass','m(e)',Conversion(_kg,me),
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
    Unit('proton mass','m(p)',Conversion(_kg,mp),
         order=0,add_symbol=True,html_symbol='<i>m</i><sub>p</sub>',latex_symbol='\t\tm_{\\mathrm{p}}')

    Unit('Planck length','l(P)',Conversion(_m,sqrt(hbar*G*c**-3)),order=0,
         add_symbol=True,html_symbol='<i>l</i><sub>P</sub>',
         latex_symbol='\t\tl_{\\mathrm{P}}',description='natural unit of length')
    Unit('Planck mass','m(P)',Conversion(_kg,sqrt(hbar*G**-1*_c)),order=0,
         add_symbol=True,html_symbol='<i>m</i><sub>P</sub>',
         latex_symbol='\t\tm_{\\mathrm{P}}',description='natural unit of mass')
    Unit('Planck time','t(P)',Conversion(_s,sqrt(hbar*G*c**-5)),order=0,
         add_symbol=True,html_symbol='<i>t</i><sub>P</sub>',
         latex_symbol='\t\tt_{\\mathrm{P}}',description='natural unit of time')
    Unit('Planck charge','q(P)',Conversion(_C,e/sqrt(alph)),order=0,
         add_symbol=True,html_symbol='<i>q</i><sub>P</sub>',
         latex_symbol='\t\tq_{\\mathrm{P}}',description='natural unit of charge')
    Unit('Planck temperature','T(P)',Conversion(_K,sqrt(hbar*G**-1*c**5*k**-2)),order=0,
         add_symbol=True,html_symbol='<i>T</i><sub>P</sub>',
         latex_symbol='\t\tT_{\\mathrm{P}}',description='natural unit of temperature')
    
    PrefixedUnit('volt ampere reactive','var',Conversion(_V*_A,1),add_symbol=True,
                 order=0,description='volt-ampere reactive power')
    PrefixedUnit('volt ampere','VA',Conversion(_V*_A,1),add_symbol=True,
                 order=0,description='volt-ampere apparent power')
    
    
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
    
    _nat = Unit('natural unit of information','nat',Conversion(_bit,1/log(2)),
         add_symbol=True,order=0,description='unit of information')
    Unit.alias('nit',_nat)
    Unit.alias('nepit',_nat)
    
    _Hart = Unit('hartley','Hart',Conversion(_bit,log2(10)),add_symbol=True,
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


    # 1990 conventional unit definitions
    _V90 = PrefixedUnit('volt 90','V(90)',Conversion(_V,KJ90/KJ),add_symbol=True,
                        html_symbol='V<sub>90</sub>',latex_symbol='V_{90}',
                        description='1990 conventional definition of voltage, electrical potential')
    _ohm90 = PrefixedUnit('ohm 90','\u03A9(90)',Conversion(_ohm,RK/RK90),add_symbol=True,
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

