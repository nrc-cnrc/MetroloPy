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
This module is loaded by the unit module and is not intended be be
imported directly.  

The most of the units here are from the SI Brochure, 9th edition.
"""

from numpy import pi
from .ummy import ummy,MFraction
from .unit import Conversion,Unit
from .prefixedunit import PrefixedUnit
from .offsetunit import OffsetUnit,OffsetConversion


with Unit._builtin():
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
                Conversion(_kg,ummy(1,1.2e-8)),add_symbol=True,
                html_symbol='<i>m</i>(&#x1d4a6;)',latex_symbol='m({\\mathcal {K}})',
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
    Unit.alias('Ohm',_ohm)
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
               latex_symbol='^{\\circ}C',ascii_symbol = 'degC',add_symbol=True,order=0,
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
        
    #Non-SI units accepted for use with the International System of Units
    _min = Unit('minute','min',Conversion(_s,60),add_symbol=True,order=0,
                description='unit of time')
    _h = Unit('hour','h',Conversion(_min,60),add_symbol=True,order=0,
              description='unit of time')
    _d = Unit('day','d',Conversion(_h,24),add_symbol=True,order=0,
              description='unit of time')
    Unit.alias('D',_d)
    _deg = Unit('degree','\t\u00B0',Conversion(_rad,pi/180),add_symbol=True,order=0,
                 latex_symbol='^{\\circ}',ascii_symbol='deg',
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