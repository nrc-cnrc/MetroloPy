# -*- coding: utf-8 -*-

# module usunits

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
imported directly.  The gummy.siunits module must be loaded before loading 
this module.

Most units here are from:
    
NIST Special Publication 1038, "The International System of Units (SI) – 
Conversion Factors for General Use", May 2006.
"""

from .unit import _BuiltinLib,Unit,Conversion
from .offsetunit import OffsetUnit,OffsetConversion
from .prefixedunit import PrefixedUnit
from .ummy import MFraction

with _BuiltinLib():
    _in = PrefixedUnit('inch','in',Conversion('m',MFraction('0.0254')),prefixes=['micro'],
                       add_symbol=True,description='unit of length')
    Unit('hand','hand',Conversion(_in,4),add_symbol=True,description='unit of length')
    _ft = Unit('foot','ft',Conversion(_in,12),add_symbol=True,description='unit of length')
    _yd = Unit('yard','yd',Conversion(_ft,3),add_symbol=True,description='unit of length')
    Unit('mile','mi',Conversion(_yd,1760),add_symbol=True,description='unit of length')
    _pica = Unit('pica','P/',Conversion(_in,MFraction(1,6)),add_symbol=True,description='unit of length')
    Unit('point','p',Conversion(_pica,MFraction(1,12)),add_symbol=True,description='unit of length')
    Unit('link','li',Conversion(_ft,MFraction(33,50)),add_symbol=True,description='unit of length')
    Unit('survey foot','ft',Conversion('m',MFraction(1200,3937)),add_symbol=False,description='unit of length')
    _rod = Unit('rod','rd',Conversion(_ft,25),add_symbol=True,description='unit of length')
    _ch = Unit('chain','ch',Conversion(_rod,4),add_symbol=True,description='unit of length')
    _fur = Unit('furlong','fur',Conversion(_ch,10),add_symbol=True,description='unit of length')
    _mi = Unit('survey mile','mi',Conversion(_fur,8),add_symbol=False,description='unit of length')
    Unit.alias('statute mile',_mi)
    Unit('league','lea',Conversion(_mi,3),add_symbol=True,description='unit of length')
    _ftm = Unit('fathom','ftm',Conversion(_yd,2),add_symbol=True,description='unit of length')
    Unit('cable','cb',Conversion(_ftm,120),add_symbol=True,description='unit of length')
    _mil = Unit('thousandth of an inch','mil',Conversion(_in,MFraction(1,1000)),add_symbol=True,description='unit of length')
    Unit.alias('thou',_mil)
    Unit.alias('thousandth',_mil)
          
    _acre = Unit('acre','acre',Conversion(_ch**2,10),add_symbol=True,description='unit of area')
    _section = Unit('section','section',Conversion(_acre,640),add_symbol=True,description='unit of area')
    Unit('survey township','twp',Conversion(_section,36),add_symbol=True,description='unit of area')
    
    Unit('cubic inch','cu in',Conversion(_in**3,1),add_symbol=True,
         short_name='cu-in',description='unit of volume')
    Unit('cubic foot','cu ft',Conversion(_ft**3,1),add_symbol=True,
         short_name='cu-ft',description='unit of volume')
    Unit('cubic yard','cu yd',Conversion(_yd**3,1),add_symbol=True,
         short_name='cu-yd',description='unit of volume')
    Unit('acre-foot','acre ft',Conversion(_ft**3,43560),add_symbol=True,
         short_name='acre-ft',description='unit of volume')
    
    _gal = Unit('gallon','gal',Conversion(_in**3,231),add_symbol=True,description='unit of volume')
    Unit.alias('liquid gallon',_gal)
    Unit.alias('liquid gal',_gal)
    
    _qt = Unit('quart','qt',Conversion(_gal,MFraction(1,4)),add_symbol=True,description='unit of volume')
    Unit.alias('liquid quart',_qt)
    Unit.alias('liquid qt',_qt)
    
    _pt = Unit('pint','pt',Conversion(_qt,MFraction(1,2)),add_symbol=True,description='unit of volume')
    Unit.alias('liquid pint',_pt)
    Unit.alias('liquid pt',_pt)
    
    _cp = Unit('cup','cp',Conversion(_pt,MFraction(1,2)),add_symbol=True,description='unit of volume')
    _gi = Unit('gill','gi',Conversion(_cp,MFraction(1,2)),add_symbol=True,description='unit of volume')
    _floz = Unit('fluid ounce','fl oz',Conversion(_gi,MFraction(1,4)),add_symbol=True,
                 short_name='fl-oz',description='unit of volume')
    _tbsp = Unit('tablespoon','Tbsp',Conversion(_floz,MFraction(1,2)),add_symbol=True,description='unit of volume')
    _tsp = Unit('teaspoon','tsp',Conversion(_tbsp,MFraction(1,3)),add_symbol=True,description='unit of volume')
    _minum = Unit('minim','min',Conversion(_tsp,MFraction('0.125')),add_symbol=False,description='unit of volume')
    Unit('fluid dram','fl dr',Conversion(_minum,60),add_symbol=True,
         short_name='fl-dr',description='unit of volume')
    Unit('shot','jig',Conversion(_tbsp,3),add_symbol=True,description='unit of volume')
    
    _bbl = Unit('barrel','bbl',Conversion(_gal,MFraction('31.5')),add_symbol=True,description='unit of volume')
    Unit.alias('liquid barrel',_bbl)
    Unit.alias('liquid bbl',_bbl)
    
    Unit('oil barrel','bbl',Conversion(_gal,42),add_symbol=False,description='unit of volume')
    Unit('hogshead','hogshead',Conversion(_gal,65),add_symbol=True,description='unit of volume')
    
    _dgal = Unit('dry gallon','gal',Conversion(_in**3,MFraction('268.8025')),add_symbol=False,
         short_name='dry-gal',description='unit of volume')
    _dqt = Unit('dry quart','qt',Conversion(_dgal,MFraction(1,4)),add_symbol=False,
         short_name='dry-qt',description='unit of volume')
    Unit('dry pint','pt',Conversion(_dqt,MFraction(1,2)),add_symbol=False,
         short_name='dry-pt',description='unit of volume')
    _pk = Unit('peck','pk',Conversion(_dgal,2),add_symbol=True,description='unit of volume')
    Unit('bushel','bu',Conversion(_pk,4),add_symbol=True,description='unit of volume')
    Unit('dry barrel','bu',Conversion(_in**3,7056),add_symbol=False,
         short_name='dry-bbl',description='unit of volume')
         
    _lb = PrefixedUnit('pound','lb',Conversion('kg',MFraction('0.45359237')),add_symbol=True,
                       prefixes=['micro'],description='unit of mass')
    Unit.alias('avoirdupois pound',_lb)
    Unit.alias('avdp lb',_lb)
    Unit.alias('lbm',_lb)
    Unit.alias('pound mass',_lb)
    Unit.alias('ulbm','ulb')
    Unit.alias('micropound mass','ulb')
    
    _oz = Unit('ounce','oz',Conversion(_lb,MFraction('0.0625')),add_symbol=True,description='unit of mass')
    Unit.alias('avoirdupois ounce',_oz)
    Unit.alias('avdp oz',_oz)
    
    Unit('dram','dr',Conversion(_oz,MFraction('0.0625')),add_symbol=True,description='unit of mass')
    
    _gr = Unit('grain','gr',Conversion(_lb,MFraction(1,7000)),add_symbol=True,description='unit of mass')
    
    Unit('hundredweight','cwt',Conversion(_lb,100),add_symbol=True,description='unit of mass')
    Unit('long hundredweight','long cwt',Conversion(_lb,112),add_symbol=True,description='unit of mass')
    Unit('short ton','tn',Conversion(_lb,2000),add_symbol=True,description='unit of mass')
    Unit('long ton','long ton',Conversion(_lb,2240),add_symbol=False,description='unit of mass')
    
    _dwt = Unit('pennyweight','dwt',Conversion(_gr,24),add_symbol=True,description='unit of mass')
    _toz = Unit('troy ounce','oz t',Conversion(_dwt,20),add_symbol=True,
         short_name='oz-t',description='unit of mass')
    Unit('troy pound','lb t',Conversion(_toz,12),add_symbol=True,
         short_name='lb-t',description='unit of mass')
         
    _lbf = Unit('pound force','lbf',Conversion('lb m s**-2',MFraction('9.80665')),add_symbol=True,description='unit of force')
    Unit('slug','slug',Conversion('lbf s**2 ft**-1',1),add_symbol=True,description='unit of force')
    
    _degR = Unit('degree Rankine','\u00B0R',Conversion('K',MFraction(5,9)),
                 latex_symbol='^{\circ}R',ascii_symbol='degR',add_symbol=True,
                 description='unit of temperature')
    Unit.alias('degree R',_degR)
    Unit.alias('deg R',_degR)
    
    _degF = OffsetUnit('degree Fahrenheit','\u00B0F',OffsetConversion('degR',MFraction('459.67')),
                latex_symbol='^{\circ}F',ascii_symbol='degF',add_symbol=True,
                description='unit of temperature')
    Unit.alias('degree F',_degF)
    Unit.alias('deg F',_degF)
               
    Unit('pound per square inch','psi',Conversion(_lbf*_in**-2),add_symbol=True,
         description='unit of pressure')
               
    Unit('board-foot','board-foot',Conversion(_ft**2*_in,1),add_symbol=True,
         short_name='board-ft',description='unit of volume')
    
    _cal = Unit('calorie','cal',Conversion('J',MFraction('4.184')),add_symbol=True,
                description='unit of energy')
    Unit.alias('thermochemical calorie',_cal)
    Unit.alias('calth',_cal)
    
    _Cal = Unit('large calorie','Cal',Conversion(_cal,1000),add_symbol=True,
                description='unit of energy')
    Unit.alias('kilocalorie',_Cal)
    Unit.alias('kcal',_Cal)
    Unit.alias('dietary calorie',_Cal)
    Unit.alias('Calorie',_Cal)
    
    _btu = Unit('IT British thermal unit','BTU',Conversion('J',MFraction('1055.05585262')),
                add_symbol=True,description='unit of energy')
    Unit.alias('Btu',_btu)
   
    Unit('horsepower','hp',Conversion('lbf ft s**-1',550), add_symbol=True,
         description='unit of power')
    
    _ftlb = Unit('foot-pound','ft\u00B7lb',Conversion(_ft*_lbf,1),add_symbol=True,
        ascii_symbol='ft-lb',description='unit of work or energy')
    Unit.alias('ft\u00B7lbf',_ftlb)
    Unit.alias('ft-lbf',_ftlb)
    
    _lbft = Unit('pound-foot','lb\u00B7ft',Conversion(_ft*_lbf,1),add_symbol=True,
        ascii_symbol='lb-ft',description='unit of torque')
    Unit.alias('lbf\u00B7ft',_lbft)
    Unit.alias('lbf-lft',_lbft)
    
    Unit('rack unit','U',Conversion(_in,MFraction('1.75')),add_symbol=False,
         description='unit of length')
    Unit('square','square',Conversion(_ft**2,100),add_symbol=False,
         description='unit of area used in construction')
    Unit('gasoline gallon equivalent','gasoline-gallon-equivalent',
         Conversion('kW h',33.7),add_symbol=True,description='unit of energy')
    _ttnt = PrefixedUnit('tons of TNT equivalent','t(TNT)',Conversion('GJ',4.184),add_symbol=True,
         description='unit of energy in an explosion',prefixes=['kilo','mega','giga'])
    Unit.alias('tons of TNT',_ttnt)
    Unit.alias('tons of tnt',_ttnt)
