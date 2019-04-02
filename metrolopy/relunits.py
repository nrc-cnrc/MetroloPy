# -*- coding: utf-8 -*-

# module relunits

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
imported directly.  Dimensionless units are defined here.
"""
from .unit import _BuiltinLib,Unit,Conversion,one
from .ummy import MFraction

class RatioUnit(Unit):
    """RatioUnit is used for dimensionless units like % where powers, e.g. %**2,
    are not desired.
    """
    def _cpow(self,v):
        if v == 1 or v == -1:
            return ([[self,v]],1)
        c = self.convert(1,one)
        if int(v) == v:
            if v > 0:
                return ([[self,1]],c**(v - 1))
            else:
                return ([[self,1]],c**(v + 1))
        return ([[one,1]],c**v)

with _BuiltinLib():
    # The \t at the beggining of the symbol causes the space between the numerical
    # value and the unit symbol per the Chicago Manual of Style (but not the 
    # SI brochure).
    RatioUnit('percent','\t%',Conversion(one,MFraction('0.01')),short_name='%',
         latex_symbol='\t\\%',add_symbol=False)
    
    RatioUnit('part per million','ppm',Conversion(one,MFraction('1e-6')),add_symbol=True)
    RatioUnit('part per billion','ppb',Conversion(one,MFraction('1e-9')),add_symbol=True)
    RatioUnit('part per trillion','ppt',Conversion(one,MFraction('1e-12')),add_symbol=True)
    RatioUnit('part per quadrillion','ppq',Conversion(one,MFraction('1e-15')),add_symbol=True)
