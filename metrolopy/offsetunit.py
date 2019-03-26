# -*- coding: utf-8 -*-

# module offsetunit

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
The OffsetUnit class was created to handle Celsius and Fahrenheit units.
"""

from .unit import Unit,Conversion
from .nonlinearunit import NonlinearUnit,NonlinearConversion
from .exceptions import IncompatibleUnitsError

class OffsetConversion(NonlinearConversion):
    """
    Represents a conversion of the form x -> x + offset.
    """
    def __init__(self,unit,offset):
        self._unit = unit
        self.offset = offset
        
    def _to(self,g):
        return g + self.offset
    
    def _frm(self,g):
        return g - self.offset
        
    def copy(self):
        r = OffsetConversion(self._unit,self.offset)
        r.parent = self.parent
        return r
        

class _IntervalUnit(Unit):
    
    def __init__(self,parent,*params,**kwds):
        conv = Conversion(parent.conversion.unit,1)
        if 'short_name' in kwds:
            kwds['short_name'] += '-i'
        elif kwds.get('add_symbol'):
            if 'ascii_symbol' in kwds:
                kwds['short_name'] = kwds['ascii_symbol'] + '-i'
            else:
                kwds['short_name'] = params[1] + '-i'
        kwds['add_symbol'] = False
        
        kwds['parent'] = parent
        Unit.__init__(self,params[0] + ' interval',params[1],conversion=conv,**kwds)
    
    def _getme(self,ul,e):
        return self.parent
    
class OffsetUnit(NonlinearUnit):
    """
    This class was created to handle units such as the degree Celsius and the
    degree Fahrenheit.  This class takes the same parameters as the Unit class,
    but actually creates two unit instances...
    """
    def __init__(self,*params,**kwds):            
        Unit.__init__(self,*params,**kwds)
        
        self.interval_unit = _IntervalUnit(self,*params,**kwds)
        
    def zero(self):
        return -self.conversion.offset
        
    def _add_alias(self,alias):
        self._aliases.add(alias)
        Unit.alias(alias + '-i',self.interval_unit)
        
    def _add(self,a,bunit,b,aconv):
        if bunit is not self.interval_unit:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may only be added to its interval unit counterpart')
        return (self._nummy_add(a,b),self)
    
    def _radd(self,a,bunit,b,aconv):
        if bunit is not self.interval_unit:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may only be added to its interval unit counterpart')
        return (self._nummy_radd(a,b),self)
    
    def _sub(self,a,bunit,b,aconv):
        if bunit is self:
            return (self._nummy_sub(a,b),self.interval_unit)
        if bunit is self.interval_unit:
            return (self._nummy_sub(a,b),self)
        raise IncompatibleUnitsError('a quantity with unit ' + bunit.tostring() + ' may not be subtracted from a quantity with unit ' + self.tostring() + '; automatic conversion is disabled with offset unit instances')
    
    def _rsub(self,a,bunit,b,aconv):
        if bunit is self:
            return (self._nummy_rsub(a,b),self.interval_unit)
        raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be subtracted from a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with offset unit instances')

        