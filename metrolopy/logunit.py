# -*- coding: utf-8 -*-

# module logunit

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
Unit and Conversion sub-classes are defined here to implement logarithmic 
units.
"""

import numpy as np
from .unit import one,_CompositeUnit,Quantity
from .nonlinearunit import NonlinearConversion,NonlinearUnit
from .exceptions import IncompatibleUnitsError

class LogConversion(NonlinearConversion):   
    def __init__(self,reference,multiplier,log_base,log_func,offset=0):
        if isinstance(reference,Quantity):
            self._unit = reference.unit
            self._rf = reference.value
        else:
            self._unit = one
            self._rf = reference
    
        self.reference = reference
        
        self.multiplier = multiplier           
        self.log_base = log_base
        self.log_func = lambda x: log_func(x) if x > 0 else -float('inf')
        self._lnbase = np.log(log_base)
        self.offset = offset
        
    def to(self,g):
        g = (g - self.offset)/self.multiplier
        def f(x):
            if x == -float('inf'):
                return 0.0
            try:
                return self.log_base**x
            except OverflowError:
                return float('inf')
        return self._rf*f(g)
        
    def frm(self,g):
        return self.multiplier*self.log_func(g/self._rf) + self.offset

    def copy(self):
        r = LogConversion(self.reference,self.multiplier,self.log_base,self.log_func)
        r.parent = self.parent
        return r
    
    
class LogUnit(NonlinearUnit):
    
    _ufunc_dict = {}
    
    def __init__(self,*p,**kwds):
        super().__init__(*p,**kwds)
        self.reference = self.conversion.reference
        self.multiplier = self.conversion.multiplier
        self.log_base = self.conversion.log_base
        
    def get_composite(self,ul):
        n = 0
        for u,e in ul:
            if not u.linear:
                n += 1
        if n > 1:
            raise IncompatibleUnitsError('unit ' + self.tostring() + ' may not be combined with other nonlinear units')
            
        return _LogCompositeUnit
        
    def zero(self):
        return -float('inf')
    
    def _add(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (a + b,self)
    
    def _radd(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (b + a,self)
    
    def _sub(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + bunit.tostring() + ' may not be subtracted from a quantity with unit ' + self.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (a - b,self)
    
    def _rsub(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be subtracted from a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (a + b,self)
    
    def _mul(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return ((a*b)*c,un)
    
    def _rmul(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return ((b*a)*c,un)
    
    def _truediv(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._div_cancel(bunit)
        else:
            un,c = bunit._rdiv_cancel(self)
        
        return ((a/b)*c,un)
    
    def _neg(self,a):
        return (-a,self)
    
    def _pos(self,a):
        return (+a,self)
    
    def _abs(self,a):
        return (abs(a),self)
        
    def __mul__(self,a):
        if not a.linear:
            raise IncompatibleUnitsError('only linear units may multiply or divide unit ' + self.tostring())
        return _CompositeUnit(self.units + a.units)
    
    def __truediv__(self,a):
        if not a.linear:
            raise IncompatibleUnitsError('only linear units may multiply or divide unit ' + self.tostring())
        vi = [(e[0],-e[1]) for e in a.units]
        return _CompositeUnit(self.units + vi)


class _LogCompositeUnit(LogUnit,_CompositeUnit):
    pass