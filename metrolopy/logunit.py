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
from .gummy import gummy
from .unit import one,_CompositeUnit
from .nonlinearunit import NonlinearConversion,NonlinearUnit
from .exceptions import IncompatibleUnitsError

class LogConversion(NonlinearConversion):   
    def __init__(self,reference,multiplier,log_base,log_func,offset=0):
        if isinstance(reference,gummy):
            self._unit = reference._unit
            self._rf = gummy(reference.x,reference.u)
            self._rf._ref = reference._ref
            self._rf._refs = reference._refs
        else:
            self._unit = one
            self._rf = reference
    
        self.reference = reference
        
        self.multiplier = multiplier           
        self.log_base = log_base
        self.log_func = lambda x: log_func(x) if x > 0 else -float('inf')
        self._lnbase = np.log(float(log_base))
        self.offset = offset
        
    def _to(self,g):
        g = (g - self.offset)/self.multiplier
        def f(x):
            if x == -float('inf'):
                return 0.0
            try:
                return self.log_base**x
            except OverflowError:
                return float('inf')
        def d(x):
            try:
                return self._lnbase*self.log_base**x
            except OverflowError:
                return 1
        ret = self._rf*gummy.apply(f,d,g)
        return ret
        
        
    def _frm(self,g):
        g = g/self._rf
        def f(x):
            return self.log_func(x)
        def d(x):
            if x != 0:
                return 1/(x*self._lnbase)
            return 1
        ret = self.multiplier*gummy.apply(f,d,g) + self.offset
        return ret
        
    def copy(self):
        r = LogConversion(self.reference,self.multiplier,self.log_base,self.log_func)
        r.parent = self.parent
        return r
    
    
class LogUnit(NonlinearUnit):   
    def __init__(self,*p,**kwds):
        super().__init__(*p,**kwds)
        self.reference = self.conversion.reference
        if isinstance(self.reference,gummy):
            self.referencef = self.reference.x
        else:
            self.referencef = self.reference
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
        return (self._nummy_add(a,b),self)
    
    def _radd(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (self._nummy_radd(a,b),self)
    
    def _sub(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + bunit.tostring() + ' may not be subtracted from a quantity with unit ' + self.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (self._nummy_sub(a,b),self)
    
    def _rsub(self,a,bunit,b,aconv):
        if bunit is not self:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be subtracted from a quantity with unit ' + bunit.tostring() + '; automatic conversion is disabled with LogUnit instances')
        return (self._nummy_rsub(a,b),self)
    
    def _mul(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return (self._nummy_mul(self._nummy_mul(a,b),c),un)
    
    def _rmul(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return (self._nummy_mul(self._nummy_rmul(a,b),c),un)
    
    def _truediv(self,a,bunit,b,aconv):
        if not bunit.linear:
            raise IncompatibleUnitsError('only quantities with linear units may multiply or divide a quanity with unit ' + self.tostring())
        
        if aconv:
            un,c = self._div_cancel(bunit)
        else:
            un,c = bunit._rdiv_cancel(self)
        
        return (self._nummy_mul(self._nummy_truediv(a,b),c),un)
    
    def _neg(self,a):
        return (self._nummy_neg(a),self)
    
    def _pos(self,a):
        return (self._nummy_pos(a),self)
    
    def _abs(self,a):
        return (self._nummy_abs(a),self)
        
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