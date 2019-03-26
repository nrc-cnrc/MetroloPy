# -*- coding: utf-8 -*-

# module nonlinearunit

# Copyright (C) 2019 National Research Council Canada
# Author:  Harold Parks

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
This module defines the Nonlinear Unit and Nonlinear Conversion abstract 
classes that are the base class for the LogUnit and OffsetUnit classes.
"""

from .unit import Unit, Conversion
from .exceptions import IncompatibleUnitsError

class NonlinearConversion(Conversion):
    """
    Base class for non-linear conversions.
    """
    # self._unit must be set in __init__ and must be linear

    # self._parent is the owning unit and must be set in the owning unit's __init__
    
    linear = False
        
    def chain(self,c):
        return _ChainedConversion(self,c)
    
    def rchain(self,c):
        return _ChainedConversion(c,self)
        
    def copy(self):
        raise NotImplementedError('this is an "abstract" method and must be overriden in a subclass')
        
    def pow(self,e):
        if e != 1:
            raise ValueError()
        return self
        
    
class _ChainedConversion(NonlinearConversion):
    linear = False
    
    def __init__(self,c1,c2):
        self.parent = c1.parent
        self._unit = c2.unit
        
        if isinstance(c1,_ChainedConversion):
            frmlst = c1._frmlst
            tolst = c1._tolst
        elif not c1.linear:
            frmlst = [c1._frm]
            tolst = [c1._to]
        else:
            frmlst = []
            tolst = []
        
        if isinstance(c2,_ChainedConversion):
            frmlst = c2._frmlst + frmlst
            tolst += c2._tolst
        elif not c2.linear:
            frmlst = [c1._frm] + frmlst
            tolst += [c1._to]
            
        self._frmlst = frmlst
        self._tolst = tolst
        
        factor = 1
        if c1.linear or isinstance(c1,_ChainedConversion):
            factor *= c1.factor
        if c2.linear or isinstance(c2,_ChainedConversion):
            factor *= c2.factor
        self.factor = factor
        
        self._c1 = c1
        self._c2 = c2
                    
    def to(self,g):
        g = super().to(g)
        return g*self.factor
        
    def frm(self,g):
        g = g/self.factor
        return super().frm(g)
        
    def _to(self,g):
        for f in self._tolst:
            g = f(g)
        return g
        
    def _frm(self,g):
        for f in self._frmlst:
            g = f(g)
        return g
        
    def chain(self,c):
        return _ChainedConversion(self,c)
    
    def copy(self):
        return _ChainedConversion(self._c1,self._c2)
    

class NonlinearUnit(Unit):
    """Base class of non-linear units.
    """
    linear = False
        
    #def _getme(self,unit_list,power):
        # add this method to return a different instance depending on whether
        # the unit is in a composite unit or not. e.g. an OffsetUnit where
        # 'degC/m' actually means 'degC-i/m'.
        # This us called only from Unit.unit when returning a composite unit
        # from a parsed string.  unit_list is the raw unit list returned from 
        # the parser an power is the exponent for this unit.
        
    #def _format_xu(g,fmt,style,norm,nsig,solidus=solidus,mulsep=mulsep):
        # add this method to override the gummy display formatting
        # g is a gummy and fmt is a gummy format
        # return a tuple: ('override','put the string to be displayed here')
        
    def zero(self):
        raise NotImplementedError('this is an "abstract" method and must be overriden in a subclass')
        
    def get_composite(self,ul):
        # If this method raises a NotImplementedError then the unit cannot be
        # combined with other units.  Override this method to return a subclass
        # of _CompositeUnit to allow combinations.
        raise NotImplementedError()
        
    def from_uunit(self,u,unit):
        raise IncompatibleUnitsError('uunit may not be used with unit ' + self.tostring())
        
    def to_uunit(self,u,unit):
        raise IncompatibleUnitsError('uunit may not be used with unit ' + self.tostring())
    
    # override any of the methods below to allow the operations
    def _add(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the + operation is not defined with unit ' + self.tostring())
        
    def _radd(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the + operation is not defined with unit ' + self.tostring())
        
    def _sub(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the - operation is not defined with unit ' + self.tostring())
        
    def _rsub(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the - operation is not defined with unit ' + self.tostring())
        
    def _mul(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the * operation is not defined with unit ' + self.tostring())
        
    def _rmul(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the * operation is not defined with unit ' + self.tostring())
        
    def _truediv(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the / operation is not defined with unit ' + self.tostring())
        
    def _rtruediv(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the / operation is not defined with unit ' + self.tostring())
        
    def _pow(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the ** operation is not defined with unit ' + self.tostring())
        
    def _rpow(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the ** operation is not defined with unit ' + self.tostring())
        
    def _mod(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the mod operation is not defined with unit ' + self.tostring())
        
    def _rmod(self,a,bunit,b,aconv):
        raise IncompatibleUnitsError('the mod operation is not defined with unit ' + self.tostring())
        
    def _neg(self,a):
        raise IncompatibleUnitsError('the - operation is not defined with unit ' + self.tostring())
    
    def _pos(self,a):
        raise IncompatibleUnitsError('the + operation is not defined with unit ' + self.tostring())
    
    def _abs(self,a):
        raise IncompatibleUnitsError('the abs operation is not defined with unit ' + self.tostring())
        
    def __pow__(self,a):
        raise IncompatibleUnitsError('the pow operation is not defined with unit ' + self.tostring())
        
    def __mul__(self,a):
        raise IncompatibleUnitsError('the * operation is not defined with unit ' + self.tostring())
    
    def __truediv__(self,a):
        raise IncompatibleUnitsError('the / operation is not defined with unit ' + self.tostring())
    
    def __rtruediv__(self,a):
        raise IncompatibleUnitsError('the reciprocal of unit ' + self.tostring() + ' is not allowed')
        
class ReciprocalConversion(NonlinearConversion):
    def __init__(self,conversion):
        self._conversion = conversion
        if hasattr(conversion,'factor'):
            self.factor = 1/conversion.factor
        self._unit = 1/conversion._unit
        self.parent = 1/conversion.parent
        
    def to(self,g):
        return self._conversion.frm(g)
        
    def frm(self,g):
        return self._conversion.to(g)
        
    def chain(self,c):
        return _ChainedConversion(self,c)
    
    def copy(self):
        return ReciprocalConversion(self._conversion)