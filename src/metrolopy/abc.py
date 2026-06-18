# -*- coding: utf-8 -*-

# module abc

# Copyright (C) 2026 National Research Council Canada
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

from abc import abstractmethod
from numbers import Number

class UncertainValue(Number):
    @property
    @abstractmethod
    def x(self):
        pass
    
    @property
    @abstractmethod
    def u(self):
        pass
    
    @abstractmethod
    def tostring(fmt=None,**kwds):
        pass
    
class UncertainComplexValue(UncertainValue):
    pass
    
class AbcQuantity(Number):
    autoconvert = False
    
    @property
    @abstractmethod
    def value(self):
        pass
    
    @property
    @abstractmethod
    def unit(self):
        pass
    
    @property
    @abstractmethod
    def unit_is_one(self):
        pass
    
    @property
    @abstractmethod
    def c(self):
        pass
    
    @abstractmethod
    def convert(self,unit):
        pass
    
    @abstractmethod
    def tostring(fmt=None,**kwds):
        pass
    
    @abstractmethod
    def copy(self,unit):
        pass
    
    @abstractmethod
    def tofloat(self,unit):
        pass
    
    @abstractmethod
    def tobaseunit(self,unit):
        pass
    
    @abstractmethod
    def totuple(self,unit):
        pass
    
    @abstractmethod
    def splonk(self,unit):
        pass
    
    @abstractmethod
    def _cmp(self,v):
        pass
    
class AbcQuantityArray(AbcQuantity):
    pass
    
class AbcUnit(Number):
    @property
    @abstractmethod
    def conversion(self):
        pass
    
    @property
    @abstractmethod
    def is_dimensionless(self):
        pass
    
    @property
    @abstractmethod
    def linear(self):
        pass
    
    @property
    @abstractmethod
    def base(self):
        pass
    
    @property
    @abstractmethod
    def units(self):
        pass
    
    @abstractmethod
    def convert(self,unit):
        pass
    
    @abstractmethod
    def _mul(self,v):
        pass
    
    @abstractmethod
    def _rmul(self,v):
        pass
    
    @abstractmethod
    def _truediv(self,v):
        pass
    
    @abstractmethod
    def _rtruediv(self,v):
        pass
    
    @abstractmethod
    def _mul_cancel(self,v):
        pass
    
    @abstractmethod
    def tostring(fmt=None,**kwds):
        pass