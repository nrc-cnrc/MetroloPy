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
    
class UncertainComplexValue(UncertainValue):
    pass
    
class AbcQuantity(Number):
    @property
    @abstractmethod
    def value(self):
        pass
    
    @property
    @abstractmethod
    def unit(self):
        pass
    
class AbcQuantityArray(AbcQuantity):
    pass
    
class AbcUnit(Number):
    @property
    @abstractmethod
    def conversion(self):
        pass