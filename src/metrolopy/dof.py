# -*- coding: utf-8 -*-

# module dof

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

from math import isnan

class DoF:
    __slots__ = 'value'

    def __init__(self,value):
        if not isinstance(value,int):
            value = float(value)
        if value <= 0:
            raise ValueError('dof value ' + str(value) + ' is less than or equal to zero')
        if isnan(value):
            raise ValueError('dof value is NaN')
        self.value = value

    def __repr__(self):
        return 'DoF(' + str(self.value) + ')'
    
_DoF_inf = DoF(float('inf'))