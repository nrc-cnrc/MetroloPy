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
<<<<<<< HEAD
    """
    Represents the degrees of freedom for an uncertainty.  ummy or gummy
    instances that share the same DoF instance are assumed to have uncertainties
    based on the samples from the same underlying normal distribution.
    
    For example if we define:
        
         a = gummy(1.2,u=3.3,dof=5)
         b = gummy(0.7,u=3.3,dof=5)
         
    we are assuming that, even though both a and b have the number of degrees 
    of freedom, their uncertainites are esimated from samples from two 
    different distributions that may not have exaclty the same population 
    variance. However if we define:
        
        c = gummy(1.2,u=3.3,dof=DoF(5))
        d = gummy(0.7,u=2.4,dof=Dof(5))
        
    we are now assuming that the uncertanties are both estimated using the 
    same samples from the same distribution.  E.g. c and d may represent
    the fitted parameters of a least squares regression where the variance
    of the same data points is used to estimate the uncertainty of both c and
    d.
    """
    
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    __slots__ = '_value'

    def __init__(self,value):
        if not isinstance(value,int):
            value = float(value)
        if value <= 0:
            raise ValueError('dof value ' + str(value) + ' cannot be less than or equal to zero')
        if isnan(value):
            raise ValueError('dof value cannot be NaN')
        self._value = value
        
    @property
    def value(self):
        return self._value

    def __repr__(self):
        return 'DoF(' + str(self._value) + ')'

    
_DoF_inf = DoF(float('inf'))