# -*- coding: utf-8 -*-

# module mfraction

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

from decimal import Decimal
from fractions import Fraction


class MFraction(Fraction):
    """
    A fraction.Fraction sub-class that works with Decimal and mpmath.mpf objects
    """
    
    @classmethod
    def fromnum(cls,x):
        return MFraction(Fraction(x))
    
    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and not (isinstance(args[0],str) 
                                   or isinstance(args[0],Fraction)):
            return args[0]
        ret = super(MFraction, cls).__new__(cls, *args, **kwargs)
        if ret.denominator == 1:
            return ret.numerator
        return ret
    
    def _mpmath_(self,p,r):
        from mpmath import rational
        return rational.mpq(self.numerator,self.denominator)
    
    def todecimal(self):
        return Decimal(self.numerator)/Decimal(self.denominator)
    
    def __add__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__add__(v)
        return MFraction(super().__add__(v))
    
    def __radd__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__radd__(v)
        return MFraction(super().__radd__(v))
    
    def __sub__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__sub__(v)
        return MFraction(super().__sub__(v))
    
    def __rsub__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rsub__(v)
        return MFraction(super().__rsub__(v))
    
    def __mul__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__mul__(v)
        return MFraction(super().__mul__(v))
    
    def __rmul__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rmul__(v)
        return MFraction(super().__rmul__(v))
    
    def __truediv__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__truediv__(v)
        return MFraction(super().__truediv__(v))
    
    def __rtruediv__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rtruediv__(v)
        return MFraction(super().__rtruediv__(v))
        
    def __pow__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__pow__(v)
        return MFraction(super().__pow__(v))
    
    def __rpow__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rpow__(v)
        if isinstance(v,Fraction):
            return MFraction(v).__pow__(self)
        return MFraction(super().__rpow__(v))
        
    def __floordiv__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__floordiv__(v)
        return MFraction(super().__floordiv__(v))
        
    def __rfloordiv__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rfloordiv__(v)
        return MFraction(super().__rfloordiv__(v))
        
    def __mod__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__mod__(v)
        return MFraction(super().__mod__(v))
    
    def __rmod__(self,v):
        if isinstance(v,Decimal):
            return self.todecimal().__rmod__(v)
        return MFraction(super().__rmod__(v))
    
    def __abs__(self):
        return MFraction(super().__abs__())
    
    def __neg__(self):
        return MFraction(super().__neg__())
    
    def __pos__(self):
        return MFraction(super().__pos__())

