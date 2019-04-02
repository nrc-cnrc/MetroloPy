# -*- coding: utf-8 -*-

# module functions

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
A number of mathematical functions are defined here that can be used with
gummys.
"""

import numpy as np
from .gummy import gummy,jummy
from .ummy import ummy
from numbers import Complex
from .dfunc import _call

def _callg(f,*args):
    # decide whether to call jummy.apply, gummy.apply or ummy.apply
    g = None
    u = None
    c = False
    for a in args:
        if isinstance(a,jummy):
            return a.__array_ufunc__(f,'__call__',*args)
        if isinstance(a,gummy):
            g = a
        if isinstance(a,ummy):
            u = a
        if isinstance(a,Complex):
            c = True
    if g is not None:
        if c:
            return jummy(g).__array_ufunc__(f,'__call__',*args)
        return g.__array_ufunc__(f,'__call__',*args)
    if u is not None:
        return u.__array_ufunc__(f,'__call__',*args)
    return _call(f,*args)

def _bcallg(f,*args):
    # broadcast f across any arrays in args.
    bargs = np.broadcast(*args)
    if bargs.shape == ():
        return _callg(f,*args)
    ret = np.array([_callg(f,*a) for a in bargs])
    ret = ret.reshape(bargs.shape)
    return ret


def sin(x):
    """
    Returns the sine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.sin,x)
    
def cos(x):
    """
    Returns the cosine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.cos,x)

def tan(x):
    """
    Returns the tangent of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.tan,x)
 
def arcsin(x):
    """
    Returns the inverse sine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.arcsin,x)

def arccos(x):
    """
    Returns the inverse cosine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.arccos,x)

def arctan(x):
    """
    Returns the inverse tangent of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.arctan,x)

def arctan2(x1,x2):
    """
    Returns the inverse tangent of x1/x2 choosing the quadrant correctly where 
    x1 and x2 can be float, complex, gummy or jummy.  
    """
    return _bcallg(np.arctan2,x1,x2)
    
def sinh(x):
    """
    Returns the hyperbolic sine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.sinh,x)

def cosh(x):
    """
    Returns the hyperbolic cosine of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.cosh,x)

def tanh(x):
    """
    Returns the hyperbolic tangent of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.tanh,x)

def arcsinh(x):
    """
    Returns the inverse hyperbolic sine of x where x can be float, complex, gummy 
    or jummy.
    """
    return _bcallg(np.arcsinh,x)

def arccosh(x):
    """
    Returns the inverse hyperbolic cosine of x where x can be float, complex, gummy 
    or jummy.
    """
    return _bcallg(np.arccosh,x)

def arctanh(x):
    """
    Returns the inverse hyperbolic tangent of x where x can be float, complex, 
    gummy or jummy.
    """
    return _bcallg(np.arctanh,x)

def exp(x):
    """
    Returns e to the power of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.exp,x)

def exp2(x):
    """
    Returns 2 to the power of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.exp2,x)

def expm1(x):
    """
    Returns exp(x) - 1 where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.expm1,x)

def log(x):
    """
    Returns the natural logrithm of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.log,x)
    
def log2(x):
    """
    Returns the log base 2 of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.log2,x)

def log10(x):
    """
    Returns the log base 10 of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.log10,x)

def log1p(x):
    """
    Returns the natural logrithm of x plus 1 where x can be float, complex, gummy 
    or jummy.
    """
    return _bcallg(np.log1p,x)

def logaddexp(x1,x2):
    """
    Returns the log(exp(x1) + exp(x2)) where x1 and x2 can be float, complex, gummy 
    or jummy and log is the natural logrithm.
    """
    return _bcallg(np.logaddexp,x1,x2)

def logaddexp2(x1,x2):
    """
    Returns the log2(2**x1 + 2**x2) where x1 and x2 can be float, complex, gummy 
    or jummy and log2 is the logrithm to base 2.
    """
    return _bcallg(np.logaddexp2,x1,x2)

def sqrt(x):
    """
    Returns the square root of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.sqrt,x)

def cbrt(x):
    """
    Returns the cube root of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.cbrt,x)

def square(x):
    """
    Returns the square of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.square,x) 

def add(x1,x2):
    """
    returns x1 + x2
    """
    return _bcallg(np.add,x1,x2)

def subtract(x1,x2):
    """
    returns x1 - x2
    """
    return _bcallg(np.subtract,x1,x2)

def negative(x):
    """
    returns -x
    """
    return _bcallg(np.negative,x)

def multiply(x1,x2):
    """
    returns x1 * x2
    """
    return _bcallg(np.multiply,x1,x2)

def divide(x1,x2):
    """
    returns x1 / x2
    """
    return _bcallg(np.divide,x1,x2)

def true_divide(x1,x2):
    """returns x1 / x2
    """
    return _bcallg(np.true_divide,x1,x2)

def floor_divide(x1,x2):
    """
    returns x1 // x2
    """
    return _bcallg(np.floor_divide,x1,x2)

def reciprocal(x):
    """
    returns 1/x
    """
    return _bcallg(np.reciprocal,x)

def power(x1,x2):
    """
    returns x1**x2
    """
    return _bcallg(np.power,x1,x2)

def absolute(x):
    """
    Returns the absoulte value of x where x can be float, complex, gummy or jummy.
    This is equivalent to abs(x).
    """
    return _bcallg(np.absolute,x)

def mod(x1,x2):
    """
    returns x1 % x2
    """
    return _bcallg(np.mod,x1,x2)

def remainder(x1,x2):
    """
    returns x1 % x2
    """
    return _bcallg(np.remainder,x1,x2)

def divmod(x1,x2):
    """
    returns (x1 // x2, x1 % x2)
    """
    return _bcallg(np.divmod,x1,x2)

def modf(x1,x2):
    """
    returns (x1 % 1, x1 // 1), a tuple of integer and fractional parts
    """
    return _bcallg(np.modf,x1,x2)

def angle(x):
    """
    Returns the complex argument of x, where x can be float, complex, gummy or 
    jummy and the return value is a float if x is float or complex and gummy
    if x is gummy or jummy.
    """
    return _bcallg(np.angle,x)

def real(x):
    """
    returns x.real
    """
    return _bcallg(np.real,x)

def imag(x):
    """
    returns x.imag
    """
    return _bcallg(np.imag,x)

def conj(x):
    """
    returns x.conjugate(), the complex conjugate of x
    """
    return _bcallg(np.conj,x)

def around(x,n=0):
    """
    Returns x rounded to n digits where x can be float, complex, gummy or jummy
    but n must be int.
    """
    return _bcallg(np.around,x,n)

def rint(x):
    """
    Returns x rounded to the nearest integer value where x can be float, complex, 
    gummy or jummy.
    """
    return _bcallg(np.rint,x)

def fix(x):
    """
    Returns x rounded towards zero where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.fix,x)

def floor(x):
    """
    Returns the floor of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.floor,x)

def ceil(x):
    """
    Returns the ceiling of x where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.ceil,x)

def trunc(x):
    """
    Returns x rounded towards zero where x can be float, complex, gummy or jummy.
    """
    return _bcallg(np.trunc,x)

def heaviside(x,h0):
    """
    Heavyside function of x, h0 is the value at x = 0
    """
    return _bcallg(np.heavyside,x,h0)

def sign(x):
    """
    sign of x
    """
    return _bcallg(np.sign,x)

def sum(*args,**kwds):
    """
    Alias for numpy.sum
    """
    return np.sum(*args,**kwds)

def prod(*args,**kwds):
    """
    Alias for numpy.prod
    """
    return np.prod(*args,**kwds)

def cumsum(*args,**kwds):
    """
    Alias for numpy.cumsum
    """
    return np.cumsum(*args,**kwds)

def cumprod(*args,**kwds):
    """
    Alias for numpy.cumprod
    """
    return np.cumprod(*args,**kwds)

def diff(*args,**kwds):
    """
    Alias for numpy.diff
    """
    return np.diff(*args,**kwds)

def ediff1d(*args,**kwds):
    """
    Alias for numpy.ediff1d
    """
    return np.ediff1d(*args,**kwds)

def gradient(*args,**kwds):
    """
    Alias for numpy.gradient
    """
    return np.gradient(*args,**kwds)

def cross(*args,**kwds):
    """
    Alias for numpy.cross
    """
    return np.cross(*args,**kwds)
    

