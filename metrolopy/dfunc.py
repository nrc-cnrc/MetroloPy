# -*- coding: utf-8 -*-

# module dfunc

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
Class Dfunc is an abstract base class inherited by gummy and jummy to provides 
some support for numpy broadcasting for functions and operators.
"""

import numpy as np
from numbers import Real,Integral

def _f_darctan2(x1,x2):
    return (x2/(x1**2 + x2**2),-x1/(x1**2 + x2**2))

def _f_dlogaddexp(x1,x2):
    return (np.exp(x1)/(np.exp(x1) + np.exp(x2)),
            np.exp(x2)/(np.exp(x1) + np.exp(x2)))
    
def _f_dlogaddexp2(x1,x2):
    return (2**x1*np.log(2)/(2**x1 + 2**x2),
           2**x2*np.log(2)/(2**x1 + 2**x2))

def _f_heaviside(x,h0):
    if not isinstance(h0,Real):
        raise TypeError('h0 must be a real number')
    return x._apply(np.heaviside,0,x,h0)

def _f_around(x,n=0):
    if not isinstance(n,Integral):
        raise TypeError('n must be an int')
    return x._nprnd(lambda y: np.around(y,n))


def _f_add(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._add(x2)
    else:
        return x2._radd(x1)
    
def _f_sub(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._sub(x2)
    else:
        return x2._rsub(x1)

def _f_mul(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._mul(x2)
    else:
        return x2._rmul(x1)
    
def _f_div(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._truediv(x2)
    else:
        return x2._rtruediv(x1)
    
def _f_fdiv(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._floordiv(x2)
    else:
        return x2._rfloordiv(x1)
    
def _f_pow(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._pow(x2)
    else:
        return x2._rpow(x1)
    
def _f_mod(x1,x2):
    if isinstance(x1,Dfunc):
        return x1._mod(x2)
    else:
        return x2._rmod(x1)


ddict = {np.sin: np.cos,
         np.cos: lambda x: -np.sin(x),
         np.tan:  lambda x: 1/np.cos(x)**2,
         np.arctan2: _f_darctan2,
         np.arcsin: lambda x: 1/np.sqrt(1 - x**2),
         np.arccos: lambda x: -1/np.sqrt(1 + x**2),
         np.arctan: lambda x: 1/(1 + x**2),
         np.sinh:  np.cosh,
         np.cosh:  np.sinh,
         np.tanh:  lambda x: -np.tanh(x)**2 + 1,
         np.arcsinh:  lambda x: 1/np.sqrt(1 + x**2),
         np.arccosh:  lambda x: 1/np.sqrt(1 - x**2),
         np.arctanh:  lambda x: 1/(1 - x**2),
         np.exp:  np.exp,
         np.exp2:  lambda x: np.log(2)*2**x,
         np.expm1:  np.exp,
         np.log:  lambda x: 1/x,
         np.log2:  lambda x: 1/(x*np.log(2)),
         np.log10:  lambda x: 1/(x*np.log(10)),
         np.log1p:  lambda x: 1/(1+x),
         np.logaddexp:  _f_dlogaddexp,
         np.logaddexp2:  _f_dlogaddexp2,
         np.heaviside: lambda x: 0, 
         np.sign:  lambda x: 0, 
        }

fdict = {np.angle:  lambda x: x.angle(),
         np.around:  _f_around,
         np.round_:  _f_around,
         np.heaviside:  _f_heaviside,
         np.absolute:  lambda x: abs(x),
         np.add:  _f_add,
         np.subtract:  _f_sub,
         np.negative:  lambda x: -x,
         np.multiply:  _f_mul,
         np.divide:  _f_div,
         np.true_divide:  _f_div,
         np.floor_divide:  _f_fdiv,
         np.reciprocal:  lambda x: 1/x,
         np.power:  _f_pow,
         np.mod:  _f_mod,
         np.remainder: _f_mod,
         np.divmod:  lambda x1,x2: (_f_fdiv(x1,x2),_f_mod(x1,x2)),
         np.modf:  lambda x: (x%1,x//1),
         np.sqrt:  lambda x: x**0.5,
         np.square:  lambda x: x**2,
         np.cbrt:  lambda x: x**(1/3),
         np.real:  lambda x: x.real,
         np.imag:  lambda x: x.imag,
         np.conj:  lambda x: x.conjugate(),
         np.rint:  lambda x: x._nprnd(np.rint),
         np.fix:  lambda x: x._nprnd(np.fix),
         np.floor:  lambda x: x._nprnd(np.floor),
         np.ceil:  lambda x: x._nprnd(np.ceil),
         np.trunc:  lambda x: x._nprnd(np.trunc),
         np.isnan:  lambda x: np.isnan(x.x),
         np.isinf:  lambda x: np.isinf(x.x),
         np.isfinite:  lambda x: np.isfinite(x.x),
         np.isneginf:  lambda x: np.isneginf(x.x),
         np.isposinf:  lambda x: np.isposinf(x.x)
        }


try_fconvert = True
def _call(f,*x):
    if try_fconvert:
        try:
            return f(*x)
        except:
            x = [a.tofloat() if isinstance(a,Dfunc) else float(a) for a in x]
            return f(*x)

    return f(*x)
   
        
def _broadcast(f,x):
    # used for the binary operations __add__, ...
    if isinstance(x,np.ndarray):
        bx = np.broadcast(x)
        ret = np.array([f(*b) for b in bx])
        if bx.shape == ():
            ret = ret.item()
        else:
            ret = ret.reshape(bx.shape)
        return ret
    return _call(f,x)


class Dfunc:
    """
    Class `Dfunc` is an abstract base class that provides some support for numpy
    broadcasting for functions and operators.  An inheriting class must implement
    the `_apply(self,function,derivative,*args)`, `_napply(self,function,*args)`,
    and `tofloat(self)` methods, as well as `_add(x)`, `_radd(x)`, `_sub(x)`, ...
    """
    
    def tofloat(self):
        # this should return a copy of self with the x and u properties 
        # converted to float values
        raise NotImplementedError()
        
    @classmethod
    def apply(cls,function,derivative,*args):
        """
        A classmethod that applies a function to one or more gummy or jummy 
        objects propagating the uncertainty.
        
        Parameters
        ----------
        function: `function`
              The the function to be applied. For `gummy.apply`, 'function'
              should take one or more float arguments and return a float value 
              or float array.  For `jummy.apply`, 'function' may also take and
              return complex values.

        derivative:  `function`
              The name of a second function that gives the derivatives
              with respect to the arguments of `function`.  `derivative` should
              take an equal number of arguments as `function`.  If `function`
              takes one argument `derivative` should return a float and if
              `function` takes more than one argument then `derivative` should
              return a tuple, list or array of floats that contains the derivatives
              with respect to each argument.  In the case of `jummy.apply`, the
              derivatives with respect to each argument may be real or complex
              values, in which case `function` is assumed to be holomorphic.  Or
              the derivative may be a 2 x 2 matrix of the form:

                              [[ du/dx, du/dy ],
                               [ dv/dx, dv/dy ]]

             where function(x + j*y) = u + j*v.

        *args:  `gummy`, `jummy`, or `float`
              One or more arguments to which `function` will be applied.  These
              arguments need not all be `Dfunc` objects; arguments  such as
              floats will be taken to be constants with no uncertainty.
              They may also be numpy ndarrays in which case the usual numpy
              broadcasting rules apply.
              
        Returns
        -------
        `gummy`, `jummy`:
            If none of the arguments are `gummy` or `jummy`
            then the return value is the same type as the return value of `function`.
            Otherwise `gummy.apply` returns a `gummy` and `jummy.apply` returns either a
            `gummy` or a `jummy` depending on whether `function` has a float or
            a complex return value.
            
        
        Examples
        --------
            
        >>> import numpy as np
        >>> x = gummy(0.678,u=0.077)
        >>> gummy.apply(np.sin,np.cos,x)
        0.627 +/- 0.060
        
        >>> x = gummy(1.22,u=0.44)
        >>> y = gummy(3.44,u=0.67)
        >>> def dhypot(x,y):
        ...     return (x1/sqrt(x1**2 + x2**2),x2/np.sqrt(x1**2 + x2**2))
        >>> gummy.apply(np.hypot,dhypot,x,y)
        3.65 +/- 0.65
        """
        
        bargs = np.broadcast(*args)
        if bargs.shape == ():
            return _call(lambda *x: cls._apply(function,derivative,*x), *args) 
        
        ret = np.array([_call(lambda *x: cls._apply(function,derivative,*x), *a) for a in bargs])
        ret = ret.reshape(bargs.shape)
        return ret
    
    @classmethod
    def napply(cls,function,*args):
        """
        gummy.napply(function, arg1, arg2, ...) and
        jummy.napply(function, arg1, arg2, ...)
        
        A classmethod that applies a function to one or more gummy or jummy 
        objects propagating the uncertainty.  This method is similar to apply 
        except that the derivatives are computed numerically so a derivative 
        function does not need to be supplied.
        
        Parameters
        ----------
        function: `function`
            The the function to be applied. For `gummy.apply`, 'function'
            should take one or more float arguments and return a float value
            r float array.  For `jummy.apply`, 'function' may also take and
            return complex values.

        *args:  `gummy`, `jummy`, or `float`
              One or more arguments to which `function` will be applied.  These
              arguments need not all be `Dfunc` objects; arguments  such as
              floats will be taken to be constants with no uncertainty.
              They may also be numpy ndarrays in which case the usual numpy
              broadcasting rules apply.

        Returns
        -------
        `gummy`, `jummy`:
            If none of the arguments are `gummy` or `jummy`
            then the return value is the same type as the return value of `function`.
            Otherwise `gummy.apply` returns a `gummy` and `jummy.apply` returns either a
            `gummy` or a `jummy` depending on whether `function` has a float or
            a complex return value.
            
        
        Examples
        --------
            
        >>> import numpy as np
        >>> x = gummy(0.678,u=0.077)
        >>> gummy.napply(np.sin,x)
        0.627 +/- 0.060
        
        >>> x = gummy(1.22,u=0.44)
        >>> y = gummy(3.44,u=0.67)
        >>> gummy.napply(np.hypot,x,y)
        3.65 +/- 0.65
        """
        
        bargs = np.broadcast(*args)
        if bargs.shape == ():
            return _call(lambda *x: cls._napply(function,*x), *args)
        
        ret = np.array([_call(lambda *x: cls._napply(function,*x), *a) for a in bargs])
        ret = ret.reshape(bargs.shape)
        return ret
       
    def __array_ufunc__(self,ufunc,method,*args,**kwds):
        if method != '__call__':
            return None
        
        try:
            return _call(lambda *x: self._apply(ufunc,ddict[ufunc],*x), *args)
        except KeyError:
            try:
                return _call(fdict[ufunc],*args)
            except KeyError:
                return None
            
    def __add__(self,b):
        return _broadcast(self._add,b)
    
    def __radd__(self,b):
        return _broadcast(self._radd,b)
    
    def __sub__(self,b):
        return _broadcast(self._sub,b)
    
    def __rsub__(self,b):
        return _broadcast(self._rsub,b)
    
    def __mul__(self,b):
        return _broadcast(self._mul,b)
        
    def __rmul__(self,b):
        return _broadcast(self._rmul,b)
        
    def __truediv__(self,b):
        return _broadcast(self._truediv,b)
    
    def __rtruediv__(self,b):
        return _broadcast(self._rtruediv,b)
    
    def __floordiv__(self,b):
        return _broadcast(self._floordiv,b)
    
    def __rfloordiv__(self,b):
        return _broadcast(self._rfloordiv,b)
    
    def __mod__(self,b):
        return _broadcast(self._mod,b)
    
    def __rmod__(self,b):
        return _broadcast(self._rmod,b)
    
    def __pow__(self,b):
        return _broadcast(self._pow,b)
            
    def __rpow__(self,b):
        return _broadcast(self._rpow,b)
