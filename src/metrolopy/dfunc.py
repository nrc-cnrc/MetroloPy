# -*- coding: utf-8 -*-

# module dfunc

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

"""
Class Dfunc is an abstract base class inherited by ummy and immy to provide 
support for numpy ufunc functions.
"""

import numpy as np
from numbers import Real,Integral
from .util import _isscalar,_replq,_replu


def _f_darctan2(x1,x2):
    return (x2/(x1**2 + x2**2),-x1/(x1**2 + x2**2))

def _f_dlogaddexp(x1,x2):
    return (np.exp(x1)/(np.exp(x1) + np.exp(x2)),
            np.exp(x2)/(np.exp(x1) + np.exp(x2)))
    
def _f_dlogaddexp2(x1,x2):
<<<<<<< HEAD
    return (2**x1/(2**x1 + 2**x2),
            2**x2/(2**x1 + 2**x2))
=======
    return (2**x1*np.log(2)/(2**x1 + 2**x2),
           2**x2*np.log(2)/(2**x1 + 2**x2))
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5

def _f_heaviside(x,h0):
    if not isinstance(h0,Real):
        raise TypeError('h0 must be a real number')
<<<<<<< HEAD
    return x.apply(np.heaviside,lambda *x:[0,0],x,h0)
=======
    return x._apply(np.heaviside,0,x,h0)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5

def _f_around(x,n=0):
    if not isinstance(n,Integral):
        raise TypeError('n must be an int')
<<<<<<< HEAD
    return x._nprnd(lambda y: np.around(y,decimals=n))
=======
    return x._nprnd(lambda y: np.around(y,n))
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5


def _f_add(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__add__(x2)
    else:
        return x2.__radd__(x1)
    
def _f_sub(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__sub__(x2)
    else:
        return x2.__rsub__(x1)

def _f_mul(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__mul__(x2)
    else:
        return x2.__rmul__(x1)
    
def _f_div(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__truediv__(x2)
    else:
        return x2.__rtruediv__(x1)
    
def _f_fdiv(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__floordiv__(x2)
    else:
        return x2.__rfloordiv__(x1)
    
def _f_pow(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__pow__(x2)
    else:
        return x2.__rpow__(x1)
    
def _f_mod(x1,x2):
    if isinstance(x1,Dfunc):
        return x1.__mod__(x2)
    else:
        return x2.__rmod__(x1)
<<<<<<< HEAD
    
def _f_fmod(x1,x2):
    if isinstance(x1,Dfunc):
        return np.sign(x1)*np.abs(x1.__mod__(x2))
    else:
        return np.sign(x1)*np.abs(x2.__rmod__(x1))
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5


ddict = {np.sin: np.cos,
         np.cos: lambda x: -np.sin(x),
         np.tan:  lambda x: 1/np.cos(x)**2,
         np.arctan2: _f_darctan2,
         np.arcsin: lambda x: 1/np.sqrt(1 - x**2),
<<<<<<< HEAD
         np.arccos: lambda x: -1/np.sqrt(1 - x**2),
=======
         np.arccos: lambda x: -1/np.sqrt(1 + x**2),
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
         np.arctan: lambda x: 1/(1 + x**2),
         np.sinh:  np.cosh,
         np.cosh:  np.sinh,
         np.tanh:  lambda x: -np.tanh(x)**2 + 1,
         np.arcsinh:  lambda x: 1/np.sqrt(1 + x**2),
<<<<<<< HEAD
         np.arccosh:  lambda x: 1/np.sqrt(x**2 - 1),
=======
         np.arccosh:  lambda x: 1/np.sqrt(1 - x**2),
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
<<<<<<< HEAD
=======
         np.heaviside: lambda x: 0, 
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
         np.sign:  lambda x: 0, 
        }

fdict = {np.angle:  lambda x: x.angle(),
         np.around:  _f_around,
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
<<<<<<< HEAD
         np.fmod: _f_fmod,
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
         np.round: _f_around,
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
            x = [a.tofloat() if hasattr(a,'tofloat') else float(a) for a in x]
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
    def apply(cls,function,derivative,*args,**kwds):
        """
        A classmethod that applies a function to one or more gummy or jummy 
        objects propagating the uncertainty.
        
        Parameters
        ----------
        function: `function`
              The the function to be applied. For `gummy.apply`, 'function'
              should take one or more float arguments and return a float value 
              or float array.  For `jummy.apply`, 'function' may also take and
              return complex values.  If the function returns an array like 
              value, it must be convertable to a numpy homogeneous array.

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
        `gummy`, `jummy` or a `numpy.ndarray` of `gummy` or `jummy`:
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
        
<<<<<<< HEAD

        b = np.broadcast(*args)
        if b.shape == ():
            return cls._iapply(function,derivative,*args,**kwds)

        ret = np.array([cls._iapply(function,derivative,*a,**kwds) for a in b])
=======
        if np.all([_isscalar(a) for a in args]):
            return cls._iapply(function,derivative,*args,**kwds)
        
        b = np.broadcast(*args)
        ret = np.empty(b.shape,dtype=object)
        ret = [cls._iapply(function,derivative,*a,**kwds) for a in b]
        ret = np.array(ret)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        
        shape = b.shape
        if isinstance(ret[0],np.ndarray):
            shape += ret[0].shape
<<<<<<< HEAD
            
        return np.reshape(ret,shape)
=======
        np.reshape(ret,shape)
        
        return ret
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    
    @classmethod
    def _iapply(cls,function,derivative,*args,**kwds):
        args = [_replq(a) for a in args]
        x = [_replu(a) for a in args]
        
        fx = function(*x,**kwds)
        d = derivative(*x,**kwds)
        if _isscalar(fx) and len(args) > 1:
           d = [i[0] if not _isscalar(i) and len(i) == 1 else i for i in d]
        
        if _isscalar(fx):
            r = cls._apply(function,derivative,fx,d,*args,**kwds)
        else:
            fx = np.asarray(fx)
            d = np.asarray(d)
            r = np.empty(fx.shape,dtype=object)
            with np.nditer(fx,flags=['multi_index']) as it:
                for i in it:
                    idx = tuple(it.multi_index)
                    r[idx] = cls._apply(function,derivative,fx[idx],d[idx],
                                        *args,**kwds)
        return r
    
    @classmethod
    def _apply(cls,function,derivative,fx,d,*a,**kwds):
        raise NotImplementedError()
    
    @classmethod
    def napply(cls,function,*args,**kwds):
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
            return complex values.  If the function returns an array like 
            value, it must be convertable to a numpy homogeneous array.

        *args:  `gummy`, `jummy`, or `float`
              One or more arguments to which `function` will be applied.  These
              arguments need not all be `Dfunc` objects; arguments  such as
              floats will be taken to be constants with no uncertainty.
              They may also be numpy ndarrays in which case the usual numpy
              broadcasting rules apply.

        Returns
        -------
        `gummy`, `jummy` or a `numpy.ndarray` of `gummy` or `jummy`:
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
        
<<<<<<< HEAD
        b = np.broadcast(*args)
        if b.shape == ():
            return cls._inapply(function,*args,**kwds)
        
        ret = np.array([cls._inapply(function,*a,**kwds) for a in b])
=======
        if np.all([_isscalar(a) for a in args]):
            return cls._inapply(function,*args,**kwds)
        
        b = np.broadcast(*args)
        ret = np.empty(b.shape,dtype=object)
        ret = [cls._inapply(function,*a,**kwds) for a in b]
        ret = np.array(ret)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        
        shape = b.shape
        if isinstance(ret[0],np.ndarray):
            shape += ret[0].shape
<<<<<<< HEAD
            
        return np.reshape(ret,shape)
=======
        np.reshape(ret,shape)
        
        return ret
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    
    @classmethod
    def _inapply(cls,function,*args,**kwds):
        args = [_replq(a) for a in args]
        x = [_replu(a) for a in args]
        
        fxx = function(*x,**kwds)
        
        if _isscalar(fxx):
            r = cls._napply(function,fxx,*args,**kwds)
        else:
            fxx = np.asarray(fxx)
            r = np.empty(fxx.shape,dtype=object)
            with np.nditer(fxx,flags=['multi_index']) as it:
                for i in it:
                    idx = tuple(it.multi_index)
                    r[idx] = cls._napply(lambda *x:np.asarray(function(*x))[idx],
                                         fxx[idx],*args)
        return r
        
    @classmethod
    def _napply(cls,function,fxx,*args,**kwds):
        raise NotImplementedError()
    def _ufunc(self,ufunc,*args,**kwds):
        if any(isinstance(a,np.ndarray) for a in args):
            args = [np.array(a) if isinstance(a,Dfunc) else a for a in args]
        try:
            return _call(lambda *x: self.apply(ufunc,ddict[ufunc],*x), *args)
        except KeyError:
            try:
                return _call(fdict[ufunc],*args)
            except KeyError:
                return self.napply(ufunc,*args)
       
    def __array_ufunc__(self,ufunc,method,*args,**kwds):
        if method != '__call__':
            return None
        
        return self._ufunc(ufunc,*args,**kwds)
            
<<<<<<< HEAD
    def __array_function__(self,func,method,*args,**kwds): 
        #return self._ufunc(func,*args,**kwds)
        return self._ufunc(func,*args[0],**args[1])
=======
    def __array_function__(self,func,method,*args,**kwds):        
        return self._ufunc(func,*args,**kwds)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
