# -*- coding: utf-8 -*-

# module ummy

## Copyright (C) 2019 National Research Council Canada
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
The ummy object defined here was created primarily to separate the code for
first-order uncertainty propagation from the other methods in the gummy object.  
But the ummy object can be used as a stripped down verison of gummy without the
Monte-Carlo uncertainty propagation, unit, and pretty printing functionality
of gummy.
"""

import numpy as np
import weakref
from .dfunc import Dfunc
from .exceptions import UncertiantyPrecisionWarning
from collections import namedtuple
from warnings import warn
from math import isnan,isinf,isfinite,sqrt,log
from fractions import Fraction
from numbers import Rational,Integral

try:
    from mpmath import mp,mpf,rational
except:
    mp = mpf = None

_FInfo = namedtuple('_FIinfo','rel_u precision warn_u ')
_iinfo = _FInfo(rel_u=0,precision=0,warn_u=0)
_finfodict = {}
def _getfinfo(x):
    if isinstance(x,Rational):
        return (_iinfo,x)
    
    if mpf is not None and isinstance(x,mpf):
        n = mp.prec
        x = mpf(x)
        try:
            return (_finfodict[n],x)
        except:
            f = _FInfo(rel_u=0.4222/2**n,
                       precision=mp.dps,
                       warn_u = 0)
            _finfodict[n] = f
            return (f,x)
        
    t = type(x)
    try:
        return (_finfodict[t],x)
    except:
        try:
            f = np.finfo(t)
            f = _FInfo(rel_u=f.eps*0.2111,
                       precision=f.precision,
                       warn_u=100*f.tiny*f.eps)
            _finfodict[t] = f
        except:
            _finfodict[t] = _iinfo
            return _iinfo
        return (f,x)

def _combu(a,b,c):
    # function for combining uncertainties,
    # computes sqrt(a**2 + 2*c*a*b + b**2)
    # with 1 >= c >= -1,
    # this method avoids overflows or underflows when taking the squares
    
    if isnan(a) or isnan(b):
        try:
            return type(a)('nan')
        except:
            return float('nan')
    
    if abs(b) > abs(a):
        a,b = b,a
    if b == 0:
        return abs(a)
        
    r = b/a
    x = 1 + r**2 + 2*c*r
    
    if x < 0:
        if x < -1e6:
            # something horrible happened, should only get here due to a bug
            raise FloatingPointError('u is sqrt of ' + str(a**2*x))
        # maybe possibly get here due to a floating point rounding error somewhere
        return 0
    
    return abs(a)*x**0.5

def _isscalar(x):
    try:
        len(x)
        return False
    except:
        return True
    
def _floor(x):
    if x < 0:
        return int(x) - 1
    return int(x)

class ummy(Dfunc):
    max_dof = 10000 # any larger dof will be rounded to float('inf')
    nsig = 2 # The number of digits to quote the uncertainty to
    thousand_spaces = True # If true spaces will be placed between groups of three digits.
    
    # Set sci_notation to None to automatically switch to scientific notation
    # when the exponent is greater than sci_notation_high or less than 
    # sci_notation_low.  Or set sci_notation to True (or False) to force 
    # scientific notation to (not) be displayed if possible.
    sci_notation = None
    sci_notation_high = 7
    sci_notation_low = -3
    
    # if rounding_u is true then uncertianty will be added for floating point 
    #   rounding errors
    rounding_u = False
    
    max_digits = 15
    
    def __init__(self,x,u=0,dof=float('inf'),utype=None):
        if isinstance(x,ummy):
             self._copy(x,self,formatting=False)
             return
         
        try:
            float(x)
        except TypeError:
            raise TypeError('x must be a real number')
            
        try:
            float(u)
            if not isnan(u) and u < 0:
                raise TypeError()
        except TypeError:
            raise TypeError('u must be a real number >= 0')
                
        if self.rounding_u and not isinstance(x,Rational):
            finfo,x = _getfinfo(x)
            if finfo is _iinfo:
                warn('numpy.finfo cannot get the floating point accuracy for x\nno uncertainty will be included to account for floating point rounding errors',UncertiantyPrecisionWarning)
            else:
                uc = _combu(u,float(abs(x)*finfo.rel_u),0)
            
            uinfo = _getfinfo(u)[0]
            if uc < uinfo.warn_u and u is not 0:
                warn('a gummy/ummy was defined with an uncertainty too small to be accurately represented by a float\nuncertainty values may contain significant numerical errors',UncertiantyPrecisionWarning)
            u = uc
            self._finfo = finfo
        else:
            self._finfo = None
            
        if isinstance(x,Fraction):
            x = MFraction(x)
        if isinstance(u,Fraction):
            u = MFraction(u)
            
        self._x = x
        self._u = u
        
        try:
            if isinstance(dof,Integral):
                dof = int(dof)
            else:
                dof = float(dof)
            if dof <= 0 or isnan(dof):
                raise TypeError()
        except TypeError:
            raise TypeError('dof must be a real number > 0')
        if dof > self.max_dof:
            self._dof = float('inf')
        else:
            self._dof = dof
                    
        self._refs = 1
        if not isinf(u) and not isnan(u) and u > 0:
            self._ref = _GummyRef()
            if utype is not None:
                GummyTag.set_tag(self,utype)
        else:
            self._ref = None
        
    @property
    def x(self):
        return self._x
        
    @property
    def u(self):
        return self._u
        
    @property
    def dof(self):
        """
        `float`, read-only

        Returns the number or degrees of freedom that the uncertainty of the 
        gummy is based on.  If the gummy was created as the result of an 
        operation between two or more other gummys, then the dof is the effective
        number of degrees of freedom calculated using a version of the Welch-
        Satterthwaite approximation that takes into account correlations, see
        [R. Willink, Metrologia, 44, 340 (2007)].
        """
        if self._dof < 1:
            # Occasionally correlations can result in a _dof less than 1;
            # see the _get_dof method.
            return 1
        return self._dof
        
    @property
    def const(self):
        """`bool`, read-only
        
        Returns `True` if the gummy represents an exact value with no uncertainty.
        """
        return (self._u == 0)
    
    @property
    def finfo(self):
        if self._finfo is None:
            self._finfo = _getfinfo(self._x)
        return self._finfo
        
    @staticmethod
    def _copy(s,r,formatting=True,tofloat=False):
        # copies attributes of s to r
        if tofloat:
            r._x = float(s._x)
            r._u = float(s._u)
        else:
            r._x = s._x
            r._u = s._u
        r._dof = s._dof
        r._ref = s._ref
        r._refs = s._refs
        r._finfo = s._finfo
        if formatting:
            if s.nsig != ummy.nsig:
                r.nsig = s.nsig
    
    def copy(self,formatting=True,tofloat=False):
        """
        Returns a copy of the gummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new gummy.  The default for `formatting` is `True`.  If tofloat
        is true the x and u properties will be converted to float values
        before copying.
        """
        r = type(self)(self._x,u=self._u,dof=self._dof)
        self._copy(self,r,formatting=formatting,tofloat=tofloat)
        return r
        
    def tofloat(self):
        """
        Returns a copy of the gummy with x an u converted to floats.
        """
        return self.copy(formatting=False,tofloat=True)
    
    def correlation(self, g):
        """
        Returns the correlation coefficient between `self` and `g`."""
        if not isinstance(g,ummy):
            return 0
        if self._ref is None or g._ref is None:
            return 0
        return self._refs*g._refs*self._ref.cor(g._ref)
        
    def covariance(self,g):
        """
        Returns the covariance between `self` and `g`.
        """
        if not isinstance(g,ummy):
            return 0
        if self._ref is None or g._ref is None:
            return 0
        return self._u*g._u*self._refs*g._refs*self._ref.cor(g._ref)
        
    @staticmethod
    def correlation_matrix(gummys):
        """
        Returns the correlation matrix of a list or array of gummys.  The return
        value is `numpy.ndarray`.
        """
        return np.array([[b.correlation(a) if isinstance(b,ummy) else 0 
                            for b in gummys] for a in gummys])
        
    @staticmethod
    def covariance_matrix(gummys):
        """
        Returns the variance-covariance matrix of a list or array of gummys.  The return
        value is `numpy.ndarray`.
        """
        return np.array([[b.covariance(a) if isinstance(b,ummy) else 0 for b in gummys] 
                            for a in gummys])
                                
    @staticmethod
    def _set_correlation_matrix(gummys, matrix):
        n = len(gummys)
        m = np.asarray(matrix)
        
        if not np.allclose(m.T,m,atol = 10*_GummyRef._cortol):
            raise ValueError('the correlation matrix is not symmetric')
        e,v = np.linalg.eigh(m)
        if any(e <= -_GummyRef._cortol):
            raise ValueError('the correlation matrix is not positive semi-definate')
            
        if m.shape != (n,n):
            raise ValueError('matrix must have shape len(gummys) x len(gummys)')
        
        # Check that the diagonal elements are 1 (leaving room for rounding errors).
        for i in range(n):
            if abs(m[i][i] - 1.0) > 0.01:
                raise ValueError('correlation matrix diagonal elements must be 1')
            
        # Set the correlations using the off diagonal elements.
        for i in range(0,n):
            for j in range ((i+1),n):
                if m[i][j] != 0:
                    if gummys[i]._ref is not None and gummys[j]._ref is not None:
                        gummys[i]._ref.set_cor(gummys[j],m[i][j])
                    else:
                        raise ValueError('constants cannot have a non-zero correlation')
                    
    @staticmethod
    def _set_covariance_matrix(gummys, matrix):
        n = len(gummys)
        m = np.array(matrix)
        if m.shape != (n,n):
            raise ValueError('matrix must have shape len(gummys) x len(gummys)')
        
        # Set the u of each gummy using the diagonal of the matrix.
        for i in range(n):
            if m[i][i] == 0:
                gummys[i]._u = 0
                gummys[i]._ref = None
            else:
                if gummys[i]._ref is None:
                    gummys[i]._ref = _GummyRef()  
                e = m[i][i]**0.5
                gummys[i]._u = e
        
        # Set the correlations
        for i in range(0,n):
            for j in range(0,n):
                m[i][j] /= gummys[i]._u*gummys[j]._u
        ummy._set_correlation_matrix(gummys,m)
        
    @classmethod
    def create(cls,x,u=None,dof=None,correlation_matrix=None,
               covariance_matrix=None):
        """
        A class method that creates a list of (possibly) correlated ummys.
        
        Parameters
        ----------
        x:
            A list of floats corresponding to the x-value of each ummy.
        u, dof, k, loc, utype:  optional
            Lists that correspond to the
            parameters in the ummy initializer (with the i-th value in each
            list passed to the initializer for the i-th ummy).  These may also be a single
            value with this same value is to passed to each initializer.
        correlation_matrix:  optional
            A list or array to be used as the correlation
            matrix of the ummys. This is optional and must be set to the
            default value of None if the covariance_matrix is specified.
            If both the correlation_matrix and the covariance_matrix are
            None (or omitted) then the ummys will be uncorrelated.
        covariance_matrix:   optional
            A list or array to be used as the variance-
            covariance matrix of the ummys.  If the covariance matrix is
            specified the u parameter will be ignored  This parameter is
            optional and must be set to the default value of None if the
            correlation_matrix is specified.  If both the correlation_matrix
            and the covariance_matrix are None (or omitted) then the ummys
            will be uncorrelated.
                
        Returns
        -------
        a list of ummys
        """
       
            
        if covariance_matrix is not None:
            if correlation_matrix is not None:
                raise TypeError('correlation_matrix and covariance_matrix cannot both be specified')
                
        n = len(x)
        if u is None:
            u = [0]*n
        elif _isscalar(u):
            u = [u]*n
        if dof is None:
            dof = [float('inf')]*n
        elif _isscalar(dof):
            dof = [dof]*n
            
        ret = [cls(x[i],u=u[i],dof=dof[i]) for i in range(n)]
            
        if correlation_matrix is not None:
            cls._set_correlation_matrix(ret, correlation_matrix)
        
        if covariance_matrix is not None:
            cls._set_covariance_matrix(ret, covariance_matrix)
            
        return ret
        
    def __str__(self):
        x = float(self._x)
        try:
            xabs = abs(x)
           
            if xabs != 0:
                xexp = _floor(np.log10(xabs))
            else:
                xexp = None
                
            u = float(self.u)
            
            if self._u == 0 or isnan(self._u) or isinf(self._u):
                if xexp is None:
                    return '0'
                else:
                    if (((self.sci_notation is None and (xexp > self.sci_notation_high or xexp < self.sci_notation_low)) 
                        or self.sci_notation) ):
                        xtxt = self._format_mantissa('unicode',x*10**(-xexp),
                                                     self.finfo.precision)
                        etxt = 'e'
                        if xexp >= 0:
                            etxt += '+'
                        etxt += str(xexp)
                        return xtxt + etxt
                    else:
                        return self._format_mantissa('unicode',x,
                                                     self.finfo.precision)
                    
            uabs = abs(u)
            if uabs != 0:
                uexp = _floor(np.log10(uabs))
                
            if xexp is None:
                xp = uexp
            else:
                xp = max(xexp,uexp)
                
            psn = (uexp - self.nsig + 1 > 0)
            
            if (((self.sci_notation is None and (xp > self.sci_notation_high or xp < self.sci_notation_low)) 
                        or self.sci_notation) or psn or (xexp is not None and (uexp > xexp and xexp == 0))):
                xtxt = self._format_mantissa('unicode',x*10**(-xp),-uexp+xp+self.nsig-1)
                etxt = 'e'
                if xp >= 0:
                    etxt += '+'
                etxt += str(xp)
            else:
                xtxt = self._format_mantissa('unicode',x,-uexp+self.nsig-1)
                etxt = ''
                
            if self.nsig <= 0:
                return xtxt + etxt
    
            utxt = self._format_mantissa('unicode',u*10**(-uexp),self.nsig-1,parenth=True)
            return xtxt + '(' + utxt + ')' + etxt
        except:
            try:
                return(str(self.x) + '{' + str(self.u) + '}' + '??')
            except:
                try:
                    return(str(self.x) + '{??}')
                except:
                    return('??')
                
    def __repr__(self):
        return self.__str__()
          
    def _format_mantissa(self,fmt,x,sig,parenth=False):            
        try:
            if isnan(x):
                return 'nan'
            if np.isposinf(x):
                if fmt == 'html':
                    return '&infin;'
                if fmt == 'latex':
                    return r'\infty'
                if fmt == 'ascii':
                    return 'inf'
                return '\u221E'
            if np.isneginf(x):
                if fmt == 'html':
                    return '-&infin;'
                if fmt == 'latex':
                    return r'-\infty'
                if fmt == 'ascii':
                    return '-inf'
                return '-\u221E'
        except:
            # for large int values the above np functions can generate an error
            pass
        
        ellipses = ''
        if isinstance(x,Rational) and sig is None:
            if x == int(x):
                s = str(int(x))
            else:
                xf = float(x)
                e = 15 - _floor(np.log10(abs(xf)))
                if e < 0:
                    e = 0
                n = 10**e*x
                if int(n) != n:
                    ellipses = '...'
                s = str(round(xf,e))
        elif sig is None:
            s = str(x)
        else:
            if x != 0 and self.max_digits is not None:
                if mpf is not None and isinstance(x,mpf):
                    em = _floor(mp.log10(abs(x)))
                else:    
                    em = _floor(np.log10(abs(float(x))))
                e = self.max_digits - em
                if e < 0:
                    e = 0
                if sig > e and sig > 0:
                    ellipses = '...'
                    sig = e
            if sig < 0:
                sig = 0
            if mpf is not None and isinstance(x,mpf):
                s = mp.nstr(x,n=sig+1,strip_zeros=False)
            else:
                s = round(float(x),sig)
                s = ('{:.' + str(sig) + 'f}').format(s)
        
        if parenth:
            s = s.replace('.','')
            s = s.lstrip('0')
            if len(s) == 0:
                s = '0'
            return s
        elif self.thousand_spaces:
            # Now insert spaces every three digits
            if fmt == 'html':
                sp = '&thinsp;'
            elif fmt == 'latex':
                sp = r'\,'
            else:
                sp = ' '
                
            d = s.find('.')
            if d != -1:
                i = d + 1
                while i < len(s) and s[i].isdigit():
                    i += 1
                nr = i - d - 1
            else:
                nr = 0
                d = len(s)
                
            i = d - 1
            while i > 0 and s[i].isdigit():
                i -= 1
            nl = d - i - 1
            
            if nr <= 4 and nl <= 4:
                return s + ellipses
                
            s = list(s)    
            if nr > 4:
                for i in range(int(nr/3)):
                    s.insert(d + i + (i+1)*3 + 1, sp)
            if nl > 4:
                for i in range(int(nl/3)):
                    s.insert(d - (i+1)*3, sp)
                    
            s = ''.join(s)
            if s.endswith(sp):
                s = s[:-(len(sp))]
                
        return s + ellipses
        
    @classmethod
    def _get_dof(cls,dof1, dof2, u, du1, du2,c):
        # Use the Welch-Satterhwaite approximation to combine dof1 and dof2.
        # See [R. Willink, Metrologia, 44, 340 (2007)] concerning the use
        # of the Welch-Satterwaite approximation with correlated values.
        # Note that we allow the return value to be less than one,
        # but when the dof is retrieved using the dof property get method, 
        # values less than one will be rounded up to one.
        
        if (isinf(dof1) and isinf(dof2)) or u == 0:
            return float('inf')
            
        du1 = du1/u
        du2 = du2/u
            
        # Willink's formula, [R. Willink, Metrologia, 44, 340 (2007)], for
        # the W-S approximations for correlated uncertainties is:
        # d = (du1**2 + c*du1*du2)**2/dof1 + (du2**2 + c*du1*du2)**2/dof2
        # Willink's formula is used in the _apply method but not here.
        
        d = du1**4/dof1 + du2**4/dof2
        d += 2*c**2*(((du1 + du2)**4 - du1**4 - du2**4)/
                     ((1 - 2*c + 2*c**2)*(dof1+dof2)))

        if d == 0:
            return float('inf')
        
        r = 1/d
        if r > cls.max_dof:
            return float('inf')
        
        if r < -1e-6:
            raise ValueError('dof is negative')
        if r < 0:
            r == 1
        
        return r
    
    @classmethod
    def _apply(cls,function,derivative,*args,fxdx=None):
        n = len(args)
        if n == 0:
            return cls(function(),0)
        
        if fxdx is None:
            x = [a.x if isinstance(a,ummy) else a for a in args]
            fx = function(*x)
            d = derivative(*x)
        else:
            fx,d,x = fxdx
            
        if n == 1:
            arg = args[0]
            if not isinstance(arg,ummy):
                return function(arg)
            u = abs(d*arg.u)
            if d == 0:
                dof = float('inf')
            else:
                dof = arg.dof
            r = cls(fx,u,dof=dof)
            if r._ref is not None:
                r._ref = arg._ref
                r._refs = np.sign(d)
            return r
        
        c = ummy.correlation_matrix(args)
        du = np.array([a.u*p if isinstance(a,ummy) else 0 for a,p in zip(args,d)])
        mu = np.max([abs(a) for a in du])
        if mu != 0:
            dun = du/mu
        else:
            dun = du
        u = c.dot(dun).dot(dun)
        if not isnan(u):
            if u < 0:
                if u < -1e6:
                    # something horrible happened, should only get here due to a bug
                    raise FloatingPointError('u is sqrt of ' + str(mu**2*u))
                u = 0
                
            u = u**0.5
            u = mu*u
            
        r = cls(fx,u)
        if u == 0 or isnan(u):
            return r
            
        dm = 0

        for i,a in enumerate(args):
            if isinstance(a,ummy):
                if not isinf(a.dof):
                    # See [R. Willink, Metrologia, 44, 340 (2007)] for this
                    # extension of the W-S approximation to correlated values.
                    dm += np.sum(((c.dot(du)[i]/u)*(du[i]/u))**2)/a.dof
                if a._ref is not None:
                    a._ref.combl(r,a._refs*c[i].dot(du)/u,a._refs*du[i]/u,args)
        if r._ref is not None:
            _GummyRef.check_cor(r)
        if dm > 0:
            r._dof = 1/dm
        return r
        
    @classmethod
    def _napply(cls,function,*args,fxx=None):
        n = len(args)
        if n == 0:
            return cls(function(),0)
        
        d = _der(function,*args)
        
        if fxx is not None:
            fxx = (fxx[0],d,fxx[1])
        
        if n == 1:
            return cls._apply(function,lambda x: d,args[0],fxdx=fxx)

        return cls._apply(function,lambda *x: d,*args,fxdx=fxx)
    
    def _add(self,b):
        if not isinstance(b,ummy):
            r = type(self)(self._x + b, self._u, dof=self._dof)
            if r._ref is not None:
                r._ref = self._ref
                r._refs = self._refs
            return r
    
        c = self.correlation(b)
        x = self._x + b._x
        u = _combu(self._u,b._u,c)
 
        r = type(b)(x,u)
    
        if r._ref is None:
            return r
            
        if self._ref is None:
            r._ref = b._ref
            r._refs = b._refs
            r._dof = b._dof
            return r
        if b._ref is None:
            r._ref = self._ref
            r._refs = self._refs
            r._dof = self._dof
            return r
            
        self._ref.comb(r,self._refs*(b._u*c + self._u)/r._u,self._refs*self._u/r._u)
        b._ref.comb(r,b._refs*(self._u*c + b._u)/r._u,b._refs*b._u/r._u,self._ref)

        r._dof = self._get_dof(self._dof, b._dof, r._u, self._u, b._u,c)
        
        return r
    
    def _radd(self,b):
        r = type(self)(b + self._x, self._u, dof=self._dof)
        if r._ref is not None:
            r._ref = self._ref
            r._refs = self._refs
        return r
    
    def _sub(self,b):           
        if not isinstance(b,ummy):
            r = type(self)(self._x - b, self._u, dof=self._dof)
            if r._ref is not None:
                r._ref = self._ref
                r._refs = self._refs
            return r
            
        c = self.correlation(b)
        x = self._x - b._x
        u = _combu(self._u,-b._u,c)
        
        r = type(b)(x,u)
        if r._ref is None:
            return r
            
        if self._ref is None:
            r._ref = b._ref
            r._refs = -b._refs
            r._dof = b._dof
            return r
        if b._ref is None:
            r._ref = self._ref
            r._refs = self._refs
            r._dof = self._dof
            return r
        
        self._ref.comb(r,self._refs*(-b._u*c + self._u)/r._u,self._refs*self._u/r._u)
        b._ref.comb(r,b._refs*(self._u*c - b._u)/r._u,-b._refs*b._u/r._u,self._ref)

        r._dof = self._get_dof(self._dof, b._dof, r._u, self._u, -b._u, c)
        return r
    
    def _rsub(self,b):
        r = type(self)(b - self._x, self._u, dof=self._dof)
        if r._ref is not None:
            r._ref = self._ref
            r._refs = -self._refs
        return r
    
    def _mul(self,b):
        if not isinstance(b,ummy):
            r = type(self)(self._x*b, abs(self._u*b), dof=self._dof)
            if r._ref is not None:
                r._ref = self._ref
                r._refs = np.sign(b)*self._refs
            return r
                                        
        c = self.correlation(b)
        x = self._x*b._x
        u = _combu(b._x*self._u,self._x*b._u,c)
            
        r = type(b)(x,u)
        
        if r._ref is None:
            return r
        if self._ref is None:
            r._ref = b._ref
            r._refs = np.sign(self._x)*b._refs
            r._dof = b._dof
            return r
        if b._ref is None:
            r._ref = self._ref
            r._refs = np.sign(b._x)*self._refs
            r._dof = self._dof
            return r
        
        self._ref.comb(r,self._refs*(self._x*b._u*c + b._x*self._u)/r._u,
                           self._refs*b._x*self._u/r._u)
        b._ref.comb(r,b._refs*(b._x*self._u*c + self._x*b._u)/r._u,
                        b._refs*self._x*b._u/r._u,self._ref)

        r._dof = self._get_dof(self._dof,b._dof,r._u,b._x*self._u,self._x*b._u,c)
        return r
    
    def _rmul(self,b):
        r = type(self)(b*self._x, abs(b*self._u), dof=self._dof)
        if r._ref is not None:
            r._ref = self._ref
            r._refs = np.sign(b)*self._refs
        return r
    
    def _truediv(self,b):   
        if not isinstance(b,ummy):
            if b == 0:
                raise ZeroDivisionError('division by zero')
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral)):
                x = MFraction(self._x,b)
            else:
                x = self._x/b
            r = type(self)(x, abs(self._u/b), dof=self._dof)
            if r._ref is not None:
                r._ref = self._ref
                r._refs = np.sign(b)*self._refs
            return r
                
        c = self.correlation(b)
        if b._x == 0:
            raise ZeroDivisionError('division by zero')
        else:
            if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral)):
                x = MFraction(self._x,b._x)
            else:
                x = self._x/b._x
            
        u = _combu(self._u/b._x,-self._x*b._u/b._x**2,c)

        r = type(b)(x,u)
        
        if r._ref is None:
            return r
        if self._ref is None:
            r._ref = b._ref
            r._refs = -np.sign(self._x)*b._refs
            r._dof = b._dof
            return r
        if b._ref is None:
            r._ref = self._ref
            r._refs = np.sign(b._x)*self._refs
            r._dof = self._dof
            return r
        
        self._ref.comb(r,self._refs*(-self._x*b._u*c/b._x**2 + self._u/b._x)/r._u,
                           self._refs*self._u/(b._x*r._u))
        b._ref.comb(r,b._refs*(self._u*c/b._x - self._x*b._u/b._x**2)/r._u,
                        b._refs*(-self._x*b._u/(b._x**2*r._u)),self._ref)

        r._dof = self._get_dof(self._dof, b._dof, r._u, self._u/b._x, 
                               -b._u*self._x/b._x**2, c)
        return r
    
    def _rtruediv(self,b):
        if self._x == 0:
            raise ZeroDivisionError('division by zero')
        else:
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral)):
                x = MFraction(b,self._x)
            else:
                x = b/self._x
        u = abs(b*self._u/self._x**2)
        r = type(self)(x,u,dof=self._dof)
        if r._ref is not None:
            r._ref = self._ref
            r._refs = -np.sign(b)*self._refs
        return r
    
    def _pow(self,b):
        if not isinstance(b,ummy):
            if self._x <= 0 and np.modf(b)[0] != 0:
                raise ValueError('a negative or zero value cannot be raised to a non-integer power')
            if b == 0:
                    return type(self)(1)
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral) and b < 0):
                x = MFraction(1,self._x**-b)
            else:
                x = self._x**b
            u = abs(b*self._x**(b-1)*self._u)
            r = type(self)(x,u,dof=self._dof)
            if r._ref is not None:
                r._ref = self._ref
                if self._x < 0:
                    sgn = -((-1)**b)
                else:
                    sgn = 1
                if b < 0:
                    sgn = -sgn
                r._refs = sgn*self._refs
                    
            return r

        if self._x <= 0:
            raise ValueError('a negative or zero value cannot raised to a power which has an uncertainty')            
                
        c = self.correlation(b)
        if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral) and b._x < 0):
            x = MFraction(1,self._x**-b._x)
        else:
            x = self._x**b._x
        da = b._x*self._x**(b._x-1)
        lgx = log(self._x)
        try:
            lgx = type(self._x)(lgx)
        except:
            pass
        db = lgx*self._x**b._x
        u = _combu(da*self._u,db*b._u,c)
        
        r = type(b)(x,u)
        
        if r._ref is None:
            return r
        if self._ref is None:
            r._ref = b._ref
            r._refs = np.sign(self._x)*b._refs
            r._dof = b._dof
            return r
        if b._ref is None:
            r._ref = self._ref 
            if self._x < 0:
                sgn = -((-1)**b._x)
            else:
                sgn = 1
            if b < 0:
                sgn = -sgn
            r._refs = sgn*self._refs
            r._dof = self._dof
            return r
        
        self._ref.comb(r,self._refs*(db*b._u*c + da*self._u)/r._u,
                           self._refs*da*self._u/r._u)
        b._ref.comb(r,b._refs*(da*self._u*c + db*b._u)/r._u,b._refs*db*b._u/r._u,
                        self._ref)

        r._dof = self._get_dof(self._dof, b._dof, r._u, da*self._u, db*b._u, c)
        return r
    
    def _rpow(self,b):
        if b == 0:
            return type(self)(0)
        if (isinstance(self._x,Integral) and 
            isinstance(b,Integral) and self._x < 0):
            x = MFraction(1,b**-self._x)
        else:
            x = b**self._x
        lgb = log(b)
        try:
            lgb = type(b)(lgb)
        except:
            pass
        u = abs(b**self._x*lgb*self._u)
        r = type(self)(x,u,dof=self._dof)
        if r._ref is not None:
            r._ref = self._ref
            r._refs = np.sign(b)*self._refs
        return r
    
    def _nprnd(self,f):
        return type(self)(f(self._x))
        
    def _floordiv(self,b):
        return self._truediv(b)._nprnd(np.floor)
        
    def _rfloordiv(self,b):
        return self._rtruediv(b)._nprnd(np.floor)
        
    def _mod(self,b):
        ret = ummy._apply(np.mod,lambda x1,x2: (1, np.sign(x2)*abs(x1//x2)),self,b)
        return type(self)(ret)
        
    def _rmod(self,b):
        ret = ummy._apply(np.mod,lambda x1,x2: (1, np.sign(x2)*abs(x1//x2)),b,self)
        return type(self)(ret)
        
    def __neg__(self):
        r = self.copy(formatting=False)
        r._x = -self._x
        if self._ref is not None:
            r._refs = -self._refs
        return r
        
    def __pos__(self):
        return self.copy(formatting=False)
        
    def __abs__(self):
        r = self.copy(formatting=False)
        if self._x < 0:
            r._x = -self._x
            if self._ref is not None:
                r._refs = -self._refs
        return r
    
    def __eq__(self,v):
        if self is v:
            return True
        if isinstance(v,ummy):
            if self._u == 0 and v._u == 0:
                return self._x == v._x
            if self._ref is v._ref and self._refs == v._refs:
                if self._x == v._x and self._u == v._u:
                    return True
            return False
        elif self._u == 0:
            return self._x == v
        else:
            return False
        
    def __float__(self):
        return float(self.x)
        
    def __int__(self):
        return int(self.x)
        
    def __complex__(self):
        return complex(self.x)
    
    @property
    def real(self):
        return self
    
    @property
    def imag(self):
        return type(self)(0)
    
    def conjugate(self):
        return self.copy(formatting=False)
    
    def angle(self):
        if self._x >= 0:
            return type(self)(0)
        else:
            return type(self)(np.pi)
        
    @property
    def utype(self):
        """
        `str` or `None`

        An arbitrary string value labeling the uncertainty type.
        """
        if self._ref is None or self._ref._tag is None:
            return None
        return self._ref._tag.name
        
    @staticmethod
    def _toummylist(x):
        if isinstance(x,ummy):
           return [x]
        elif isinstance(x,GummyTag):
            return x.get_values()
        elif isinstance(x,str):
            t = GummyTag.tags.get(x)
            if t is None:
                return []
            x = t.get_values()
            return x
        
        xl = []
        for e in x:
            if isinstance(e,ummy):
               xl.append(e)
            if isinstance(e,GummyTag):
                xl += e.get_values()
            if isinstance(e,str):
                t = GummyTag.tags.get(e)
                if t is not None:
                    xl += t.get_values()
        return xl
        
        
    def ufrom(self,x):
        """
        Gets the standard uncertainty contributed from particular gummys
        or utypes if all other free variables are held fixed.
        
        Parameters
        ----------
        x:  `gummy`, `str`, or array_like
            A gummy, a string referencing a utype or a list containing
            gummys and strings.
            
        Returns
        -------
        `float`
        
        Example
        -------
        >>>  a = gummy(1.2,0.2,utype='A')
        >>>  b = gummy(3.2,0.5,utype='A')
        >>>  c = gummy(0.9,0.2,utype='B')
        >>>  d = a + b + c
        >>>  d.ufrom('A')
        0.53851648071345048
        """
        x = ummy._toummylist(x)
        x = [i for i in x if self.correlation(i) != 0]
        if len(x) == 0:
            return 0

        v = ummy.correlation_matrix(x)
            
        b = [self.correlation(z) for z in x]
        s = np.linalg.lstsq(v,b)[0]
        u = 0
        
        d = [i*self._u/j._u for i,j in zip(s,x)]
        for i in range(len(x)):
            for j in range(len(x)):
                u += d[i]*d[j]*x[i].correlation(x[j])*x[i]._u*x[j]._u
                
        return sqrt(u)
        
    def doffrom(self,x):
        """
        Gets the degrees of freedom contributed from particular gummys or
        utypes if all other free variables are held fixed.
        
        Parameters
        ----------
        x:  `gummy`, `str`, or array_like
            A gummy, a string referencing a utype or a list containing
            gummys and strings.

        Returns
        -------
        `float`

        Example
        -------
        >>>  a = gummy(1.2,0.2,dof=5,utype='A')
        >>>  b = gummy(3.2,0.5,dof=7,utype='A')
        >>>  c = gummy(0.9,0.2,utype='B')
        >>>  d = a + b + c
        >>>  d.doffrom('A')
        9.0932962619709627
        """
        x = ummy._toummylist(x)
        x = [i for i in x if self.correlation(i) != 0]
        if len(x) == 0:
            return float('inf')
            
        v = ummy.correlation_matrix(x)
        b = [self.correlation(z) for z in x]
        s = np.linalg.lstsq(v,b)[0]
        d = [i*self._u/j._u for i,j in zip(s,x)]
        usq = 0
        dm = 0
        for i in range(len(x)):
            for j in range(len(x)):
                usqi = d[i]*d[j]*x[i].correlation(x[j])*x[i]._u*x[j]._u
                usq += usqi
                dm += usqi**2/x[i]._dof
                
        if dm == 0:
            return float('inf')
        dof = usq**2/dm
        if dof > ummy.max_dof:
            return float('inf')
        if dof < 1:
            return 1
        return dof
        
        
class _GummyRef:
    # Instances of this class are hashable unlike gummys, so we can use this in 
    # a gummy's ._cor weak key dictionary and a GummyTag's .values weak set.

    _cortol = 1e-8 # correlations smaller than this are rounded to zero
    _cortolp = 1 - 1e-8 # correlations larger than this are rounded to 1
    _cortoln = -1 + 1e-8 # correlations smaller than this are rounded to -1
    
    def __init__(self):
        self._cor = weakref.WeakKeyDictionary()
        self._tag = None
        self._tag_refs = []
        self._tag_deps = []
        
    def cor(self,ref):
        if ref is self:
            return 1
        c = self._cor.get(ref)
        if c is None:
            return 0
        else:
            return c
            
    def set_cor(self,g,c):
        # We allow input values slighly greater than 1 without raising an 
        # exception to give room for rounding errors (but values greater than
        # 1 will be rounded to 1 below).
        
        if abs(c) > 1.01:
            raise ValueError('abs(c) > 1')
            
        if abs(c) < _GummyRef._cortol:
            if g._ref in self._cor:
                del self._cor[g._ref]
                del g._ref._cor[self]
            return
        
        if c > _GummyRef._cortolp:
            self.copyto(g,1)
            return
            
        if c < _GummyRef._cortoln:
            self.copyto(g,-1)
            return
        
        # We can't let an ummy with a utype die and be collected with the garbage
        # if it is correlated with a living ummy in case ufrom() or doffrom() is 
        # called with the utype as an argument.  So we add a strong reference 
        # pointing to the utyped ummy to the GummyRef of all other ummys that 
        # are correlated with it.
        if self._tag is not None:
            g._ref._tag_deps.append(self.get_tag_ref())
        elif len(self._tag_deps) > 0:
            g._ref._tag_deps += self._tag_deps
        
        self._cor[g._ref] = c
        g._ref._cor[self] = c
            
    def add_cor(self,g,c):
        if c == 0:
            return
        v = self._cor.get(g._ref)
        if v is None:
            v = c
        else:
            v += c
        self._cor[g._ref] = v
        g._ref._cor[self] = v
        
    @staticmethod
    def check_cor(r):
        a = list(r._ref._cor.items())
        for k,v in a:
            if k is not None:
                if abs(v) > 1.01:
                    raise ValueError('abs(correlation) > 1')
                if abs(v) < _GummyRef._cortol:
                    del k._cor[r._ref]
                    del r._ref._cor[k]
                    return
                if v > _GummyRef._cortolp:
                    k.copyto(r,1)
                    return
                if v < _GummyRef._cortoln:
                    k.copyto(r,-1)
                    return
                
    def copyto(self,g,refs=1):
        g._ref = self
        g._refs = refs
        if self._tag is not None:
            self._tag_refs.append(weakref.ref(g))
                            
    def comb(self,g,c1,c2,r2=None):
        self.set_cor(g,c1)
        if g._ref is not self:
            a = list(self._cor.items())
            for k,v in a:
                if k is not None and k is not g._ref and k is not r2:
                    k.add_cor(g,c2*v)
                    
            # comb is always called twice, the second time with r2 not None.
            # We must hold off the check till the end of the second call since
            # there may be cancellations which assure that the correlation
            # coefficients are less than 1.
            if r2 is not None:
                _GummyRef.check_cor(g)
        
    def combl(self,g,c1,c2,gl):
        # This is a generalization of comb() (above) called by gummy.apply() for
        # functions of more than two variables.
        self.set_cor(g,c1)
        if g._ref is not self:
            a = list(self._cor.items())
            for k,v in a:
                if k is not None and k is not g._ref:
                    ad = True
                    for s in gl:
                        if isinstance(s,ummy) and s._ref is k:
                            ad = False
                    if ad:
                        k.add_cor(g,c2*v)
                        
    def get_tag_ref(self):
        try:
            ret = self._tag_refs[-1]()
            while ret is None:
                # The GummyRef may be shared by several copys of the original
                # ummy.  Find one that is still alive.
                self._tag_refs.pop()
                ret = self._tag_refs[-1]()
            return ret
        except IndexError:
            return None
            

class GummyTag:
    # Do not create GummyTag objects directly.  A GummyTag instance is created
    # each time a gummy instance with a new and unique utype parameter is created.

    tags = weakref.WeakValueDictionary()
        
    @staticmethod
    def set_tag(g,tag):
        if isinstance(tag,str):
            tag = tag.strip()
            if tag in GummyTag.tags:
                tag = GummyTag.tags[tag]
            else:
                tag = GummyTag(tag)

        elif not isinstance(tag,GummyTag):
            raise ValueError('utype must be either None or a str')
            
        tag.values.add(g._ref)
        g._ref._tag = tag
        g._ref._tag_refs.append(weakref.ref(g))
        
    def __init__(self,tag_name):
        GummyTag.tags[tag_name] = self
        
        self.values = weakref.WeakSet()
        
        self.name = tag_name
        
    def get_values(self):
        vals = list(self.values)

        # Can't put this in a list comprehension since v.ref is a weakref and
        # it could die in the middle of the loop.
        ret = []
        for v in vals:
            if v is not None:
                r = v.get_tag_ref()
                if r is not None:
                    ret.append(r)
                
        return ret
    
                        
def _der(function,*args):
    # computes the numerical derivative of function with respect to args.
    # puts zeros for any of args that are not an ummy
    n = len(args)
    if n == 1:
        arg = args[0]
        if not isinstance(arg,ummy) or arg._u == 0:
            return 0
                
        df = None
        a = np.empty([8,8])
        x = arg.x
        s = 2*arg.u
        a[0,0] = (function(x+s)-function(x-s))/(2*s)
        for k in range(1,8):
            s = s/1.4
            a[0,k] = (function(x+s)-function(x-s))/(2*s)
            w = 1.96
            for j in range(1,k):
                a[j,k] = (a[j-1,k]*w-a[j-1,k-1])/(w-1)
                w = 1.96*w
                dff = max(abs(a[j,k]-a[j-1,k]),abs(a[j,k]-a[j-1,k-1]))
                if df is None or dff <= df:
                    df = dff
                    d = a[j,k]
            if df is not None and abs(a[k,k]-a[k-1,k-1]) > 2*df:
                break
            
        return d
        
    v = np.empty(n)
    for i,a in enumerate(args):
        if isinstance(a,ummy):
            v[i] = a.x
        else:
            v[i] = a
    d = np.zeros(n)
    for i,p in enumerate(args):
        if isinstance(p,ummy) and p._u > 0:          
            df = None
            s = 2*p.u
            x1 = np.array(v)
            x1i = x1[i]
            x1[i] = x1i+s
            x2 = np.array(v)
            x2i = x2[i]
            x2[i] = x2i-s
            a = np.empty([8,8])
            a[0,0] = (function(*x1)-function(*x2))/(2*s)
            for k in range(1,8):
                s = s/1.4
                x1[i] = x1i+s
                x2[i] = x2i-s
                a[0,k] = (function(*x1)-function(*x2))/(2*s)
                w = 1.96
                for j in range(1,k):
                    a[j,k] = (a[j-1,k]*w-a[j-1,k-1])/(w-1)
                    w = 1.96*w
                    dff = max(abs(a[j,k]-a[j-1,k]),abs(a[j,k]-a[j-1,k-1]))
                    if df is None or dff <= df:
                        df = dff
                        di = a[j,k]
                if df is not None and abs(a[k,k]-a[k-1,k-1]) > 2*df:
                    break
            d[i] = di
            
    return d
        
class MFraction(Fraction):
    """
    A fraction.Fraction sub-class that works with mpmath.mpf objects
    """
    def _mpmath_(self,p,r):
        return rational.mpq(self.numerator,self.denominator)
