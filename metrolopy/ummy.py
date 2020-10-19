# -*- coding: utf-8 -*-

# module ummy

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
The ummy object defined here was created primarily to separate the code for
first-order uncertainty propagation from the other methods in the gummy object.  
But the ummy object can be used as a stripped down verison of gummy without the
Monte-Carlo uncertainty propagation, unit, and pretty printing functionality
of gummy.
"""

import numpy as np
import weakref
from collections import namedtuple
from warnings import warn
from math import isnan,isinf,log,log10,pi
from fractions import Fraction
from numbers import Rational,Integral,Real,Complex
from .printing import PrettyPrinter,MetaPrettyPrinter
from .dfunc import Dfunc
from .exceptions import UncertiantyPrecisionWarning
from .unit import Unit,one,Quantity

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
            return (_iinfo,x)
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

def _sign(x):
    if x < 0:
        return -1
    return 1

def _format_exp(fmt,xp):
    if fmt == 'html':
        ex = ' &times; 10<sup>' + str(xp) + '</sup>'
    elif fmt == 'latex':
        ex = r' \times\,10^{' + str(xp) + '}'
    else:
        ex = 'e'
        if xp > 0:
            ex += '+'
        ex += str(xp)
    return ex


class ummy(Dfunc,PrettyPrinter):
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
    # rounding errors
    rounding_u = False
    
    max_digits = 15
    
    def __init__(self,x,u=0,dof=float('inf'),utype=None):
        if isinstance(x,Quantity):
            x = x.convert(one).value
            
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
            
        if u == 0:
            dof = float('inf')
            utype = None
            
        if isinstance(x,Fraction):
            x = MFraction(x)
        if isinstance(u,Fraction):
            u = MFraction(u)
                
        if self.rounding_u and not isinstance(x,Rational):
            finfo,x = _getfinfo(x)
            if finfo is _iinfo:
                warn('numpy.finfo cannot get the floating point accuracy for x\nno uncertainty will be included to account for floating point rounding errors',UncertiantyPrecisionWarning)
            else:
                u = _combu(u,float(abs(x)*finfo.rel_u),0)
            
            #uinfo = _getfinfo(u)[0]
            #if uc < uinfo.warn_u and u is not 0:
                #warn('a gummy/ummy was defined with an uncertainty too small to be accurately represented by a float\nuncertainty values may contain significant numerical errors',UncertiantyPrecisionWarning)
            
            self._finfo = finfo
        else:
            self._finfo = None

            
        self._x = x
        self._u = u
        
        self._refs = 1
        if not isinf(u) and not isnan(u) and u > 0:
            if isinstance(dof,_GummyRef):
                self._ref = dof
                self._refs = utype
            else:
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
                    dof = float('inf')
                    
                self._ref = _GummyRef(dof)
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
        ummy is based on.  If the ummy was created as the result of an 
        operation between two or more other gummys, then the dof is the effective
        number of degrees of freedom calculated using the Welch-Satterthwaite 
        approximation.  Caution:  A variation of the the  Welch-Satterthwaite
        approximation is used that takes into account correlations, see
        [R. Willink, Metrologia, 44, 340 (2007)].  However correlations are
        not handled perfectly.  So if accurate dof calculations are need, care
        should be taken to ensure that correlations are not generated in
        intermediate calculations.
        """
        if self._ref is None:
            return float('inf')
        if self._ref.dof < 1:
            # Occasionally correlations can result in a _dof less than 1;
            # see the _get_dof method.
            return 1
        return self._ref.dof
    
    @property
    def _dof(self):
        if self._ref is None:
            return float('inf')
        return self._ref.dof
    @_dof.setter
    def _dof(self,v):
        if self._ref is not None:
            self._ref.dof = v
    
    @property
    def finfo(self):
        if self._finfo is None:
            self._finfo = _getfinfo(self._x)[0]
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
        r._ref = s._ref
        r._refs = s._refs
        r._finfo = s._finfo

        if formatting:
            if s.nsig != type(r).nsig:
                r.nsig = s.nsig
            if s.thousand_spaces != type(r).thousand_spaces:
                r.thousand_spaces = s.thousand_spaces
    
    def copy(self,formatting=True,tofloat=False):
        """
        Returns a copy of the ummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new gummy.  The default for `formatting` is `True`.  If tofloat
        is true the x and u properties will be converted to float values
        before copying.
        """
        r = type(self).__new__(type(self))
        self._copy(self,r,formatting=formatting,tofloat=tofloat)
        return r
        
    def tofloat(self):
        """
        Returns a copy of the gummy with x an u converted to floats.
        """
        return self.copy(formatting=False,tofloat=True)
    
    def splonk(self):
        """
        returns self.x if u == 0 else returns self
        """
        if self.u == 0:
            return self.x
        return self
    
    def correlation(self, g):
        """
        Returns the correlation coefficient between `self` and `g`.
        """
        if not isinstance(g,ummy):
            return 0
        if self._ref is None or g._ref is None:
            return 0
        return (self._refs*g._refs)*self._ref.cor(g._ref)
        
    def covariance(self,g):
        """
        Returns the covariance between `self` and `g`.
        """
        if not isinstance(g,ummy):
            return 0
        if self._ref is None or g._ref is None:
            return 0
        return self._u*g._u*(self._refs*g._refs)*self._ref.cor(g._ref)
        
    @staticmethod
    def correlation_matrix(gummys):
        """
        Returns the correlation matrix of a list or array of ummys.
        """
        return [[b.correlation(a) if isinstance(b,ummy) else 0 
                for b in gummys] for a in gummys]
        
    @staticmethod
    def covariance_matrix(gummys):
        """
        Returns the variance-covariance matrix of a list or array of ummys.
        """
        return [[b.covariance(a) if isinstance(b,ummy) else 0 for b in gummys] 
                for a in gummys]
                                
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
        try:
            assert len(matrix) == n
            for r in matrix:
                assert len(r) == n
        except:
            raise ValueError('the covariance matrix must be a square matrix of the same length as gummys')
        
        # Set the u of each gummy using the diagonal of the matrix.
        for i in range(n):
            if matrix[i][i] == 0:
                gummys[i]._u = 0
                gummys[i]._ref = None
            else:
                if gummys[i]._ref is None:
                    gummys[i]._ref = _GummyRef()
                gummys[i]._u = matrix[i][i]**0.5
        
        # Set the correlations
        m = [[e/(gummys[i]._u*gummys[j]._u) for j,e in enumerate(r)]
             for i,r in enumerate(matrix)]
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
        
    def tostring(self,fmt='unicode',nsig=None,**kwds):
        if nsig is None:
            nsig = self.nsig
            
        x = float(self._x)
        try:
            xabs = abs(x)
           
            if xabs != 0:
                xexp = _floor(log10(xabs))
            else:
                xexp = None
                
            u = float(self.u)
            
            if self._u == 0 or isnan(self._u) or isinf(self._u):
                if xexp is None:
                    return '0'
                else:
                    if (((self.sci_notation is None and (xexp > self.sci_notation_high or xexp < self.sci_notation_low)) 
                        or self.sci_notation) ):
                        xtxt = self._format_mantissa(fmt,x*10**(-xexp),
                                                     self.finfo.precision)
                        etxt = _format_exp(fmt,xexp)
                        return xtxt + etxt
                    else:
                        return self._format_mantissa(fmt,x,
                                                     self.finfo.precision)
                    
            uabs = abs(u)
            if uabs != 0:
                uexp = _floor(log10(uabs))
                
            if xexp is None:
                xp = uexp
            else:
                xp = max(xexp,uexp)
                
            psn = (uexp - nsig + 1 > 0)
            
            if (((self.sci_notation is None and (xp > self.sci_notation_high or xp < self.sci_notation_low)) 
                        or self.sci_notation) or psn or (xexp is not None and (uexp > xexp and xexp == 0))):
                xtxt = self._format_mantissa(fmt,x*10**(-xp),-uexp+xp+nsig-1)
                etxt = _format_exp(fmt,xp)
            else:
                xtxt = self._format_mantissa(fmt,x,-uexp+nsig-1)
                etxt = ''
                
            if nsig <= 0:
                return xtxt + etxt
    
            utxt = self._format_mantissa(fmt,u*10**(-uexp),nsig-1,parenth=True)
            return xtxt + '(' + utxt + ')' + etxt
        except:
            try:
                return(str(self.x) + '{' + str(self.u) + '}' + '??')
            except:
                try:
                    return(str(self.x) + '{??}')
                except:
                    return('??')
          
    def _format_mantissa(self,fmt,x,sig,parenth=False):            
        try:
            if isnan(x):
                return 'nan'
            if isinf(x) and x > 0:
                if fmt == 'html':
                    return '&infin;'
                if fmt == 'latex':
                    return r'\infty'
                if fmt == 'ascii':
                    return 'inf'
                return '\u221E'
            if isinf(x) and x < 0:
                if fmt == 'html':
                    return '-&infin;'
                if fmt == 'latex':
                    return r'-\infty'
                if fmt == 'ascii':
                    return '-inf'
                return '-\u221E'
        except:
            pass
        
        ellipses = ''
        if isinstance(x,Rational) and sig is None:
            if x == int(x):
                s = str(int(x))
            else:
                xf = float(x)
                e = 15 - _floor(log10(abs(xf)))
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
                    em = _floor(log10(abs(float(x))))
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
    def _get_dof(cls,dof1,dof2,du1,du2,c):
        # Use the Welch-Satterhwaite approximation to combine dof1 and dof2.
        # See [R. Willink, Metrologia, 44, 340 (2007)] concerning the use
        # of the Welch-Satterwaite approximation with correlated values.
        # Note that we allow the return value to be less than one,
        # but when the dof is retrieved using the dof property get method, 
        # values less than one will be rounded up to one.
            
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
        if fxdx is None:
            x = [a.x if isinstance(a,ummy) else a for a in args]
            fx = function(*x)
            d = derivative(*x)
        else:
            fx,d,x = fxdx
        if len(args) == 1:
            d = [d]
            
        ad = [(a,a.u*p) for a,p in zip(args,d) 
              if isinstance(a,ummy) and a._u != 0]
        if len(ad) == 0:
            return cls(fx)
        args,du = list(map(list, zip(*ad)))
            
        if len(args) == 1:
            refs = -args[0]._refs if du[0] < 0 else args[0]._refs
            r = cls(fx,u=abs(du[0]),dof=args[0]._ref,utype=refs)
            return r
        
        maxdu = max(abs(a) for a in du)
        if maxdu == 0:
            return cls(fx)
        dun = [d/maxdu for d in du]
        
        u = sum(d**2 for d in dun)
        u += 2*sum(args[i].correlation(args[j])*dun[i]*dun[j] 
                 for i in range(len(args)) for j in range(i+1,len(args)))
        
        if not isnan(u):
            if u < 0:
                if u < -1e6:
                    # something horrible happened, should only get here due to a bug
                    raise FloatingPointError('u is sqrt of ' + str(maxdu**2*u))
                u = 0
                
            u = u**0.5
            u = maxdu*u
            
        if u == 0 or isnan(u):
            return cls(fx)
        
        du = [d/u for d in du]
        dof = float('inf')
        if any(not isinf(a._dof) for a in args):
            dm = sum(sum(args[i].correlation(args[j])*du[i]*du[j] 
                         for i in range(len(du)))**2/args[j]._dof 
                         for j in range(len(du)))
            if dm > 0:
                dof = 1/dm
                
        r = cls(fx,u=u,dof=dof)
        
        rl = {a._ref for a in args}
        for i,a in enumerate(args):
            s = sum(a.correlation(args[j])*du[j] for j in range(len(du)))
            a._ref.combl(r,float(a._refs*s),float(a._refs*du[i]),rl)
            if r._ref in rl:
                break
        r._check_cor()
            
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
    
    def _check_cor(self):
        if self._ref is None:
            return
        
        rl = list(self._ref._cor.items())
        for k,v in rl:
            if k is not None:
                if abs(v) > 1.01:
                    raise ValueError('abs(correlation) > 1')
                if abs(v) < _GummyRef._cortol:
                    del k._cor[self._ref]
                    del self._ref._cor[k]
                    return
                if v > _GummyRef._cortolp:
                    k.copyto(self,1)
                    return
                if v < _GummyRef._cortoln:
                    k.copyto(self,-1)
                    return
                
    def _combc(self,a,b,dua,dub,c):
        dua = float(dua)
        dub = float(dub)

        a._ref.combl(self,c*dub*a._refs + dua*a._refs,dua*a._refs)
        if self._ref is not a._ref:
            b._ref.combl(self,c*dua*b._refs + dub*b._refs,dub*b._refs,[a._ref])
            if self._ref is not b._ref:
                self._check_cor()   
                
    def __add__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) + b
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__radd__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) + b
        
        if not isinstance(b,ummy):
            return type(self)(self._x + b,self._u,dof=self._ref,utype=self._refs)
    
        c = self.correlation(b)
        x = self._x + b._x
        u = _combu(self._u,b._u,c)
        
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            return type(b)(x,u,dof=b._ref,utype=b._refs)
        
        if b._u == 0:
            return type(b)(x,u,dof=self._ref,utype=self._refs)
        
        dua = self._u/u
        dub = b._u/u
        if not isinf(self._ref.dof) or not isinf(b._ref.dof):
            dof = self._get_dof(self._dof,b._dof,dua,dub,c)
        else:
            dof = float('inf')
            
        r = type(b)(x,u=u,dof=dof)
    
        r._combc(self,b,dua,dub,c)
        
        return r
    
    def __radd__(self,b):
        if isinstance(b,np.ndarray):
            return b + np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b + immy(self)
        
        if isinstance(b,ummy):
            return b.__add__(self)
        
        return type(self)(self._x + b,self._u,dof=self._ref,utype=self._refs)
    
    def __sub__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) - b
          
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rsub__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) - b
        
        if not isinstance(b,ummy):
            return type(self)(self._x - b,self._u,dof=self._ref,utype=self._refs)
            
        c = self.correlation(b)
        x = self._x - b._x
        u = _combu(self._u,-b._u,c)
            
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            return type(b)(x,u,dof=b._ref,utype=-b._refs)
        
        if b._u == 0:
            return type(b)(x,u,dof=self._ref,utype=self._refs)
        
        dua = self._u/u
        dub = -b._u/u
        if not isinf(self._ref.dof) or not isinf(b._ref.dof):
            dof = self._get_dof(self._dof,b._dof,dua,dub,c)
        else:
            dof = float('inf')
        
        r = type(b)(x,u=u,dof=dof)
        
        r._combc(self,b,dua,dub,c)

        return r
    
    def __rsub__(self,b):
        if isinstance(b,np.ndarray):
            return b - np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b - immy(self)
        
        if isinstance(b,ummy):
            return b.__sub__(self)
        
        r = type(self)(b - self._x,self._u,dof=self._ref,utype=-self._refs)
        return r
    
    def __mul__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self)*b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmul__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self)*b
        
        if not isinstance(b,ummy):
            refs = -self._refs if b < 0 else self._refs
            return type(self)(self._x*b,abs(self._u*b),dof=self._ref,utype=refs)
                                        
        c = self.correlation(b)
        x = self._x*b._x
        dua = b._x*self._u
        dub = self._x*b._u
        u = _combu(dua,dub,c)
            
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            refs = -b._refs if self._x < 0 else b._refs
            return type(b)(x,u,dof=b._ref,utype=refs)
        
        if b._u == 0:
            refs = -self._refs if b._x < 0 else self._refs
            return type(b)(x,u,dof=self._ref,utype=refs)
        
        dua = dua/u
        dub = dub/u
        if not isinf(self._ref.dof) or not isinf(b._ref.dof):
            dof = self._get_dof(self._dof,b._dof,dua,dub,c)
        else:
            dof = float('inf')
            
        r = type(b)(x,u=u,dof=dof)
        
        r._combc(self,b,dua,dub,c)

        return r
    
    def __rmul__(self,b):
        if isinstance(b,np.ndarray):
            return b*np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b*immy(self)
        
        if isinstance(b,ummy):
            return b.__mul__(self)
        
        refs = -self._refs if b < 0 else self._refs
        return type(self)(self._x*b,abs(self._u*b),dof=self._ref,utype=refs)
    
    def __truediv__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self)/b
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rtruediv__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self)/b
        
        if not isinstance(b,ummy):
            if b == 0:
                raise ZeroDivisionError('division by zero')
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral)):
                x = MFraction(self._x,b)
            else:
                x = self._x/b
            refs = -self._refs if b < 0 else self._refs
            return type(self)(x,abs(self._u/b),dof=self._ref,utype=refs)
                
        c = self.correlation(b)
        if b._x == 0:
            raise ZeroDivisionError('division by zero')
        else:
            if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral)):
                x = MFraction(self._x,b._x)
            else:
                x = self._x/b._x
            
        dua = self._u/b._x
        dub = -(self._x*b._u/b._x)/b._x
        
        u = _combu(dua,dub,c)
        
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            refs = b._refs if self._x < 0 else -b._refs
            return type(b)(x,u,dof=b._ref,utype=refs)
        
        if b._u == 0:
            refs = -self._refs if b._x < 0 else self._refs
            return type(b)(x,u,dof=self._ref,utype=refs)
        
        dua = dua/u
        dub = dub/u
        if not isinf(self._ref.dof) or not isinf(b._ref.dof):
            dof = self._get_dof(self._dof,b._dof,dua,dub,c)
        else:
            dof = float('inf')
            
        r = type(b)(x,u=u,dof=dof)
        
        r._combc(self,b,dua,dub,c)

        return r
    
    def __rtruediv__(self,b):
        if isinstance(b,np.ndarray):
            return b/np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b/immy(self)
        
        if isinstance(b,ummy):
            return b.__div__(self)
        
        if self._x == 0:
            raise ZeroDivisionError('division by zero')
        else:
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral)):
                x = MFraction(b,self._x)
            else:
                x = b/self._x
        u = abs(b*self._u/self._x**2)
        refs = self._refs if b < 0 else -self._refs
        return type(self)(x,u,dof=self._ref,utype=refs)
    
    def __pow__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self)**b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rpow__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self)**b
        
        if isinstance(b,ummy) and b._u == 0:
            b = b._x
        if not isinstance(b,ummy):
            if b == 0:
                    return type(self)(1)
            if (isinstance(self._x,Integral) and 
                isinstance(b,Integral) and b < 0):
                x = MFraction(1,self._x**-b)
            else:
                x = self._x**b
            u = abs(b*self._x**(b-1)*self._u)
            if self._x < 0:
                sgn = int(-((-1)**b))
            else:
                sgn = 1
            if b < 0:
                sgn = -sgn
            return type(self)(x,u,dof=self._ref,utype=sgn*self._refs)

        if self._x <= 0:
            raise ValueError('a negative or zero value cannot raised to a power which has an uncertainty')            
                
        c = self.correlation(b)
        if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral) and b._x < 0):
            x = MFraction(1,self._x**-b._x)
        else:
            x = self._x**b._x
        dua = b._x*self._x**(b._x-1)*self._u
        lgx = log(self._x)
        
        # don't exactly remember why I did this,
        # but it causes problems if x is int.
        #try:
            #lgx = type(self._x)(lgx)
        #except:
            #pass
            
        dub = lgx*self._x**b._x*b._u
        u = _combu(dua,dub,c)
        
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            refs = -b._refs if self._x < 0 else b._refs
            return type(b)(x,u,dof=b._ref,utype=refs)
        
        dua = dua/u
        dub = dub/u
        
        if not isinf(self._ref.dof) or not isinf(b._ref.dof):
            dof = self._get_dof(self._dof,b._dof,dua,dub,c)
        else:
            dof = float('inf')
            
        r = type(b)(x,u=u,dof=dof)
        
        r._combc(self,b,dua,dub,c)

        return r
    
    def __rpow__(self,b):
        if isinstance(b,np.ndarray):
            return b**np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b**immy(self)
        
        if isinstance(b,ummy):
            return b.__pow__(self)
        
        if b == 0:
            return type(self)(0)
        if (isinstance(self._x,Integral) and 
            isinstance(b,Integral) and self._x < 0):
            x = MFraction(1,b**-self._x)
        else:
            x = b**self._x
            
        if self._u == 0:
            type(self)(x)
            
        if b < 0:
            raise ValueError('a negative number cannot raised to a power which has an uncertainty')
            
        lgb = log(b)
        
        # don't exactly remember why I did this,
        # but it causes problems if x is int.
        #try:
            #lgb = type(b)(lgb)
        #except:
            #pass
            
        u = abs(b**self._x*lgb*self._u)
        refs = -self._refs if b < 0 else self._refs
        return type(self)(x,u,dof=self._ref,utype=refs)
    
    def _nprnd(self,f):
        return type(self)(f(self._x))
        
    def __floordiv__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) // b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmfloordiv__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) // b
        
        return self._truediv(b)._nprnd(_floor)
        
    def __rfloordiv__(self,b):
        if isinstance(b,np.ndarray):
            return b // np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b // immy(self)
        
        if isinstance(b,ummy):
            return b.__floordiv__(self)
        
        return self._rtruediv(b)._nprnd(_floor)
        
    def __mod__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) % b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmod__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) % b
        
        ret = ummy._apply(lambda x1,x2: x1%x2,
                          lambda x1,x2: (1, _sign(x2)*abs(x1//x2)),self,b)
        return type(self)(ret)
        
    def __rmod__(self,b):
        if isinstance(b,np.ndarray):
            return b % np.array(self)
        
        if isinstance(b,ummy):
            return b.__mod__(self)
        
        ret = ummy._apply(lambda x1,x2: x1%x2,
                          lambda x1,x2: (1, _sign(x2)*abs(x1//x2)),b,self)
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
            if self.u == 0 and v.u == 0:
                return self._x == v._x
            if self._ref is v._ref and self._refs == v._refs:
                return self.x == v.x and self.u == v.u
            return False
        
        if self.u == 0:
            return self.x == v
        
        return False
        
    #def __float__(self):
        #return float(self.x)
        
    #def __int__(self):
        #return int(self.x)
        
    #def __complex__(self):
        #return complex(self.x)
    
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
            return type(self)(pi)
        
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
        Gets the standard uncertainty contributed from particular ummys
        or utypes if all other free variables are held fixed.
        
        Parameters
        ----------
        x:  `ummy`, `str`, or array_like
            A ummy, a string referencing a utype or a list containing
            ummys and strings.
            
        Returns
        -------
        `float`
        
        Example
        -------
        >>>  a = ummy(1.2,0.2,utype='A')
        >>>  b = ummy(3.2,0.5,utype='A')
        >>>  c = ummy(0.9,0.2,utype='B')
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
        s = np.linalg.lstsq(v,b,rcond=None)[0]
        u = 0
        
        d = [i*self.u/j.u for i,j in zip(s,x)]
        for i in range(len(x)):
            for j in range(len(x)):
                u += d[i]*d[j]*x[i].correlation(x[j])*x[i].u*x[j].u
                
        return u**0.5
        
    def doffrom(self,x):
        """
        Gets the degrees of freedom contributed from particular ummys or
        utypes if all other free variables are held fixed.  Caution:  any
        correlations in the calculations can cause errors in dof calculations.
        
        Parameters
        ----------
        x:  `ummy`, `str`, or array_like
            A ummy, a string referencing a utype or a list containing
            ummys and strings.

        Returns
        -------
        `float`

        Example
        -------
        >>>  a = ummy(1.2,0.2,dof=5,utype='A')
        >>>  b = ummy(3.2,0.5,dof=7,utype='A')
        >>>  c = ummy(0.9,0.2,utype='B')
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
        s = np.linalg.lstsq(v,b,rcond=None)[0]
        d = [i*self.u/j.u for i,j in zip(s,x)]
        usq = 0
        dm = 0
        for i in range(len(x)):
            for j in range(len(x)):
                usqi = d[i]*d[j]*x[i].correlation(x[j])*x[i].u*x[j].u
                usq += usqi
                dm += usqi**2/x[i].dof
                
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

    _cortol = 1e-15 # correlations smaller than this are rounded to zero
    _cortolp = 1 - 1e-15 # correlations larger than this are rounded to 1
    _cortoln = -1 + 1e-15 # correlations smaller than this are rounded to -1
    
    def __init__(self,dof=float('inf')):
        self._cor = weakref.WeakKeyDictionary()
        self._tag = None
        self._tag_refs = []
        self._tag_deps = []
        self.dof = dof
        
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
            
    def combl(self,g,c1,c2,rl=None):
        self.set_cor(g,c1)
        if g._ref is not self:
            a = list(self._cor.items())
            for k,v in a:
                if k is not None and k is not g._ref:
                    if rl is None or not any(i is k for i in rl):
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


class MetaImmy(MetaPrettyPrinter):
    @property
    def style(cls):
        """
        `str` in {'polar','cartesian'}
        
        Sets whether to print the value using a cartesian or polar 
        representaion.  This property can be set at the class or instance
        level.  The default is 'cartesian'.
        """
        return cls._style
    @style.setter
    def style(cls,value):
        value = value.lower().strip()
        if value in ['cartesian','ri','real-imag']:
            cls._style = 'cartesian'
        elif value in ['polar','r-phi']:
            cls._style = 'polar'
        else:
            raise ValueError('the style must be either "cartesian" or "polar"')
        
    @property
    def imag_symbol(cls):
        """
        `str`
        
        The symbol for the unit imaginary number.  The default is 'j'.
        """
        return immy._imag_symbol
    @imag_symbol.setter
    def imag_symbol(cls,value):
        immy._imag_symbol = str(value)
    
class immy(PrettyPrinter,Dfunc,metaclass=MetaImmy):
    
    _style = 'cartesian'
    _imag_symbol = 'j'
    _element_type = ummy # in the jummy subclass this is set to gummy
    
    def __init__(self,real=None,imag=None,r=None,phi=None,cov=None):
        """
        An immy object represents a complex valued quantity with `ummy`
        real and imaginary components.
        
        Parameters
        ----------
        real, imag, r, phi:  `float` or `ummy`
            The value may be specified in  either cartesian coordinates
            using `real` and `imag` or polar coordinates with `r` and `phi`.
            The pair `real`, `imag` or `r`, `phi` may both be `ummy` or
            both be `float`.  If they are `float` then `cov` may
            also be specified.
                
        cov:  2 x 2 array_like of `float`, optional
            The variance-covariance matrix for either the pair `real`,
            `imag` or the pair `r`, `phi`.
        """
        if isinstance(real,immy):
            if real._ridef:
                self._real = self._element_type(real._real)
                self._imag = self._element_type(real._imag)
                self._r = None
                self._phi = None
                self._ridef = True
                return
            else:
                self._r = self._element_type(real._r)
                self._phi = self._element_type(real._phi)
                self._real = None
                self._imag = None
                self._ridef = False
                return
        
        if cov is not None:
            cov = np.asarray(cov)
        
        if real is not None:
            self._ridef = True
            if r is not None or phi is not None:
                raise ValueError('r and phi may not be specified if real is specified')
            
            if imag is None:
                imag = real.imag
                real = real.real
                    
            if cov is not None and (isinstance(real,(ummy,self._element_type)) or 
                                    isinstance(imag,(ummy,self._element_type))):
                raise ValueError('cov may not be specified if real or imag has an uncertainty')
    
            if cov is None:
                self._real = self._element_type(real)
                self._imag = self._element_type(imag)
            else:
                if (isinstance(real,(ummy,self._element_type)) or 
                    isinstance(imag,(ummy,self._element_type))):
                    raise ValueError('cov may not be specified if real or imag ihas an uncertainty')
            
            try:
                self._real,self._imag = self._element_type.create([real,imag],
                                                                  covariance_matrix=cov)
            except ValueError as e:
                if str(e).startswith('matrix must have shape'):
                    raise ValueError('cov must be a 2 x 2 matrix')
                raise
            self._r = None
            self._phi = None
        else:
            self._ridef = False
            if phi is None or r is None:
                raise ValueError('if real is not specified then both r and phi must be specified')
            
            if (isinstance(r,(ummy,self._element_type)) or 
                isinstance(phi,(ummy,self._element_type))):
                if cov is not None:
                    raise ValueError('cov may not be specified if r or phi is a gummy')
                self._r = self._element_type(r)
                self._phi = self._element_type(phi)
            else: 
                try:
                    self._r,self._phi = self._element_type.create([r,phi],
                                                                  covariance_matrix=cov)
                except ValueError as e:
                    if e[0].startswith('matrix must have shape'):
                        raise ValueError('cov must be a 2 x 2 matrix')
                    raise
            self._real = None 
            self._imag = None 
            
    @property
    def x(self):
        """
        Returns `complex(self.real.x,self.imag.x)`, read-only
        """
        return complex(self.real.x,self.imag.x)
    
    @property
    def cov(self):
        """
        Returns the variance-covariance matrix between the real and imaginary
        parts of the value, read-only.
        """
        return ummy.covariance_matrix([self.real,self.imag])
    
    @property
    def real(self):
        """
        read-only
        Returns real part of the value.
        """
        if self._real is None: 
            self._real = self._r*np.cos(self._phi)
        return self._real
    
    @property
    def imag(self):
        """
        read-only
        Returns the imaginary part of the value.
        """
        if self._imag is None:
            self._imag = self._r*np.sin(self._phi)
        return self._imag
    
    @property
    def r(self):
        """
        read-only
        Returns the magnitude of the value.
        """
        if self._r is None:
            self._r = abs(self)
        return self._r
    
    @property
    def phi(self):
        """
        read-only
        Returns the polar angle of the value.
        """
        if self._phi is None:
            self._phi = self.angle()
        return self._phi
    
    def conjugate(self):
        """
        Returns the complex conjugate.
        """
        return type(self)(real=self.real,imag=-self.imag)
    
    def angle(self):
        """
        Returns the polar angle in radians.
        """
        return np.arctan2(self.imag,self.real)
    
    def copy(self,formatting=True,tofloat=False):
        """
        Returns a copy of the jummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new jummy.  The default for `formatting` is `True`.  If the
        tofloats parameter is True x and u for both the real and
        imaginary components will be converted to floats.
        """
        if self._ridef:
            r = self.real.copy(formatting=formatting,tofloat=tofloat)
            i = self.imag.copy(formatting=formatting,tofloat=tofloat)
            return type(self)(real=r,imag=i)
        else:
            r = self.r.copy(formatting=formatting,tofloat=tofloat)
            phi = self.phi.copy(formatting=formatting,tofloat=tofloat)
            return type(self)(r=r,phi=phi)
    
    def tofloat(self):
        """
        Returns a copy of the gummy with x an u (for both the real and
        imaginary components) converted to floats.
        """
        return self.copy(formatting=False,tofloat=True)
    
    def splonk(self):
        if self.real.u == 0 and self.imag.u == 0:
            return self.x
        if self.imag == 0:
            return self.real
        return self
            
    @classmethod
    def _apply(cls,function,derivative,*args,fxx=None,rjd=None):
        
        n = len(args)
        if fxx is None:
            rargs = [a.real for a in args]
            jargs = [a.imag for a in args]
            args = rargs + jargs
            args = [a.convert(one).value if isinstance(a,Quantity) else a for 
                    a in args]
            x = [a.x if isinstance(a,ummy) else a for a in args]
            func = lambda *a: function(*[complex(r,j) for r,j in zip(a[:n],a[n:])])
            der = lambda *a: derivative(*[complex(r,j) for r,j in zip(a[:n],a[n:])])
            fx = func(*x)
        else:
            fx,x = fxx
            
        if not _isscalar(fx):
            return [cls._apply(lambda *y: func(*y)[i],
                               lambda *y: der(*y)[i],
                               *args,fxx=(fx[i],x)) 
                    for i in range(len(fx))]
        
        d = der(*x)
        
        if n == 1:
            d = [d]

        rda = []
        rdb = []
        jda = []
        jdb = []
        for i,p in enumerate(d):
            try:
                if len(p) == 2 and len(p[0]) == 2 and len(p[1]) == 2:
                    rda.append(p[0][0])
                    rdb.append(p[0][1])
                    jda.append(p[1][0])
                    jdb.append(p[1][1])
            except:
                if p is None:
                     p = 0
                rda.append(p.real)
                jdb.append(p.real)
                jda.append(-p.imag)
                rdb.append(p.imag)
        rd = rda + rdb
        jd = jda + jdb

        if len(rd) == 1:
            rd = rd[0]
            jd = jd[0]
            
        if not isinstance(fx,Real) and isinstance(fx,Complex):
            r = cls._element_type._apply(lambda *a: func(*a).real,None,*args,
                                         fxdx=(fx.real,rd,x))
            j = cls._element_type._apply(lambda *a: func(*a).imag,None,*args,
                                         fxdx=(fx.imag,jd,x))
            if (isinstance(r,(ummy,cls._element_type)) or 
                isinstance(j,(ummy,cls._element_type))):
                return cls(real=r,imag=j)
            return complex(r,j)
        
        return cls._element_type._apply(function,der,*args,fxdx=(fx,rd,x))
    
    @classmethod
    def _napply(cls,function,*args,fxx=None):
        n = len(args)
        if fxx is None:
            rargs = [a.real for a in args]
            jargs = [a.imag for a in args]
            args = rargs + jargs
            args = [a.convert(one).value if isinstance(a,Quantity) else a for 
                    a in args]
            x = [a.x if isinstance(a,ummy) else a for a in args]
            func = lambda *a: function(*[complex(r,j) for r,j in zip(a[:n],a[n:])])
            fx = func(*x)
        else:
            fx,x = fxx
        
        if not _isscalar(fx):
            return [cls._napply(lambda *y: func(*y)[i],*args,fxx=(fx[i],x)) 
                    for i in range(len(fx))]
            
        if not isinstance(fx,Real) and isinstance(fx,Complex):
            r = cls._element_type._napply(lambda *a: func(*a).real,*args,
                                          fxx=(fx.real,x))
            j = cls._element_type._napply(lambda *a: func(*a).imag,*args,
                                          fxx=(fx.imag,x))
            return cls(real=r,imag=j)
        
        return cls._element_type._napply(func,*args,fxx=(fx,x))
            
    @property
    def style(self):
        """
        `str` in {'polar','cartesian'}
        
        Sets whether to print the value using a cartesian or polar 
        representaion.  This property can be set at the class or instance
        level.  The default is 'cartesian'.
        """
        return self._style
    @style.setter
    def style(self,value):
        value = value.lower().strip()
        if value in ['cartesian','ri','real-imag']:
            self._style = 'cartesian'
        elif value in ['polar','r-phi']:
            self._style = 'polar'
        else:
            raise ValueError('the style must be either "cartesian" or "polar"')
    
    def tostring(self,fmt='unicode',norm=None,nsig=None,style=None):
        if self.imag == 0 and self.real == 0:
            return '0'
        
        if style is None:
            style = self.style
        
        if style == 'polar':
            r = self.r.tostring(fmt=fmt,nsig=nsig)
            
            if self.phi.x < 0:
                i = abs(self.phi).tostring(fmt=fmt,nsig=nsig)
                sign = '-'
            else:
                i = self.phi.tostring(fmt=fmt,nsig=nsig)
                sign = ''
            
            if fmt == 'html':
                if sign == '':
                    sign ='&thinsp;'
                ret = r + '&thinsp;&middot;&thinsp;<i>e</i><sup>' + sign
                ret += '<i>' + self._imag_symbol + '</i>&thinsp;'
                ret += i + '</sup>'
            elif fmt == 'latex':
                ret = r + '\\,\\cdot\\,e^{' + sign + self._imag_symbol + i + '}'
            else:
                ret = r + ' exp(' + sign + self._imag_symbol + i + ')'
            return ret
        
        if self.real == 0:
            r = ''
        else:
            r = self.real.tostring(fmt=fmt,nsig=nsig)
        if self.imag == 0:
            return r
        if self.imag.x < 0:
            i = abs(self.imag).tostring(fmt=fmt,nsig=nsig)
            if self.real == 0:
                sign = '-'
            else:
                sign = ' - '
        else:
            i = self.imag.tostring(fmt=fmt,nsig=nsig)
            if self.real == 0:
                sign = ''
            else:
                sign = ' + '
            
        if fmt == 'html':
            i = '<i>' + self._imag_symbol + '</i>&thinsp;' + i
        elif fmt == 'latex':
            i = self._imag_symbol + '\\,' + i
        else:
            i = self._imag_symbol + i

        return r + sign + i
            
    def __add__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) + v
        
        r = self.real + v.real
        i = self.imag + v.imag
        return type(self)(real=r,imag=i)
    
    def __radd__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) + v
        
        r = v.real + self.real
        i = v.imag + self.imag
        return type(self)(real=r,imag=i)
    
    def __sub__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) - v
        
        r = self.real - v.real
        i = self.imag - v.imag
        return type(self)(real=r,imag=i)
    
    def __rsub__(self,v):
        if isinstance(v,np.ndarray):
            return v - np.array(self)
        
        r = v.real - self.real
        i = v.imag - self.imag
        return type(self)(real=r,imag=i)
    
    def __mul__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)*v
        
        r = self.real*v.real - self.imag*v.imag
        i = self.imag*v.real + self.real*v.imag
        return type(self)(real=r,imag=i)
    
    def __rmul__(self,v):
        if isinstance(v,np.ndarray):
            return v*np.array(self)
        
        r = v.real*self.real - v.imag*self.imag
        i = v.imag*self.real + v.real*self.imag
        return type(self)(real=r,imag=i)
    
    def __truediv__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)/v
        
        h2 = v.real*v.real + v.imag*v.imag
        r = (self.real*v.real + self.imag*v.imag)/h2
        i = (self.imag*v.real - self.real*v.imag)/h2
        return type(self)(real=r,imag=i)
    
    def __rtruediv__(self,v):
        if isinstance(v,np.ndarray):
            return v/np.array(self)
        
        h2 = self.real*self.real + self.imag*self.imag
        r = (v.real*self.real + v.imag*self.imag)/h2
        i = (v.imag*self.real - v.real*self.imag)/h2
        return type(self)(real=r,imag=i)
        
    def __pow__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)**v
        
        if self.real == 0 and self.imag == 0:
            return self.copy(formatting=False)
        
        h2 = self.real*self.real + self.imag*self.imag
        a = np.arctan2(self.imag,self.real)
        c = h2**(v.real/2)
        t = v.real*a
        if v.imag != 0:
            t += 0.5*v.imag*np.log(h2)
            c *= np.exp(-v.imag*a)
        r = c*np.cos(t)
        i = c*np.sin(t)
        return type(self)(real=r,imag=i)
    
    def __rpow__(self,v):
        if isinstance(v,np.ndarray):
            return v**np.array(self)
        
        if v.real == 0 and v.imag == 0:
            return type(self)(real=0,imag=0)
        
        h2 = v.real*v.real + v.imag*v.imag
        a = np.arctan2(v.imag,v.real)
        c = h2**(self.real/2)
        t = self.real*a
        if self.imag != 0:
            t += 0.5*self.imag*np.log(h2)
            c *= np.exp(-self.imag*a)
        r = c*np.cos(t)
        i = c*np.sin(t)
        return type(self)(real=r,imag=i)
    
    def _nprnd(self,f):
        self._real = self.real._nprnd(f)
        self._imag = self.imag._nprnd(f)
        
    def __floordiv__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) // v
        
        h2 = v.real*v.real + v.imag*v.real
        r = (self.real*v.real + self.imag*v.imag)//h2
        i = (self.imag*v.real - self.real*v.imag)//h2
        return type(self)(real=r,imag=i)
        
    def __rfloordiv__(self,v):
        if isinstance(v,np.ndarray):
            return v // np.array(self)
        
        h2 = self.real*self.real + self.imag*self.imag
        r = (v.real*self.real + v.imag*self.imag)//h2
        i = (v.imag*self.real - v.real*self.imag)//h2
        return type(self)(real=r,imag=i)
        
    def __mod__(self,v):
        raise TypeError("can't mod immy")
    
    def __rmod__(self,v):
        raise TypeError("can't mod immy")
    
    def __abs__(self):
        return (self.real**2 + self.imag**2)**0.5
    
    def __neg__(self):
        return type(self)(real=-self.real.copy(formatting=False),
                          imag=-self.imag.copy(formatting=False))
    
    def __pos__(self):
        return type(self)(real=self.real.copy(formatting=False),
                          imag=self.imag.copy(formatting=False))
    
    #def __complex__(self):
        #return complex(self.real.x,self.imag.x)
    
    #def __float__(self):
        #raise TypeError("can't convert immy to float")
        
    #def __int__(self):
        #raise TypeError("can't convert immy to int")
        
    def __eq__(self,v):
        return self.real == v.real and self.imag == v.imag

        
class MFraction(Fraction):
    """
    A fraction.Fraction sub-class that works with mpmath.mpf objects
    """
    
    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and not (isinstance(args[0],str) 
                                   or isinstance(args[0],Fraction)):
            return args[0]
        ret = super(MFraction, cls).__new__(cls, *args, **kwargs)
        if ret.denominator == 1:
            return ret.numerator
        return ret
    
    def _mpmath_(self,p,r):
        return rational.mpq(self.numerator,self.denominator)
    
    def __add__(self,v):
        return MFraction(super().__add__(v))
    
    def __radd__(self,v):
        return MFraction(super().__radd__(v))
    
    def __sub__(self,v):
        return MFraction(super().__sub__(v))
    
    def __rsub__(self,v):
        return MFraction(super().__rsub__(v))
    
    def __mul__(self,v):
        return MFraction(super().__mul__(v))
    
    def __rmul__(self,v):
        return MFraction(super().__rmul__(v))
    
    def __truediv__(self,v):
        return MFraction(super().__truediv__(v))
    
    def __rtruediv__(self,v):
        return MFraction(super().__rtruediv__(v))
        
    def __pow__(self,v):
        return MFraction(super().__pow__(v))
    
    def __rpow__(self,v):
        if isinstance(v,Fraction):
            return MFraction(v).__pow__(self)
        return MFraction(super().__rpow__(v))
        
    def __floordiv__(self,v):
        return MFraction(super().__floordiv__(v))
        
    def __rfloordiv__(self,v):
        return MFraction(super().__rfloordiv__(v))
        
    def __mod__(self,v):
        return MFraction(super().__mod__(v))
    
    def __rmod__(self,v):
        return MFraction(super().__rmod__(v))
    
    def __abs__(self):
        return MFraction(super().__abs__())
    
    def __neg__(self):
        return MFraction(super().__neg__())
    
    def __pos__(self):
        return MFraction(super().__pos__())
