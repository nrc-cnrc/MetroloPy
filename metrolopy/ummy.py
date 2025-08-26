# -*- coding: utf-8 -*-

# module ummy

# Copyright (C) 2025 National Research Council Canada
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
from collections import namedtuple
from warnings import warn
from math import isnan,isinf,log,pi,sqrt,copysign,floor
from fractions import Fraction
from numbers import Rational,Integral,Real,Complex,Number
from .printing import PrettyPrinter,MetaPrettyPrinter
from .dfunc import Dfunc
from .exceptions import UncertiantyPrecisionWarning
from .unit import Unit,one,Quantity,MFraction
from decimal import Decimal,localcontext,InvalidOperation
from abc import ABCMeta

try:
    from mpmath import mp,mpf
except:
    mp = mpf = None


_FInfo = namedtuple('_FIinfo','rel_u precision warn_u ')
_iinfo = _FInfo(rel_u=0,precision=0,warn_u=0)
_finfodict = {}
def _getfinfo(x):
    if isinstance(x,Rational) or isinstance(x,Decimal):
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
    
    if b == 0:
        return a
    if a == 0:
        return b
    
    if isnan(a) or isnan(b):
        try:
            return type(a)('nan')
        except:
            return float('nan')
    
    if abs(b) > abs(a):
        a,b = b,a
        
    r = b/a
    x = 1 + r**2 + 2*c*r
    
    if x < 0:
        if x < -1e6:
            # something horrible happened, a bug must have generated a value of
            # for c that allows this
            raise FloatingPointError('u is sqrt of ' + str(a**2*x))
        # maybe possibly get here due to a floating point rounding error
        return 0
    
    if hasattr(x,'sqrt'):
        return abs(a)*x.sqrt()
    return abs(a)*x**0.5   

def _isscalar(x):
    if isinstance(x,str):
        return True
    try:
        len(x)
        return False
    except:
        return True

def sign(x):
    return copysign(1,x)

def _format_exp(fmt,xp):
    if xp == 0:
        return ''
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

def _to_decimal(x,max_digits=None):
    # called from the ummy and gummy tostring methods to convert x and u values.
    # The Decimal data type has some convinent methods for rounding numbers to
    # show a given number of significant figures.
    
    if mpf is not None and isinstance(x,mpf):
         xd = Decimal(mp.nstr(x,n=mp.dps))
    elif isinstance(x,Rational):
        if max_digits is None:
            max_digits = ummy.max_digits
        with localcontext(prec=max_digits):
            try:
                xd = Decimal(x.numerator)/Decimal(x.denominator)
            except:
                try:
                    xd = Decimal(str(x.numerator))/Decimal(str(x.denominator))
                except:
                    xd = Decimal(float(x.numerator))/Decimal(float(x.denominator))
    else:
        try:
            xd = Decimal(x)
        except:
            try:
                xd = Decimal(str(x))
            except:
                xd = Decimal(float(x))
    return xd

def _decimal_str(x,fmt='unicode',rtrim=False,dplace=0,dalign=None,th_spaces=True):
    # string formats a Decimal number, putting digits in groups of three
    # separated by spaces if necessary
    
    if x.is_nan():
        return 'nan'
    if x.is_infinite() and x > 0:
        if fmt == 'html':
            return '&infin;'
        if fmt == 'latex':
            return r'\infty'
        if fmt == 'ascii':
            return 'inf'
        return '\u221E'
    if x.is_infinite() and x < 0:
        if fmt == 'html':
            return '-&infin;'
        if fmt == 'latex':
            return r'-\infty'
        if fmt == 'ascii':
            return '-inf'
        return '-\u221E'
    
    s,d,e = x.as_tuple()
    
    if dalign is not None:
        dplace = len(d) + e - 1 - dalign
        
    if rtrim:
        if dplace is not None and dplace > 0:
            dp = dplace
        else:
            dp = 0
            
        i = len(d)
        while i > dp + 1 and d[i-1] == 0:
            i -= 1
        if i != len(d):
            d = d[:i]
            
    if dplace is None:
        th_spaces = False
    else:
        th_spaces &= (dplace > 3 or len(d) - dplace > 5)
        if dplace < 0:
            d = (0,)*(-dplace) + d
            dplace = 0
        elif x != 0 and dplace > len(d) - 1:
            d = d + (0,)*(dplace - len(d) + 1)
        
    if th_spaces:
        if fmt == 'html':
            sp = '&thinsp;'
        elif fmt == 'latex':
            sp = r'\,'
        else:
            sp = ' '
        
    if s and x != 0:
        txt = '-'
    else:
        txt = ''
    
    for i in range(len(d)):
        if i > 0 and dplace is not None:
            if i == dplace + 1:
                txt += '.'
            elif th_spaces and (i - dplace - 1) % 3 == 0:
                txt += sp
        txt += chr(48 + d[i])
    
    return txt
                
def _xtype(x):
    # this is used when type conversions are needed, e.g. to convert float 
    # correlations or uncertainties to Decimal when dealing with ummy's
    # with Decimal x values.
    
    if isinstance(x,np.ndarray):
        tp = x.dtype.type
    else:
        tp = type(x)
    if not issubclass(tp,Number):
        return float
    if issubclass(tp,Integral):
        return MFraction.fromnum
    return tp

def _combrd(xl,du):
    ret = _udict(xl[0]._ref)
    for r in ret:
        #ret[r] *= sign(xl[0]._u)*du[0]
        ret[r] *= du[0]

    for i in range(len(xl) - 1):
        c = sign(xl[i+1]._u)*du[i+1]
        for k,v in xl[i+1]._ref.items():
            f = v*c
            if k in ret:
                ret[k] += f
            else:
                ret[k] = f
    
    n = sum([i**2 for i in ret.values()])
    if abs(n-1) > ummy.correlation_tolerance:
        n = sqrt(n)
        for k,v in ret.items():
            ret[k] /= n
        
    for k,v in list(ret.items()):
        x = ret[k]
        if x < -1 + ummy.correlation_tolerance:
            return _udict(((k,-1),))
        elif x > 1 - ummy.correlation_tolerance:
            return _udict(((k,1),))
        elif abs(ret[k]) < ummy.correlation_tolerance:
            del ret[k]
                
    return ret


# Correlations between ummy's are stored in a _udict instance.  
# The keys in the _udict are _UmmyRef instances that reprensent independent 
# variables and store the dof and utype of the variable.  The values are the 
# correlation of the ummy with those independent variables (to within a sign) 
# and are also equal to the deriviative of the ummy with respect to each of 
# those variables divided by the combined standard uncertainty of the ummy.
class _udict(dict):
    pass

class _UmmyRef:
    __slots__ = 'dof','utype'
    def __init__(self,dof=None,utype=None,dist=None):
        self.dof = dof
        self.utype = utype


class MetaUmmy(MetaPrettyPrinter,ABCMeta):
    pass

class ummy(Dfunc,PrettyPrinter,Number,metaclass=MetaUmmy):
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
    
    max_digits = 20
    
    # correlations with absolute values smaller that this are rounded to
    # zero and correlations samaller than -1 + correlation_tolerance or
    # bigger than 1 - correlation_tolerance are rounded to -1 or 1 respectively.
    correlation_tolerance = 1e-14
    
    def __init__(self,x,u=0,dof=float('inf'),utype=None):
        if isinstance(x,Quantity):
            x = x.convert(one).value
            
        if isinstance(x,ummy):
             self._copy(x,self,formatting=False)
             return
         
        # MFractions convert implicitly to Decimal or mpf
        if isinstance(x,Fraction):
            x = MFraction(x) 
        if isinstance(u,Fraction):
            u = MFraction(u)
            
        # sometimes numpy.array of shape () appear and cause problems
        if isinstance(x,np.ndarray):
            x = x.item()
        if isinstance(u,np.ndarray):
            u = u.item()
            
        try:
            float(x)
            hash(x)
        except TypeError:
            raise TypeError('x must be a real number') from None
        except ValueError:
            raise ValueError('x must be a real number') from None
                
        if self.rounding_u and not isinstance(x,Rational):
            finfo,x = _getfinfo(x)
            if finfo is _iinfo:
                warn('numpy.finfo cannot get the floating point accuracy for x\nno uncertainty will be included to account for floating point rounding errors',UncertiantyPrecisionWarning)
            else:
                u = _combu(u,float(abs(x)*finfo.rel_u),0)
            
            self._finfo = finfo
        else:
            self._finfo = None

            
        self._x = x
        self._u = u
        
        if not isinf(u) and not isnan(u) and u != 0:
            if isinstance(dof,_udict):
                # If the ummy has been created from a mathematical operation
                # between other ummy's then correlations with the independent 
                # variables are stored in the _udict.  The keys in the _udict
                # are _UmmyRef instances that represent independent 
                # variables and store the dof and utype of the variable.  The
                # values are the correlation of the ummy with those independent
                # variables (to within a sign) and are also equal to the
                # deriviative of the ummy with respect to each of those variables
                # divided by self._u.  The actual value of the
                # dof property will be computed with the W-S apporoximation.
                # u may be negative indicating that it is sharing a _ref 
                # _udict with another ummy and has correlation = -1 with that 
                # ummy: e.g. if the ummy is the result of the operation 1 - x
                # we don't bother to create a new _udict for the ummy, we just
                # use x._ref with self._u having the opposite sign from x._u.
                self._ref = dof
                
                # _dof and _utype are set to None for now and will be computed
                # if and when the dof and utype properties are called for the
                # first time.
                self._dof = None
                self._utype = None
            else:
                # The ummy represents an independent variable, that is u > 0
                # and it is not (at this moment) correlated with any other
                # ummy's.  
                    
                try:
                    float(u)
                    hash(u)
                    if not isnan(u) and u < 0:
                        raise ValueError()
                except TypeError:
                    raise TypeError('u must be a real number >= 0') from None
                except ValueError:
                    raise ValueError('u must be a real number >= 0') from None
                    
                try:
                    if dof > self.max_dof:
                        dof = float('inf')
                    else:
                        if isinstance(dof,Integral):
                            dof = int(dof)
                        else:
                            dof = float(dof)
                        if dof <= 0 or isnan(dof):
                            raise ValueError()
                except ValueError:
                    raise ValueError('dof must be a real number > 0') from None
                except TypeError:
                    raise TypeError('dof must be a real number > 0') from None
                if utype is not None and not isinstance(utype,str):
                    raise TypeError('utype must be str or None')
                
                self._dof = dof
                if utype is None:
                    self._utype = [None]
                else:
                    self._utype = [utype]
                    
                # create an _UmmyRef to represent the actual independent 
                # variable and create a _udict with one entry defining 
                # the correlation with that _UmmyRef to be 1
                self._ref = _udict(((_UmmyRef(dof,utype),1),))
        else:
            # if u == 0, then the ummy will have no correlation with any
            # independent variables, no utype and the dof property will return
            # float('inf')
            
            self._ref = None
            self._dof = None
            self._utype = [None]
        
    @property
    def x(self):
        return self._x
        
    @property
    def u(self):
        return abs(self._u)
        
    @property
    def dof(self):
        """
        `float`, read-only

        Returns the number or degrees of freedom that the uncertainty of the 
        ummy is based on.  If the ummy was created as the result of an 
        operation between two or more other gummys, then the dof is the effective
        number of degrees of freedom calculated using the Welch-Satterthwaite 
        approximation.
        """
        if self._dof is not None:
            return self._dof
        
        if self._ref is None:
            self._dof = float('inf')
        else:
            # if _dof has not been set yet, calculate it using the 
            # Welch-Satterthwaite approximation and store the result in 
            # self._dof.  self._ref is a dictionary and the keys are _UmmyRef
            # instances representing the indedpendent variables and store the
            # dof for those variables.  The values are the correlations of 
            # self with those independent variables and these correlations
            # are equal to the derivative of self with respect to each 
            # independent variable divided by self.u.
            
            rdof = 0
            for k,v in self._ref.items():
                if k.dof is not None:
                    rdof += v**4/k.dof
            if rdof > 0:
                rdof = 1/rdof
                if rdof > self.max_dof:
                    rdof = float('inf')
                elif rdof < 1:
                    rdof = 1
                self._dof = rdof
            else:
                self._dof = float('inf')

        return self._dof
    
    @property
    def utype(self):
        """
        `str`, `None` or a list containing strings and possibly `None`
    
        An arbitrary string value labeling the uncertainty type or or a 
        list of types if the gummy was constructed from independent
        variables with different utypes.
        """
        if self._utype is None:
            # if _utype has not been set yet, go through all the independent 
            # variables in self._ref and put all the utypes found in a set.
            self._utype = []
            for k in self._ref:
                if k.utype not in self._utype:
                    self._utype.append(k.utype)
                
        if len(self._utype) == 1:
            return self._utype[0]
        return self._utype

    @property
    def isindependent(self):
        """
        `bool`, read-only

        Returns `True` if the ummy is an independent variable.  That is the 
        ummy has u > 0 and was not correlated with any other ummy's when it was 
        created or is perfectly correlated or anti-correlated (correlation 
        coefficeint 1 or -1) with such an ummy.'
        """
        return (self._ref is not None) and len(self._ref) == 1
    
    @property
    def finfo(self):
        if self._finfo is None:
            self._finfo = _getfinfo(self._x)[0]
        return self._finfo
        
    @staticmethod
    def _copy(s,r,formatting=True,totype=None):
        # copies attributes of s to r
        if totype is not None:
            r._x = totype(s._x)
            r._u = totype(s._u)
        else:
            r._x = s._x
            r._u = s._u
        r._ref = s._ref
        r._dof = s._dof
        r._utype = s._utype
        r._finfo = s._finfo

        if formatting:
            if s.nsig != type(r).nsig:
                r.nsig = s.nsig
            if s.thousand_spaces != type(r).thousand_spaces:
                r.thousand_spaces = s.thousand_spaces
    
    def copy(self,formatting=True,totype=None):
        """
        Returns a copy of the ummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new gummy.  The default for `formatting` is `True`.  If totype
        is defined the x and u properties will be converted to type `totype`
        before copying.
        """
        r = type(self).__new__(type(self))
        self._copy(self,r,formatting=formatting,totype=totype)
        return r
        
    def tofloat(self):
        """
        Returns a copy of the gummy with x an u converted to `float`.
        """
        return self.copy(formatting=False,totype=float)
    
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
        if not isinstance(g,ummy) or self._ref is None or g._ref is None:
            return 0
        
        c = sign(self._u)*sign(g._u)*sum(g._ref[k]*v for k,v in self._ref.items() if k in g._ref)
                
        if abs(c) < ummy.correlation_tolerance:
            return 0
        if c < -1 + ummy.correlation_tolerance:
            return -1
        if c > 1 - ummy.correlation_tolerance:
            return 1
            
        return float(c)
        
    def covariance(self,g):
        """
        Returns the covariance between `self` and `g`.
        """
        if isinstance(g,ummy):
            return self.u*g.u*self.correlation(g)
        return 0
        
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
        
    @classmethod
    def create(cls,x,u=None,dof=None,utype=None,correlation_matrix=None,
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
       
        n = len(x)
        if covariance_matrix is not None:
            if correlation_matrix is not None:
                raise TypeError('correlation_matrix and covariance_matrix cannot both be specified')
            try:
                u = [sqrt(covariance_matrix[i][i]) for i in range(n)]
            except ValueError:
                raise ValueError('the diagonal elements of the covariance matrix must be positive real numbers')
            except IndexError:
                raise ValueError('covariance must have shape len(gummys) x len(gummys)')
            correlation_matrix = [[covariance_matrix[i][j]/(u[i]*u[j]) for i in range(n)] for j in range(n)]
        else:
            if u is None:
                u = [0]*n
            elif _isscalar(u):
                u = [u]*n
            
        if correlation_matrix is not None:
            if not _isscalar(dof):
                if any(i != dof[0] for i in dof):
                    raise TypeError('dof cannot be set individually of a correlation of covariance matrix is specified')
                dof = dof[0]
            if not _isscalar(utype):
                if any(i != utype[0] for i in utype):
                    raise TypeError('utype cannot be set individually of a correlation of covariance matrix is specified')
                utype = utype[0]
                
            m = np.asarray(correlation_matrix)
            
            if not np.allclose(m.T,m,atol = sqrt(ummy.correlation_tolerance)):
                raise ValueError('the correlation matrix is not symmetric')
                
            if m.shape != (n,n):
                raise ValueError('the correlation matrix must have shape len(gummys) x len(gummys)')
            
            for i in range(n):
                if abs(m[i][i] - 1.0) > sqrt(ummy.correlation_tolerance):
                    raise ValueError('correlation matrix diagonal elements must be 1')
                    
            # The correlated values can be constructed from a linear combination
            # of indepedent variables with coefficients equal to the matrix
            # square root of the correlation matrix.
            val,vec = np.linalg.eig(m)
            
            # Round eigenvalues very close to zero to zero, allowing for small
            # negative values that may appear due to the finite numerical 
            # precision.
            for i in range(n):
                if val[i] < -ummy.correlation_tolerance:
                    raise ValueError('the correlation matrix is not positive semi-definate')
                if val[i] < ummy.correlation_tolerance:
                    val[i] = 0
            
            # matrix sqaure root of m        
            sqrtm = vec*np.sqrt(val)@np.linalg.inv(vec) 
            
            ind = [_UmmyRef(dof,utype) for _ in range(n)]
            idof = [_udict(((ind[j],float(sqrtm[i][j])) for j in range(n))) 
                    for i in range(n)]
                        
            ret = [cls(x[i],u=u[i],dof=idof[i]) for i in range(n)]
            
        else:
            if dof is None:
                dof = [float('inf')]*n
            elif _isscalar(dof):
                dof = [dof]*n
                
            if utype is None or isinstance(utype,str):
                utype = [utype]*n
        
            ret = [cls(x[i],u=u[i],dof=dof[i],utype=utype[i]) for i in range(n)]
            
        return ret
        
    def tostring(self,fmt='unicode',nsig=None,**kwds):
        if nsig is None:
            nsig = self.nsig
        if nsig > self.max_digits - 2:
            nsig = self.max_digits - 2
        if nsig < 1:
            nsig = 1
            
        try:
            uexp = 0
            xexp = 0
            dp = 0
            xd = _to_decimal(self._x,max_digits=self.max_digits)
            ud = _to_decimal(self.u)
            with localcontext(prec=self.max_digits):
                if ud == 0 or not ud.is_finite():
                    if xd.is_finite():
                        nm = False
                        nd = self.max_digits
                        if self.finfo is not None and self.finfo.precision > 0 and self.finfo.precision <= nd:
                            nd = self.finfo.precision
                            nm = True
                        xd = round(xd,nd-xd.adjusted()-1)
                        if xd == 0:
                            xexp = 0
                        else:
                            xexp = xd.adjusted()
                            
                        if self._x == xd or nm:
                            xd = xd.normalize()
                else:
                    ud = round(ud,nsig-ud.adjusted()-1)
                    ud = round(ud,nsig-ud.adjusted()-1) # round again in case 9 round up to 10 in the last line
                    uexp = ud.adjusted()
        
                    if xd.is_finite():
                        try:
                            xd = xd.quantize(ud)
                        except InvalidOperation:
                            xd = round(xd,self.max_digits-xd.adjusted()-1)
                            xexp = xd.adjusted()
                                
                        if abs(xd) < abs(ud):
                            dp = xd.adjusted() - uexp
                            xexp = uexp
                        else:
                            xexp = xd.adjusted()
                    else:
                        xexp = uexp
                        
            scin = False
            if self.sci_notation is None:
                if xexp > self.sci_notation_high or xexp < self.sci_notation_low:
                    scin = True
            else:
                scin = self.sci_notation
            if xexp > self.max_digits - 1:
                scin = True
            elif xexp < 0 and len(xd.as_tuple().digits) + self.sci_notation_low > self.max_digits - 1:
                scin = True
            elif ud.is_finite() and ud != 0 and xd.is_finite() and (xd.as_tuple()[2] > 0 or len(xd.as_tuple()[1]) < nsig):
                scin = True

            if not scin:
                xexp = 0
                
            if xexp == 0:
                txt =_decimal_str(xd,dalign=0,fmt=fmt)
            else:
                txt = _decimal_str(xd,dplace=dp,fmt=fmt)
                
            if xd.is_finite():
                if ud != 0:
                    txt += '(' + _decimal_str(ud,dplace=None,fmt=fmt) + ')'
                txt += _format_exp(fmt,xexp)
            else:
                if ud != 0:
                    txt += '(' + _decimal_str(ud,dplace=None,fmt=fmt)
                    txt += _format_exp(fmt,xexp) + ')'
            return txt
       
        except:
            try:
                return(str(self.x) + '{' + str(self.u) + '}' + '??')
            except:
                try:
                    return(str(self.x) + '{??}')
                except:
                    return('??')
        
    #@classmethod
    #def _get_dof(cls,dof1,dof2,du1,du2,c):
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
        
        #d = du1**4/dof1 + du2**4/dof2
        #d += 2*c**2*(((du1 + du2)**4 - du1**4 - du2**4)/
                     #((1 - 2*c + 2*c**2)*(dof1+dof2)))

        #if d == 0:
            #return float('inf')
        
        #r = 1/d
        #if r > cls.max_dof:
            #return float('inf')
        
        #if r < 0:
            #r = 1
        
        #return r
    
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
            
        xtype = _xtype(fx)
            
        try:
            ad = [(a,a._u*p) for a,p in zip(args,d) 
                  if isinstance(a,ummy) and a._u != 0]
        except TypeError:
            ad = [(a,xtype(a._u)*p) for a,p in zip(args,d) 
                  if isinstance(a,ummy) and a._u != 0]
        if len(ad) == 0:
            return cls(fx)
        args,du = list(map(list, zip(*ad)))
            
        if len(args) == 1:
            r = cls(fx,u=du[0],dof=args[0]._ref)
            return r
        
        maxdu = max(abs(a) for a in du)
        if maxdu == 0:
            return cls(fx)
        try:
            dun = [d/maxdu for d in du]
        except TypeError:
            dun = [xtype(d)/xtype(maxdu) for d in du]
        
        u = sum(d**2 for d in dun)
        try:
            u += 2*sum(sum(args[i]._ref[k]*v for k,v in args[j]._ref.items() 
                           if k in args[i]._ref)*dun[i]*dun[j] 
                       for i in range(len(args)) for j in range(i+1,len(args)))
        except TypeError:
            u += 2*sum(xtype(sum(args[i]._ref[k]*v for k,v in args[j]._ref.items() 
                                 if k in args[i]._ref))*dun[i]*dun[j] 
                       for i in range(len(args)) for j in range(i+1,len(args)))
        
        if not isnan(u):
            if u < 0:
                if u < -1e6:
                    # something horrible and unexpected happened
                    raise FloatingPointError('u is sqrt of ' + str(maxdu**2*u))
                u = 0
                
            if hasattr(u,'sqrt'):
                u = u.sqrt()
            else:
                u = u**0.5
            u = maxdu*u
            
        if u == 0 or isnan(u):
            return cls(fx)
        
        try:
            du = [float(d/u) for d in du]
        except TypeError:
            du = [float(xtype(d)/xtype(u)) for d in du]
            
        r = cls(fx,u=u,dof=_combrd(args,du))
            
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
    
    def _aop(self,b,f,d1,d2):
        # computes the result of a binary operation f(self,b).  d1 and d2 are 
        # the derivatives with respect to self and b respectively.
        
        # correlation between self and b (/sign(self._u)*sign(b._u))
        if not isinstance(b,ummy) or b._ref is None or self._ref is None:
            c = 0
        else:
            c = sum(b._ref[k]*v for k,v in self._ref.items() if k in b._ref)
            
        x = f(self._x,b._x)
        try:
            dua = d1(self._x,b._x)*self._u
            dub = d2(self._x,b._x)*b._u
            u = _combu(dua,dub,c) # combined standard uncertainty
        except TypeError:
            # e.g. if self.x and b.x are Decimal we may need to convert the
            # correlation c which is float to Decimal
            xtype = _xtype(x)
            dua = d1(self._x,b._x)*xtype(self._u)
            dub = d2(self._x,b._x)*xtype(b._u)
            u = _combu(dua,dub,xtype(c))
            
        # some special cases where we don't need to bother to work out all the
        # correlations with the independent variables
        if u == 0:
            return type(b)(x)
        
        if self._u == 0:
            return type(b)(x,u=u,dof=b._ref)
        
        if b._u == 0:
            return type(b)(x,u=u,dof=self._ref)
        
        if self._ref is b._ref:
            s = sign(dua + dub)
            return type(b)(x,u=s*u,dof=self._ref)
        
        # not a special case, so combine the correlation dictionaries of self 
        # and b
        
        dua = float(dua/u)
        dub = float(dub/u)
       
        ref = _udict(((k,dua*v) for k,v in self._ref.items()))
        for k,v in b._ref.items():
            if k in ref:
                ref[k] += v*dub
            else:
                ref[k] = v*dub
            
        # check that sum squared correlations add to 1 in case rounding errors
        # have accumulated
        n = sum([i**2 for i in ref.values()])
        if abs(n-1) > ummy.correlation_tolerance:
            n = sqrt(n)
            for k,v in ref.items():
                ref[k] /= n
            
        # round correlations to 1, -1 or 0 if within correlation_tolerance of 
        # these values
        for k,v in list(ref.items()):
            c = ref[k]
            if c < -1 + ummy.correlation_tolerance:
                ref = _udict(((k,1),))
                u = -u
                break
            elif c > 1 - ummy.correlation_tolerance:
                ref = _udict(((k,1),))
                break
            elif abs(c) < ummy.correlation_tolerance:
                del ref[k]
        
        return type(b)(x,u=u,dof=ref)
                
    def __add__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) + b
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__radd__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) + b
        
        if not isinstance(b,ummy):
            return type(self)(self._x + b,u=self._u,dof=self._ref)
        
        return self._aop(b,
                         lambda x1,x2:x1 + x2,
                         lambda x1,x2:1,
                         lambda x1,x2:1)
    
    def __radd__(self,b):
        if isinstance(b,np.ndarray):
            return b + np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b + immy(self)
        
        if isinstance(b,ummy):
            return b.__add__(self)
        
        return type(self)(self._x + b,u=self._u,dof=self._ref)
    
    def __sub__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) - b
          
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rsub__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) - b
        
        if not isinstance(b,ummy):
            return type(self)(self._x - b,u=self._u,dof=self._ref)
            
        return self._aop(b,
                         lambda x1,x2:x1 - x2,
                         lambda x1,x2:1,
                         lambda x1,x2:-1)
    
    def __rsub__(self,b):
        if isinstance(b,np.ndarray):
            return b - np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b - immy(self)
        
        if isinstance(b,ummy):
            return b.__sub__(self)
        
        r = type(self)(b - self._x,u=-self._u,dof=self._ref)
        return r
    
    def __mul__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self)*b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmul__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self)*b
        
        if not isinstance(b,ummy):
            x = self._x*b
            try:
                ub = self._u*b
            except TypeError:
                if isinstance(x,Integral):
                    ub = MFraction.fromnum(self._u)*b
                else:
                    ub = type(x)(self._u)*b
            return type(self)(x,u=ub,dof=self._ref)
              
        return self._aop(b,
                         lambda x1,x2:x1*x2,
                         lambda x1,x2:x2,
                         lambda x1,x2:x1)
    
    def __rmul__(self,b):
        if isinstance(b,np.ndarray):
            return b*np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b*immy(self)
        
        if isinstance(b,ummy):
            return b.__mul__(self)
        
        x = self._x*b
        try:
            ub = self._u*b
        except TypeError:
            ub = _xtype(x)(self._u)*b
        return type(self)(x,u=ub,dof=self._ref)
    
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
            try:
                ub = self._u/b
            except TypeError:
                ub = _xtype(x)(self._u)/b
            return type(self)(x,ub,dof=self._ref)
                
        if b._x == 0:
            raise ZeroDivisionError('division by zero')
        else:
            if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral)):
                x = MFraction(self._x,b._x)
            else:
                x = self._x/b._x
            
        return self._aop(b,
                         lambda x1,x2:x1/x2,
                         lambda x1,x2:1/x2,
                         lambda x1,x2:-x1/x2**2)
    
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
        try:
            ub = -b*self._u/self._x**2
        except TypeError:
            ub = -b*_xtype(x)(self._u)/self._x**2
        return type(self)(x,ub,dof=self._ref)
    
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
            try:
                ub = b*self._x**(b-1)*self._u
            except TypeError:
                ub = b*self._x**(b-1)*_xtype(x)(self._u)
            return type(self)(x,ub,dof=self._ref)

        if self._x <= 0:
            raise ValueError('a negative or zero value cannot raised to a power which has an uncertainty')            
                
        if (isinstance(self._x,Integral) and 
                isinstance(b._x,Integral) and b._x < 0):
            f = lambda x1,x2:MFraction(1,x1**-x2)
        else:
            f = lambda x1,x2:x1**x2
        
        if hasattr(self._x,'ln'):
            lgx = self._x.ln()   
        else:
            lgx = log(self._x)
            
        return self._aop(b,
                         f,
                         lambda x1,x2:x2*x1**(x2 - 1),
                         lambda x1,x2:lgx*x1**x2)
    
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
            
        if hasattr(b,'ln'):
            lgb = b.ln()
        else:
            lgb = log(b)
            
        try:
            ub = b**self._x*lgb*self._u
        except TypeError:
            ub = b**self._x*lgb*_xtype(x)(self._u)
        return type(self)(x,ub,dof=self._ref)
    
    def _nprnd(self,f):
        return type(self)(f(self._x))
        
    def __floordiv__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) // b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmfloordiv__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) // b
        
        return self.__truediv__(b)._nprnd(floor)
        
    def __rfloordiv__(self,b):
        if isinstance(b,np.ndarray):
            return b // np.array(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return b // immy(self)
        
        if isinstance(b,ummy):
            return b.__floordiv__(self)
        
        return self.__rtruediv__(b)._nprnd(floor)
        
    def __mod__(self,b):
        if isinstance(b,np.ndarray):
            return np.array(self) % b
        
        if isinstance(b,(immy,Unit,Quantity)):
            return b.__rmod__(self)
        
        if isinstance(b,Complex) and not isinstance(b,Real):
            return immy(self) % b
        
        ret = ummy._apply(lambda x1,x2: x1%x2,
                          lambda x1,x2: (1,copysign(abs(x1//x2),x2)),self,b)
        return type(self)(ret)
        
    def __rmod__(self,b):
        if isinstance(b,np.ndarray):
            return b % np.array(self)
        
        if isinstance(b,ummy):
            return b.__mod__(self)
        
        ret = ummy._apply(lambda x1,x2: x1%x2,
                          lambda x1,x2: (1,copysign(abs(x1//x2),x2)),b,self)
        return type(self)(ret)
    
    def __divmod__(self,b):
        return (self//b,self%b)
    
    def __rdivmod__(self,b):
        return (b//self,b%self)
        
    def __neg__(self):
        r = self.copy(formatting=False)
        r._x = -self._x
        if self._ref is not None:
            r._u = -self._u
        return r
        
    def __pos__(self):
        return self.copy(formatting=False)
        
    def __abs__(self):
        r = self.copy(formatting=False)
        if self._x < 0:
            r._x = -self._x
            if self._ref is not None:
                r._u = -self._u
        return r
    
    def __eq__(self,v):
        if self is v:
            return True
        
        if isinstance(v,ummy):
            if self.u == 0 and v.u == 0:
                return self._x == v._x
            if self.correlation(v) == 1:
                return self.x == v.x and self.u == v.u
            return False
        
        if self.u == 0:
            return self.x == v
        
        return False
    
    def __ne__(self, v):
        return not self == v
        
    def __lt__(self, v):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be ordered')
        return self.x < v
    
    def __le__(self, v):
       if self.u != 0:
           raise TypeError('an ummy with non-zero uncertainty cannot be ordered')
       return self.x <= v
        
    def __gt__(self, v):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be ordered')
        return self.x > v
        
    def __ge__(self, v):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be ordered')
        return self.x >= v
        
    def __float__(self):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be converted to float')
        return float(self.x)
    
    def __int__(self):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be converted to int')
        return int(self.x)
        
    def __complex__(self):
        if self.u != 0:
            raise TypeError('an ummy with non-zero uncertainty cannot be converted to complex')
        return complex(self.x)
    
    def __bool__(self):
        return self != 0
    
    #def __trunc__(self):
        #try:
            #return self.x.__trunc__()
        #except:
            #return float(self.x).__trunc__()
        
    #def __floor__(self):
        #try:
            #return self.x.__floor__()
        #except:
            #return float(self.x).__floor__()
        
    #def __ceil__(self):
        #try:
            #return self.x.__ceil__()
        #except:
            #return float(self.x).__ceil__()
        
    #def __round__(self,ndigits=None):
        #try:
            #ret = round(self.x,ndigits)
        #except:
            #ret = round(float(self.x),ndigits)
        
        #if ndigits is None:
            #return ret
        #else:
            #return type(self)(ret)
        
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
        
    def __hash__(self):
        if self._ref is None:
            return hash(self.x)
        s = sign(self._u)
        return hash((tuple((k,s*v) for k,v in self._ref.items()),self.x,self.u))
    
    def _ufromrefs(self,x):
        if self._ref is None:
            return dict()
        
        if isinstance(x,ummy) or isinstance(x,str) or isinstance(x,Quantity):
            x = [x]
        
        x = [i.value if isinstance(i,Quantity) else i for i in x]
            
        d = dict()
        for g in x:
            if isinstance(g,str):
                for k in self._ref:
                    if k.utype == g:
                        if k in d:
                            d[k] += 1
                        else:
                            d[k] = 1
            elif isinstance(g,ummy):
                for k,v in g._ref.items():
                    if k in d:
                        d[k] += v**2
                    else:
                        d[k] = v**2
                    
        for k,v in list(d.items()):
            if k in self._ref and v > 0:
                if v > 1:
                    d[k] = self._ref[k]**2
                else:
                    d[k] *= self._ref[k]**2
            else:
                del d[k]
                    
        return d
            
        
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
       
        #x = ummy._toummylist(x)
        #x = [i for i in x if self.correlation(i) != 0]
        #if len(x) == 0:
            #return 0

        #v = ummy.correlation_matrix(x)
            
        #b = [self.correlation(z) for z in x]
        #s = np.linalg.lstsq(v,b,rcond=None)[0]
        #u = 0
        
        #d = [i*self.u/j.u for i,j in zip(s,x)]
        #for i in range(len(x)):
            #for j in range(len(x)):
                #u += d[i]*d[j]*x[i].correlation(x[j])*x[i].u*x[j].u
                
        #return u**0.5
        
        return float(self.u*np.sqrt(sum(v for v in self._ufromrefs(x).values())))
        
    def doffrom(self,x):
        """
        Gets the degrees of freedom contributed from particular ummys or
        utypes if all other free variables are held fixed.
        
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
        #x = ummy._toummylist(x)
        #x = [i for i in x if self.correlation(i) != 0]
        #if len(x) == 0:
            #return float('inf')
            
        #v = ummy.correlation_matrix(x)
        #b = [self.correlation(z) for z in x]
        #s = np.linalg.lstsq(v,b,rcond=None)[0]
        #d = [i*self.u/j.u for i,j in zip(s,x)]
        #usq = 0
        #dm = 0
        #for i in range(len(x)):
            #for j in range(len(x)):
                #usqi = d[i]*d[j]*x[i].correlation(x[j])*x[i].u*x[j].u
                #usq += usqi
                #dm += usqi**2/x[i].dof
                
        #if dm == 0:
           # return float('inf')
        #dof = usq**2/dm
        #if dof > ummy.max_dof:
            #return float('inf')
        #if dof < 1:
            #return 1
        #return dof
        
        dof = sum(v**2/k.dof for k,v in self._ufromrefs(x).items() if k.dof is not None)
        if dof > 0:
            dof = sum(v for v in self._ufromrefs(x).values())**2/dof
        else:
            return float('inf')
        if dof > ummy.max_dof:
            return float('inf')
        if dof < 1:
            return 1
        return float(dof)
    
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
        if isinstance(p,ummy) and p._u != 0:          
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


class MetaImmy(MetaPrettyPrinter,ABCMeta):
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
    
class immy(PrettyPrinter,Dfunc,Number,metaclass=MetaImmy):
    
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
        Returns the magnitude of the value (`abs(self)`).
        """
        if self._r is None:
            self._r = (self.real**2 + self.imag**2)**0.5
        return self._r
    
    @property
    def phi(self):
        """
        read-only
        Returns the polar angle of the value (`self.angle()`)
        """
        if self._phi is None:
            self._phi = np.arctan2(self.imag,self.real)
        return self._phi
    
    def conjugate(self):
        """
        Returns the complex conjugate.
        """
        return type(self)(real=self.real,imag=-self.imag)
    
    def angle(self):
        """
        Returns the polar angle in radians (which is also the phi read-only 
        property value).
        """
        return self.phi
    
    def copy(self,formatting=True,totype=None):
        """
        Returns a copy of the jummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new jummy.  The default for `formatting` is `True`.  If the
        tofloats parameter is True x and u for both the real and
        imaginary components will be converted to floats.
        """
        if self._ridef:
            r = self.real.copy(formatting=formatting,totype=totype)
            i = self.imag.copy(formatting=formatting,totype=totype)
            return type(self)(real=r,imag=i)
        
        else:
            r = self.r.copy(formatting=formatting,totype=totype)
            phi = self.phi.copy(formatting=formatting,totype=totype)
            return type(self)(r=r,phi=phi)
    
    def tofloat(self):
        """
        Returns a copy of the gummy with x an u (for both the real and
        imaginary components) converted to floats.
        """
        return self.copy(formatting=False,totype=float)
    
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
        return self.r
    
    def __neg__(self):
        return type(self)(real=-self.real.copy(formatting=False),
                          imag=-self.imag.copy(formatting=False))
    
    def __pos__(self):
        return type(self)(real=self.real.copy(formatting=False),
                          imag=self.imag.copy(formatting=False))
    
    def __eq__(self,v):
        return self.real == v.real and self.imag == v.imag
    
    def __hash__(self):
        return hash((self.real,self.imag))
    
    def __complex__(self):
        return complex(float(self.real),float(self.imag))
    
    def __bool__(self):
        return self != 0
    
    #def __float__(self):
        #raise TypeError("can't convert immy to float")
        
    #def __int__(self):
        #raise TypeError("can't convert immy to int")