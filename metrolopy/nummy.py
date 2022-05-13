# -*- coding: utf-8 -*-

# module nummy

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
The nummy object defined here integrates the Monte-Carlo uncertainty 
propagation code in the distributions module with the ummy object.  The nummy 
class is not intended to be used directly; rather it is utilized by the gummy
class.
"""
import numpy as np
from .ummy import ummy
from .distributions import (Distribution,TDist,NormalDist,MultivariateElement,
                            MultivariateDistribution,MultiNormalDist,MultiTDist)
from .exceptions import NoSimulatedDataError
from math import isinf,isfinite,isnan,sqrt

def _bop(f,npf,s,b):    
    if isinstance(b,nummy):
        f._dist = Distribution.apply(npf,s._dist,b._dist)
    elif isinstance(b,ummy):
        f._dist = Distribution.apply(npf,s._dist,nummy(b)._dist)
    else:
        f._dist = Distribution.apply(npf,s._dist,b)
    return f

def _rbop(f,npf,s,b):
    f._dist = Distribution.apply(npf,b,s._dist)
    return f

def _uop(f,npf,s):        
    f._dist = Distribution.apply(npf,s._dist)
    return f

    
class nummy(ummy):
    #  This class is not intended to be used directly and was created to contain
    #  the code that the gummy object uses for Monte-Carlo uncertainty propagation.

    _cimethod = 'shortest'
    _bayesian = False # see the gummy bayesian property
    _fp = None
    _nsim = None
    
    def __init__(self,x,u=0,dof=float('inf'),utype=None,name=None):
        self._bayesian = nummy._bayesian
        
        if isinstance(x,ummy):
            self._copy(x,self,formatting=False)
            return
        
        self.name = name
        if isinstance(x,Distribution):
            if x._used or not x.isindependent:
                raise ValueError('Distribution instances may only be used as the x parameter of a new gummy only if they\nrepresent independant variables and have not previously been used with another gummy')                
            x._used = True
            if isinstance(x,MultivariateElement):
                raise TypeError('a MultivariateElement may not be used in a gummy initializer\nuse the create static method with a MultivariateDistribution')
            
            self._dist = x
            
            if hasattr(x,'dof'):
                if nummy._bayesian:
                    u = float(x.u())*np.sqrt(x.dof/(x.dof-2))
                    dof = float('inf')
                else:
                    u = float(x.u())
                    dof = x.dof
            else:
                dof = float('inf')
                u = float(x.u())
                
            super().__init__(x.x(),u=u,dof=dof,utype=utype)
            
        else:
            super().__init__(x,u=u,dof=dof,utype=utype)
            x = float(x)
            u = float(self._u)
            if self._u == 0 or isnan(x) or isinf(x) or isnan(u):
                self._dist = x
                return
            dof = self._dof
            if isinf(dof):
                self._dist = NormalDist(x,u)
            else:
                if nummy._bayesian:
                    self._dist = TDist(x,u*np.sqrt((dof-2)/dof),dof)
                    self._dof = float('inf')
                else:
                    self._dist = TDist(x,u,dof)
            
    @property
    def distribution(self):
        """
        read-only

        Returns ths `Distribution` instance associated with the gummy.
        """
        return self._dist
    
    @property
    def dof(self):
        """
        float, read-only

        Returns the number or degrees of freedom that the uncertainty of the 
        gummy is based on.  If `gummy.bayesian` is set to `False`, then the Welch-
        Satterthwaite approximation is used to calculate the effective number
        of degrees of freedom for gummys that result from an operation between
        two or more other gummys.  A version of the Welch-Satterthwaite 
        approximation that takes into account correlations is used here, see
        [R. Willink, Metrologia, 44, 340 (2007)].  If `gummy.bayesian` is `True`
        then gummys that are the result from an opertaion between other gummys
        will always have dof = float('inf').
        """
        if isinf(self._dof) and hasattr(self._dist,'dof'):
            dof = self._dist.dof
        else:
            dof = self._dof
        if dof < 1:
            # Occasionally correlations can result in a dof less than 1;
            # see the ummy._get_dof method.
            return 1
        return dof
    
    @property
    def name(self):
        if self._name is None:
            return None
        if isinstance(self._name,str):
            return self._name
        return self._name[0]
    @name.setter
    def name(self,v):
        if v is None:
            self._name = None
            return
        elif isinstance(v,str):
            self._name = v.strip()
            return
        
        try:
            if len(v) != 4:
                raise ValueError('the name must be a string or a length 4 tuple or str')
        except TypeError:
            raise ValueError('the name must be a string or a length 4 tuple of str')
            
        try:
            n = v[0].strip()
            self._name = tuple([n if e is None else e.strip() for e in v])
        except AttributeError:
            raise ValueError('the name must be a string or a length 4 tuple of str')
            
    def get_name(self,fmt='unicode',norm=None):
        fmt = fmt.strip().lower()
        if fmt not in {'unicode','html','latex','ascii'}:
            raise ValueError('fmt "' + str(fmt) + '" is not recognized')
            
        if self._name is None:
            return None
        
        if isinstance(self._name,str):
            name = str(self._name).strip()
            if fmt == 'html' and len(name) == 1:
                return '<i>' + name + '</i>'
            if fmt == 'latex' and len(name) > 1:
                if norm is None:
                    norm = type(self).latex_norm
                return norm(self.name)
            return self._name
        
        if fmt == 'unicode':
            return self._name[0]
        if fmt == 'html':
            return self._name[1]
        if fmt == 'latex':
            return self._name[2]
        if fmt == 'ascii':
            return self._name[3]
    
    @property
    def bayesian(self):
        """
        `bool`

        Read/write at the class level, but read-only at the instance level.
        The default value is `False`; this should only be changed once at the
        beginning of the session.  This property affects how the level of
        confidence `p` (sometimes called coverage probability) of an expanded
        uncertainty is related to the coverage factor `k` for a gummy based on
        data with finite degrees of freedom.

        Standard uncertainties are often based on the standard deviation of a set
        of measurements (and the assumption that these measurements are drawn
        from a normally distributed population).  Traditionally (e.g. the GUM
        2008 edition) the standard uncertainty is taken to be the standard
        deviation of the mean (s/sqrt(n), where s is the sample standard deviation
        and n is the number of measurements).  However there is some "extra
        uncertainty" because the sample standard devation not exactly equal to
        the population standard deviation.  This is taken into account by using
        a Student's t distribution to calculate the expanded uncertainty.  However
        it has been pointed out, by those who advocate a Bayesian point of view,
        that the probability distribution for the measurand here is best described
        by a shifted and scaled Student's t distribution.  So the standard
        uncertainty should be the standard deviation of this distribution which
        is s*sqrt{(n-1)/[n*(n-3)]}.  Thus

        u(bayesian) = sqrt[dof/(dof - 2)]*u(traditional)

        where dof = n - 1 and the "extra uncertainty" is incorporated directly
        into the standard uncertainty.

        Example
        -------
        >>> gummy.bayesian = True
        >>> g = gummy(1,0.03,dof=5)
        >>> g.bayesian
        True
        """
        return self._bayesian
        
    @staticmethod
    def simulate(nummys,n=100000,ufrom=None):
        if ufrom is not None:
            ufrom = [g._dist for g in ufrom]
        if isinstance(nummys,nummy):
            nummys = [nummys]
        Distribution.simulate([g.distribution for g in nummys],n,ufrom)
        nummy._nsim = n
        
    @property
    def simdata(self):
        """
        `numpy.ndarray`, read-only

        Returns an array containing the Monte-Carlo simulation data.  A 
        `NoSimulatedDataError` is raised if no Monte-Carlo data is available.
        """
        if not isinstance(self._dist,Distribution):
            if nummy._nsim is None:
                raise NoSimulatedDataError()
            return np.full(nummy._nsim,self._dist)
        return self.distribution.simdata
        
    @property
    def simsorted(self):
        """
        `numpy.ndarray`, read-only

        Returns a sorted array containing the Monte-Carlo simulation data.  A 
        `NoSimulatedDataError` is raised if no Monte-Carlo data is available.
        """
        if not isinstance(self._dist,Distribution):
            if nummy._nsim is None:
                raise NoSimulatedDataError
            return np.full(nummy._nsim,self._dist)
        return self.distribution.simsorted
        
    @property
    def xsim(self):
        if not isinstance(self._dist,Distribution):
            return self._dist
        return self.distribution.mean
        
    @property
    def usim(self):
        if not isinstance(self._dist,Distribution):
            return 0
        return self.distribution.stdev
        
    @property
    def cisim(self):
        if not isinstance(self._dist,Distribution):
            return [self._dist,self._dist]
        if self._cimethod == 'shortest':
            return self.distribution.ci(self.p)
        else:
            return self.distribution.cisym(self.p)
              
    @property
    def Usim(self):
        if not isinstance(self._dist,Distribution):
            return 0
        x = self.distribution.mean
        
        if self._cimethod == 'shortest':
            ci = self.distribution.ci(self.p)
        else:
            ci = self.distribution.cisym(self.p)
            
        return (ci[1]-x,x-ci[0])
        
    @property
    def ksim(self):
        """
        read-only

        Returns ``0.5*(gummy.Usim[0] + gummy.Usim[1])/gummy.usim``
        """
        if self.usim == 0 or not isinstance(self._dist,Distribution):
            return float('inf')
        return 0.5*(self.Usim[0] + self.Usim[1])/self.usim
        
    @property
    def independent(self):
        """
        `bool`, read-only

        Returns `False` if the owning gummy was created from a operation involving
        other gummys and `True` otherwise.
        """
        if isinstance(self._dist,Distribution):
            return self._dist.isindependent
        else:
            return False
        
    @property
    def cimethod(self):
        """
        `str` in {'shortest', 'symmetric'}), default is 'shortest'
        
        Get or set the method for calculating the confidence interval from 
        Monte-Carlo data.  If this property is set at the class level, it will
        change the default `cimethod` value for new gummys but will no affect
        gummys that have already been created.
        
        Can be set either to the string 'shortest' or the string 'symmetric'.
        This property gets or sets the method for calculating confidence
        intervals from Monte-Carlo data.  
        
        If it is set to 'shortest', the confidence interval will be taken to be 
        the shortest interval that includes the desired fraction of the probability 
        distribution.  
        
        If it is set to 'symmetric', then the confidence interval will be set so 
        that, for n Monte-Carlo samples and a coverage probability of `p`, then
        `n`*(1-`p`)/2 samples lie below the lower limit of the confidence interval
        and the same number of samples lie above the upper limit of the confidence 
        interval.
        """
        
        return self._cimethod
    @cimethod.setter
    def cimethod(self,value):
        value = value.lower().strip()
        if value not in ['shortest','symmetric']:
            raise ValueError('cimethod ' + str(value) + ' is not recognized')
        self._cimethod = value
        
    @property
    def p(self):
        if self._fp is None:
            return 0.68268949213708585
        return self._fp()
        
    @staticmethod
    def set_seed(seed):
        """
        Sets the seed for the numpy.random.RandomState object shared by all 
        `Distribution` instances.
        """
        
        Distribution.set_seed(seed)
        
    def toummy(self):
        """
        returns an ummy representaion of the nummy
        """
        r = ummy(self.x,u=self.u,dof=self.dof)
        r._ref = self._ref
        r._refs = self._refs
        return r
    
    def splonk(self):
        """
        splonks the nummy
        """
        return self.toummy().splonk()
    
        
    @staticmethod
    def _copy(s,r,formatting=True,tofloat=False):
        # copies attributes of s to r, called from ummy.copy()
        super(nummy,nummy)._copy(s,r,formatting=formatting,tofloat=tofloat)
        if isinstance(s,nummy):
            r._dist = s._dist
            if formatting:
                r.name = s.name
                r._cimethod = s._cimethod
            else:
                r.name = None
        else:
            r.name = None
            if s._u == 0:
                r._dist = r._x
                return
            dof = s._dof
            if isinf(dof):
                r._dist = NormalDist(s._x,s._u)
            else:
                if r._bayesian:
                    r._dist = TDist(s._x,s._u*np.sqrt((dof-2)/dof),dof)
                    r._dof = float('inf')
                else:
                    r._dist = TDist(s._x,s._u,dof)
        
    @classmethod
    def _apply(cls,function,derivative,*args,fxdx=None):
        # called from ummpy.apply()

        a = list(args)
        for i,e in enumerate(args):
            if isinstance(e,nummy):
                a[i] = e._dist
            elif isinstance(e,ummy):
                a[i] = nummy(e)._dist
        r = super(nummy,cls)._apply(function,derivative,*args,fxdx=fxdx)
        if isinstance(r,nummy):
            r._dist = Distribution.apply(function,*a)
        return r
        
    @classmethod
    def _napply(cls,function,*args,fxx=None):
        # called from ummpy.napply()
            
        a = list(args)
        for i,e in enumerate(args):
            if isinstance(e,nummy):
                a[i] = e._dist
            elif isinstance(e,ummy):
                a[i] = nummy(e)._dist
        r = super(nummy,cls)._napply(function,*args,fxx=fxx)
        if isinstance(r,nummy):
            r._dist = Distribution.apply(function,*a)
        return r
        
    @staticmethod
    def _set_correlation_matrix(gummys, matrix):
        super(nummy,nummy)._set_correlation_matrix(gummys,matrix)
        cov = nummy.covariance_matrix(gummys)
        mean = [g.x for g in gummys]
        if isfinite(gummys[0].dof) and all([g.dof==gummys[0].dof for g in gummys]):
            mvdist = MultiTDist(mean,cov,gummys[0].dof)
        else:
            mvdist = MultiNormalDist(mean,cov)
        for i,g in enumerate(gummys):
            g._dist = mvdist[i]
            
    @staticmethod
    def _set_covariance_matrix(gummys, matrix):
        super(nummy,nummy)._set_covariance_matrix(gummys,matrix)
        mean = [g.x for g in gummys]
        if isfinite(gummys[0].dof) and all([g.dof==gummys[0].dof for g in gummys]):
            mvdist = MultiTDist(mean,matrix,gummys[0].dof)
        else:
            mvdist = MultiNormalDist(mean,matrix)
        for i,g in enumerate(gummys):
            g._dist = mvdist[i]
        
    @classmethod
    def create(cls,x,u=0,dof=float('inf'),name=None,correlation_matrix=None,
               covariance_matrix=None):
        if name is None:
            name = [None]*len(x)
                   
        if isinstance(x,MultivariateDistribution):

            # Do this the hard way because we want to raise an exception if
            # a MultivariateElement is passed the the nummy initializer in 
            # case it is tried outside create and the .ref is not set up
            # correctly.
            u = x.u()
            nd = len(u)
            if hasattr(x,'dof'):
                if nummy._bayesian:
                    dof = [float('inf')]*nd
                    for i,uu in u:
                        u[i] *= np.sqrt(x.dof[i]/(x.dof[i]-2))
                else:
                    dof = x.dof
            else:
                dof = [float('inf')]*nd
            
            ret = super(nummy,cls).create(x.x(),u=u,dof=dof,covariance_matrix=x.cov)
                   
            for i,r in enumerate(ret):
                r._dist = x[i]
                r.name = name[i]
                
            return ret
                
        if any([isinstance(v,Distribution) for v in x]):
            if correlation_matrix is not None or covariance_matrix is not None:
                raise TypeError('Distribtuion instances may not be used in x if a correlation_matrix nor a covariance_matrix is defined')
            
        if dof is None:
            dof = [float('inf')]*len(x)
        if u is None:
            u = [0]*len(x)
        d = [None]*len(x)
        for i,v in enumerate(x):
            if isinstance(v,Distribution):
                d[i] = v
                x[i] = v.x()
                u[i] = v.u()
                if hasattr(v,'dof'):
                    dof[i] = v.dof
                    if nummy._bayesian:
                        u[i] = u[i]*np.sqrt(v.dof/(v.dof-2))
                        dof[i] = float('inf')
            
        ret = super(nummy,cls).create(x,u,dof,correlation_matrix,covariance_matrix)
        for i,r in enumerate(ret):
            r.name = name[i]
            if d[i] is None:
                if u[i] == 0:
                    r._dist = x[i]
                elif isinf(dof[i]):
                    r._dist = NormalDist(x[i],u[i])
                else:
                    if nummy._bayesian:
                        r._dist = TDist(x,u[i]*np.sqrt((dof[i]-2)/dof[i]),dof[i])
                        r._dof = float('inf')
                    else:
                        r._dist = TDist(x,u[i],dof[i])
            else:
                r._dist = d[i]
                    
        return ret

    def hist(self,xlabel=None,title=None,hold=False,**kwds):
        if xlabel is None and self.name is not None:
            xlabel = str(self.name)
            if len(xlabel) == 1:
                xlabel = '$ ' + xlabel + ' $'
        self.distribution.hist(xlabel=xlabel,title=title,hold=hold,**kwds)
        
    def covariance_sim(self,g):
        """
        Returns the covariance, calculated from Monte-Carlo data, between the 
        owning gummy and the gummy `g.`
        
        See the method `gummy.covariance(g)` for the corresponding result based
        on first order error propagation.
        """
        #if self._u == 0 or g._u == 0:
            #return 0
        return self._dist.covsim(g._dist)
        
    def correlation_sim(self,g):
        """
        Returns the correlation coefficient, calculated from Monte-Carlo data, 
        between the owning gummy and the gummy `g`.
        
        See the method `gummy.correlation(g)` for the corresponding result based
        on first order error propagation.
        """
        #if self._u == 0 or g._u == 0:
            #return 0
        return self.covariance_sim(g)/(self.usim*g.usim)
        
    @staticmethod
    def covariance_matrix_sim(gummys):
        """
        The staticmethod takes a list of gummys an returns the variance-covariance
        matrix calculated from Monte-Carlo data.  The return value is numpy
        ndarray.
        
        See the method gummy.covariance_matrix(gummys) for the corresponding
        result based on first order error propagation.
        """
        d = [v._dist for v in gummys]
        return Distribution.covsim_matrix(*d)
        
    @staticmethod
    def correlation_matrix_sim(gummys):
        """
        The staticmethod takes a list of gummys an returns the correlation
        matrix calculated from Monte-Carlo data.  The return value is numpy 
        ndarray.
        
        See the method `gummy.correlation_matrix(gummys)` for the corresponding
        result based on first order error propagation.
        """
        cov = nummy.covariance_matrix_sim(gummys)
        return [[cov[i][j]/sqrt(cov[i][i]*cov[j][j]) 
                for i in range(len(gummys))] for j in range(len(gummys))]
                    
    @staticmethod
    def covplot(x,y,title=None,xlabel=None,ylabel=None,hold=False,**kwds):
        if xlabel is None and x.name is not None:
            xlabel = str(x.name)
            if len(xlabel) == 1:
                xlabel = '$ ' + xlabel + ' $'
        if ylabel is None and y.name is not None:
            ylabel = str(y.name)
            if len(ylabel) == 1:
                ylabel = '$ ' + ylabel + ' $'

        Distribution.covplot(x.distribution,y.distribution,xlabel=xlabel,ylabel=ylabel,
                             title=title,hold=hold,**kwds)
        
    def ufrom(self,x,sim=False):
        """
        Gets the standard uncertainty contributed from particular gummys
        or utypes if all other free variables are held fixed.

        Parameters
        ---------_
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
        if not sim:
            return super().ufrom(x)
        
        x = [g for g in ummy._toummylist(x) if self.correlation(g) != 0]
            
        v = nummy.correlation_matrix_sim(x)
            
        b = [self.correlation(z) for z in x]
        s = np.linalg.lstsq(v,b,rcond=-1)[0]
        u = 0
        
        d = [i*self.usim/j.usim for i,j in zip(s,x)]
        for i in range(len(x)):
            for j in range(len(x)):
                u += d[i]*d[j]*x[i].correlation(x[j])*x[i].usim*x[j].usim
                
        return sqrt(u)
        
    def __add__(self,b):
        return _bop(super().__add__(b),np.add,self,b)
        
    def __radd__(self,b):
        return _rbop(super().__radd__(b),np.add,self,b)
        
    def __sub__(self,b):
        return _bop(super().__sub__(b),np.subtract,self,b)
        
    def __rsub__(self,b):
        return _rbop(super().__rsub__(b),np.subtract,self,b)
        
    def __mul__(self,b):
        return _bop(super().__mul__(b),np.multiply,self,b)
        
    def __rmul__(self,b):
        return _rbop(super().__rmul__(b),np.multiply,self,b)
        
    def __truediv__(self,b):
        return _bop(super().__truediv__(b),np.divide,self,b)
        
    def __rtruediv__(self,b):
        return _rbop(super().__rtruediv__(b),np.divide,self,b)
    
    def __floordiv__(self,b):
        return _bop(super().__floordiv__(b),np.floor_divide,self,b)
        
    def __rfloordiv__(self,b):
        return _rbop(super().__rfloordiv__(b),np.floor_divide,self,b)
        
    def __pow__(self,b):
        return _bop(super().__pow__(b),np.power,self,b)
        
    def __rpow__(self,b):
        return _rbop(super().__rpow__(b),np.power,self,b)
    
    def _nprnd(self,f):
        ret = super()._nprnd(f)
        ret._dist = Distribution.apply(f,self._dist)
        return ret
        
    def __mod__(self,b):
        return _bop(super().__mod__(b),np.mod,self,b)
    
    def __rmod__(self,b):
        return _rbop(super().__rmod__(b),np.mod,self,b)
        
    def __neg__(self):
        return _uop(super().__neg__(),np.negative,self)
        
    def __pos__(self):
        return _uop(super().__pos__(),np.positive,self)
        
    def __abs__(self):
        return _uop(super().__abs__(),np.absolute,self)