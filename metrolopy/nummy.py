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
from .ummy import ummy,_udict,_isscalar
from .distributions import (Distribution,TDist,NormalDist,MultivariateElement,
                            MultivariateDistribution,Convolution)
from .exceptions import NoSimulatedDataError
from math import isinf,isnan,sqrt
from html import escape

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

def _getnrefs(x,fr):
    if isinstance(x,ummy):
        x = [x]
    ret = set()
    for i in x:
        for r in x._getrefs(fr):
            if r.dist is not None:
                ret.add(r.dist)
    return ret

def get_name(name,fmt,norm):
    fmt = fmt.strip().lower()
    if fmt not in {'unicode','html','latex','ascii'}:
        raise ValueError('fmt "' + str(fmt) + '" is not recognized')
        
    if name is None:
        return None
    
    if isinstance(name,str):
        name = name.strip()
        if fmt == 'html':
            name = escape(name)
            if len(name) == 1:
                name = '<i>' + name + '</i>'
        elif fmt == 'latex':
            sc = {
                '&': r'\&',
                '%': r'\%',
                '$': r'\$',
                '#': r'\#',
                '_': r'\_',
                '{': r'\{',
                '}': r'\}',
                '~': r'\_',
                '^': r'\_',
                '\\': r'\_'
                }
            name = ''.join(sc[c] if c in sc else c for c in name)
            if len(name) > 1:
                name = norm(name)
        elif fmt == 'ascii':
            name = ''.join(i if ord(i) < 128 else '_' for i in name)
        return name
    
    if fmt == 'unicode':
        return name[0]
    if fmt == 'html':
        return name[1]
    if fmt == 'latex':
        return name[2]
    if fmt == 'ascii':
        return name[3]

    
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
            
            if utype is None:
                utype = x.utype
            else:
                x.utype = utype
            self._dist = x
            
            if hasattr(x,'dof'):
                ndof = x.dof
                if nummy._bayesian:
                    u = float(x.u())*np.sqrt(x.dof/(x.dof-2))
                else:
                    u = float(x.u())
            else:
                ndof = float('inf')
                u = float(x.u())
                
            super().__init__(x.x(),u=u,dof=ndof,utype=utype)
            
        else:
            super().__init__(x,u=u,dof=dof,utype=utype)
            x = float(x)
            u = float(self.u)
            if u == 0 or isnan(x) or isinf(x) or isnan(u) or isinf(u):
                self._dist = x
                return
            if isinstance(dof,_udict):
                self._dist = None
            else:
                if dof is None or dof > self.max_dof:
                    self._dist = NormalDist(x,u)
                else:
                    if nummy._bayesian:
                        self._dist = TDist(x,u*np.sqrt((dof-2)/dof),dof)
                    else:
                        self._dist = TDist(x,u,dof)
                self._dist.utype = utype
            
    @property
    def distribution(self):
        """
        read-only

        Returns ths `Distribution` instance associated with the gummy.
        """
        return self._dist
    
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
        if norm is None:
            norm = type(self).latex_norm
        return get_name(self._name,fmt,norm)
    
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
            if isinstance(ufrom,(nummy,str,Distribution)):
                ufrom = [ufrom]
            ufrom = [g.distribution if isinstance(g,nummy) else g for g in ufrom]
        if isinstance(nummys,nummy):
            nummys = [nummys]
        Distribution.simulate([g.distribution for g in nummys if isinstance(g,nummy)],n=n,ufrom=ufrom)
        nummy._nsim = n
        
    def sim(self,n=100000,ufrom=None):
        return nummy.simulate([self],n=n,ufrom=ufrom)
    
    def clear(self):
        if isinstance(self._dist,Distribution):
            self._dist.clear()
        
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
    def isindependent(self):
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
        r = ummy(self.x,u=self.u)
        r._ref = self._ref
        r._refs = self._refs
        return r
    
    def splonk(self):
        """
        splonks the nummy
        """
        return self.toummy().splonk()
    
        
    @staticmethod
    def _copy(s,r,formatting=True,totype=None):
        # copies attributes of s to r, called from ummy.copy()
        super(nummy,nummy)._copy(s,r,formatting=formatting,totype=totype)
        if isinstance(s,nummy):
            r._dist = s._dist
            if formatting:
                r.name = s.name
                r._cimethod = s._cimethod
            else:
                r.name = None
        else:
            r.name = None
            if s.u == 0:
                r._dist = r.x
                return
            dof = s.dof
            if isinf(dof):
                r._dist = NormalDist(s.x,s.u)
            else:
                if r._bayesian:
                    r._dist = TDist(s.x,s.u*np.sqrt((dof-2)/dof),dof)
                    r.dof = float('inf')
                else:
                    r._dist = TDist(s.x,s.u,dof)
        
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
        
    @classmethod
    def create(cls,x,u=0,dof=float('inf'),name=None,utype=None,
               correlation_matrix=None,covariance_matrix=None):
        if _isscalar(name):
            name = [name]*len(x)
                   
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
                dof = float('inf')
            
            ret = super(nummy,cls).create(x.x(),u=u,dof=dof,utype=utype,
                                          covariance_matrix=x.cov)
                   
            for i,r in enumerate(ret):
                r._dist = x[i]
                r.name = name[i]
                
            return ret
                
        if any([isinstance(v,Distribution) for v in x]):
            if correlation_matrix is not None or covariance_matrix is not None:
                raise TypeError('Distribtuion instances may not be used in x if a correlation_matrix nor a covariance_matrix is defined')
            
        if dof is None:
            dof = [float('inf')]*len(x)
        elif _isscalar(dof):
            dof = [dof]*len(x)
        else:
            dof = [float('inf') if i is None else i for i in dof]
            
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
            
        ret = super(nummy,cls).create(x,u=u,dof=dof,utype=utype,
                                      correlation_matrix=correlation_matrix,
                                      covariance_matrix=covariance_matrix)
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
                    else:
                        r._dist = TDist(x,u[i],dof[i])
            else:
                r._dist = d[i]
                    
        return ret

    def hist(self,**kwds):
        if not isinstance(self._dist,Distribution):
            raise TypeError('hist may not be called from a constant nummy')
        self.distribution.hist(**kwds)
        
    def covariance_sim(self,g):
        """
        Returns the covariance, calculated from Monte-Carlo data, between the 
        owning gummy and the gummy `g.`
        
        See the method `gummy.covariance(g)` for the corresponding result based
        on first order error propagation.
        """
        return self._dist.covsim(g._dist)
        
    def correlation_sim(self,g):
        """
        Returns the correlation coefficient, calculated from Monte-Carlo data, 
        between the owning gummy and the gummy `g`.
        
        See the method `gummy.correlation(g)` for the corresponding result based
        on first order error propagation.
        """
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
            xlabel = x.get_name(fmt='latex')
        if ylabel is None and y.name is not None:
            ylabel = y.get_name(fmt='latex')

        Distribution.covplot(x.distribution,y.distribution,xlabel=xlabel,ylabel=ylabel,
                             title=title,hold=hold,**kwds)
        
    #def ufrom(self,x,sim=False):

        #if not sim:
            #return super().ufrom(x)
            
        #v = nummy.correlation_matrix_sim(x)
            
        #b = [self.correlation(z) for z in x]
        #s = np.linalg.lstsq(v,b,rcond=-1)[0]
        #u = 0
        
        #d = [i*self.usim/j.usim for i,j in zip(s,x)]
        #for i in range(len(x)):
            #for j in range(len(x)):
               # u += d[i]*d[j]*x[i].correlation(x[j])*x[i].usim*x[j].usim
                
        #return sqrt(u)
    
    def ufromsim(self,x):
        """
        Gets the standard deviation of the Monte-Carlo data only allowing the
        independent variables in `x` to vary.  Independent istributions not in 
        `x` are held fixed.  `sim` or `simulate` must be called to generate
        Monte-Carlo data before calling this method.
        
        Parameters
        ----------
        x:  `nummy`, `str`, or array_like
            A nummy, a string referencing a utype or a list containing
            nummys and strings.
            
        Returns
        -------
        `float`
        """
        return float(np.std(self.datafrom(x,save=False),ddof=1))
    
    def datafrom(self,x,save=True):
        """
        Recomputes the convolution with only the varaibles in `x` allowed to
        vary.  `sim` or `simulate` must be called to generate
        Monte-Carlo data before calling this method.  This method cannot be 
        called with save == `True` from from a nummy representing an independent 
        variable (that is from a nummy not created by by mathematical operations
        between two or more other nummy's).
        
        Parameters
        ----------
        x:  list containing `nummy` or `str`
            all independent Distributions not in the list or having a utype
            not in the list are held fixed at their `.x()` value
            
        save:  If `save` is `True` the recomputed data is stored in the `simdata`
            attribute and `None` is returned.  If `save` is `False` then the
            recomputed data is returned and the `simdata` attribute is not
            overwritten.
        
        Returns
        -------
        'numpy.array' if `save` is `False`, otherwise returns `None`

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate`.
        `RuntimeError`:  
            if this method is called from an independent `nummy` 
        """
        if not isinstance(self._dist,Distribution):
            if nummy._nsim is None:
                raise NoSimulatedDataError
            return np.full(nummy._nsim,self._dist)
        
        if not isinstance(self.distribution,Convolution) and save == 'True':
            raise RuntimeError('datafrom may not be called from an independent nummy with save == True')
        
        if self.distribution.simdata is None:
            raise NoSimulatedDataError('no simulated data')
            
        if isinstance(x,(nummy,str,Distribution)):
            x = [x]
        x = [i.distribution if isinstance(i,nummy) else i for i in x]
        
        if not isinstance(self.distribution,Convolution):
            if self.distribution in x or self.distribution.utype in x:
                ret = self.distribution.simdata
            else:
                ret = np.full(self.distribution.simdata,self.distribution.x())
        else:
            ret = self.distribution.datafrom(x,save=save)
            
        if save:
            self.distribution.simdata = ret
        else:
            return ret
            
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