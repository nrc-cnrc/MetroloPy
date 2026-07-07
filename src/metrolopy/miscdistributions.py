# -*- coding: utf-8 -*-
# module miscdistributions

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
This module contains the classes that represent distributions for Monte Carlo
calculations.
"""
import numpy as np
from .distributions import (Distribution,MultivariateDistribution,
                            MultivariateElement,ScipyStatsDist)
from .dof import DoF

class MultiNormalDist(MultivariateDistribution):
    def __init__(self,mean,cov):
        """
        Multivariate normal distribution.  The parameter `mean` is a list of mean
        values for each dimension and `cov` is the variance-covariance matrix.
        """
        self.mean = np.asarray(mean)
        self._cov = np.asarray(cov)
        nd = len(self.mean)
        if self._cov.shape != (nd,nd):
            ValueError('cov.shape != (len(mean),len(mean))')
        super().__init__(nd)
            
    def _simulate(self,n):
        self.simdata = Distribution.random_rng().multivariate_normal(self.mean,self.cov,n).T
        
    def x(self):
        return self.mean
        
    def u(self):
        return self.cov.diagonal()
        
    @property
    def cov(self):
        return self._cov
    
class MultiTElement(MultivariateElement):
    def __init__(self,parent,index,dof):
        """
        Represents one dimension of a `MultiTDist`.  Do not create
        instances of this class directly.  Instances are created by the parent
        `MultiTDist __init__` method.
        """
        self.parent = parent
        self.index = index
        self.DoF = dof
        
class MultiTDist(MultivariateDistribution):
    def __init__(self,mean,cov,dof):
        """
        Multivariate shifted and scaled Student's t distribution.  The parameter 
        mean is a list of mean values for each dimension and cov is the 
        variance-covariance matrix.  The parameter dof is the number of
        degrees of freedom and must be scalar; all dimensions must have the same 
        number of degrees of freedom.
        """
        self.mean = np.asarray(mean)
        self._cov = np.asarray(cov)
        nd = len(self.mean)
        
        if not isinstance(dof,DoF):
            dof = DoF(dof)
        
        if self._cov.shape != (nd,nd):
            ValueError('cov.shape != (len(mean),len(mean))')
        self._elements = [MultiTElement(self,i,dof) for i in range(nd)]
        self._len = nd
        self.DoF = dof
    
    def _simulate(self,n):
        gamma = np.tile(np.random.gamma(self.DoF.value/2,2/self.DoF.value,n),(len(self),1)).T
        Z = np.random.multivariate_normal(np.zeros(len(self)),self.cov,n)
        self.simdata = (self.mean + Z/np.sqrt(gamma)).T
        
    def x(self):
        return self.mean
        
    def u(self):
        return self.cov.diagonal()
        
    @property
    def cov(self):
        return self._cov
        
class UniformDist(Distribution):
    _param_names = ('center','half_width')
    
    def __init__(self,center=None,half_width=None,lower_limit=None,upper_limit=None):
        """
        A uniform distribution.  Two and only two of the parameters `center`,
        `half_width`, `lower_limit` and `upper_limit` must be specified.
        """
        k = 0
        if center is not None:
            k += 1
            self.center = center
            if half_width is not None:
                self.lower_limit = center - half_width
                self.upper_limit = center + half_width
            elif lower_limit is not None:
                self.half_width = center - lower_limit
                self.upper_limit = center + self.half_width
            elif upper_limit is not None:
                self.half_width = upper_limit - center
                self.lower_limit = center - self.half_width
        if half_width is not None:
            k += 1
            if half_width <= 0:
                raise ValueError('half_width <= 0')
            self.half_width = half_width
            if lower_limit is not None:
                self.center = half_width + lower_limit
                self.upper_limit = self.center + half_width
            elif upper_limit is not None:
                self.center = upper_limit - half_width
                self.lower_limit = self.center - half_width 
        if lower_limit is not None:
            k += 1
            self.lower_limit = lower_limit
            if upper_limit is not None:
                if lower_limit >= upper_limit:
                    raise ValueError('lower_limit >= upper_limit')
                self.center = (upper_limit + lower_limit)/2
                self.half_width = (upper_limit - lower_limit)/2
        if upper_limit is not None:
            k += 1
            self.upper_limit = upper_limit
        if k != 2:
            raise ValueError('Two and only two of the parameters may be specified.')
        self._u = self.half_width/np.sqrt(3)
    
    def random(self,n=None):
        return Distribution.random_rng().uniform(self.lower_limit,self.upper_limit,n)
    
    def x(self):
        return self.center
        
    def u(self):
        return self._u
    
    @staticmethod
    def cdf(z,center,half_width):
        from scipy.stats import uniform
        return uniform.cdf(z,loc=center-half_width,scale=2*half_width)
    
    @staticmethod
    def _est_params(z):
        a = np.min(z)
        b = np.max(z)
        return [(a + b)/2,(b - a)/2]
    
class GammaDist(Distribution):
    _param_names = ('shape','scale')
    
    def __init__(self,shape,scale):
        """
        Gamma distribution with the `shape` and `scale` parameters.
        """
        if shape <= 0:
            raise ValueError('shape < 0')
        if scale <= 0:
            raise ValueError('scale < 0')
        self.shape = shape
        self.scale = scale
        
    def random(self,n=None):
        return Distribution.random_rng().gamma(self.shape,self.scale,n)
        
    def x(self):
        return self.shape*self.scale
        
    def u(self):
        return self.shape*np.sqrt(self.scale)
    
    @staticmethod
    def cdf(z,shape,scale):
        from scipy.stats import gamma
        return gamma.cdf(z,shape,scale=scale)
    
    @staticmethod
    def _est_params(z):
        from scipy.stats import gamma
        p = gamma.fit(z)
        return [p[0],p[2]]
        
class LaplaceDist(Distribution):
    _param_names = ('x','scale')
    
    def __init__(self,x,scale):
        """
        Laplace distribution with location parameter `x` and `scale` parameter.
        """
        if scale <= 0:
            raise ValueError('scale < 0')
        self._x = x
        self.scale = scale
        
    def random(self,n=None):
        return Distribution.random_rng().laplace(self.x(),self.scale,n)
        
    def x(self):
        return self._x
        
    def u(self):
        return self.scale*np.sqrt(2)
    
    @staticmethod
    def cdf(z,x,scale):
        from scipy.stats import laplace
        return laplace.cdf(z,loc=x,scale=scale)
    
    @staticmethod
    def _est_params(z):
        from scipy.stats import laplace
        return laplace.fit(z)
        
class TriangularDist(Distribution):
    _param_names = ('mode','left_width','right_width')
    
    def __init__(self,mode,left_width=None,right_width=None,half_width=None,lower_limit=None,upper_limit=None):
        """
        Triangular distribution.  `mode` is required, then, for a symmetric distribution
        specify `half_width`,  and for non-symmetric distributions specify two and only
        two of the parameters `left_width`, `right_width`, `lower_limit`, `upper_limit`.
        The x value is taken to be the mode.
        """
        self.mode = mode
        if half_width is not None:
            if half_width <= 0:
                raise ValueError('half_width <= 0')
            self.half_width = half_width
            self.left_width = half_width
            self.right_width = half_width
            if left_width is not None:
                raise ValueError('left_width and half_width may not both be specified')
            if right_width is not None:
                raise ValueError('right_width and half_width may not both be specified')
            if mode is not None:
                self.mode = mode
                if lower_limit is not None or upper_limit is not None:
                    raise ValueError('The distribution parameters are over specified.')
                self.lower_limit = mode-half_width
                self.upper_limit = mode+half_width
        else:
            if left_width is not None:
                if left_width <= 0:
                    raise ValueError('left_width <= 0')
                if lower_limit is not None:
                    raise ValueError('The distribution parameters are over specified.')
                self.lower_limit = mode - left_width
                self.left_width = left_width
            else:
                if lower_limit is None:
                    raise ValueError('The distribution parameters are under specified.')
                if upper_limit is not None and upper_limit <= lower_limit:
                    raise ValueError('upper_limit <= lower_limit')
                self.left_width = mode - lower_limit
                self.lower_limit = lower_limit
                
            if right_width is not None:
                if right_width <= 0:
                    raise ValueError('right_width <= 0')
                if upper_limit is not None:
                    raise ValueError('The distribution parameters are over specified.')
                self.upper_limit = mode + right_width
                self.right_width = right_width
            else:
                if upper_limit is None:
                    raise ValueError('The distribution parameters are under specified.')
                self.right_width = upper_limit - mode
                self.upper_limit = upper_limit
            
    def random(self,n=None):
        return Distribution.random_rng().triangular(self.lower_limit,self.mode,self.upper_limit,n)
        
    def x(self):
        return self.mode
        
    def u(self):
        return np.sqrt((self.left_width**2+self.right_width**2+self.left_width*self.right_width)/18)
    
    @staticmethod
    def cdf(z,mode,left_width,right_width):
        from scipy.stats import triang
        return triang.cdf(z,left_width/(left_width+right_width),loc=mode-left_width,scale=left_width+right_width)
    
    @staticmethod
    def _est_params(z):
        m = np.mean(z)
        return [m,np.abs(m - np.min(z)),np.abs(np.max(z) - m)]
        
class ExponentialDist(Distribution):
    _param_names = ('scale',)
    
    def __init__(self,scale=None,rate=None):
        """
        Exponential distribution:
        
        f(x;scale) = exp(-scale*x)/scale
        
        Specify either `scale` or `rate` (`scale` = 1/`rate`).
        """
        if scale is not None and rate is not None:
            raise ValueError('Both scale and rate may not both be specified.')
        if scale is not None:
            if scale <= 0:
                raise ValueError('scale <= 0')
            self.scale = scale
            self.rate = 1/scale
        if rate is not None:
            if rate <= 0:
                raise ValueError('rate <= 0')
            self.rate = rate
            self.scale = 1/rate
            
    def random(self,n=None):
        return Distribution.random_rng().exponential(self.scale,n)
        
    def x(self):
        return self.scale
        
    def u(self):
        return self.scale
    
    @staticmethod
    def cdf(z,scale):
        from scipy.stats import expon
        return expon.cdf(z,scale=scale)
    
    @staticmethod
    def _est_params(z):
        return [np.mean(z)]
        
class PoissonDist(Distribution):
    discrete = True
    
    def __init__(self,lam):
        """
        Poisson distribution with rate parameter `lam`.
        """
        if lam <= 0:
            raise ValueError('lam <= 0')
        self.lam = lam
            
    def random(self,n=None):
        return Distribution.random_rng().poisson(self.lam,n)
        
    def x(self):
        return self.lam
        
    def u(self):
        return np.sqrt(self.lam)
        
        
class BinomialDist(Distribution):
    discrete = True
    
    def __init__(self,n,p):
        """
        Binomial distribution with number of trials `n` and success probability `p`.
        """
        if n <= 0:
            raise ValueError('n <= 0')
        if p <= 0:
            raise ValueError('p <= 0')
        self.n = n
        self.p = p
            
    def random(self,n=None):
        return Distribution.random_rng().binomial(self.n,self.p,n)
        
    def x(self):
        return self.n*self.p
        
    def u(self):
        return np.sqrt(self.n*self.p*(1-self.p))
        

class CurvlinearTrapDist(Distribution):
    def __init__(self,center=None,half_width=None,limit_half_range = None,
                 lower_limit=None,upper_limit=None):
        """
        Curvlinear trapezoidal distribution, `limit_half_range` is required.  Also
        either `center` and `half_width` or `lower_limit` and `upper_limit` are required.
        
        This is intended to represent a variable that follows a uniform distribution
        but where the upper and lower limits are not exactly known and may vary by
        up to the `limit_half_range` from the given lower and upper limit values.
        """
        if limit_half_range is None:
            raise ValueError('The limit_half_range must be specified.')
        self.limit_half_range = limit_half_range
        k = 0
        if center is not None:
            k += 1
            self.center = center
            if half_width is not None:
                self.lower_limit = center - half_width
                self.upper_limit = center + half_width
            elif lower_limit is not None:
                self.half_width = center - lower_limit
                self.upper_limit = center + self.half_width
            elif upper_limit is not None:
                self.half_width = upper_limit - center
                self.lower_limit = center - self.half_width
        if half_width is not None:
            k += 1
            if half_width <= 0:
                raise ValueError('half_width <= 0')
            self.half_width = half_width
            if lower_limit is not None:
                self.center = half_width + lower_limit
                self.upper_limit = self.center + half_width
            elif upper_limit is not None:
                self.center = upper_limit - half_width
                self.lower_limit = self.center - half_width 
        if lower_limit is not None:
            k += 1
            self.lower_limit = upper_limit
            if upper_limit is not None:
                if lower_limit >= upper_limit:
                    raise ValueError('lower_limit >= upper_limit')
                self.center = (upper_limit + lower_limit)/2
                self.half_width = upper_limit - lower_limit
        if upper_limit is not None:
            k += 1
            self.upper_limit = upper_limit
        if k != 2:
            raise ValueError('two and only two of the parameters may be specified')
    
    def random(self,n=None):
        r1 = Distribution.random_rng().uniform(0,1,n)
        r2 = Distribution.random_rng().uniform(0,1,n)
        sa = (self.lower_limit - self.limit_half_range) + 2*self.limit_half_range*r1
        sb = (self.lower_limit + self.upper_limit) - sa
        return sa+(sb - sa)*r2

    def x(self):
        return self.center
        
    def u(self):
        return np.sqrt(self.half_width**2/3 + self.limit_half_range**2/9)
        
class TrapezoidalDist(Distribution):
    _param_names = ('lower_limit','upper_limit','top_to_base_ratio')
    
    def __init__(self,lower_limit,upper_limit,top_to_base_ratio):
        """
        Symmetric trapezoidal distribution
        """
        if top_to_base_ratio > 1:
            raise ValueError('top_to_base_ratio > 1')
        if lower_limit >= upper_limit:
            raise ValueError('lower_limit >= upper_limit')
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.top_to_base_ratio = top_to_base_ratio
        
    def random(self,n=None):
        r1 = Distribution.random_rng().uniform(0,1,n)
        r2 = Distribution.random_rng().uniform(0,1,n)
        return (self.lower_limit + ((self.upper_limit - self.lower_limit)/2)*
                ((1+self.top_to_base_ratio)*r1 + (1-self.top_to_base_ratio)*r2))
    
    def x(self):
        return (self.upper_limit + self.lower_limit)/2
        
    def u(self):
        return (self.upper_limit - self.lower_limit)*np.sqrt((1+self.top_to_base_ratio**2)/24)
    
    @staticmethod
    def cdf(z,lower_limit,upper_limit,top_to_base_ratio):
        from scipy.stats import trapezoid
        c = (1 - top_to_base_ratio)/2
        d = (1 + top_to_base_ratio)/2
        loc = lower_limit
        scale = upper_limit - lower_limit
        return trapezoid.cdf(z,c,d,loc=loc,scale=scale)
    
    @staticmethod
    def _est_params(z):
        return [np.min(z),np.max(z),0.5]
        
class ArcSinDist(Distribution):
    def __init__(self,center=None,half_width=None,lower_limit=None,upper_limit=None):
        """
        Arcsin distribution.  Specify either `center` and `half_width` or `lower_limit`
        and `upper_limit`.
        """
        k = 0
        if center is not None:
            k += 1
            self.center = center
            if half_width is not None:
                self.lower_limit = center - half_width
                self.upper_limit = center + half_width
            elif lower_limit is not None:
                self.half_width = center - lower_limit
                self.upper_limit = center + self.half_width
            elif upper_limit is not None:
                self.half_width = upper_limit - center
                self.lower_limit = center - self.half_width
        if half_width is not None:
            k += 1
            if half_width <= 0:
                raise ValueError('half_width <= 0')
            self.half_width = half_width
            if lower_limit is not None:
                self.center = half_width + lower_limit
                self.upper_limit = self.center + half_width
            elif upper_limit is not None:
                self.center = upper_limit - half_width
                self.lower_limit = self.center - half_width 
        if lower_limit is not None:
            k += 1
            self.lower_limit = upper_limit
            if upper_limit is not None:
                if lower_limit >= upper_limit:
                    raise ValueError('lower_limit >= upper_limit')
                self.center = (upper_limit + lower_limit)/2
                self.half_width = upper_limit - lower_limit
        if upper_limit is not None:
            k += 1
            self.upper_limit = upper_limit
        if k != 2:
            raise ValueError('two and only two of the parameters may be specified')
            
    def random(self,n=None):
        return (self.center + 
            self.half_width*np.sin(2*np.pi*Distribution.random_rng().uniform(0,1,n)))
            
    def x(self):
        return self.center
        
    def u(self):
        return self.half_width/(2*np.sqrt(2))
        
class LogNormalDist(Distribution):
    _param_names = ('mu','sigma')
    
    def __init__(self,mu,sigma):
        """
        Log-normal distribution, `mu` and `sigma` are the mean and standard deviation
        of the logarithm of the random variable.
        """
        self.mu = mu
        self.sigma = sigma
            
    def random(self,n=None):
        return Distribution.random_rng().lognormal(self.mu,self.sigma,n)
        
    def x(self):
        return np.exp(self.mu + self.sigma**2/2)
        
    def u(self):
       return np.sqrt((np.exp(self.sigma**2) - 1)*np.exp(2*self.mu + self.sigma**2))
   
    @staticmethod
    def cdf(z,mu,sigma):
        from scipy.stats import lognorm
        return lognorm.cdf(z,sigma,loc=mu)
    
    @staticmethod
    def _est_params(z):
        from scipy.stats import lognorm
        p = lognorm.fit(z)
        return [p[1],p[0]]
        
class WeibullDist(Distribution):
    _param_names = ('shape','scale')
    
    def __init__(self,shape,scale):
        """
        Weibull distribution with `shape` and `scale` parameters.
        """
        from scipy.special import gamma
        self._gamma = gamma
        if shape < 0:
            raise ValueError('shape < 0')
        if scale <= 0:
            raise ValueError('scale <= 0')
        self.shape = shape
        self.scale = scale
            
    def random(self,n=None):
        return self.scale*Distribution.random_rng().weibull(self.shape,n)
        
    def x(self):
        return self.scale*self._gamma(1 + 1/self.shape)
        
    def u(self):
        return self.scale*np.sqrt(self._gamma(1+2/self.shape) - self._gamma(1+1/self.shape)**2)
    
    @staticmethod
    def cdf(z,shape,scale):
        from scipy.stats import weibull_min
        return weibull_min.cdf(z,shape,scale=scale)
    
    @staticmethod
    def _est_params(z):
        from scipy.stats import weibull_min
        p = weibull_min.fit(z)
        return [p[0],p[2]]
        
class AveragedDist(Distribution):
    def __init__(self,distribution,nsamples):
        """
        Each sample from is the average of `nsamples` drawn from `distribution`.
        """
        if not isinstance(distribution,Distribution):
            from scipy.stats.distributions import rv_frozen
            if not isinstance(distribution,rv_frozen):
                raise TypeError('distribution is not a Distribution or rv_frozen instance')
            distribution = ScipyStatsDist(distribution)
            
        self._dist = distribution
        self.nsamples = nsamples
        
    def random(self,n=None):
        ret = self.dist.random(n*self.nsamples)
        return np.mean(np.reshape(ret,(self.nsamples,n)),axis=0)
    
    def x(self):
        return self._dist.x()
    
    def u(self):
        return self._dist.u()/np.sqrt(self.nsamples)
    
    @property
    def dist(self):
        return self._dist
    

class AveragedErrDist(Distribution):
    def __init__(self,distribution,nsamples):
        """
        Each sample from is the average of `nsamples` drawn from `distribution`
        minus the distribution mean and divided by the standard deviation of the
        nsamples.
        """
        if not isinstance(distribution,Distribution):
            from scipy.stats.distributions import rv_frozen
            if not isinstance(distribution,rv_frozen):
                raise TypeError('distribution is not a Distribution or rv_frozen instance')
            distribution = ScipyStatsDist(distribution)
            
        self._dist = distribution
        self.nsamples = nsamples
        
    def random(self,n=None):
        ret = self.dist.random(n*self.nsamples)
        z = np.reshape(ret,(self.nsamples,n))
        return (np.mean(z,axis=0) - self.dist.x())/(np.std(z,axis=0,ddof=1)/np.sqrt(self.nsamples))
    
    def x(self):
        return 0
    
    def u(self):
        return 1
    
    @property
    def dist(self):
        return self._dist