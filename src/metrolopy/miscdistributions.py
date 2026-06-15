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
from .distributions import Distribution,MultivariateDistribution,MultivariateElement

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
        self.dof = dof
        
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
        
        try:
            if len(dof) == nd:
                dof = np.asarray(dof)
                for d in dof:
                    if d != dof[0]:
                        raise ValueError('all dimensions must have the same dof')
            else:
                raise TypeError('dof must be a scalar or be the same length as mean')
        except TypeError:
            self.dof = np.array(nd*[dof])
        
        if self._cov.shape != (nd,nd):
            ValueError('cov.shape != (len(mean),len(mean))')
        self._elements = [MultiTElement(self,i,dof) for i in range(nd)]
        self._len = nd
    
    def _simulate(self,n):
        gamma = np.tile(np.random.gamma(self.dof[0]/2,2/self.dof[0],n),(len(self),1)).T
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
    
class GammaDist(Distribution):
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
        
class LaplaceDist(Distribution):
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
        
class TriangularDist(Distribution):
    def __init__(self,mode,half_width=None,left_width=None,right_width=None,lower_limit=None,upper_limit=None):
        """
        Triangular distribution.  `mode` is required, then, for a symmetric distribution
        specify `half_width`,  and for non-symmetric distributions specify two and only
        two of the parameters `left_width`, `right_width`, `lower_limit`, `upper_limit`.
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
        
class ExponentialDist(Distribution):
    def __init__(self,scale=None,rate=None):
        """
        Exponential distribution:
        
        f(x;rate) = rate*exp(-rate*x)
        
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
        
class PoissonDist(Distribution):
    def __init__(self,lam):
        """
        Poisson distribution with rate parameter `lam`.
        """
        if lam <= 0:
            raise ValueError('lam <= 0')
        if np.modf(lam)[1] != 0:
            raise ValueError('lam is not an integer')
        self.lam = lam
            
    def random(self,n=None):
        return Distribution.random_rng().poisson(self.lam,n)
        
    def x(self):
        return self.lam
        
    def u(self):
        return np.sqrt(self.lam)
        
        
class BinomialDist(Distribution):
    def __init__(self,n,p):
        """
        Binomial distribution with number of trials `n` and success probability `p`.
        """
        if n <= 0:
            raise ValueError('n <= 0')
        if np.modf(n)[1] != 0:
            raise ValueError('n is not an integer')
        if p <= 0:
            raise ValueError('p <= 0')
        if np.modf(p)[1] != 0:
            raise ValueError('p is not an integer')
        self.n = n
        self.p = p
            
    def random(self,n=None):
        return Distribution.random_rng().binomial(self.n,self.p,n)
        
    def x(self):
        return self.n*self.p
        
    def u(self):
        return np.sqrt(self.n*self.p(1-self.p))
        

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
        sa = (self.lower_limt - self.limit_half_range) + 2*self.limit_half_range*r1
        sb = (self.lower_limit + self.upper_limit) - sa
        return sa+(sb - sa)*r2

    def x(self):
        return self.center
        
    def u(self):
        return np.sqrt(self.half_width**2/3 + self.limit_half_range**2/9)
        
class TrapezoidalDist(Distribution):
    def __init__(self,lower_limit,upper_limit,top_to_base_ratio):
        """
        Trapezoidal distribution
        """
        if top_to_base_ratio > 1:
            raise ValueError('top_to_base_ratio > 1')
        if lower_limit >= upper_limit:
            raise ValueError('lower_limit >= upper_limit')
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        
    def random(self,n=None):
        r1 = Distribution.random_rng().uniform(0,1,n)
        r2 = Distribution.random_rng().uniform(0,1,n)
        return (self.lower_limit + ((self.upper_limit - self.lower_limit)/2)*
                ((1+self.top_to_base_ratio)*r1 + (1-self.top_to_base_ratio)*r2))
    
    def x(self):
        return (self.upper_limit - self.lower_limit)/2
        
    def u(self):
        return (self.upper_limit - self.lower_limit)*np.sqrt((1+self.top_to_base_ratio**2)/24)
        
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
        return np.exp(self.mu)
        
    def u(self):
       return self.sigma*np.exp(self.mu)
        
class WeibullDist(Distribution):
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
        
class AveragedFrom(Distribution):
    def __init__(self,distribution,nsamples):
        """
        Each sample from is the average of `nsamples` drawn from `distribution`.
        """
        self._dist = distribution
        self.nsamples = nsamples
        self.dof = nsamples - 1
        
    def random(self,n=None):
        ret = self._dist.random(n*self.nsamples)
        return np.mean(np.reshape(ret,(self.nsamples,n)),axis=0)
    
    def x(self):
        return self._dist.x()
    
    def u(self):
        return self._dist.u()/np.sqrt(self.nsamples)
    
