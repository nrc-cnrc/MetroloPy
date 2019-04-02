# -*- coding: utf-8 -*-

# module distributions

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
This module contains the classed that represent distributions for use with 
gummy objects
"""

from .exceptions import NoSimulatedDataError
#from .unit import Unit
import numpy as np

class Distribution:
    """
    Abstract base class for distributions used for Monte-Carlo uncertainty
    propagation.
    
    In a derived class define the following methods:
        
    random(n=None):  Return a numpy array of n values drawn from the distribution.
        If `n` is `None` then a single scalar value should be returned.  Preferably
        use, as a random number generator, the `numpy.RandomState` object accessed
        with the `Distribution.random_state` static method.
        
    x():  A scalar "center" of the distribution.  This is used to get the `x` value
        of a gummy defined with the distribution.
        
    u():  A scalar "standard uncertainty" of the distribution (usually the standard
        deviation).  This is used to get the u value of a gummy defined with
        the distribution.

    Attributes
    ----------
    simdata: numpy.ndarray or `None`
        Simulated data following a call from the `simulate` static method.

    isindependent: `bool`
        `False` if the class is a `Convolution`, `True` otherwise
        
    Example
    -------
    >>> class ChiSquaredDist(Distribution):
    ...    def __init__(self,dof):
    ...        self.dof = dof
    ...
    ...    def random(self,n=None):
    ...        return Distribution.random_state().chisquare(self.dof,n)
    ...
    ...    def x(self):
    ...        return self.dof
    ...
    ...    def u(self):
    ...        return 2*self.dof
    
    """
    _random_state = None
    simdata = None
    
    isindependent = True
    
    _called = []
    _mean = None
    _stdev =  None
    _simsorted = None
    _ci = None
    _cisym = None
    _cov = None
    _used = False
    
    @staticmethod
    def random_state():
        """
        Returns the `numpy.random.RandomState` object shared by all distributions.
        """
        if Distribution._random_state is None:
            Distribution._random_state = np.random.RandomState()
        return Distribution._random_state
    
    @staticmethod
    def set_seed(seed):
        """
        Sets the seed of the `numpy.random.RandomState` object shared by all
        distributions.
        """
        Distribution.random_state().seed(seed)
    
    @staticmethod
    def apply(f,*d):
        """
        Applies a function `f` to distribution(s) `d1`, `d2`, ..., `f`(`d1`,`d2`,...)
        and returns a distribution resulting from the convolution.

        Parameters
        ----------
        f: function
            Any function that takes an appropriate number of scalar values
            and returns a scalar.
        
        d1, d2, ...: `Distribution` or scalar
            the distributions to be used as arguments to f.  Fixed values
            (e.g. `float`) can also be used for some or all arguments.

        Returns
        -------
        `Convolution` or `float`:
            if any of `d1`, `d2`, ... is a `Distribution`, f(d1,d2,...)
            otherwise
        """
        if not any([isinstance(v,Distribution) for v in d]):
            return f(*d)
        return Convolution(f,*d)
        
    @staticmethod
    def simulate(distributions,n=100000,ufrom=None):
        """
        Generates simulated data for the desired distributions.

        Parameters
        ----------
        distributions: `list` of `Distribution`
            Each Distribution in the list will have the simulated data
            available in the `simdata` attribute. Data will be cleared
            from all Distribution instances not in the list.
        
        n:  int > 0, optional
            The number of samples to generate for each distribution.  The
            default is 100000
            
        ufrom:  `list` of `Distributions` or `None`
            If this parameter is not `None`, then only variables represented
            by distributions in this list will be allowed to vary.  All other
            variables will be held fixed at the distribution x value.
        """
        for d in Distribution._called:
            d.clear()
        Distribution._called = []

        for d in distributions:
            if isinstance(d,Distribution):
                d._simulate(int(n),ufrom)
            
        for d in Distribution._called:
            if d not in distributions:
                d.clear()
                Distribution._called.remove(d)
         
    @property
    def mean(self):
        """
        Returns the mean of the simulated data.

        Raises
        ------
        A `NoSimulatedDataError` will be raised if no simulated data is available
        from a call to `Distribution.simulate()`.
        """
        if self._mean is not None:
            return self._mean
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        self._mean = np.mean(self.simdata)
        return self._mean
        
    @property
    def stdev(self):
        """
        Returns the standard deviation of the simulated data.

        Raises
        ------
        A NoSimulatedDataError will be raised if no simulated data is available 
        from a call to Distribution.simulate().
        """
        if self._stdev is not None:
            return self._stdev
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        self._stdev = np.std(self.simdata,ddof=1)
        
        if np.isnan(self._stdev):
            if np.any(np.isinf(self.simdata)):
                self._stdev = float('inf')
                
        return self._stdev
        
    @property
    def simsorted(self):
        """
        numpy.ndarray, read-only
        
        Returns a sorted numpy array containing the simulated data values.

        Raises
        ------
        A NoSimulatedDataError will be raised if no simulated data is available 
        from a call to Distribution.simulate().
        """
        if self._simsorted is not None:
            return self._simsorted
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        self._simsorted = np.sort(self.simdata)
        return self._simsorted
        
    def ci(self,p):
        """
        Returns the shortest interval that contains the fraction p of the simulated
        data values.

        Returns
        -------
        `tuple` of `float`:
            a tuple containing the lower and upper limits of the interval


        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate`.

        See Also
        --------
        cisym
        """
        if p > 1:
            raise ValueError('p > 1')
        if p < 0:
            raise ValueError('p < 0')
        if self._ci is None:
            self._ci = {}
        elif p in self._ci:
            return self._ci[p]
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        s = self.simsorted
        n = int((1-p)*len(self.simdata)+0.5)
        if n == 0:
            mn = np.argmin(s[-n:] - s[:n])
            lower = upper = (s[mn] + s[mn-1])/2
        else:
            mn = np.argmin(s[-n:] - s[:n])
            lower = s[mn]
            upper = s[mn-n]
        self._ci[p] = (lower, upper)
        return (lower, upper)
        
    def cisym(self,p):
        """
        Returns the interval that contains the fraction p of the simulated
        data values and with an equal number of values below the lower limit of
        the interval and above the upper limit of the interval.

        Returns
        -------
        `tuple` of `float`:
            a tuple containing the lower and upper limits of the interval

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate`.

        See Also
        --------
        ci
        """
        if p > 1:
            raise ValueError('p > 1')
        if p < 0:
            raise ValueError('p < 0')
        if self._cisym is None:
            self._cisym = {}
        elif p in self._cisym:
            return self._cisym[p]
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        s = self.simsorted
        lower = np.percentile(s,50*(1-p))
        upper = np.percentile(s,50*(1+p))
        self._cisym = (lower, upper)
        return (lower, upper)
        
    def hist(self,hold=False,xlabel=None,ylabel='$ \\mathrm{probability\\:density} $',
             title=None,**kwds):
        """
        Generates a histogram from the simulated data.
        
        Parameters
        ----------
        hold:  `bool`, optional
            If this is `False`, pyplot.show() will be called before the method
            exits.  The default is `False`
            
        xlabel:  `str` or `None`, optional
            a label for the histogram horizontal axis, the default is `None`
        
        ylabel: `str` or `None`, optional
            a label for the histogram vertical axis, the default is 'probability
            density'
            
        title: `str` or `None`
            a title for the histogram, the default is `None`
        
        kwds:
            additional key words that will be passed to the `pyplot.hist` method
            that actually created the histogram

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate`.
        """
        import matplotlib.pyplot as plt
        if self.simdata is None:
            raise NoSimulatedDataError('simulate must be called to generate some data')
        if not 'bins' in kwds:
            kwds['bins'] = 100
        if not ('normed' in kwds or 'density' in kwds):
            try:
                from matplotlib import __version__ as version
                version = [int(v) for v in version.split('.')[:2]]
                usedensity = (version[0] > 2 or 
                              (version[0] == 2 and version[1] > 0))
            except:
                usedensity = True
            if usedensity:
                kwds['density'] = True
            else:
                kwds['normed'] = True
        if not 'histtype' in kwds:
            kwds['histtype'] ='stepfilled'
        plt.hist(self.simsorted[np.abs(self.simsorted) != np.inf],**kwds)
        if ylabel is not None:
            plt.ylabel(ylabel)
        if xlabel is not None:
            plt.xlabel(xlabel)
        if title is not None:
            plt.title(title)
        if not hold:
            plt.show()
                
    def clear(self):
        """
        Clears the simulated data.
        """
        self.simdata = None
        self._mean = None
        self._stdev = None
        self._simsorted = None
        self._ci = {}
        self._cisym = {}
        self._cov = {}
        
    def covsim(self,d):
        """
        Returns the covariance between this `Distribution` instance and another
        `Distribution` `d`.

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate` for both `self` and `d`.
        """
        if self._cov is None:
            self._cov = {}
        if d in self._cov:
            return self._cov[d]
        if self.simdata is None or d.simdata is None:
            raise NoSimulatedDataError('simulated data does not exist for both distributions')
        cov = np.cov([self.simdata,d.simdata])[1][0]
        self._cov[d] = cov
        return cov
        
    @staticmethod
    def covsim_matrix(*d):
        """
        Returns the variance-covariance matrix of the Distributions `d1`, `d2`, ...

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate` for any of `d1`, `d2`, ...
        """
        if any([(v.simdata is None) for v in d]):
            raise NoSimulatedDataError('simulated data does not exist for all distributions')
        return np.cov([v.simdata for v in d])
        
    @staticmethod
    def covplot(x,y,fmt='ko',xlabel=None,ylabel=None,title=None,hold=False,**kwds):
        """
        Plots the Distribution `x` versus the Distribution `y`

        Parameters
        ----------
        x, y:  `Distribution`
            the distributions to be plotted
            
        fmt:   `str`, optional
            Format parameter passed to `pyplot.plot()`, the default is 'ko'
        
        xlabel:  `str` or `None`
         a label for the plot x-axis, the default is `None`
        
        ylabel:  `str` or `None`
         a label for the plot y-axis, the default is `None`
            
        title: `str` or `None`
            a title for the histogram, the default is `None`
        
        hold:  `bool` i
            If this is `False`, `pyplot.show()` will be called before the
            method exits.  The default is `False`
            
        kwds:
            additional key words that will be passed to `pyplot.plot()`

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate` for `x` and `y`.
        """
        import matplotlib.pyplot as plt
        
        if x.simdata is None or y.simdata is None:
            raise NoSimulatedDataError('simulated data does not exist for both distributions')
        if 'ms' not in kwds and 'markersize' not in kwds:
            n = len(x.simdata)
            if n >= 100000:
                kwds['ms'] = 0.2
            elif n >= 50000:
                kwds['ms'] = 0.2
            elif n >= 20000:
                kwds['ms'] = 0.5
            else:
                kwds['ms'] = 1
        plt.plot(x.simdata,y.simdata,fmt,**kwds)
        if title is None:
            cov = Distribution.covsim_matrix(x,y)
            cor = cov[1][0]/np.sqrt(cov[1][1]*cov[0][0])
            if abs(cor) < 0.0005:
                cor = '0.000'
            elif abs(cor) < 0.1:
                cor = '{:.3f}'.format(cor)
            else:
                cor = '{:.2f}'.format(cor)
            title = '$ \\mathrm{correlation} = ' + cor + ' $'
            plt.title(title)
        if xlabel is None:
            xlabel = '$ x $'
        plt.xlabel(xlabel)
        if ylabel is None:
            ylabel = '$ y $'
        plt.ylabel(ylabel)
        if not hold:
            plt.show()
        
    def _simulate(self,n,ufrom):
        if ufrom is not None and self not in ufrom:
            self.simdata = np.full(n,self.x())
            return
        Distribution._called.append(self)
        self.simdata = self.random(n)
        
    def __add__(self,v):
        return Convolution(np.add,self,v)
        
    def __radd__(self,v):
        return Convolution(np.add,v,self)
        
    def __sub__(self,v):
        return Convolution(np.subtract,self,v)
        
    def __rsub__(self,v):
        return Convolution(np.subtract,v,self)
        
    def __mul__(self,v):
        return Convolution(np.multiply,self,v)
        
    def __rmul__(self,v):
        return Convolution(np.multiply,v,self)
        
    def __truediv__(self,v):
        return Convolution(np.divide,self,v)
        
    def __rtruediv__(self,v):
        return Convolution(np.divide,v,self)
        
    def __pow__(self,v):
        return Convolution(np.power,self,v)
        
    def __rpow__(self,v):
        return Convolution(np.power,v,self)
        
    def __abs__(self):
        return Convolution(np.abs,self)
        
    def __neg__(self):
        return Convolution(np.negative,self)
        
    def __pos__(self):
        return self
        
    def random(self,n=None):
        """
        Override this method in a derived class

        Return a numpy array of n values drawn from the distribution.
        If n is None then a single scalar value should be returned.
        """

        raise NotImplementedError()
        
    def x(self):
        """
        Override this method in a derived class

        Return a scalar "center" of the distribution.  This is used to
        get the `x` value of a gummy defined with the distribution.
        """
        raise NotImplementedError()
        
    def u(self):
        """
        Override this method in a derived class

        Return a scalar "standard uncertainty" of the distribution (usually
         the standard deviation).  This is used to get the `u` value of a
         gummy defined with the distribution.
        """
        # Override in a derived class.
        # A
        raise NotImplementedError()
    

class Convolution(Distribution):
    isindependent = False
    
    def __init__(self,func,*args):
        """
        Represents the distribution resulting from applying func to d1, d2, ...
        func(d1,d2,...)

        Parameters
        ----------
        func:  function
            a function that takes an appropriate number of scalar arguments and
            returns a scalar.
        
        *args: `Distribution` or `float`
            arguments for `func`
        """
        self.func = func
        self.args = args
        
    def random(self,n=None):
        raise NotImplementedError('use the simulate static method to generate data from a Convolution')
        
    def x(self):
        args = [a.x() if isinstance(a,Distribution) else a for a in self.args]
        return self.func(*args)
        
    def _random(self,n,ufrom):
        for a in self.args:
            if isinstance(a,Distribution) and (a not in Distribution._called):
                a._simulate(n,ufrom)
        v = [a.simdata if isinstance(a,Distribution) else a for a in self.args]
        
        if not isinstance(self.func,np.ufunc):
            # If func is not a numpy ufunc, see if it broadcasts like one
            try:
                ntst = max(n,self.args+5)
                vtst = [a.simdata[:ntst] if isinstance(a,Distribution) else a for a in self.args]
                tst = self.func(*vtst)
                if len(tst) != ntst:
                    raise TypeError()
                    
                func = self.func
            except:
                # if it doesn't broadcast, make it broadcast with frompyfunc
                func = np.frompyfunc(self.func,len(v),1)
        else:
            func = self.func
            
        ret = func(*v)

        if ret.dtype != np.float64:
            ret = np.array(ret,dtype=np.float64)
            
        return ret
                
    def _simulate(self,n,ufrom):
        Distribution._called.append(self)
        self.simdata = self._random(n,ufrom)
        

class MultivariateDistribution:
    def __init__(self,nd):
        """
        Base class for multivarate distributions.  `nd` is the number of dimensions
        of the distribution.
        
        To create a multi-variate distribution, inherit from the 
        MultivariateDistribution and define the following methods:

        _simulate(n):  Return a numpy array of n samples drawn from the distribution. 
            Preferably use, as a random number generator, the numpy RandomState 
            object accessed with the Distribution.random_state() static method.
                
        x():  a list or array with the "center" of the distribution for each 
            dimension.  This is used to get the x value of a gummys defined with 
            the distribution.
                
        u():  A list or array with "standard uncertainty" of the distribution 
            (usually the standard deviation) for each dimension.  This is used 
            to get the u value of a gummy defined with the distribution.
        
        and the following read-only property
        
        cov:  Returns the variance-covariance matrix
        
        The __init__ function must also call the MultivariateDistribution 
        __init__ with the number of dimensions nd of the distribution 
        e.g. super().__init__(nd).
        
        Example
        -------

                class DirechletDist(MultivariateDistribution):
                    def __init__(self,alpha):
                        self.alpha = alpha

                        a0 = sum(alpha)
                        b = (a0**2*(a0 + 1))
                        self._x = [a/a0 for a in alpha]
                        self._u = [a*(a0 - a)/b for a in alpha]
                        self._cov = [[ai*(a0 - ai)/b if i == j else -ai*aj/b
                                      for i,ai in enumerate(alpha)]
                                      for j,aj in enumerate(alpha)]

                        super().__init__(len(alpha))

                    def _simulate(self,n):
                        self.simdata = Distribution.random_state().dirichlet(self.alpha,n).T

                    def x(self):
                        return self._x

                    def u(self):
                        return self._u

                    @property
                    def cov(self):
                        return self._cov
        """
        self._elements = [MultivariateElement(self,i) for i in range(nd)]
        self._len = nd
        
    def __getitem__(self,i):
        return self._elements[i]
        
    def __iter__(self):
        return iter(self.elements)
        
    def __len__(self):
        return self._len
        
    def x(self):
        raise NotImplementedError()
        
    def u(self):
         raise NotImplementedError()
         
    @property
    def cov(self):
        raise NotImplementedError()
        
    def clear(self):
        self.simdata = None
         

class MultivariateElement(Distribution):
    def __init__(self,parent,index):
        """
        Represents one dimension of a multivariate distribution.  Do not create
        instances of this class directly.  Instances are created by the parent
        `MultivariateDistribution __init__` method.
        """
        self.parent = parent
        self.index = index
        
    def random(self,n=None):
        raise NotImplementedError('with a MuliVaraiteElement, data can only be generated with the simulate static method')
        
    def _simulate(self,n,ufrom):
        if ufrom is not None and self not in ufrom:
            self.simdata = np.full(n,self.x())
            return
        if self.parent not in Distribution._called:
            Distribution._called.append(self.parent)
            self.parent._simulate(n)
        self.simdata = self.parent.simdata[self.index]
        
    def cov(self,d):
        if not isinstance(d,MultivariateElement):
            return 0.0
        if d.parent is not self.parent:
            return 0.0
        return self.parent.cov[self.index][d.index]
        
    def clear(self):
        self.parent.simdata = None
        super().clear()
        
    def x(self):
        return self.parent.x()[self.index]
        
    def u(self):
        return self.parent.cov.diagonal()[self.index]

class NormalDist(Distribution):
    def __init__(self,x,s):
        """
        Normal distribution with mean `x` and standard deviation `s`.
        """
        if s < 0:
            raise ValueError('s < 0')
        self._x = x
        self._u = s
        
    def random(self,n=None):
        if self._u == 0:
            if n == None:
                return self._x
            else:
                return np.full(n,self._x)
        return Distribution.random_state().normal(self._x,self._u,n)
        
    def x(self):
        return self._x
        
    def u(self):
        return self._u
        
        
class TDist(Distribution):
    bayesian_default = False
    
    def __init__(self,x,s,dof):
        """
        Shifted and scaled Student's t distribution with `dof` degrees of freedom,
        mean `x`, and scale factor `s`.
        """
        if dof <= 0:
            raise ValueError('dof <= 0')
        self._x = x
        self._s = s
        self.dof = dof
        
    def random(self,n=None):
        if self._s == 0:
            if n == None:
                return self._x
            else:
                return np.full(n,self._x)
        return self._x + self._s*Distribution.random_state().standard_t(self.dof,n)
        
    def x(self):
        return self._x
        
    def u(self):
        return self._s

            
class MultiNormalDist(MultivariateDistribution):
    def __init__(self,mean,cov):
        """
        Multivariate normal distribution.  The parameter `mean` is a list of mean
        values for each dimension and `cov` is the variance-covariance matrix.
        """
        self.mean = np.asarray(mean)
        self._cov = np.asarray(cov)
        nd = len(self.mean)
        if cov.shape != (nd,nd):
            ValueError('cov.shape != (len(mean),len(mean))')
        super().__init__(nd)
            
    def _simulate(self,n):
        self.simdata = Distribution.random_state().multivariate_normal(self.mean,self.cov,n).T
        
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
        
        if cov.shape != (nd,nd):
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
            raise ValueError('Two and only two of the parameters may be specified.')
        self._u = self.half_width/np.sqrt(3)
    
    def random(self,n=None):
        return Distribution.random_state().uniform(self.lower_limit,self.upper_limit,n)
    
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
        return Distribution.random_state().gamma(self.shape,self.scale,n)
        
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
        return Distribution.random_state().laplace(self.x(),self.scale,n)
        
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
                raise ValueError('left_width and half_width both specified')
            if right_width is not None:
                raise ValueError('right_width and half_width both specified')
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
                self.lower_limit = mode-left_width
            else:
                if lower_limit is None:
                    raise ValueError('The distribution parameters are under specified.')
                if upper_limit is not None and upper_limit <= lower_limit:
                    raise ValueError('upper_limit <= lower_limit')
                self.left_width = mode-upper_limit
            if right_width is not None:
                if right_width <= 0:
                    raise ValueError('right_width <= 0')
                if upper_limit is not None:
                    raise ValueError('The distribution parameters are over specified.')
                self.upper_limit = mode-right_width
            else:
                if upper_limit is None:
                    raise ValueError('The distribution parameters are under specified.')
                self.right_width = mode-upper_limit
            
    def random(self,n=None):
        return Distribution.random_state().triangular(self.mode,self.lower_limit,self.upper_limit,n)
        
    def x(self):
        return self.mode
        
    def u(self):
        return np.sqrt((self.lower_limit**2+self.upper_limit**2)/6)
        
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
        return Distribution.random_state().exponential(self.scale,n)
        
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
        return Distribution.random_state().poisson(self.lam,n)
        
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
        return Distribution.random_state().binomial(self.n,self.p,n)
        
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
        r1 = Distribution.random_state().uniform(0,1,n)
        r2 = Distribution.random_state().uniform(0,1,n)
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
        r1 = Distribution.random_state().uniform(0,1,n)
        r2 = Distribution.random_state().uniform(0,1,n)
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
            self.half_width*np.sin(2*np.pi*Distribution.random_state().uniform(0,1,n)))
            
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
        return Distribution.random_state().lognormal(self.mu,self.sigma,n)
        
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
        return self.scale*Distribution.random_state().weibull(self.shape,n)
        
    def x(self):
        return self.scale*self._gamma(1 + 1/self.shape)
        
    def u(self):
        return self.scale*np.sqrt(self._gamma(1+2/self.shape) - self._gamma(1+1/self.shape)**2)
        
        
