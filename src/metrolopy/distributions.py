# -*- coding: utf-8 -*-

# module distributions

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

from .exceptions import NoSimulatedDataError
from .dof import DoF

import numpy as np


class Distribution:
    """
    Abstract base class for distributions used for Monte-Carlo uncertainty
    propagation.
    
    In a derived class define the following methods:
        
    random(n=None):  Return a numpy array of n values drawn from the distribution.
        If `n` is `None` then a single scalar value should be returned.  Preferably
        use, as a random number generator, the `numpy.random.Generator` object
        accessed with the `Distribution.random_rng` static method.
        
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
    ...        return Distribution.random_rng().chisquare(self.dof,n)
    ...
    ...    def x(self):
    ...        return self.dof
    ...
    ...    def u(self):
    ...        return 2*self.dof
    
    """
    _random_state = None
    _random_rng = None
    simdata = None
    
    utype = None
    
    p_names = None
    p_names_html = None
    p_names_latex = None
    p_names_ascii = None
    
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
        Returns a legacy `numpy.random.RandomState` instance that may be
        used as a random number generator by derived classes.  Use of the newer
        `numpy.random.Generator` class is now preferred.  An instance of the
        `numpy.random.Generator` class shared by all distributions is returned 
        with the `Distribution.random_rng` static method.
        """
        if Distribution._random_state is None:
            Distribution._random_state = np.random.RandomState()
        return Distribution._random_state
        
    @staticmethod
    def random_rng():
        """
        Returns the `numpy.random.Generator` object shared by all
        distributions.
        """
        if Distribution._random_rng is None:
            Distribution._random_rng = np.random.default_rng()
        return Distribution._random_rng
    
    @staticmethod
    def set_seed(seed):
        """
        Reinitalized the `numpy.random.Generator` object shared by all
        distributions with `seed`.
        """
        Distribution._random_rng = np.random.default_rng(seed=seed)
        Distribution._random_state = np.random.RandomState(seed=seed)
    
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
    def clear_all():
        """
        Clears the Monte-Carlo data from all existing Distribution instances.
        """
        for d in Distribution._called:
            d.clear()
        Distribution._called = []
        
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
        Distribution.clear_all()
        
        if ufrom is not None and any(isinstance(i,Convolution) for i in ufrom):
            raise TypeError('Distributions in ufrom must be independent')
        
        for d in distributions:
            if isinstance(d,Distribution):
                d._simulate(int(n),ufrom)
            
        #for d in Distribution._called:
            #if d not in distributions:
                #d.clear()
                #Distribution._called.remove(d)
                
    def sim(self,n=100000,ufrom=None):
        self.simulate([self],n=n,ufrom=ufrom)
                
    @property
    def isindependent(self):
        return True
         
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
        
    def hist(self,hold=False,xlabel='$ \\mathrm{value} $',ylabel='$ \\mathrm{probability\\:density} $',
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
        if not 'density' in kwds:
            kwds['density'] = True
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
        if isinstance(d,Distribution):
            dt = d.simdata
        if self.simdata is None or dt is None:
            raise NoSimulatedDataError('simulated data does not exist for both distributions')
        cov = np.cov([self.simdata,dt])[1][0]
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
        d = [i.simdata if hasattr(i,'simdata') else i for i in d]
        if any([(v is None) for v in d]):
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
        if ufrom is not None and self not in ufrom and self.utype not in ufrom:
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
    
    def __floordiv__(self,v):
        return Convolution(np.floor_divide,self,v)
        
    def __rfloordiv__(self,v):
        return Convolution(np.floor_divide,v,self)
        
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

        raise NotImplementedError()
        
        
    def cdf(self,*p):
        raise NotImplementedError()
    

def _callarg(a,n,ufrom):
    if isinstance(a,Convolution):
        return a.func(*(_callarg(i,n,ufrom) for i in a.args))
    
    if isinstance(a,Distribution):
        if n is None and a.simdata is None:
            raise NoSimulatedDataError('no simulated data is available')
            
        if ufrom is not None and a not in ufrom and a.utype not in ufrom:
            if n is None:
                n = len(a.simdata)
            return np.full(n,a.x())
        
        if a.simdata is None:
            a._simulate(n,None)
        return a.simdata
    
    return a
        
    
class Convolution(Distribution):
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
        self.args = args
        
        if isinstance(func,np.ufunc):
            self.func = func
        else:
            self.func = np.frompyfunc(func,len(self.args),1)
        
    @property
    def isindependent(self):
        return False
    
    def random(self,n=None):
        raise NotImplementedError('use the simulate static method to generate data from a Convolution')
        
    def x(self):
        args = [a.x() if isinstance(a,Distribution) else a for a in self.args]
        return self.func(*args)
        
    def _random(self,n,ufrom):
        ret = self.func(*(_callarg(i,n,ufrom) for i in self.args))
        
        if ret.dtype != np.float64:
            ret = np.array(ret,dtype=np.float64)
        if n is not None and ret.shape == ():
            ret = np.full(n,ret)
            
        return ret
                
    def _simulate(self,n,ufrom):
        Distribution._called.append(self)
        self.simdata = self._random(n,ufrom)
        
    def datafrom(self,ufrom,save=True):
        """
        Recomputes the convolution with only Distributions in `ufrom` allowed to
        vary.  `sim` or `simulate` must be called to generate
        Monte-Carlo data before calling this method.
        
        Parameters
        ----------
        ufrom:  list containing `Distribution` (not `Convolution`) or `str`
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
        """
        if any(isinstance(i,Convolution) for i in ufrom):
            raise TypeError('Distributions in ufrom must be independent')
            
        ret = self._random(None,ufrom)
        if save:
            self.simdata = ret
        else:
            return ret
        

class MultivariateDistribution:
    def __init__(self,nd):
        """
        Base class for multivarate distributions.  `nd` is the number of dimensions
        of the distribution.
        
        To create a multi-variate distribution, inherit from the 
        MultivariateDistribution and define the following methods:

        _simulate(n):  Return a numpy array of n samples drawn from the distribution. 
            Preferably use, as a random number generator, the numpy.random.Generator 
            object accessed with the Distribution.random_rng() static method.
                
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
                        self.simdata = Distribution.random_rng().dirichlet(self.alpha,n).T

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
        return Distribution.random_rng().normal(self._x,self._u,n)
        
    def x(self):
        return self._x
        
    def u(self):
        return self._u
    
    def cdf(self,z,x,s):
        from scipy.stats import norm
        return norm.cdf(z,loc=x,scale=s)
        
        
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
        if not isinstance(dof,DoF):
            dof = DoF(dof)
        self.DoF = dof
        
    def random(self,n=None):
        if self._s == 0:
            if n == None:
                return self._x
            else:
                return np.full(n,self._x)
        return self._x + self._s*Distribution.random_rng().standard_t(self.DoF.value,n)
        
    def x(self):
        return self._x
        
    def u(self):
        return self._s
    
    @staticmethod
    def cdf(z,x,s,dof):
        from scipy.stats import t as tdist
        return tdist.cdf(z,dof,loc=x,scale=s)
    

class ScipyStatsDist(Distribution):
    def __init__(self,distribution):
        self.distribution = distribution

    def random(self,n=None):
        return self.distribution.rvs(size=n,random_state=self.random_rng())
    
    def x(self):
        return self.distribution.mean()
    
    def u(self):
        return self.distribution.std()
    
    def cdf(self,z,*params):
        return self.distribution.dist.cdf(z,*params)