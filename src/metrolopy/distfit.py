# -*- coding: utf-8 -*-

# module distfit

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

import numpy as np
from inspect import isfunction,ismethod

from .fit import Fit
from .distributions import Distribution,ScipyStatsDist

class DistFit(Fit):
    """
    Fits a one-dimensional continuous distribution.
    Parameters
    ----------
    x: array of `float`
        samples from the distribution to be fitted
        
    cdf: `Distribution`, `scipy.stats.rv_continuous` distribution or 
        function or method
        
        This can be the distribution to be fit as represented by a 
        `gummy.Distribution` (either the class or an instance of the class) or 
        a `scipy.stats.rv_continuous` distribution (either the class or an 
        instance of the class).
        
        It can also be a function or a method that gives the cumlative
        distribution function of the distribtion to be fit.
        
    fix: array of `bool`
        A mask for the fit parameters.  For any element in `fix` that is
        `True`, the corresponding fit parameter will be held constant at
        its initial value.
        
    p0: array of `float`, optional
        Inital values for the fit parameters.
        
        This is required if cdf is a function or method, but optional for most
        `Distribution` or `scipy.stats` distributions.

    Attributes
    ----------
    p:  `numpy.array` of `gummy`
        The fitted values for the fit function parameters as gummys
        including uncertainties and units.
        
    pf:  `numpy.array` of `float`
        The fitted values for the fit function parameters as floats
        
    Properies
    ---------
    dist:  Returns a `gummy.Distribution` instance representing the fitted
    distribution.
        
    Methods
    -------
    plot(...):
        plots the histogram of the data and the fitted distribution
        
    cdf_plot(...):
        plots cumulative distribution function of the fitted distribution
        along with the data points.
    """
    def __init__(self,x,cdf,p0=None,fix=None,rmin=None,p_names=None,**kwds):
        x = np.sort(x)
        
        if isinstance(cdf,Distribution):
            if isinstance(cdf,ScipyStatsDist):
                cdf = cdf.distribution
            else:
                cdf = type(cdf)
            
        if isinstance(cdf,type) and issubclass(cdf,Distribution):
            if cdf.discrete:
                raise TypeError('DistFit does not yet support discrete distributions')
                
            if p_names is None:
                p_names = cdf.param_names()
                
            self.cdf = lambda *args:cdf.cdf(*args)
            self._disttp = cdf
            self._wrap = False
        elif not isfunction(cdf) and not ismethod(cdf):
            from scipy.stats.distributions import rv_frozen,rv_continuous,rv_discrete
            if isinstance(cdf,rv_frozen):
                dist = cdf.dist
                gen = type(cdf.dist)
            elif isinstance(cdf,rv_continuous):
                dist = cdf
                gen = type(cdf)
            elif isinstance(cdf,type) and issubclass(cdf,rv_continuous):
                dist = cdf
                gen = cdf
            else:
                raise TypeError('cdf must be a function, method, or a continuous distrbution')
                
            if isinstance(dist,rv_discrete) or (isinstance(dist,type) and issubclass(dist,rv_discrete)):
                raise TypeError('DistFit does not yet support discrete distributions')
                
            if p_names is None:
                shapes = dist.shapes
                if shapes is not None and shapes.strip() != '':
                    shapes = [s.strip() for s in dist.shapes.split(',')]
                else:
                    shapes = []
                p_names = shapes + ['loc','scale']
                 
            self.cdf = dist.cdf
            self._disttp = gen
            self._wrap = True
        else:
            self.cdf = cdf
            self.p_names = None
            self._disttp = None
            self._wrap = False
            
        self._dist = None
        if rmin is None:
            rmin = 1e-17
        self.rmin = rmin
        self.rmax = -np.log(rmin)
        super().__init__(x,p0=p0,fix=fix,p_names=p_names,**kwds)
        
    @property
    def dist(self):
        """
        Returns a `gummy.Distribution` instance representing the fitted
        distribution.
        """
        if self._dist is None:
            dist = self._disttp(*self.pf)
            if self._wrap:
                dist = ScipyStatsDist(dist)
            self._dist = dist
        return self._dist
    
    def get_p0(self):
        if self._disttp is None:
            raise NotImplementedError()
            
        if self._wrap:
            try:
                return self._disttp().fit(self.xf)
            except:
                try:
                    return self._disttp().fit(self.xf,method='MM')
                except:
                    raise NotImplementedError()
                    
        return self._disttp._est_params(self.xf)
        
    def f(self,x,*p,**kwds):
        n = len(x)
        z = np.empty(n + 1,dtype=float)
        z[:-1] = n*np.diff(self.cdf(x,*p,**kwds),prepend=0)
        z[-1] = n*(1 - self.cdf([x[-1]],*p,**kwds)[0])
        z[z < self.rmin] = self.rmin
        z[z > self.rmax] = self.rmax
        return (1.7325*np.log(z) + z)/3.0662 #np.log(z) + z - 1
    
    def plot(self,ylabel=None,xlabel=None,title=None,hold=False,
             plot_points=None,show_fit=True,show_data=True,xmin=None,xmax=None,fit_options={},
             fig_options={},subplot_options={},hist_options={}):
        """
        plots the histogram of the data and the fitted distribution

        Parameters
        ----------
        ylabel : str, optional
            label for the y-axis
        xlabel : str, optional
            label for the x-axis
        title : str, optional
            plot title
        hold : bool, optional
            If this is True, pyplot.show is called before this method returns.
            The default is False.
        plot_points : int, optional
            The number of points in the fitted curve. The default is 500.
        show_fit : bool, optional
            Whether or not to plot the fitted curve The default is True.
        show_data : bool, optional
            Whether or not to plot a histogra. The default is True.
        xmin : float, optional
            Minimum value for the x-axis. The default is None.
        xmax : float, optional
            Maximum value for the x-axis. The default is None.
        fit_options : dict, optional
            keywords to be passed to the `pyplot.plot` called for the fitted
            curve.
        fig_options : dict, optional
            keywords to be passed when the `pyplot.figure` is created.
        subplot_options : dict, optional
            keywords to be passed when the `pyplot.figure.add_subplot` is called.
        hist_options : dict, optional
            keywords to be passed to `pyplot.hist`

        Returns
        -------
        Figure,Axes
        """
        import matplotlib.pyplot as plt
        
        if plot_points is None:
            plot_points = self.plot_points
            
        if 'ylabel' not in subplot_options and ylabel is not None:
            subplot_options['ylabel'] = ylabel
            
        if 'xlabel' not in subplot_options and xlabel is not None:
            subplot_options['xlabel'] = xlabel
            
        if 'title' not in subplot_options and title is not None:
            subplot_options['title'] = title
            
        fig = plt.figure(**fig_options)
        ax = fig.add_subplot(**subplot_options)
            
        if show_fit:
            if xmin is None:
                p1 = self.xf[0]
            else:
                p1 = xmin
            if xmax is None:
                p2 = self.xf[-1]
            else:
                p2 = xmax
            r = p2 - p1
            if xmin is None:
                st = p1 - self.over_plot*r
            else:
                st = p1
            if xmax is None:
                en = self.xf[-1] + self.over_plot*r
            else:
                en = p2
            fx = np.linspace(st,en,plot_points)
            dx = (en - st)/plot_points
            fy = np.diff(self.cdf(fx,*self.pf))/dx
            fx = fx[1:]
            ax.plot(fx,fy,**fit_options)
        
        if show_data:
            if not 'bins' in hist_options:
                hist_options['bins'] = 100
            if not 'density' in hist_options:
                hist_options['density'] = True
            if not 'histtype' in hist_options:
                hist_options['histtype'] ='stepfilled'
            ax.hist(self.xf,**hist_options)
            
        if not hold:
            plt.show()
            
        return fig,ax
            
    def plot_cdf(self,ylabel=None,xlabel=None,title=None,hold=False,
                 plot_points=None,show_fit=True,show_data=True,xmin=None,
                 xmax=None,fit_options={},fit_format='k-',
                 fig_options={},subplot_options={},data_options={},data_format='ko'):
        """
        plots the cumlative distribution function of the fitted distribution
        along with the data points.

        Parameters
        ----------
        ylabel : str, optional
            label for the y-axis
        xlabel : str, optional
            label for the x-axis
        title : str, optional
            plot title
        hold : bool, optional
            If this is True, pyplot.show is called before this method returns.
            The default is False.
        plot_points : int, optional
            The number of points in the fitted curve. The default is 500.
        show_fit : bool, optional
            Whether or not to plot the fitted curve The default is True.
        show_data : bool, optional
            Whether or not to plot a histogra. The default is True.
        xmin : float, optional
            Minimum value for the x-axis. The default is None.
        xmax : float, optional
            Maximum value for the x-axis. The default is None.
        fit_options : dict, optional
            keywords to be passed to the `pyplot.plot` called for the fitted
            curve.
        fit_format : str, optional
            Format string for the fitted curve. The default is 'k-'.
        fig_options : dict, optional
            keywords to be passed when the `pyplot.figure` is created.
        subplot_options : dict, optional
            keywords to be passed when the `pyplot.figure.add_subplot` is called
        fit_options : dict, optional
            keywords to be passed to the `pyplot.plot` called for the data
            points.
        data_format : str, optional
            Format string for the data points. The default is 'ko'.

        Returns
        -------
        Figure, Axes
        """
        import matplotlib.pyplot as plt
        
        if plot_points is None:
            plot_points = self.plot_points
            
        if 'ylabel' not in subplot_options and ylabel is not None:
            subplot_options['ylabel'] = ylabel
            
        if 'xlabel' not in subplot_options and xlabel is not None:
            subplot_options['xlabel'] = xlabel
            
        if 'title' not in subplot_options and title is not None:
            subplot_options['title'] = title
            
        fig = plt.figure(**fig_options)
        ax = fig.add_subplot(**subplot_options)
            
        if show_fit:
            if xmin is None:
                p1 = self.xf[0]
            else:
                p1 = xmin
            if xmax is None:
                p2 = self.xf[-1]
            else:
                p2 = xmax
            r = p2 - p1
            if xmin is None:
                st = p1 - self.over_plot*r
            else:
                st = p1
            if xmax is None:
                en = self.xf[-1] + self.over_plot*r
            else:
                en = p2
            fx = np.linspace(st,en,plot_points)
            fy = self.cdf(fx,*self.pf)
            if fit_format is None:
                ax.plot(fx,fy,**fit_options)
            else:
                ax.plot(fx,fy,fit_format,**fit_options)
        
        if show_data:
            if 'ms' not in data_options and 'markersize' not in data_options:
                if self.count > 100:
                    data_options['ms'] = 1
                elif self.count > 30:
                    data_options['ms'] = 2
                elif self.count > 20:
                    data_options['ms'] = 3
            if 'ls' not in data_options and 'linestyle' not in data_options:
                data_options['ls'] = 'none'
                    
            nx = len(self.xf)
            if data_format is None:
                ax.plot(self.xf,(np.arange(nx) + 1)/nx,**data_options)
            else:
                ax.plot(self.xf,(np.arange(nx) + 1)/nx,data_format,**data_options)
        
        if not hold:
            plt.show()
            
        return fig,ax
        