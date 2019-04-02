# -*- coding: utf-8 -*-

# module fit

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
This module defines some classes to facilitate curve fitting.
"""

from .ummy import ummy,_der
from .gummy import gummy,_isscalar
from .exceptions import FitWarning
from .unit import Unit,one
from .printing import PrettyPrinter

import numpy as np
from warnings import warn


class _Fit:
    # Base class for fitting.  This implements plotting and unpacks  any
    # uncertainty or unit information in the data.
    
    plot_points = 100
    over_plot = 0.05
    
    xlabel = None
    ylabel = None
    
    def __init__(self,x,y,ux=None,uy=None,xweights=None,yweights=None,
                 sigma_is_known=True,xunit=None,yunit=None):
        self.sigma_is_known = sigma_is_known
            
        if isinstance(x,np.ma.masked_array) or isinstance(y,np.ma.masked_array):
            if not isinstance(x,np.ma.masked_array):
                x = np.ma.masked_array(x)
                
            if x.ndim > 1:
                mask = x[0].mask
                for i in  range(1,x.shape[0]):
                    mask = np.ma.mask_or(mask,x[i].mask)
            else:
                mask = x.mask
                
            if y is not None:
                if not isinstance(y,np.ma.masked_array):
                    y = np.ma.masked_array(y)
                    
                if y.ndim > 1:
                    for i in  range(y.shape[0]):
                        mask = np.ma.mask_or(mask,y[i].mask)
                else:
                    mask = np.ma.mask_or(mask,y.mask)
            
            x.mask = mask
            x = x.compressed()
            
            if y is not None:
                y.mask = mask
                y = y.compressed()
            
        self.x,self.xf,self.ux,self._xunit,self.xdim,self._xgummies = _Fit._getvalues('x',x,ux,xunit)
        if y is not None and not _isscalar(y):
            self.y,self.yf,self.uy,self._yunit,self.ydim,self._ygummies = _Fit._getvalues('y',y,uy,yunit)
            self.implicit = False
            if self.xf.shape[-1] != self.yf.shape[-1]:
                raise TypeError('x and y do not have equal lengths')
        else:
            self.implicit = True
            self.y = y
            self.yf = None
            self.uy = None
            self._yunit = None
            self.ydim = None
            self._ygummies = None
        
        self.count = self.xf.shape[-1]
            
    @staticmethod
    def _getvalues(txt,x,u,units):
        try:
            x = np.asarray_chkfinite(x)
        except ValueError:
            raise ValueError('the '+txt+'-values must not contain infs or NaNs')
        
        if x.ndim > 2 or x.ndim == 0:
            raise TypeError(txt+' must be a 1-d or 2-d array')
            
        count = x.shape[-1]
        
        if u is not None:
            if _isscalar(u):
                u = float(u)
            else:
                u = np.asarray_chkfinite(u,dtype=np.float64)
                if u.shape != x.shape:
                    raise TypeError('u'+txt+' is an array and u'+txt+'.shape != '+txt+'.shape')

        if x.ndim == 1:
            if not isinstance(x[0],ummy):
                if units is None:
                    units = one
                else:
                    units = Unit.unit(units)
                return x,np.asarray_chkfinite(x,dtype=np.float64),u,units,1,False
            x = x.reshape((1,count))
            if u is not None and not _isscalar(u):
                u = u.reshape((1,count))
            dim = 1
        else:
            dim = x.shape[0]
            if not isinstance(x[0][0],ummy):
                if units is None:
                    units = [one]*x.shape[0]
                else:
                    try:
                        if len(units) != x.shape[0]:
                            raise ValueError(txt + 'unit must be a list with length equal to the dimension of the ' + txt + '-parameters')
                    except TypeError:
                        raise ValueError(txt + 'unit must be a list with length equal to the dimension of the ' + txt + '-parameters')
                    units = [Unit.unit(un) for un in units]
                return x,np.asarray_chkfinite(x,dtype=np.float64),u,units,dim,False

        if x[0][0].const:
            gummies = False
        else:
            gummies = True            
            if u is not None:
                raise ValueError('u'+txt+' may not be specified if the '+txt+'-values contain non-constant gummies')
                    
        if units is None:
            if isinstance(x[0][0],gummy):
                try:
                    units = [un.unit for un in (x.T)[0]]
                except:
                    raise TypeError('gummies and non gummies may not be mixed in the '+txt+'-values')
            else:
                units = [one]*x.shape[0]
        else:
            if dim == 1:
                units = [units]
            units = [Unit.unit(un) for un in units]
            
        rxf = np.empty(x.shape,dtype=np.float64)
        if not gummies:
            ru = u
        else:
            ru = np.empty(x.shape,dtype=np.float64)
        
        for i,e in enumerate(x):
            for j,v in enumerate(e):
                if v.const != (not gummies) or not isinstance(v,ummy):
                    raise TypeError('constant and non-constant quantities may be mixed in the '+txt+'-values')
                if isinstance(v,gummy):
                    v = v.convert(units[i])
                rxf[i][j] = v.x
                if gummies:
                    ru[i][j] = v.u
            
        if dim == 1:
            n = rxf.shape[1]
            rxf = rxf.reshape(n)
            x = x.reshape(n)
            if ru is not None and not _isscalar(ru):
                ru = ru.reshape(n)
            units = units[0]
            
        return x,rxf,ru,units,dim,gummies
            
    @property
    def xunit(self):
        """
        Read-only, returns the units on the x data array either as a gummy.Unit 
        object if the x data is 1 dimensional or a list of gummy.Unit objects
        is the x data has more than 1 dimension.
        """
        return self._xunit
        
    @property
    def yunit(self):
        """
        Read-only, returns the units on the y data array either as a gummy.Unit 
        object if the y data is 1 dimensional or a list of gummy.Unit objects
        is the y data has more than 1 dimension.
        """
        return self._yunit
            
    def plot(self,data_format='ko',data_options={},show_data=True,
             error_bars=True,error_bar_k=1,
             fit_format='k-',fit_options={},show_fit=True,
             cik=None,cip=None,ciformat='g-',cioptions={},
             clk=None,clp=None,clformat='r-',cloptions = {},
             xmin=None,xmax=None,xlabel=None,ylabel=None,hold=False,
             plot_points=None):
        """
        Plots the data points, fitted curve, as well as confidence limits and 
        control limits around the fitted curve.
        
        
        Parameters
        ----------
        data_format: `str`, optional
            The format string passed to `pyplot.plot` or `pyplot.errorbar` when
            plotting the data points.  The default is 'ko'.
        data_options: `dict`, optional
            A dictionary containing key words that are passed to `pyplot.plot` or
            `pyplot.errorbar` when plotting the data points.
        show_data: `bool`, optional
            Whether or not to plot the data points. The default is `True`.
        error_bars: `bool`, optional
            Whether or not to plot error bars on the data points (if uncertainty
            values were defined for the data).  The default is `True`.
        error_bar_k: `int` or `float`, optional
            The length of the error bars are determined by multiplying the
            uncertainty for each data point by this quantity. The default value
            is 1.
        fit_format: `str`, optional
            The format string passed to `pyplot.plot` or` pyplot.errorbar` when
            plotting the fitted curve.  The default is 'k-'.
        fit_options: `dict`, optional
            A dictionary containing key words that are passed to `pyplot.plot`
            or `pyplot.errorbar` when plotting the fitted curve.
        show_fit: `bool`, optional
            Whether or not to plot the fitted curve. The default is `True`.
        xmin and xmax: `float`, optional
            The lower and upper limits of the fitted, confidence interval
            and control limit curves.  If this is None, the limits are equal
            to x1 +/- (x2 - x1)*`Fit.over_plot` where x1 is the x value of the
            first data point, x2 is the x value of the last data  point and
            `Fit.over_plot` is an attribute of the `Fit` object with  default
            value 0.05.
        xlabel and ylabel:  `str`, optional
            Labels for the x and y axes. If units are defined for the x or y axes,
            the unit symbol will be added to the end of the labels defined here.
            If these are set to `None`, then the values of the `Fit.xlabel` and
            `Fit.ylabel` attributes will be used.  The default is `None`.
        plot_points:  `int`, optional
            The number of points to use in each curve when plotting the fit,
            confidence interval, and control limit curves.  If this is set to
            `None`, then the value of the `Fit.plot_points` attribute will be used,
            which has a default value of 100.
        hold: `bool`, optional
            If hold is `False` then ``pyplot.show()`` is executed just before this
            function returns.  The detault is `False`
        cik: `float`, `int`, or `None`', optional
            Coverage factor for the uncertainty bands  in the plot.  If `cik` and `cip`
            are `None` then uncertainty bands will not be shown.  Do not specify both
            `cik` and `cip`.
        `cip`: `float`, `int`, or `None`, optional
            Confidence level for the uncertainty bands in the plot.  If `cik`
            and `cip` are None then uncertainty bands will not be shown.  Do
            not specify both `cik and `cip`.
        ciformat: `str`, optional
            Format string passes to the `pyplot.plot` command that plots the
            uncertainty bands. The default is 'g-'.
        cioptions:  `dict`, optional
            Keywork options passed to the `pyplot.plot` command that plots the
            uncertainty bands.
        clk,clp,clformat, cloptions:  optional
            Control limit options, same as above for the uncertainty bands.  The
            control limit band if the control limit k factor multiplied by the
            RSS of the fit uncertainty and the standard deviation of the residuals.
        """
        import matplotlib.pyplot as plt
         
        if self.xdim != 1 or self.ydim != 1:
            raise NotImplementedError('plot is only available for 1-d data')

        if self.implicit:
            raise NotImplementedError('plot is only available for an explict fit')
            
        if xlabel is None:
            xlabel = self.xlabel
            
        if ylabel is None:
            ylabel = self.ylabel
            
        if plot_points is None:
            plot_points = self.plot_points
            
        xf = self.xf
            
        yf = self.yf
        
        if show_data:
            if 'ms' not in data_options and 'markersize' not in data_options:
                if self.count > 100:
                    data_options['ms'] = 1
                elif self.count > 30:
                    data_options['ms'] = 2
                elif self.count > 20:
                    data_options['ms'] = 3
            if error_bars and (self.ux is not None or self.uy is not None):
                ux = self.ux
                if ux is not None:
                    if _isscalar(ux):
                        ux = error_bar_k*ux
                    else:
                        ux = error_bar_k*np.asarray(ux)
                uy = self.uy
                if uy is not None:
                    if _isscalar(uy):
                        uy = error_bar_k*uy
                    else:
                        uy = error_bar_k*np.asarray(uy)
                if data_format is None:
                    data_format = ''
                if 'ls' not in data_options and 'linestyle' not in data_options:
                    data_options['ls'] = 'none'
                plt.errorbar(xf,yf,xerr=ux,yerr=uy,fmt=data_format,
                               **data_options)
            else:
                if data_format is None:
                    plt.plot(xf,yf,**data_options)
                else:
                    plt.plot(xf,yf,data_format,**data_options)            
            
        if show_fit or cik is not None or clk is not None or cip is not None or clp is not None:
            if xmin is None:
                p1 = xf[0]
            else:
                p1 = xmin
            if xmax is None:
                p2 = xf[-1]
            else:
                p2 = xmax
            r = p2 - p1
            if xmin is None:
                st = p1 - self.over_plot*r
            else:
                st = p1
            if xmax is None:
                en = xf[-1] + self.over_plot*r
            else:
                en = p2
            fx = np.linspace(st,en,plot_points)
            
        if not hasattr(self,'dof'):
            dof = float('inf')
        else:
            dof = self.dof
        if cip is not None:
            if cik is not None:
                raise ValueError('cip and cik may not both be specified')
            cik = gummy._p_method.fptok(cip,dof,False)
        if clp is not None:
            if clk is not None:
                raise ValueError('clp and clk may not both be specified')
            clk = gummy._p_method.fptok(clp,dof,False)
        
        if show_fit:
            fy = np.array([self.ypredf(x) for x in fx])
            if fit_format is None:
                plt.plot(fx,fy,**fit_options)
            else:
                plt.plot(fx,fy,fit_format,**fit_options)
        
        if cik is not None:
            u = np.array([self.ypred(x).u for x in fx])
            up = fy + u*cik
            un = fy - u*cik
            if ciformat is None:
                plt.plot(fx,up,**cioptions)
                plt.plot(fx,un,**cioptions)
            else:
                plt.plot(fx,up, ciformat, **cioptions)
                plt.plot(fx,un, ciformat, **cioptions)
        
        if clk is not None:
            u = np.array([self.control_limit(x, clk) for x in fx])
            upl = fy + u
            unl = fy - u
            if clformat is None:
                plt.plot(fx,upl,**cloptions)
                plt.plot(fx,unl,**cloptions)
            else:
                plt.plot(fx,upl,clformat,**cloptions)
                plt.plot(fx,unl,clformat,**cloptions)
            
        xlabel = gummy._plotlabel(xlabel,self.xunit.tostring(fmt='latex'))
        if xlabel is not None:
            plt.xlabel(xlabel)
            
        ylabel = gummy._plotlabel(ylabel,self.yunit.tostring(fmt='latex'))
        if ylabel is not None:
            plt.ylabel(ylabel)
            
        if not hold:
            plt.show()
            
    def control_limit(self,x,k=2):
        if self.sigma is None:
            s = self.s
        else:
            s = self.sigma
        return k*(np.sqrt(self.ypred(x).u + s**2))
        
    
class Fit(_Fit,PrettyPrinter):
    
    latex_math = None
    
    def __init__(self,x,y=None,p0=None,ux=None,uy=None,sigma_is_known=True,xunit=None,
                 yunit=None,solver=None,maxiter=None,nprop=False,**kw):
        """
        Performs a non-linear fit.  The function may be passed in the arguments
        or may be specified by overriding the Fit.f(...) method in a subclass.
        The fitting is performed as soon as the instance is created.
        
        Parameters
        ----------
        x: array_like
           The x-coordinates of the data.  This is a list or numpy array of
           floats or gummys (all point must be of the same type, floats and gummys 
           may not be mixed).  The x-coordinates may be one dimensional or may 
           be multi-dimensional.  For d-dimensional coordinates with (with N total 
           data points) this parameter should be of the form:
               
               [[x1[1], x1[2], ... , x1[N]],
                [x2[1], x2[2], ... , x2[N]],
                .
                .
                .
                [xd[1], xd[2], ... , xd[N]]]
               
            If gummys are given, then the must be dimensionless (unit=one) unless
            the get_puints method is implemented in a subclass.
        y:  array_like, optional
            The y-coordinates of the data (shape and type requirements are the
            same as for the x-coordinates).  This may be omitted only if the
            odr solver is used.
        f:  function
            The fit function.  For d dimensional x-coordinates and k fit parameters
            it should be of the form f(x1,x2,...,xd,p1,p2,...,pk) and return a
            float or (if y is multi-dimensional) a list or array of floats.  This 
            parameter is required unless the f method is overridden in a subclass.
        p0: array_like of `float`
            The initial values for the fit parameters.  This parameter is required
            unless the get_p0 method is overridden in a subclass.
        ux: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `x` values. This should not be specified if the `x`
            argument contains gummys.  If this is specified then only the odr
            solver may be used.  The default is `None`.
        uy: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `y` values. This should not be specified if the y
            argument contains gummys.  The default is `None`.
        sigma_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`. This
            parameter is ignored if `nprop` is `True`.
        xunits, yunits: `str` or `None`, optional
            units for the x and y coordinates. These should not be specified
            if the `x` and `y` parameters contain gummys. These may only be
            specified if the `get_punits` method is overridden in a subclass.
        solver:  {'nls','odr'}, optional
            If this is 'nls' then `scipy.optimize.leastsq` is used to perform
            the fit.  If it is 'odr' then `scipy.odr` is used.  'nls' may not
            be used if the y-coordinate is `None` or multi-dimensional or if
            there is uncertainty in the x-coordinates.  If this is `None`,
            then 'nls' will be used when possible.
        maxiter:  `int` or `None`, optional
            The maximum number of iterations that the solver may use. If this
            is `None` or omitted then the default value for the solver will
            be used.
        nprop:  `bool`, optional
            If this is `True` then uncertainties in the fit will be numerically
             calculated by varying each data point. This will not work if there
             are more than a few data points or if the it is not very stable.
             If this is False than the covariance matrix generated by the solver
             will be used to calculate the uncertainties.  The default value is
             `False`
        other keywords:  optional
            Any additional keyword parameters will be passed to the solver.
        
        Attributes
        ----------
        p:  `list` of `gummy`
            The fitted values for the fit function parameters as gummys
        pf:  `list` of `float`
            The fitted values for the fit function parameters as floats
        res:  `numpy.ndarray` of `float`
            the fit residuals
        s:  `float`
            the standard deviation (or, when there are uncertainties for
            the input data, the square root of the reduced chi-squared) of
            the residuals
        cov:  `numpy.ndarray` of `float`
            the covariance matrix generated by the solver
        fit_output:
            the raw output of the solver
        x:  `numpy.ndarray` of `float` or of `gummy`
            numpy array of the x-coordinates of the data.
        xf:  `numpy.ndarray` of `float`
            numpy array of the x-coordinates of the data as floats
        xdim:  `int`
            the number of dimensions of the x-coordinates
        ux:  `float`, `numpy.ndarray` of `floats` or `None`
            uncertainties in the x-coordinates
        y:  `numpy.ndarray` of `float` or of `gummy`
            numpy array of the y-coordinates of the data.
        yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
        ydim:  `int`
            the number of dimensions of the y-coordinates
        uy:  `float`, `numpy.ndarray` of `floats` or `None`
            uncertainties in the y-coordinates
        count:  `int`
            the number of data points
        p0:  `list` of `float`
            the initial values for the fit function parameters
        solver:  `str`
            the solver used
        punits:  `list` of `Unit`
            the units of the fit parameters
        nparam:  `int`
            the number of fit parameters

        Methods
        -------
        ypred(x1,x2,...):
            Takes `xdim` floats and returns a gummy representing the predicted
            value at that x-coordinate.
        ypredf(x1,x2,...):
            Takes `xdim` floats and returns a float giving the  predicted value
            at that x-coordinate.
        plot(...):
            plots the data (only available if x and y are one-dimensional)

        Notes
        -----
        When subclassing Fit you may override the following methods:
            
        f(self,x1,x2,...,xd,p1,p2,...,pk):
            The fit function. (the signature shown here assumes that x is
            d-dimentional and there are k fit parameters), It should return a
             float or (if y is multi-dimensional) a list of floats.

        jac(self,x1,x2,...,xk,p1,p2,...,pk):
            The Jacobian.  This may optionally be overridden in a derived
            class.  If this method throws a `NotImplementedError` the
            derivatives will be calculated numerically.
        
            It must have the same signature as the f method and return a list
            of derivatives of the form:
            
            [df/dx1,df/dx2,...,df/dp1,df/dp2,...] 
        
            if f returns a scalar or:
        
            [[df1/dx1,df1/dx2,...,df1/dp1,df1/dp2,...],
            [df2/dx1,df2/dx2,...,df2/dp1,df2/dp2,...],...]
        
            if f returns a 1-d array [f1,f2,...].

        get_p0(self):
            Returns an initial guess for the fit parameters based on the input x
            and y data.  This is not required, but if it is not implemented then
            the `p0` parameter is a required parameter for the `__init__` method.

        get_punits(self):
            Returns a list units for the fit parameters.  This is not required,
            but if it is not implemented then only float values or dimensionless
            gummys may be as the `x `and `y` parameters and the `xunit` and `yunit`
            parameters to the `__init__` method may not be used.

        funicode, flatex, fhtml:
            Returns a `str` containing unicode, latex, and html representations of
            the fit function.
        """
        super().__init__(x,y,ux=ux,uy=uy,sigma_is_known=sigma_is_known,xunit=xunit,yunit=yunit)
        
        self.nprop = nprop
        
        if type(self) is Fit:
            if kw.get('f') is not None: 
                self.f = kw['f']
                del kw['f']
            
        if solver is not None:
            solver = solver.strip().lower()
        self.solver = solver
        self.solver,self._solver = self._get_solver(solver,**kw)
        
        if p0 is None:
            try:
                self.p0 = self.get_p0()
            except NotImplementedError:
                raise ValueError('p0 must be specified')
        else:
            self.p0 = p0
            
        getpu = False
        if self.xdim == 1:
            getpu = getpu or (self.xunit is not one)
        else:
            for u in self.xunit:
                getpu = getpu or (u is not one)
        if self.ydim == 1:
            getpu = getpu or (self.yunit is not None and self.yunit is not one)
        else:
            for u in self.yunit:
                getpu = getpu or (u is not None and u is not one)
        if getpu:
            try:
                self.punits = self.get_punits()
            except NotImplementedError:
                raise ValueError('only dimensionless quanities are accepted')
        else:
            self.punits = [one]*len(self.p0)
            
        self.maxiter = maxiter
            
        self.nparam = len(self.p0)
        
        try:
            # See if f will broadcast properly across the xf array...
            if self.xdim == 1:
                r = self.f(self.xf,*self.p0)
            else:
                r = self.f(*(list(self.xf)+self.p0))
            if self.ydim == 1:
                if not r.shape != (self.count,):
                    raise TypeError()
            else:
                if len(r[0]) != self.count:
                    raise TypeError()
        except NotImplementedError:
            pass
        except:
            # ...if not vectorize it so that it does.
            ydim = self.ydim
            if ydim is None:
                ydim = 1
            ot = [np.float64]*ydim
            self.f = np.vectorize(self.f,ot,excluded=list(range(self.xdim,self.xdim+self.nparam)))
        
        self._solver(**kw)
        
        if self.nprop and (self._xgummies or self._ygummies):
            self._solver(**kw)
            p0 = self.p0
            self.p0 = self.pf
            xf = self.xf
            yf = self.yf
            res = self.res
            sigma = self.sigma
            s = self.s
            cov = self.cov
            fit_output = self.fit_output
            
            def f(*a):
                self.xf = np.array(a[:self.count*self.xdim],dtype=np.float64)
                if self.xdim > 1:
                    self.xf = self.xf.reshape((self.count,self.xdim))
                self.yf = np.array(a[self.count*self.xdim:],dtype=np.float64)
                if self.ydim > 1:
                    self.yf = self.yf.reshape((self.count,self.ydim))
                self._solver(**kw)
                return self.pf
            self.p = gummy._napply(f,*list(np.concatenate((self.x.flatten(),self.y.flatten()))),
                                   fxx=(self.pf,self.xf))
            self.p0 = p0
            self.xf = xf
            self.yf = yf
            self.res = res
            self.sigma = sigma
            self.s = s
            self.cov = cov
            self.fit_output = fit_output
        else:
            for p in self.p:
                if not p.unit.linear:
                    p._set_U()
            
    def _get_solver(self,solver,**kw):
        if self.ux is not None or self.implicit or self.ydim > 1:
            if solver is not None and solver != 'odr':
                raise ValueError('if there are uncertainties on the x-values, ydim >1, or the model is impicit\nthen only the odr solver can be used')
            return 'odr',self._odr
            
        if solver == 'nls' or solver is None:
            return 'nls',self._leastsq
        elif solver == 'odr':
            return 'odr',self._odr
        elif solver == 'gls':
            raise ValueError('solver gls is only available with PolyFit')
        else:
            raise ValueError('solver ' + str(solver) + ' is not recognized')
            
    def _get_corr(self,**kw):
        from scipy.linalg import inv,cholesky
        if 'ignore_correlations' in kw:
            ignore_corr = kw['ignore_correlations']
        else:
            ignore_corr = False
            
        v = None
        ic = None
        w = None
        vsc = 1
        
        if 'ycov' in kw:
            v = np.asarray(kw['ycov'])
            if v.shape != (self.count,self.count):
                raise ValueError('ycov mustbe a len(x) by len(x) matrix')
            self.uy = np.diagonal(v)
            if ignore_corr:
                v = None
        elif self._ygummies and not ignore_corr:
            v = gummy.covariance_matrix(self.y)
            
        if v is not None and np.count_nonzero(v - np.diag(np.diagonal(v))) == 0:
            v = None
            
        if v is not None:
            try:
                vsc = np.diagonal(v).sum()/self.count
                v /= vsc
                ic = cholesky(inv(v))
            except:
                ignore_corr = True
                warn('unable to invert the covariance matrix for y;\ncorrelations will be ignored and correlations between the fit parameters and y will not be set',FitWarning)
                v = None
                
        if v is None and self.uy is not None and not _isscalar(self.uy):
            vsc = np.sum(self.uy**2)/self.count
            w = 1/self.uy
            w *=  np.sqrt(vsc)
            
        return ic,w,vsc,ignore_corr
            

    def _leastsq(self,**kw):
        from scipy.optimize import leastsq
        from scipy.linalg import inv
        
        ic,w,vsc,ignore_corr = self._get_corr(**kw)
        # If there are correlations between the y values and ignore_corr is False,
        # then ic is the cholesky decomposition of the inverse of the covariance 
        # matrix of self.y. ignore_correlations may be True
        # if it was passed in kw or if there was an error when calculating ic.
        # If there are no correlations, but self.uy is an array,then w is an array 
        # of weights for each data point.  vsc is a scalar normalization factor
        # for ic or w.
        
        args = (self.xf,self.yf,self.f)
        if ic is not None:
            if self.xdim == 1:
                def func(params,x,y,f):
                    return ic@(f(x,*params) - y)      
            else:
                def func(params,x,y,f):
                    return ic@(f(*(list(x)+list(params))) - y)
                
            try:
                if self.xdim == 1:
                    self.jac(self.xf,*self.p0)
                    def dfun(*a):
                        return self.jac(a[1],*(a[0]))[1:]@ic.T
                else:
                    self.jac(*(list(self.xf)+list(self.p0)))
                    def dfun(*a):
                        return self.jac(*(list(a[1])+list(a[0])))[self.xdim:]@ic.T
            except NotImplementedError:
                dfun = None
        elif w is not None:
            if self.xdim == 1:
                def func(params,x,y,f):
                    return w*(f(x,*params) - y)    
            else:
                def func(params,x,y,f):
                    return w*(f(*(list(x)+list(params))) - y)
                
            try:
                if self.xdim == 1:
                    self.jac(self.xf,*self.p0)
                    def dfun(*a):
                        return w*self.jac(a[1],*(a[0]))[1:]
                else:
                    self.jac(*(list(self.xf)+list(self.p0)))
                    def dfun(*a):
                        return w*self.jac(*(list(a[1])+list(a[0])))[self.xdim:]
            except NotImplementedError:
                dfun = None      
        else:
            if self.xdim == 1:
                def func(params,x,y,f):
                    return f(x,*params) - y
            else:
                def func(params,x,y,f):
                    return f(*(list(x)+list(params))) - y
            try:
                if self.xdim == 1:
                    self.jac(self.xf,*self.p0)
                    def dfun(*a):
                        return self.jac(a[1],*(a[0]))[1:]
                else:
                    self.jac(*(list(self.xf)+self.p0))
                    def dfun(*a):
                        return self.jac(*(list(a[1])+list(a[0])))[self.xdim:]
            except NotImplementedError:
                dfun = None
        
        kw['full_output'] = 1
        if self.maxiter is not None:
            kw['maxfev'] = self.maxiter
            
        self.fit_output = leastsq(func,self.p0,Dfun=dfun,args=args,col_deriv=1,**kw)
        self.pf,cov,infodict,errmsg,ier = self.fit_output
             
        if ier == 5:
            raise RuntimeError('the fit did not converge: max iterations reached (maxiter = '+str(infodict['nfev'])+')\nyou may wish to try again either increasing maxiter or with a different p0')
        elif ier not in [1, 2, 3, 4]:
            raise RuntimeError('the fit did not converge: ' + errmsg)
        
        if self.xdim == 1:
            self.res = self.yf - self.ypredf(self.xf)
        else:
            self.res = self.yf - self.ypredf(*self.xf)
        
        if cov is not None:
            cov *= vsc

        if ic is not None:
            sigmasq = np.sum((ic@self.uy)**2)/self.count
            self.sigma = np.sqrt(sigmasq)
            var = np.sum((ic@self.res)**2)/(self.count - self.nparam)
            if not self.sigma_is_known:
                cov *= var/sigmasq
        elif w is not None:
            sigmasq = ((w*self.uy)**2).sum()/self.count
            self.sigma = np.sqrt(sigmasq)
            var = ((w*self.res)**2).sum()/(self.count - self.nparam)
            if not self.sigma_is_known:
                cov *= var/sigmasq
        elif self.uy is not None:
            if _isscalar(self.uy):
                self.sigma = self.uy
                sigmasq = self.uy**2
            else:
                sigmasq = np.sum(self.uy**2)/self.count
                self.sigma = np.sqrt(sigmasq)
            var = (self.res**2).sum()/(self.count - self.nparam)
            if cov is not None:
                if self.sigma_is_known:
                    cov *= sigmasq
                else:
                    cov *= var
        else:
            sigmasq = None
            self.sigma = None
            var = (self.res**2).sum()/(self.count-self.nparam)
            if cov is not None:
                cov *= var
                
        self.p = None
        self.dof = None
        if self._ygummies and not ignore_corr and self.sigma_is_known:
            # set the correlations between self.p and self.y
            try:
                if dfun is None:
                    gf = [ummy(p,np.sqrt(u)) for p,u in zip(self.pf,np.diagonal(cov))]
                    if self.xdim == 1:
                        jac = np.array([_der(self.f,x,*gf)[1:] for x in self.x],dtype=float)
                    else:
                        jac = np.array([_der(self.f,*(list(x)+list(gf)))[self.xdim:] for x in self.x.T],dtype=float)
                else:
                    if self.xdim == 1:
                        jac = np.array([self.jac(x,*self.pf)[1:] for x in self.x],dtype=float)
                    else:
                        jac = np.array([self.jac(*(list(x)+list(self.pf)))[self.xdim:] for x in self.x.T],dtype=float)
                sc = np.sqrt((jac*jac).sum(axis=0))
                jac /= sc
                y = [g.graft(one) for g in self.y]
                if ic is not None:
                    jac = ic@jac
                    y = [np.sum(y*c) for c in ic]
                elif w is not None:
                    jac *= w[:,np.newaxis]
                    y = y*w
                else:
                    y = y
                pcov = inv(jac.T@jac)
                if cov is None:
                    cov = pcov*vsc
                    cov /= np.outer(sc,sc)
                b = pcov@jac.T
                self.p = [np.sum(y*bi) for bi in b]/sc
                for g,gf in zip(self.p,self.pf):
                    g._x = gf
            except:
                self.p = None
                if cov is not None:
                    warn('unable to set the correlations between the y-values and the fit parameters',FitWarning)
                
        self.s = np.sqrt(var)
            
        if self.dof is None:
            self.dof = self.count - len(self.pf)
    
        self.cov = cov
        
        if self.p is None:
            if cov is not None and self.count > self.nparam:
                if self.sigma is None or not self.sigma_is_known:
                    dof = self.dof
                    if dof < 1:
                        dof = 1
                else:
                    dof = float('inf')
                
                try:
                    self.p = gummy.create(self.pf,unit=self.punits,dof=dof,
                                          covariance_matrix = self.cov)
                except ValueError:
                    self.p = gummy.create(self.pf,unit=self.punits)
                    warn('the covariance matrix is not positive semidefinate; uncertainties cannot be calculated',FitWarning)
                
            else:
                self.cov = None
                self.p = [gummy(v,unit=self.punits[i]) for i,v in enumerate(self.pf)]
                warn('unable get a covariance matrix of full rank; uncertainties cannot be calculated',FitWarning)
        else:
            self.p = [g.graft(u) for g,u in zip(self.p,self.punits)]
            
    @staticmethod
    def _getw(x,u,cov,w,dim,gummies):
        from scipy.linalg import inv
        
        if w is not None:
            return np.asarray(w)
            
        if gummies and dim > 1:
            return [inv(gummy.covariance_matrix(v)) for v in (x.T)]
            
        if cov is not None:
            cov = np.asarray(cov)
            if cov.ndim < 2 or cov.ndim > 3:
                raise TypeError('covx.ndim must be 2 or 3')
            if cov.ndim == 2:
                return inv(cov)
            else:
                return [inv(c) for c in cov]
                
        if u is not None:
            u = np.asarray(u)
            return 1/u**2
            
        return None
        
    def _odr(self,**kw):
        from scipy import odr
        
        covx = kw.get('covx')
        covy = kw.get('covy')
        wd = kw.get('wd')
        we = kw.get('we')
        
        if (self.ux is not None) + (covx is not None) + (wd is not None) > 1:
            raise ValueError('only one of (ux, covx, and odrwd) may be specified')
            
        if (self.uy is not None) + (covy is not None) + (we is not None) > 1:
            raise ValueError('only one of (uy, covy, and odrwe) may be specified')
            
        wd = self._getw(self.x,self.ux,covx,wd,self.xdim,self._xgummies)
        we = self._getw(self.y,self.uy,covy,we,self.ydim,self._ygummies)
                                     
        if 'fix' in kw:
            del kw['fix']
        
        try:  
            if self.xdim == 1:
                self.jac(self.xf,*self.p0)
                def fjacb(p,x):
                    return np.array(self.jac(x,*p)[1:])
                def fjacd(p,x):
                    return np.array(self.jac(x,*p)[0])
            else:
                self.jac(*(list(self.xf)+self.p0))
                def fjacb(p,x):
                    return np.array(self.jac(*(list(x)+list(p)))[self.xdim:])
                def fjacd(p,x):
                    return np.array(self.jac(*(list(x)+list(p)))[:self.xdim])
                
        except NotImplementedError:
            fjacb = None
            fjacd = None
            
        if self.xdim == 1:
            if self.ydim == 1:
                def f(p,x):
                    return self.f(x,*p)
            else:
                def f(p,x):
                    return np.asarray(self.f(*x,*p))
        else:
            if self.implicit or self.ydim == 1:
                def f(p,x):
                    return self.f(*(list(x)+list(p)))
            else:
                def f(p,x):
                    return np.asarray(self.f(*(list(x)+list(p))))
                
        if self.implicit:
            if self.y is None:
                yf = f(self.p0,self.xf)
                if yf.ndim == 1:
                    yf = 1
                else:
                    yf = yf.shape(0)
            else:
                yf = self.y
                if self.ydim == 1:
                    yf = np.array([yf])
        else:
            yf = self.yf
        self.odr_data = odr.Data(self.xf,y=yf,wd=wd,we=we,
                                     fix=kw.get('fix'))

                
        mkw = {}
        if 'estimate' in kw:
            mkw['estimate'] = kw['estimate']
            del kw['estimate']
        if 'extra_args' in kw:
            mkw['extra_args'] = kw['extra_args']
            del kw['extra_args']
            
        self.odr_model = odr.Model(f,fjacb,fjacd,implicit=self.implicit,**mkw)
        
        if self.maxiter is not None:
            kw['maxit'] = self.maxiter
        
        if 'fit_type' in kw:
            fit_type = kw['fit_type']
            del kw['fit_type']
            if isinstance(fit_type,str):
                fit_type = fit_type.strip().lower()
                if fit_type == 'leastsq':
                    fit_type = 2
                elif fit_type == 'odr':
                    fit_type = 0
                elif fit_type == 'implicit':
                    fit_type = 1
                    
            if self.implicit:
                if fit_type != 1:
                    raise ValueError('fit_type must be 1 (implicit) if there is no y-value data')
            else:
                if fit_type == 1:
                    raise ValueError('fit_type 1 (implicit) is not allowed if there is y-value data')
            
        else:
            fit_type = None
            
        if 'deriv' in kw:
            deriv = kw['deriv']
            del kw['deriv']
        else:
            deriv = None
            
        if 'var_calc' in kw:
            var_calc = kw['var_calc']
            del kw['var_calc']
        else:
            var_calc = None
            
        if 'del_init' in kw:
            del_init = 'del_init'
            del kw['del_init']
        else:
            del_init = None

        self.odr = odr.ODR(self.odr_data,self.odr_model,beta0=self.p0,**kw)
        
        if fit_type is not None or deriv is not None or var_calc is not None or del_init is not None:
            self.odr.set_job(fit_type=fit_type,deriv=deriv,var_calc=var_calc,
                             del_init=del_init)
        
        self.fit_output = self.odr.run()
        
        if self.fit_output.info > 5:
            raise odr.odr_error('Error encountered during fit:\n' + str(self.fit_output.stopreason) + '\nODR info  = ' + str(self.fit_output.info))
        if self.fit_output.info == 4:
            warn('the iteration limit was reached')

        self.pf = self.fit_output.beta

        if self.implicit:
            self.res = f(self.pf,self.xf)
        else:
            if self.xdim == 1:
                self.res = self.yf - self.ypredf(self.xf)
            else:
                self.res = self.yf - self.ypredf(*self.xf)
        
        var = self.fit_output.res_var
        
        self.s = np.sqrt(var)
        
        if self.uy is None and self.ux is not None:
            u = self.ux
        elif self.uy is not None and self.ux is None:
            u = self.uy
        elif self.uy is not None and self.ux is not None:
            u = np.sqrt(self.ux**2 + self.uy**2)
        else:
            u = None
            
        if u is not None:
            if _isscalar(u):
                self.sigma = u
            else:
                if self.ydim == 1:
                    self.sigma = np.sqrt((u**2).sum()/self.count)
                else:
                    self.sigma = 1
        else:
            self.sigma = None
       
        self.dof = self.count - len(self.pf)
        
        if self.sigma is None or not self.sigma_is_known:
            self.cov = var*self.fit_output.cov_beta
        elif self.sigma_is_known:
            self.cov = (self.sigma**2)*self.fit_output.cov_beta
        
        try:
            if self.sigma is None or not self.sigma_is_known:
                dof = self.dof
            else:
                dof = float('inf')
            self.p = gummy.create(self.pf,unit=self.punits,dof=dof,
                                  covariance_matrix = self.cov)
        except ValueError:
            self.p = gummy.create(self.pf,unit=self.punits)
            warn('the covariance matrix found by odr is not positive semidefinate; uncertainties cannot be calculated',FitWarning)
        
    def ypredf(self,*x):
        """
        returns a float representing the value predicted by the fit at x
        """
        if self.xdim == 1:
            x = x[0]
            if isinstance(x,gummy):
                x = x.convert(self._xunit)
                x = x.x
            if self.ydim == 1:
                return self.f(x,*self.pf)
            else:
                return np.asarray(self.f(x,*self.pf))
        if len(x) != self.xdim:
            raise TypeError('the length of the x parameter does not match the dimensinality of the x-value array')
        a = list(x)+list(self.pf)
        if self.ydim == 1:
            ret = self.f(*a)
        else:
            ret = np.asarray(self.f(*a))
        
        if isinstance(ret,np.ndarray) and ret.shape == ():
            return ret.item()
        else:
            return ret
            
    def ypred(self,*x):
        """
        returns a gummy representing the value predicted by the fit at x
        """
        func = np.frompyfunc(self._ypred,self.xdim,self.ydim)
        ret = func(*x)
        if hasattr(ret,'shape') and ret.shape == ():
            return ret.item()
        return ret
        
    def _ypred(self,*x):
        p = list(self.p)
        for i,e in enumerate(p):
            p[i] = p[i].graft(one)
        if self.xdim == 1:
            x = x[0]
            if isinstance(x,gummy):
                x = x.convert(self._xunit)
                x = x.graft(one)
            a = [x] + p
        else:
            if len(x) != self.xdim:
                raise TypeError('the length of the x parameter does not match the dimensionality of the x-value array')
            x = list(x)
            for i,v in enumerate(x):
                if isinstance(v,gummy):
                    x[i] = v.convert(self._xunit[i])
                    x[i] = x.graft(one)
            a = x + p
            
        if self.ydim == 1:
            try:
                r = gummy.apply(self.f,self.jac,*a)
            except NotImplementedError:
                r = gummy.napply(self.f,*a)
            r = r.graft(self._yunit)
            return r
        else:
            r = np.zeros(self.ydim,dtype=np.object)
            try:
                for i in self.ydim:
                    def f(*v):
                        return self.f(*v)[i]
                    r[i] = gummy.apply(f,self.jac,*a)
                    r[i] = r[i].graft(self._yunit)
            except NotImplementedError:
                for i in self.ydim:
                    def f(*v):
                        return self.f(*v)[i]
                    r[i] = gummy.napply(f,*a)
                    r[i] = r[i].graft(self._yunit)
            return r
                
    def ptostring(self,fmt='unicode'):
        """
        Returns a string that displays the fit parameters.
        
        Parameters
        ----------
        fmt:  {'unicode', 'ascii', 'latex' or 'html'}
            format for the output
        """
        fmt = fmt.strip().lower()
        
        if fmt == 'latex':
            txt = '\\begin{align}'
        else:
            txt = ''
        for i,p in enumerate(self.p):
            if fmt == 'unicode' or fmt == 'utf-8':
                txt += 'p('+str(i+1)+') = ' + p.tostring(fmt='unicode') + '\n'
            elif fmt == 'ascii' or fmt == 'utf-8':
                txt += 'p('+str(i+1)+') = ' + p.tostring(fmt='ascii') + '\n'
            elif fmt == 'latex':
                txt += 'p_{'+str(i+1)+'} &= ' + p.tostring(fmt='latex').strip('$') + '\\\\'
            elif fmt == 'html':
                txt += '<i>p</i><sub>'+str(i+1)+'</sub> = ' + p.tostring(fmt='html') + '<br>'
            else:
                raise ValueError('format ' + str(fmt) + ' is not understood')
                
        if fmt == 'latex':
            txt += '\\end{align}'
        
        return txt
                
    def tostring(self,fmt='unicode'):
        """
        Returns a string to display the fit function and fit parameters.
        
        Parameters
        ----------
        fmt:  {'unicode', 'ascii', 'latex' or 'html'}
            format for the output
        """
        fmt = fmt.strip().lower()

        txt = ''
        try:
            if fmt == 'unicode' or fmt == 'ascii':
                txt += self.funicode() + '\n\n'
            elif fmt == 'latex':
                txt += self.flatex().strip('$') + '\\\\[10pt]'
            elif fmt == 'html':
                txt += self.fhtml() + '<br><br>'
            else:
                raise ValueError('format ' + str(fmt) + ' is not understood')
        except NotImplementedError:
            pass

        if fmt == 'unicode' or fmt == 'ascii':
            txt += 'best fit parameters:\n'
        elif fmt == 'latex':
            txt += type(self).latex_norm('best fit parameters:') + '\\\\'
        else:
            txt += 'best fit parameters:<br>'
        txt += self.ptostring(fmt)
        
        return txt
        
        
    def tohtml(self):
        """
        Returns a string containing an html representation of the fit function
        and the best fit parameters.
        
        Equivalent to ``Fit.tostring('html')``.
        """
        return super().tohtml()
        
    def tolatex(self):
        """
        Returns a string containing a latex representation of the fit function
        and the best fit parameters.
        
        Equivalent to ``Fit.tostring('latex')``.
        """
        return super().tolatex()
        
    def toascii(self):
        """
        Returns a string containing an ascii representation of the fit function
        and the best fit parameters.
        
        Equivalent to ``Fit.tostring('ascii')``.
        """
        return super().toascii()
        
    def latex(self):
        """
        If called from an IPython or Jupyter console, prints the fit function
        and best fit parameters with the latex interpreter.
        """
        return super().latex()
    
    def html(self):
        """
        If called from an IPython or Jupyter console, prints the fit function
        and best fit parameters with the html interpreter.
        """
        return super().html()
        
    def f(self,*a):
        """
        Not implemented, implemented this in a derived class.
        
        The function to fit.

        It must either have signature:
        
        f(self,x,p1,p2,...,pn) 
        
        where there are p1 to pn are the n fit parameters and the independent
        variable x has one dimension, or:
        
        f(self,x1,x2,...,xm,p1,p2,...,pn)
        
        where the independent variable x has m dimensions at each observation.

        f should return either a float or a 1-d array of floats depending on the
        dimension of the response variable.
        """
        raise NotImplementedError()
    
    def jac(self,*a):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        The Jacobian of the fit function.
        
        If this method throws a `NotImplementedError` the derivatives will be
        calculated numerically.
        
        It must have the same signature as the f method and return a list
        of derivatives of the form:
            
        [df/dx1,df/dx2,...,df/dp1,df/dp2,...] 
        
        if f returns a scalar or:
        
        [[df1/dx1,df1/dx2,...,df1/dp1,df1/dp2,...],
         [df2/dx1,df2/dx2,...,df2/dp1,df2/dp2,...],...]
        
        if f returns a 1-d array [f1,f2,...].
        """
        raise NotImplementedError()
        
    def get_p0(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns an initial guess for the the fit parameters [p1,p2,...] based
        on the input x and y data.  
        
        If this method is not implemented then the inital values must be passed 
        in the p0 parameter when the instance is created.
        """
        raise NotImplementedError()
        
    def get_punits(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns a list of units for each fit parameter [p1,p2,...] based on the
        units of the input data.
        
        If this is not implemented then only dimensionless data (with unit one)
        may be fit.
        """
        raise NotImplementedError()
        
    def funicode(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns a string containing a unicode representation of the fit function.
        """
        raise NotImplementedError()
    
    def flatex(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns a string containing a latex representation of the fit function.
        """
        raise NotImplementedError()
    
    def fhtml(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns a string containing an html representation of the fit function.
        """
        raise NotImplementedError()
        
class PolyFit(Fit):
    def __init__(self,x,y,deg=1,ux=None,uy=None,sigma_is_known=True,p0=None,
                 xunit=None,yunit=None,solver=None,maxiter=None,
                 nprop=False,**kw):
        """
        Fits the x,y data to a polynomial

        Parameters
        ----------
        x: array_like
           The x-coordinates of the data.  This is a list or numpy array of
           floats or gummys (all point must be of the same type, floats and gummys
           may not be mixed).  The x-coordinates may be one dimensional or may
           be multi-dimensional.  For d-dimensional coordinates with (with N total
           data points) this parameter should be of the form:

               [[x1[1], x1[2], ... , x1[N]],
                [x2[1], x2[2], ... , x2[N]],
                .
                .
                .
                [xd[1], xd[2], ... , xd[N]]]

            If gummys are given, then the must be dimensionless (unit=one) unless
            the get_puints method is implemented in a subclass.
        y:  array_like
            The y-coordinates of the data (shape and type requirements are the
            same as for the x-coordinates).
        deg:  `int`
            the degree of the polynomial
        p0: array_like of `float`, optional
            The initial values for the fit parameters.  These are only used if solver
            is 'odr' or 'nls'.
        ux: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `x` values. This should not be specified if the `x`
            argument contains gummys.  If this is specified then only the odr
            solver may be used.  The default is `None`.
        uy: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `y` values. This should not be specified if the y
            argument contains gummys.  The default is `None`.
        sigma_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`. This
            parameter is ignored if `nprop` is `True`.
        xunits, yunits: `str` or `None`, optional
            units for the x and y coordinates. These should not be specified
            if the `x` and `y` parameters contain gummys. These may only be
            specified if the `get_punits` method is overridden in a subclass.
        solver:  {'ols','nls','odr', `None`}, optional
            By default a linear fit algorithm, ols, will be used if `x` and `y` are
            one dimensional and there is no uncertainty in the x-values.  The odr
            solver must be used if there is uncertainty in the x-values or if the
            y-coordinates are multi-dimensional.  By default the nonlinear least
            squares solver, 'nls' will be used if `x` is multi-dimensional.
        maxiter:  `int` or `None`, optional
            The maximum number of iterations that the solver may use. If this
            is `None` or omitted then the default value for the solver will
            be used.
        nprop:  `bool`, optional
            If this is `True` then uncertainties in the fit will be numerically
             calculated by varying each data point. This will not work if there
             are more than a few data points or if the it is not very stable.
             If this is False than the covariance matrix generated by the solver
             will be used to calculate the uncertainties.  The default value is
             `False`
        other keywords:  optional
            Any additional keyword parameters will be passed to the solver.

        Attributes
        ----------
        p:  `list` of `gummy`
            The fitted values for the fit function parameters as gummys
        pf:  `list` of `float`
            The fitted values for the fit function parameters as floats
        res:  `numpy.ndarray` of `float`
            the fit residuals
        s:  `float`
            the standard deviation (or, when there are uncertainties for
            the input data, the square root of the reduced chi-squared) of
            the residuals
        cov:  `numpy.ndarray` of `float`
            the covariance matrix generated by the solver
        fit_output:
            the raw output of the solver
        x:  `numpy.ndarray` of `float` or of `gummy`
            numpy array of the x-coordinates of the data.
        xf:  `numpy.ndarray` of `float`
            numpy array of the x-coordinates of the data as floats
        xdim:  `int`
            the number of dimensions of the x-coordinates
        ux:  `float`, `numpy.ndarray` of `floats` or `None`
            uncertainties in the x-coordinates
        y:  `numpy.ndarray` of `float` or of `gummy`
            numpy array of the y-coordinates of the data.
        yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
        ydim:  `int`
            the number of dimensions of the y-coordinates
        uy:  `float`, `numpy.ndarray` of `floats` or `None`
            uncertainties in the y-coordinates
        count:  `int`
            the number of data points
        p0:  `list` of `float`
            the initial values for the fit function parameters
        solver:  `str`
            the solver used
        punits:  `list` of `Unit`
            the units of the fit parameters
        nparam:  `int`
            the number of fit parameters
        deg:  `int`
            the degree of the polynomial

        Methods
        -------
        ypred(x1,x2,...):
            Takes xdim floats and returns a gummy representing the predicted
            value at that x-coordinate.
        ypredf(x1,x2,...):
            Takes xdim floats and returns a float giving the  predicted value
            at that x-coordinate.
        plot(...):
            plots the data (only available if x and y are one-dimensional)
        """
        
        x = np.asarray(x)
        if _isscalar(deg):
            if x.ndim == 2:
                raise TypeError('if x is muilt-dimensional, then deg must be an list, tuple or array of length xdim')
            self.deg = int(deg)
            self._order = deg+1
        else:
            if x.ndim == 1:
                raise TypeError('if x is 1-d, then deg must be a scalar')
            try:
                self.deg = np.array(deg)
                for i,d in enumerate(self.deg):
                    self.deg[i] = int(d)
                self._order = self.deg + 1
            except:
                raise TypeError('deg must be a scalar value or a 1-d array')

        super().__init__(x,y=y,p0=p0,ux=ux,uy=uy,sigma_is_known=sigma_is_known,xunit=xunit,
                 yunit=yunit,solver=solver,maxiter=maxiter,
                 nprop=nprop,**kw)

    def _get_solver(self,solver,**kw):
        if self.xdim > 1:
            if _isscalar(self.deg) or len(self.deg) != self.xdim:
                raise TypeError('len(deg) != xdim')
                
        if self.ux is not None or self.ydim > 1 or self.implicit:
            if solver is not None and solver != 'odr':
                raise ValueError('if there are uncertainties on the x-values, ydim > 1, or the model is impicit\nthen only the odr solver can be used')
            return 'odr',self._odr
        
        if solver == 'ols': solver = 'gls'
        
        if solver in ('ols','gls') and self.xdim > 1:
            raise ValueError('solvers gls and ols can only be used if xdim == 1')
            
        if solver is None:
            if self.xdim > 1:
                return 'nls',self._leastsq
            else:
                return 'gls',self._gls
            
        if solver == 'gls':
            return 'gls',self._gls
        elif solver == 'nls':
            return 'nls',self._leastsq
        elif solver == 'odr' or solver is None:
            return 'odr',self._odr
        else:
            raise ValueError('solver ' + str(solver) + ' is not recognized')
        
    def _gls(self,**kw):
        from scipy.linalg import lstsq,inv
        from scipy.special import binom
        
        if 'rcond' in kw:
            self.rcond = kw['rcond']
        else:
            self.rcond = len(self.xf)*np.finfo(self.xf.dtype).eps
            
        ic,w,vsc,ignore_corr = self._get_corr(**kw)
        # If there are correlations between the y values and ignore_corr is False,
        # then ic is the cholesky decomposition of the inverse of the covariance 
        # matrix of self.y. ignore_correlations may be True
        # if it was passed in kw or if there was an error when calculating ic.
        # If there are no correlations, but self.uy is an array,then w is an array 
        # of weights for each data point.  vsc is a scalar normalization factor
        # for ic or w.
            
        xav = np.mean(self.xf)
        xd = self.xf-xav # Bad data! Bad! (de-mean the data)
        x = np.vander(xd,self._order)
        sc = np.sqrt((x*x).sum(axis=0))
        x /= sc

        # dunno if this is the best way to do it:  solve once with a method
        # that is quite numerically stable to get the fit parameter best values, 
        # then solve again inverting x.T@x to get the covariances
        if ic is not None:
            x = ic@x
            yf = ic@self.yf
        elif w is not None:
            x *= w[:,np.newaxis]
            yf = self.yf*w
        else:
            yf = self.yf
        self.fit_output = lstsq(x,yf,self.rcond)
        pf,_,rank,_ = self.fit_output
        if rank != self._order:
            warn('the matrix is poorly conditioned',FitWarning)
        pf = (pf.T/sc).T
        pf = pf[::-1]
        
        p = None
        try:
            cov = inv(x.T@x)
            if self._ygummies and not ignore_corr and self.sigma_is_known:
                y = self.y
                if ic is not None:
                    y = [np.sum(y*c) for c in ic] # y = ic@self.y
                elif w is not None:
                    y = w*y
                b = cov@x.T
                p = [np.sum(y*bi) for bi in b]/sc # correlate p with self.y
                p = p[::-1]
            cov *= vsc
        except:
            warn('unable to find the variance-covariance matrix; uncertainties cannot be calculated',FitWarning)
            cov = None
        
        #Find the residuals:
        xm = np.array([[z**i for i in range(self._order)] for z in xd])
        res = self.yf - xm.dot(pf)
        
        #Transform the de-meaned fit parameters and covariance matrix
        #back to the un-de-meaned coordinate system:
        tr = np.array([[0 if j < i else (-xav)**(j-i)*binom(j,i) 
                        for j in range(self._order)] for i in range(self._order)])
        if p is not None:
            p = tr.dot(p)
        if pf is not None:
            pf = tr.dot(pf)
            if p is not None:
                for g,gf in zip(p,pf):
                    g._x = gf
        
        if cov is not None:
            cov /= np.outer(sc,sc)
            cov = cov[:,::-1][::-1]
            cov = tr.dot(cov).dot(tr.transpose())
            
        if ic is not None:
            sigmasq = np.sum((ic@self.uy)**2)/self.count
            self.sigma = np.sqrt(sigmasq)
            var = np.sum((ic@res)**2)/(self.count - self._order)
            if not self.sigma_is_known:
                cov *= var/sigmasq
        elif w is not None:
            sigmasq = ((w*self.uy)**2).sum()/self.count
            self.sigma = np.sqrt(sigmasq)
            var = ((w*res)**2).sum()/(self.count - self._order)
            if not self.sigma_is_known:
                cov *= var/sigmasq
        elif self.uy is not None:
            if _isscalar(self.uy):
                self.sigma = self.uy
                sigmasq = self.uy**2
            else:
                sigmasq = np.sum(self.uy**2)/self.count
                self.sigma = np.sqrt(sigmasq)
            var = (res**2).sum()/(self.count - self._order)
            if cov is not None:
                if self.sigma_is_known:
                    cov *= sigmasq
                else:
                    cov *= var
        else:
            self.sigma = None
            var = (res**2).sum()/(self.count-self._order)
            if cov is not None:
                cov *= var

        self.s = np.sqrt(var)
        self.cov = cov
        self.dof = self.count-self._order
            
        units = [self.yunit*self.xunit**-i for i in range(self._order)]
        
        if p is None:
            if self.uy is None or not self.sigma_is_known:
                dof = self.dof
            else:
                dof = float('inf')
            try:
                self.p = gummy.create(pf,unit=units,dof=dof,
                                  covariance_matrix = self.cov)
            except ValueError:
                self.p = gummy.create(pf,unit=self.punits)
                warn('the covariance matrix is not positive semidefinate; uncertainties cannot be calculated',FitWarning)
        else:
            self.p = [g.graft(units[i]) for i,g in enumerate(p)]
            if pf is None:
                pf = [g.x for g in p]
            
        self.pf = pf
        self.res = res
                     
    def get_p0(self):
        if self.solver in ('ols','gls'):   
            # Linear fit will be used, no p0 needed.
            return [0]*self._order
        # If a nonlinear fit will be used, get p0 from a linear fit:
        try:
            return PolyFit(self.xf,self.yf,deg=self.deg,solver='ols').pf
        except:
            pass
            
        return [np.mean(self.yf)] + ([0]*np.array(self.deg).sum())
        
    def get_punits(self):
        if self.xdim == 1:
            return [self._yunit*self._xunit**-i for i in range(self._order)]
            
        return [self._yunit] + [self._yunit*self._xunit[j]**-i for j in range(self.xdim) for i in range(1,self._order[j])]
        
        
    def f(self,*a):
        if self.xdim == 1:
            r = 0
            for i in range(self._order):
                r += a[i+1]*a[0]**i
            return r
            
        k = self.xdim+1
        r = a[self.xdim]
        for j in range(self.xdim):
            for i in range(1,self._order[j]):
                r += a[k]*a[j]**i
                k += 1
        return r
        
    def jac(self,*a):
        if self.xdim == 1:
            d0 = 0
            for i in range(1,self._order):
                d0 += i*a[i+1]*a[0]**(i-1)
            return [d0 if i == 0 else a[0]**(i-1) for i in range(self._order+1)]

        k = self.xdim
        d0 = [0]*self.xdim
        for j in range(self.xdim):
            for i in range(1,self._order[j]):
                d0[j] += i*a[k]*a[j]**(i-1)
                k += 1
                
        if _isscalar(a[0]):
            return d0 + [1] + [a[j]**i for j in range(self.xdim) for i in range(1,self._order[j])]
        return d0 + [[1]*len(a[0])] + [a[j]**i for j in range(self.xdim) for i in range(1,self._order[j])]
        
    def funicode(self):
        if self.xdim == 1:
            txt = 'y = p(1)'
            if self._order >= 2:
                txt += ' + p(2)*x'
            for i in range(2,self._order):
                txt += ' + p('+str(i+1)+')*x**'+str(i)
            return txt
            
        txt = 'y = p(1)'
        k = 2
        for j in range(self.xdim):
            for i in range(1,self._order[j]):
                if i == 1:
                    txt += ' + p('+str(k)+')*x('+str(j+1)+')'
                else:
                    txt += ' + p('+str(k)+')*x('+str(j+1)+')**'+str(i)
                k += 1
                
        return txt
    
    def flatex(self):
        if self.xdim == 1:
            txt = '$ y = p_{1}'
            if self._order >= 2:
                txt += ' + p_{2} x'
            for i in range(2,self._order):
                txt += ' + p_{'+str(i+1)+'} x^{'+str(i)+'}'
            return txt + ' $'
            
        txt = '$ y = p_{1}'
        k = 2
        for j in range(self.xdim):
            for i in range(1,self._order[j]):
                if i == 1:
                    txt += ' + p_{'+str(k)+'} x_{'+str(j+1)+'}'
                else:
                    txt += ' + p_{'+str(k)+'} x_{'+str(j+1)+'}^{'+str(i)+'}'
                k += 1
                
        return txt + ' $'
    
    def fhtml(self):
        if self.xdim == 1:
            txt = '<i>y</i> = <i>p</i><sub>1</sub>'
            if self._order >= 2:
                txt += ' + <i>p</i><sub>2</sub> <i>x</i>'
            for i in range(2,self._order):
                txt += ' + <i>p</i><sub>'+str(i+1)+'</sub> <i>x</i><sup>'+str(i)+'</sup>'
            return txt
            
        txt = '<i>y</i> = <i>p</i><sub>1</sub>'
        k = 2
        for j in range(self.xdim):
            for i in range(1,self._order[j]):
                if i == 1:
                    txt += ' + <i>p</i><sub>'+str(k)+'</sub> <i>x</i><sub>'+str(j+1)+'</sub>'
                else:
                    txt += ' + <i>p</i><sub>'+str(k)+'</sub> <i>x</i><sub>'+str(j+1)+'</sub><sup>'+str(i)+'</sup>'
                k += 1
                
        return txt


class DoubleExpFit(Fit):
    """
    DoubleExpFit(x,y,p0=None,ux=None,uy=None,sigma_is_known=True,xunit=None, yunit=None,
             solver=None,maxiter=None,nprop=False,**keywords)
   
    Fits the x,y data to a function of the form:
       
    p[0]*np.exp(x/p[1])+p[2]*np.exp(x/p[3])+p[4]
   
    Parameters
    ----------
    x: array_like
       The x-coordinates of the data.  This is a list or numpy array of
       floats or gummys (all point must be of the same type, floats and gummys
       may not be mixed).
    y:  array_like, optional
        The y-coordinates of the data (the type requirements are the
        same as for the x-coordinates).
    p0: array_like of `float`, optional
        The initial values for the fit parameters.
    ux: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `x` values. This should not be specified if the `x`
        argument contains gummys.  If this is specified then only the odr
        solver may be used.  The default is `None`.
    uy: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `y` values. This should not be specified if the y
        argument contains gummys.  The default is `None`.
    sigma_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`. This
        parameter is ignored if `nprop` is `True`.
    xunits, yunits: `str` or `None`, optional
        units for the x and y coordinates. These should not be specified
        if the `x` and `y` parameters contain gummys. These may only be
        specified if the `get_punits` method is overridden in a subclass.
    solver:  {'nls','odr'}, optional
        If this is 'nls' then `scipy.optimize.leastsq` is used to perform
        the fit.  If it is 'odr' then `scipy.odr` is used.  'nls' may not
        be used if the y-coordinate is `None` or multi-dimensional or if
        there is uncertainty in the x-coordinates.  If this is `None`,
        then 'nls' will be used when possible.
    maxiter:  `int` or `None`, optional
        The maximum number of iterations that the solver may use. If this
        is `None` or omitted then the default value for the solver will
        be used.
    nprop:  `bool`, optional
        If this is `True` then uncertainties in the fit will be numerically
         calculated by varying each data point. This will not work if there
         are more than a few data points or if the it is not very stable.
         If this is False than the covariance matrix generated by the solver
         will be used to calculate the uncertainties.  The default value is
         `False`
    other keywords:  optional
        Any additional keyword parameters will be passed to the solver.

    Attributes
    ----------
    p:  `list` of `gummy`
        The fitted values for the fit function parameters as gummys
    pf:  `list` of `float`
        The fitted values for the fit function parameters as floats
    res:  `numpy.ndarray` of `float`
        the fit residuals
    s:  `float`
        the standard deviation (or, when there are uncertainties for
        the input data, the square root of the reduced chi-squared) of
        the residuals
    cov:  `numpy.ndarray` of `float`
        the covariance matrix generated by the solver
    fit_output:
        the raw output of the solver
    x:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the x-coordinates of the data.
    xf:  `numpy.ndarray` of `float`
        numpy array of the x-coordinates of the data as floats
    xdim:  `int`
        the number of dimensions of the x-coordinates
    ux:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the x-coordinates
    y:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the y-coordinates of the data.
    yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
    ydim:  `int`
        the number of dimensions of the y-coordinates
    uy:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the y-coordinates
    count:  `int`
        the number of data points
    p0:  `list` of `float`
        the initial values for the fit function parameters
    solver:  `str`
        the solver used
    punits:  `list` of `Unit`
        the units of the fit parameters
    nparam:  `int`
        the number of fit parameters

    Methods
    -------
    ypred(x1,x2,...):
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
    plot(...):
        plots the data (only available if x and y are one-dimensional)
    """
    
    def get_p0(self):
        ft = ExpFit(self.xf,self.yf)
        r = np.abs(ft.pf[1])/(np.max(self.xf) - np.min(self.xf))
        if r < 0.33:
            a = int(self.count*2*r)
            b = self.count
        elif r < 1:
            a = int(2*self.count/3)
            b = self.count
        else:
            a = 0
            b = int(self.count/3)
        ft2 = ExpFit(self.xf[a:b],ft.yf[a:b])
        return[ft.pf[0],ft.pf[1],ft2.pf[0],ft2.pf[1],ft.pf[2]]
        
    def get_punits(self):
        return [self._yunit, self._xunit, self._yunit, self._xunit, self._yunit]
        
    def f(self,x,p1,p2,p3,p4,p5):
        return p1*np.exp(x/p2)+p3*np.exp(x/p4)+p5
        
    def jac(self,x,p1,p2,p3,p4,p5):
        return ((p1/p2)*np.exp(x/p2)+(p3/p4)*np.exp(x/p4), np.exp(x/p2), -x*(p1/p2**2)*np.exp(x/p2), np.exp(x/p4), 
               -x*(p3/p4**2)*np.exp(x/p4), np.ones(np.shape(x)))
               
    def funicode(self):
        return 'y = p(1)*exp(-x/p(2)) + p(3)*exp(-x/p(4)) + p(5)'
    
    def flatex(self):
        return '$ y = p_{1}\exp(-x/p_{2}) + p_{3}\exp(-x/p_{4}) + p_{5} $'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub> exp(-<i>x</i>/<i>p</i><sub>2</sub>) + <i>p</i><sub>3</sub> exp(-<i>x</i>/<i>p</i><sub>4</sub>) + <i>p</i><sub>5</sub>'
        
        
class ExpFit(Fit):
    """
    ExpFit(x,y,p0=None,ux=None,uy=None,sigma_is_known=True,xunit=None, yunit=None,
             solver=None,maxiter=None,nprop=False,**keywords)
   
    Fits the x,y data to a function of the form:
       
    p[0]*np.exp(x/p[1])+p[2]
   
    Parameters
    ----------
    x: array_like
       The x-coordinates of the data.  This is a list or numpy array of
       floats or gummys (all point must be of the same type, floats and gummys
       may not be mixed).
    y:  array_like, optional
        The y-coordinates of the data (the type requirements are the
        same as for the x-coordinates).
    p0: array_like of `float`, optional
        The initial values for the fit parameters.
    ux: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `x` values. This should not be specified if the `x`
        argument contains gummys.  If this is specified then only the odr
        solver may be used.  The default is `None`.
    uy: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `y` values. This should not be specified if the y
        argument contains gummys.  The default is `None`.
    sigma_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`. This
        parameter is ignored if `nprop` is `True`.
    xunits, yunits: `str` or `None`, optional
        units for the x and y coordinates. These should not be specified
        if the `x` and `y` parameters contain gummys. These may only be
        specified if the `get_punits` method is overridden in a subclass.
    solver:  {'nls','odr'}, optional
        If this is 'nls' then `scipy.optimize.leastsq` is used to perform
        the fit.  If it is 'odr' then `scipy.odr` is used.  'nls' may not
        be used if the y-coordinate is `None` or multi-dimensional or if
        there is uncertainty in the x-coordinates.  If this is `None`,
        then 'nls' will be used when possible.
    maxiter:  `int` or `None`, optional
        The maximum number of iterations that the solver may use. If this
        is `None` or omitted then the default value for the solver will
        be used.
    nprop:  `bool`, optional
        If this is `True` then uncertainties in the fit will be numerically
         calculated by varying each data point. This will not work if there
         are more than a few data points or if the it is not very stable.
         If this is False than the covariance matrix generated by the solver
         will be used to calculate the uncertainties.  The default value is
         `False`
    other keywords:  optional
        Any additional keyword parameters will be passed to the solver.

    Attributes
    ----------
    p:  `list` of `gummy`
        The fitted values for the fit function parameters as gummys
    pf:  `list` of `float`
        The fitted values for the fit function parameters as floats
    res:  `numpy.ndarray` of `float`
        the fit residuals
    s:  `float`
        the standard deviation (or, when there are uncertainties for
        the input data, the square root of the reduced chi-squared) of
        the residuals
    cov:  `numpy.ndarray` of `float`
        the covariance matrix generated by the solver
    fit_output:
        the raw output of the solver
    x:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the x-coordinates of the data.
    xf:  `numpy.ndarray` of `float`
        numpy array of the x-coordinates of the data as floats
    xdim:  `int`
        the number of dimensions of the x-coordinates
    ux:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the x-coordinates
    y:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the y-coordinates of the data.
    yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
    ydim:  `int`
        the number of dimensions of the y-coordinates
    uy:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the y-coordinates
    count:  `int`
        the number of data points
    p0:  `list` of `float`
        the initial values for the fit function parameters
    solver:  `str`
        the solver used
    punits:  `list` of `Unit`
        the units of the fit parameters
    nparam:  `int`
        the number of fit parameters

    Methods
    -------
    ypred(x1,x2,...):
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
    plot(...):
        plots the data (only available if x and y are one-dimensional)
    """
    
    def get_p0(self):
        # Use a second order polynomial fit and set the dervitives at the mean x
        # from the polynomial fit equal to the derivatives of the fit function.
        mnx = np.mean(self.xf)
        pft = PolyFit(self.xf,self.yf,deg=2)
        b = 0.5*(pft.pf[1]+2*pft.pf[2]*mnx)/pft.pf[2]
        a = 2*b**2*pft.pf[2]*np.exp(-mnx/b)
        c = pft.ypredf(mnx)-a*np.exp(mnx/b)
        return[a,b,c]
        
    def get_punits(self):
        return [self._yunit, self._xunit, self._yunit]
        
    def f(self,x,p1,p2,p3):
        return p1*np.exp(x/p2)+p3
        
    def jac(self,x,p1,p2,p3):
        return ((p1/p2)*np.exp(-x/p2), np.exp(x/p2), -x*(p1/p2**2)*np.exp(x/p2), np.ones(np.shape(x)))
        
    def funicode(self):
        return 'y = p(1)*exp(x/p(2)) + p(3)'
    
    def flatex(self):
        return '$ y = p_{1}\exp(x/p_{2}) + p_{3} $'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub> exp(<i>x</i>/<i>p</i><sub>2</sub>) + <i>p</i><sub>3</sub>'
        
        
        
class OneOverTFit(Fit):
    """
    DoubleExpFit(x,y,p0=None,ux=None,uy=None,sigma_is_known=True,xunit=None, yunit=None,
             solver=None,maxiter=None,nprop=False,**keywords)
   
    Fits the x,y data to a function of the form:
       
    p[0]/x + p[1]
   
    Parameters
    ----------
    x: array_like
       The x-coordinates of the data.  This is a list or numpy array of
       floats or gummys (all point must be of the same type, floats and gummys
       may not be mixed).
    y:  array_like, optional
        The y-coordinates of the data (the type requirements are the
        same as for the x-coordinates).
    p0: array_like of `float`, optional
        The initial values for the fit parameters.
    ux: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `x` values. This should not be specified if the `x`
        argument contains gummys.  If this is specified then only the odr
        solver may be used.  The default is `None`.
    uy: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `y` values. This should not be specified if the y
        argument contains gummys.  The default is `None`.
    sigma_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`. This
        parameter is ignored if `nprop` is `True`.
    xunits, yunits: `str` or `None`, optional
        units for the x and y coordinates. These should not be specified
        if the `x` and `y` parameters contain gummys. These may only be
        specified if the `get_punits` method is overridden in a subclass.
    solver:  {'nls','odr'}, optional
        If this is 'nls' then `scipy.optimize.leastsq` is used to perform
        the fit.  If it is 'odr' then `scipy.odr` is used.  'nls' may not
        be used if the y-coordinate is `None` or multi-dimensional or if
        there is uncertainty in the x-coordinates.  If this is `None`,
        then 'nls' will be used when possible.
    maxiter:  `int` or `None`, optional
        The maximum number of iterations that the solver may use. If this
        is `None` or omitted then the default value for the solver will
        be used.
    nprop:  `bool`, optional
        If this is `True` then uncertainties in the fit will be numerically
         calculated by varying each data point. This will not work if there
         are more than a few data points or if the it is not very stable.
         If this is False than the covariance matrix generated by the solver
         will be used to calculate the uncertainties.  The default value is
         `False`
    other keywords:  optional
        Any additional keyword parameters will be passed to the solver.

    Attributes
    ----------
    p:  `list` of `gummy`
        The fitted values for the fit function parameters as gummys
    pf:  `list` of `float`
        The fitted values for the fit function parameters as floats
    res:  `numpy.ndarray` of `float`
        the fit residuals
    s:  `float`
        the standard deviation (or, when there are uncertainties for
        the input data, the square root of the reduced chi-squared) of
        the residuals
    cov:  `numpy.ndarray` of `float`
        the covariance matrix generated by the solver
    fit_output:
        the raw output of the solver
    x:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the x-coordinates of the data.
    xf:  `numpy.ndarray` of `float`
        numpy array of the x-coordinates of the data as floats
    xdim:  `int`
        the number of dimensions of the x-coordinates
    ux:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the x-coordinates
    y:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the y-coordinates of the data.
    yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
    ydim:  `int`
        the number of dimensions of the y-coordinates
    uy:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the y-coordinates
    count:  `int`
        the number of data points
    p0:  `list` of `float`
        the initial values for the fit function parameters
    solver:  `str`
        the solver used
    punits:  `list` of `Unit`
        the units of the fit parameters
    nparam:  `int`
        the number of fit parameters

    Methods
    -------
    ypred(x1,x2,...):
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
    plot(...):
        plots the data (only available if x and y are one-dimensional)    """
    
    def get_p0(self):
        b = np.mean(self.yf[(self.count-2-int(self.count/10)):(self.count-1)])
        a = (self.xf[-1]-self.xf[0])/3.0
        return[a,b]
        
    def get_punits(self):
        return [1/self._xunit, self._yunit]
        
    def f(self,x,p1,p2):
        return p1/x+p2
        
    def jac(self,x,p1,p2):
        return (-p1/x**2,1/x, np.ones(np.shape(x)))
        
    def funicode(self):
        return 'y = p(1)/x + p(2)'
    
    def flatex(self):
        return '$ y = p_{1}/x + p_{2} $'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub>/x + <i>p</i><sub>2</sub>'
        
        
class SinFit(Fit):
    """
    SinFit(x,y,p0=None,ux=None,uy=None,sigma_is_known=True,xunit=None, yunit=None,
             solver=None,maxiter=None,nprop=False,**keywords)
   
    Fits the x,y data to a function of the form:
       
    p[0]*sin(p[1]*x + p[2]) + p[3]
   
    Parameters
    ----------
    x: array_like
       The x-coordinates of the data.  This is a list or numpy array of
       floats or gummys (all point must be of the same type, floats and gummys
       may not be mixed).
    y:  array_like, optional
        The y-coordinates of the data (the type requirements are the
        same as for the x-coordinates).
    p0: array_like of `float`, optional
        The initial values for the fit parameters.
    ux: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `x` values. This should not be specified if the `x`
        argument contains gummys.  If this is specified then only the odr
        solver may be used.  The default is `None`.
    uy: `float`, array_like of `float`  or `None`, optional
        Uncertainty in the `y` values. This should not be specified if the y
        argument contains gummys.  The default is `None`.
    sigma_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`. This
        parameter is ignored if `nprop` is `True`.
    xunits, yunits: `str` or `None`, optional
        units for the x and y coordinates. These should not be specified
        if the `x` and `y` parameters contain gummys. These may only be
        specified if the `get_punits` method is overridden in a subclass.
    solver:  {'nls','odr'}, optional
        If this is 'nls' then `scipy.optimize.leastsq` is used to perform
        the fit.  If it is 'odr' then `scipy.odr` is used.  'nls' may not
        be used if the y-coordinate is `None` or multi-dimensional or if
        there is uncertainty in the x-coordinates.  If this is `None`,
        then 'nls' will be used when possible.
    maxiter:  `int` or `None`, optional
        The maximum number of iterations that the solver may use. If this
        is `None` or omitted then the default value for the solver will
        be used.
    nprop:  `bool`, optional
        If this is `True` then uncertainties in the fit will be numerically
         calculated by varying each data point. This will not work if there
         are more than a few data points or if the it is not very stable.
         If this is False than the covariance matrix generated by the solver
         will be used to calculate the uncertainties.  The default value is
         `False`
    other keywords:  optional
        Any additional keyword parameters will be passed to the solver.

    Attributes
    ----------
    p:  `list` of `gummy`
        The fitted values for the fit function parameters as gummys
    pf:  `list` of `float`
        The fitted values for the fit function parameters as floats
    res:  `numpy.ndarray` of `float`
        the fit residuals
    s:  `float`
        the standard deviation (or, when there are uncertainties for
        the input data, the square root of the reduced chi-squared) of
        the residuals
    cov:  `numpy.ndarray` of `float`
        the covariance matrix generated by the solver
    fit_output:
        the raw output of the solver
    x:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the x-coordinates of the data.
    xf:  `numpy.ndarray` of `float`
        numpy array of the x-coordinates of the data as floats
    xdim:  `int`
        the number of dimensions of the x-coordinates
    ux:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the x-coordinates
    y:  `numpy.ndarray` of `float` or of `gummy`
        numpy array of the y-coordinates of the data.
    yf:  `numpy.ndarray` of `float`
            numpy array of the y-coordinates of the data as floats
    ydim:  `int`
        the number of dimensions of the y-coordinates
    uy:  `float`, `numpy.ndarray` of `floats` or `None`
        uncertainties in the y-coordinates
    count:  `int`
        the number of data points
    p0:  `list` of `float`
        the initial values for the fit function parameters
    solver:  `str`
        the solver used
    punits:  `list` of `Unit`
        the units of the fit parameters
    nparam:  `int`
        the number of fit parameters

    Methods
    -------
    ypred(x1,x2,...):
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
    plot(...):
        plots the data (only available if x and y are one-dimensional)
    """
    
    def get_p0(self):
        # Look for low frequency components with a third order polynomial fit
        # and set the deriviatives from the polynomial fit eqaul to the derivatives
        # of the fit function at the mean x.
        mnx = np.mean(self.xf)
        pft = PolyFit(self.xf,self.yf,deg=3)
        d1 = pft.pf[1] + 2*pft.pf[2]*mnx + 3*pft.pf[3]*mnx**2
        d2 = 2*pft.pf[2] + 6*pft.pf[3]*mnx
        d3 = 6*pft.pf[3]
        b = np.sqrt(abs(d3/d1))
        a = np.sqrt((d1/b)**2 + (d2/b**2)**2)
        c = np.arctan2(-d2/b,d1) - b*mnx
        d = mnx - a*np.sin(b*mnx+c)
        
        # Look for high frequency components with an fft.
        if np.all(np.diff(self.xf) > 0):
            xs = self.xf
            ys = self.yf
        else:
            s = np.argsort(xs)
            xs = self.xf[s]
            ys = self.yf[s]
            
        xi = np.linspace(xs[0],xs[-1],int((xs[-1]-xs[0])/min(np.diff(xs))+1.5))
        yi = np.interp(xi,xs,ys)
        
        fft = np.fft.rfft(yi)
        afft = abs(fft)
        am = np.argmax(afft[2:]) + 2
        aa = abs(fft[am])*2/len(yi)
        amf = am + (afft[am+1]-afft[am-1])/(afft[am-1]+afft[am]+afft[am+1])
        if amf > 1.5 and aa > a:
            a = aa
            b = amf*2*np.pi/(xi[-1]-xi[0])
            c = -np.mod(np.arctan2(fft[am].real,fft[am].imag) + 
                        4*np.pi/len(yi) + np.pi + b*xi[0],2*np.pi)
            if c > np.pi:
                c -= 2*np.pi
            elif c <= -np.pi:
                c += 2*np.pi
            d = abs(fft[0])/len(yi)
        
        return[a,b,c,d]
        
    def get_punits(self):
        return [self._yunit, 1/self._xunit, one, self._yunit]
        
    def f(self,x,p1,p2,p3,p4):
        return p1*np.sin(p2*x+p3) + p4
        
    def jac(self,x,p1,p2,p3,p4):
        return (p1*p2*np.cos(p2*x+p3), np.sin(p2*x+p3), x*p1*np.cos(p2*x+p3),
                p1*np.cos(p2*x+p3), np.ones(np.shape(x)))
        
    def funicode(self):
        return 'y = p(1)*sin(p(2)*x + p(3)) + p(4)'
    
    def flatex(self):
        return '$ y = p_{1}\\sin(p_{2} x +  p_{3}) + p_{4}$'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub> sin(<i>p</i><sub>2</sub> x + <i>p</i><sub>3</sub>) + <i>p</i><sub>4</sub>'
        
        