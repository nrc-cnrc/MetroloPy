# -*- coding: utf-8 -*-

# module fit

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
This module defines some classes to facilitate curve fitting.
"""

from .ummy import _isscalar,UncertainValue,njacobian
from .gummy import gummy
from .exceptions import FitWarning
from .unit import Unit,one
from .printing import PrettyPrinter
from .unit import Quantity

import numpy as np
from warnings import warn

def odr_jac(work,count,xdim,nparam,ydim,nwe,task):
    # get the jacobian with respect to the fit parameters and the jacobian with 
    # respect to x from the rwork array returned by odrpack.odr_fit
    a = 2*count*xdim + 3*count*ydim  + nparam*nparam + 10*nparam + 17
    b = a + count*nparam*ydim
    jp = work[a:b]
    if ydim == 1:
        jp = np.reshape(jp,(count,nparam),order='F')
    else:
        jp = np.reshape(jp,(count,nparam,ydim),order='F')

    c = b + nwe*ydim + ydim*(nparam + xdim) + 4*count*xdim + ydim*ydim
    d = c + count*xdim*ydim
    if task == 'OLS':
        jx = None
    else:
        jx = work[c:d]
        if ydim == 1:
            jx = np.reshape(jx,(count,xdim),order='F')
        else:
            jx = np.reshape(jx,(count,xdim,ydim),order='F')
    
    return (jp,jx)

def getvalues(txt,x,u,cov,w,units,ignore_corr):
    # unpack Fit x or y values
    
    from scipy.linalg import pinvh
    
    if x is None:
        return None,None,None,None,None,None,None,None,0
    elif _isscalar(x):
        if isinstance(x,UncertainValue):
            xf = float(x.x)
            u = x.u
        else:
            xf = float(x)
            u = None
        if units is None:
            units = one
        else:
            units = Unit.unit(units)
        return x,xf,u,None,w,None,None,units,0
    
    x = np.asarray(x)
        
    shape = np.shape(x)
    if np.ndim(x) == 1:
        dim = 1
        #x = np.reshape(x,(1,len(x)))
    else:
        dim = shape[0]
        
    shape = np.shape(x)
    
    if units is not None:
        if _isscalar(units):
            units = [Unit.unit(units)]*dim
        else:
            if np.shape(units) != (dim,):
                raise TypeError(txt + 'units must be either None, a single unit or a list of units with length equal to xdim')
            units = [Unit.unit(u) for u in units]
        if dim == 1:
            units = units[0]
        
    pcorr = False
    if x.dtype is np.dtype('O'):
        xx = np.empty(shape,dtype=np.dtype('O'))
        xf = np.empty(shape)
        uf = np.empty(shape)

        if dim == 1:
            units = x[0].unit if isinstance(x[0],Quantity) else one
        if units is None:
            units = [i[0].unit if isinstance(i[0],Quantity) else one for i in x]
            
        with np.nditer(x,flags=['multi_index','refs_ok']) as it:
            for i in it:
                if isinstance(i,np.ndarray) and not np.ndim(i):
                    ii = i.item()
                else:
                    ii = i
                if dim == 1:
                    iunit = units
                else:
                    iunit = units[it.multi_index[0]]
                if iunit is one:
                    ii = ii.convert(one) if isinstance(i,Quantity) else ii
                else:
                    if isinstance(ii,Quantity):
                        if ii.unit is one:
                            ii *= iunit
                        else:
                            ii = ii.convert(iunit)
                    else:
                        ii = gummy(ii,unit=iunit)
                xx[it.multi_index] = ii
                ii = ii.value if isinstance(ii,Quantity) else ii
                if isinstance(ii,UncertainValue):
                    if ii.u > 0:
                        pcorr = True
                    xf[it.multi_index] = ii.x
                    uf[it.multi_index] = ii.u
                else:
                    xf[it.multi_index] = ii
                    uf[it.multi_index] = 0
    
        if np.any(uf):
            if dim == 1:
                cov = gummy.covariance_matrix(xx)
            else:
                cov = np.array([gummy.covariance_matrix(i) for i in xx.T])
            u = uf
            x = xx
        else:
            x = np.array(list(xx))
    else:
        if units is None:
            if dim == 1:
                units = one
            else:
                units = [one]*dim
            
        if x.dtype is np.dtype('float64'):
            xf = x
        else:
            xf = x.astype(np.dtype('float64'))

    if u is None and cov is not None:
            if dim == 1:
                u = np.sqrt(np.diagonal(cov))
            else:
                u = np.sqrt(np.array([np.diagonal(i) for i in cov]).T)
           
    if u is not None and cov is not None:
        if dim == 1:
            if not np.any(cov - np.diag(u**2)):
                cov = None
        else:
            if not np.any(np.any(cov - np.diag((u.T)**2) for i in cov)):
                cov = None
                
    if ignore_corr:
        cov = None
        if w is not None and not _isscalar(w):
            if dim == 1 and len(np.shape(w)) > 1:
                w = np.diagonal(w)
            elif len(np.shape(w)) > 2 and len(w) == len(xf[0]) and len(w[0]) == dim:
                w = np.array([np.diag(np.diagonal) for i in w])
                
    var = None
    if u is not None:
        if w is None:
            if cov is not None:
                if dim == 1:
                    w = pinvh(cov)
                else:
                    w = np.array([pinvh(i) for i in cov])
            else:
                w = 1/u**2
            var = 1
        else:
            if cov is not None:
                if dim == 1:
                    if _isscalar(w):
                        ww = np.full(len(x),w)
                    else:
                        ww = w
                    var = np.sum(ww.T@cov@ww)/np.size(w)
            else:
                var = np.sum(w*u**2)/np.size(w)
        
    return x,xf,u,cov,w,var,pcorr,units,dim

def wsqrt(w):
    # used by the _least_squares and _ols solvers to get the square root of the
    # weights
    
    from scipy.linalg import cholesky
    
    if w is None:
        return 1
    if np.ndim(w) < 2:
        return np.sqrt(w)
    return cholesky(w)
    

class _Fit:
    # Base class for fitting.  This implements plotting and unpacks any
    # uncertainty or unit information in the data.
    
    plot_points = 100
    over_plot = 0.05
    
    xlabel = None
    ylabel = None
    
    def __init__(self,x,y,ux=None,uy=None,xweights=None,weights=None,xcov=None,
                 ycov=None,variance_is_known=True,xunit=None,yunit=None,ignore_correlations=False):
        self.variance_is_known = variance_is_known
            
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
                
        self.x,self.xf,self.ux,self.xcov,self.xweights,self.known_xvar,self.x_is_gummies,self._xunit,self.xdim = getvalues('x',x,ux,xcov,xweights,xunit,ignore_correlations)
        if y is not None and not _isscalar(y):
            self.y,self.yf,self.uy,self.ycov,self.weights,self.known_var,self.y_is_gummies,self._yunit,self.ydim = getvalues('y',y,uy,ycov,weights,yunit,ignore_correlations)
        else:
            self.y = y
            self.yf = None
            self.uy = None
            self._yunit = None
            self.ydim = None
            self.ypcorr = None
            self.covy = None
            self.weights = None
            self.known_var = None
            self.y_is_gummies = False
            
        self.sqrt_weights = wsqrt(self.weights)
        self.sqrt_xweights = wsqrt(self.xweights)
    
        self.implicit = y is None
            
        if self.known_var is None and self.known_xvar is None:
            self.variance_is_known = False
        
        self.count = self.xf.shape[-1]
                       
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
            fy = self.ypredf(fx)
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
        if self.known_var is None:
            s = self.s
        else:
            s = np.sqrt(self.known_var)
        return k*(np.sqrt(self.ypred(x).u**2 + s**2))
        
    
class Fit(_Fit,PrettyPrinter):
    
    latex_math = None
    
    def __init__(self,x,y=None,p0=None,ux=None,uy=None,variance_is_known=True,
                 xunit=None, yunit=None,solver=None,xweights=None,weights=None,
                 xcov=None,ycov=None,ignore_correlations=False,**kw):
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
        jac: function, optional
            The jacobian of the fit function.  
            
            If provided jac  must have the same signature as the f method and 
            return a list (or array) of derivatives of the form:
            
            [df/dx1,df/dx2,...,df/dp1,df/dp2,...] 
        
            if f returns a scalar or:
        
            [[df1/dx1,df1/dx2,...,df1/dp1,df1/dp2,...],
            [df2/dx1,df2/dx2,...,df2/dp1,df2/dp2,...],...]
        
            if f returns a 1-d array [f1,f2,...].
            
            Instead of passing the jacobian as a parameter, the jac method may 
            be overridden is a subclass.  If no jacobian is available it will 
            be calculated numerically.
        ux: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `x` values. This should not be specified if the `x`
            argument contains gummys.  If this is specified then only the odr
            solver may be used.  The default is `None`.
        uy: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `y` values. This should not be specified if the y
            argument contains gummys.  The default is `None`.
        variance_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`.
        xunits, yunits: `str` or `None`, optional
            units for the x and y coordinates. These should not be specified
            if the `x` and `y` parameters contain gummys. These may only be
            specified if the `get_punits` method is overridden in a subclass.
        solver:  {'ols','nls','odr'}, optional
            If this is 'nls' then `scipy.optimize.least_squares` is used to perform
            the fit.  If it is 'odr' then `odrpack.odr_fit` is used.  'nls' may 
            not be used if the y-coordinate is `None` or multi-dimensional or if
            there is uncertainty in the x-coordinates.  If this is `None`,
            then 'nls' will be used when possible.
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
            
        jacp(self,x1,x2,...,xk,p1,p2,...,pk):
            The Jacobian of the fit function with respect to the fit parameters.
            
            If the jac method (above) is not defined a `NotImplementedError` 
            may be raised.
            
            It has the signature:
                
            [df/dp1,df/dp2,...] 
            
            if f returns a scalar or:
            
            [[df1/dp1,df1/dp2,...],
             [df2/dp1,df2/dp2,...],...]
            
            if f returns a 1-d array [f1,f2,...].
            
        jacx(self,x1,x2,...,xk,p1,p2,...,pk):
            The derivative fit function with respect to the x-values.
            
            If the jac method (above) is not defined a `NotImplementedError` 
            may be raised.
            
            It has the signature:
                
            [df/dp1,df/dp2,...] 
            
            if f returns a scalar or:
            
            [[df1/dp1,df1/dp2,...],
             [df2/dp1,df2/dp2,...],...]
            
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
        super().__init__(x,y,ux=ux,uy=uy,variance_is_known=variance_is_known,
                         xunit=xunit,yunit=yunit, xweights=xweights,
                         weights=weights,xcov=xcov,ycov=ycov,
                         ignore_correlations=ignore_correlations)
        
        if '_ch_create_p' in kw:
            self._ch_create_p = kw.pop('_ch_create_p')
        else:
            self._ch_create_p = False
            
        if solver is not None:
            solver = solver.strip().lower()
        self.solver = solver
        self.solver,self._solver = self._get_solver(solver,**kw)
        
        self.xvar = None
        self.yvar = None
        self.xres = None
        self.yres = None
        self.njacp = None
        self.njacx = None
        
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
            
        self.nparam = len(self.p0)
            
        if 'f' in kw:
            self.f = kw['f']
            del kw['f']
            
        if 'jacp' in kw:
            self.jacp = kw['jacp']
            del kw['jacp']
            
        if 'jacx' in kw:
            self.jacx = kw['jacx']
            del kw['jacx']
                 
        try:
            # See if f will broadcast properly across the xf array...
            if self.xdim == 1:
                self.f0 = self.f(self.xf,*self.p0)
            else:
                self.f0 = self.f(*self.xf,*self.p0)
            
            if np.shape(self.f0)[-1] != self.count:
                raise TypeError
                
            if not isinstance(self.f0,np.ndarray):
                self.f = lambda *x: np.array(self.f(*x))
                self.f0 = np.array(self.f0)
                
        except NotImplementedError:
            raise TypeError('the fit function f has not been sepecified')
        except:
            self.f0 = None
            # ...if not vectorize it so that it does.
            if self.ydim is None or self.ydim == 1:
                ot = [np.float64]
            else:
                ot = [np.float64]*self.ydim
            self.f = np.vectorize(self.f,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))
            
        # _f takes two arguments x and p rather than *x,*p
        if self.xdim == 1:
            self._f = lambda x,p:self.f(x,*p)
        else:
            self._f = lambda x,p:self.f(*x,*p)
        
        if self.f0 is None:
            self.f0 = self._f(self.xf,self.p0)
            
        if self.yf is not None and np.shape(self.yf) != np.shape(self.f0):
            raise TypeError('The shape of the fit function f(x) does not match y')
            
        jok = True
        self._jacp = None
        try:
            # See if jacp will broadcast properly across the xf array...
            if self.xdim == 1:
                r = self.jacp(self.xf,*self.p0)
            else:
                r = self.jacp(*self.xf,*self.p0)
                
            if np.shape(r)[-1] != self.count:
                raise TypeError
            
            if not isinstance(r,np.ndarray):
                self.jacp = lambda *x: np.array(self.jacp(*x))
        except NotImplementedError:
            jok = False
        except:
            # ...if not vectorize it so that it does.
            ot = [np.float64]*self.nparam
            self.jacp = np.vectorize(self.jacp,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))
        finally:
            if jok:
                if self.xdim == 1:
                    self._jacp = lambda x,p:self.jacp(x,*p)
                else:
                    self._jacp = lambda x,p:self.jacp(*x,*p)
        
        jok = True
        self._jacx = None
        try:
            # See if jacp will broadcast properly across the xf array...
            if self.xdim == 1:
                r = self.jacx(self.xf,*self.p0)
            else:
                r = self.jacx(*self.xf,*self.p0)

            if np.shape(r)[-1] != self.count:
                raise TypeError
            
            if not isinstance(r,np.ndarray):
                self.jacp = lambda *x: np.array(self.jacp(*x))
        except NotImplementedError:
            jok = False
        except:
            # ...if not vectorize it so that it does.
            ot = [np.float64]*self.xdim
            self.jacx = np.vectorize(self.jacx,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))
        finally:
            if jok:
                if self.xdim == 1:
                    self._jacx = lambda x,p:self.jacx(x,*p)
                else:
                    self._jacx = lambda x,p:self.jacx(*x,*p)

        nprop = False
        if 'nprop' in kw:
            nprop = kw.pop('nprop')
            
        self._solver(**kw)
        
        if nprop and (self._xgummies or self._ygummies):
            self._solver(**kw)
            p0 = self.p0
            self.p0 = self.pf
            xf = self.xf
            yf = self.yf
            res = self.res
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
            self.s = s
            self.cov = cov
            self.fit_output = fit_output
        else:
            for p in self.p:
                if not p.unit.linear:
                    p._set_U()
            
    def _get_solver(self,solver,**kw):
        if solver is None:
            return 'nls',self._leastsq
        
        solver = solver.strip().lower()
            
        if self.ydim > 1:
            if solver is not None and solver != 'odr':
                raise ValueError(str(solver) + ' has been manually selected, but only the odr solver supports y-values with dim > 1')
            return 'odr',self._odr
        
        if self.ux is not None or self.y is None:
            if solver is None:
                return 'odr',self._odr
            
        if solver == 'nls':
            return 'nls',self._leastsq
        elif solver == 'odr':
            return 'odr',self._odr
        elif solver == 'ols' or solver == 'gls':
            return 'ols',self._ols
        else:
            raise ValueError('solver ' + str(solver) + ' is not recognized')
            
    def get_jacx(self):
        try:
            if self.xdim <= 1:
                return self.jacx(self.xf,*self.pf).T*self.sqrt_weights
            else:
                if np.ndim(self.sqrt_weights) <= 1:
                    return self.jacx(*self.xf,*self.pf)*self.sqrt_weights
                else:
                    return self.jacx(*self.xf,*self.pf)@self.sqrt_weights
        except NotImplementedError:
            # numerically calculate the jacobian of the fit function with respect to x
            f = self.f
            x = self.x
            xf = self.xf
            pf = self.pf
            w = self.sqrt_weights
            if self.xdim == 1:
                x = [x[i] if isinstance(x[i],UncertainValue) and x[i].u > 0 else 
                              gummy(xf[i],1) for i in range(len(x))]
                r = np.array([njacobian(lambda z:f(z,*pf),i) for i in x]).T
            else:
                x = [[x[j][i] if isinstance(x[i],UncertainValue) and x[i].u > 0 else 
                              gummy(xf[i],1) for i in len(x[0])] for j in range(len(x))]
                r = np.array([njacobian(lambda z:f(*z,*pf),*i) for i in x]).T
                
            if np.ndim(w) <= 1:
                return np.reshape(w*r,(len(x),1))
            else:
                if np.ndim(w) <= 1:
                    return (w*r).T
                else:
                    return (w@r).T
            
    def _create_p(self):
        # Given fit parameters pf (float array) and the covariance matrix cov for the
        # parameters, return correlated gummys representing the fit parameters. jac,
        # if provided, should be the jacobian mutiplied by any weighting of the data
        # points.
        
        if self.njacp is None:
            try:
                return gummy.create(self.pf,covariance_matrix=self.cov,dof=self.dof)
            except np.linalg.LinAlgError:
                warn('the covariance matrix is not positive semidefinate; uncertainties cannot be calculated',FitWarning)
                return gummy.create(self.pf)
            
        if self.y_is_gummies or self.x_is_gummies or self._ch_create_p:
            odr = self.solver.startswith('odr') and self.xvar > 0
            
            uy = self.uy
            if self.yvar is None:
                yvar = self.var
            else:
                yvar = self.yvar
            if self.variance_is_known:
                dof = float('inf')
            elif yvar is not None:
                uy = np.sqrt(yvar)
                dof = (self.count - 1)/self.count
            
            y = self.y
            if not self.y_is_gummies and uy is not None:
                if self.ydim == 0:
                    y = [self.y]
                y = np.empty(np.shape(self.yf),dtype=np.dtype('O'))
                scl = _isscalar(uy)
                with np.nditer(y,flags=['multi_index','refs_ok'],op_flags=['readwrite']) as it:
                    for i in it:
                        if scl:
                            iuy = uy
                        else:
                            iuy = uy[it.multi_index]
                        i[...] = gummy(0,iuy,dof=dof)
                        
            if self.weights is not None:
                if np.ndim(self.sqrt_weights) == 2:
                    y = self.sqrt_weights@y
                else:
                    y = self.sqrt_weights*y
                            
            ux = self.ux
            if self.variance_is_known:
                xdof = float('inf')
            elif odr:
                ux = np.sqrt(self.xvar)
                xdof = (self.count - 1)/self.count
            
    
            if odr:
                if self.njacx is None:
                    self.njacx = self.get_jacx()
                wr = 1/self.sqrt_xweights
                cjx = np.cos(np.arctan((self.njacx.T)[0]*wr))
                y = cjx*y
                
                
            p = self.rcov@(y@self.njacp)
    
            if ux is not None and self.x is not None:
                x = self.x
                if not self.x_is_gummies and ux is not None and x is not None:
                    if self.xdim == 0:
                        x = [x]
                    x = np.empty(np.shape(self.xf),dtype=np.dtype('O'))
                    scl = _isscalar(ux)
                    with np.nditer(x,flags=['multi_index','refs_ok'],op_flags=['readwrite']) as it:
                        for i in it:
                            if scl:
                                iux = ux
                            else:
                                iux = ux[it.multi_index]
                            i[...] = gummy(0,iux,dof=xdof)
                            
                if self.xweights is not None:
                    if np.ndim(self.sqrt_xweights) == 2:
                        x = self.sqrt_xweights@x
                    else:
                        x = self.sqrt_xweights*x
                        
                if odr:
                    x = cjx*x
                else:
                    if self.njacx is None:
                        self.njacx = self.get_jacx()
                    x = x*(self.njacx.T)[0]
    
                px = self.rcov@(x@self.njacp)
                p = p + px
                if odr and self.known_var:
                    p = p/np.sqrt(2)
                    
            for g,gf in zip(p,self.pf):
                g.value._x = gf
            return np.array([g*u if g.unit is one else g for g,u in zip(p,self.punits)])
        
        # If we have a jacobian, calculate effiective degrees of freedom for each
        # parameter as sum(weights)**2/sum(weights**2) where the weights are the
        # square of the elements of the jacobian.
        try:
            pf = np.asarray(self.pf)
            cov = np.asarray(self.cov)
            pjac = np.asarray(self.njacp)
            n = len(pf)
            
            # sqrtm below if used to transform the jacobian from a basis with 
            # correlated parameters to one where the parameters are uncorrelated
            m = np.array([[cov[i][j]/np.sqrt(cov[i][i]*cov[j][j]) if cov[i][i] != 0 and cov[j][j] != 0 else 0 for i in range(n)] for j in range(n)])
            u = np.sqrt(np.diag(cov))
            val,vec = np.linalg.eig(m)
            val = np.real(val)
            val = np.clip(val,0,None)
            vec = np.real(vec)
            sqrtm = ((vec*np.sqrt(val)@np.linalg.inv(vec)).T*u).T
            jact = sqrtm.T@pjac.T
            
            # calculate the effective degrees of freedom for the uncorrelated 
            # parameters then transform the resulting gummys back to the
            # correlated basis
            df = (np.sum(jact**2,axis=1)**2/np.sum(jact**4,axis=1))*(self.count - 1)/self.count
            df = np.array([i if i >= 1 else 1 for i in df])
            g = np.array([gummy(0,1,dof=df[i]) for i in range(n)])
            p = sqrtm@g + pf
            return np.array([g*u for g,u in zip(p,self.punits)])
                        
        except np.linalg.LinAlgError:
            warn('unable to calculate the effective degrees of freedom for the fit parameters')
    
    
    def _leastsq(self,**kw):
        # non-linear least square solver
        from scipy.optimize import least_squares
        
        from scipy.linalg import pinvh,LinAlgError
            
        w = self.sqrt_weights
        if np.ndim(w) == 2:
            def func(params,x,y,f):
                return w@(f(x,params) - y)
                
            if self._jacp is not None:
                def dfun(*a):
                    return (self._jacp(a[1],a[0])@w.T).T
            else:
                dfun = None
        else:
            def func(params,x,y,f):
                return w*(f(x,params) - y)
                    
            if self._jacp is not None:
                def dfun(*a):
                    return (w*np.asarray(self._jacp(a[1],a[0]))).T
            else:
                dfun = None
        if dfun is not None:
            kw['jac'] = dfun
            
        args = (self.xf,self.yf,self._f)
        self.fit_output = least_squares(func,self.p0,args=args,**kw)
        self.pf = self.fit_output.x
        
        if self.fit_output.status not in {1,2,3,4}:
            warn('the fit was not sucessful, status ' + str(self.fit_output.status) + ', ' + self.fit_output.message,FitWarning)

        if self.xdim == 1:
            self.res = self.yres = self.yf - self.ypredf(self.xf)
        else:
            self.res = self.yres = self.yf - self.ypredf(*self.xf)
        self.xres = None
        
        if np.ndim(w) == 2:
            var = np.sum((w@self.res)**2)/(self.count - self.nparam)
        else:
            var = ((w*self.res)**2).sum()/(self.count - self.nparam)
                    
        self.dof = self.count - len(self.pf)
        self.var = self.yvar = var
        self.xvar = None
        self.s = np.sqrt(var)
        
        try:
            cov = pinvh(self.fit_output.jac.T@self.fit_output.jac)
        except LinAlgError:
            self.cov = None
            warn('unable to calculate covariances for the fit parameters',FitWarning)
            return np.array([gummy(i) for i in self.pf])
            
        self.rcov = cov
        if self.variance_is_known and self.known_var is not None:
            self.cov = cov*self.known_var
        else:
            self.cov = cov*var
        
        self.njacp = self.fit_output.jac
        self.xweights = None
        self.sqrt_xweights = 1
        
        self.p = self._create_p()
        
    def _odr(self,**kw):
        #orthogonal distance regression solver
        from odrpack import odr_fit
        
        if np.ndim(self.weights) == 2 and len(self.weights[0]) > self.ydim:
            weight_y = np.diagonal(self.weights)
        else:
            weight_y = self.weights
            
        if np.ndim(self.xweights) == 2 and len(self.xweights[0]) > self.xdim:
            weight_x = np.diagonal(self.xweights)
        else:
            weight_x = self.xweights
                
        if self.implicit:
            ydim = 1
            yf = np.zeros(self.f0.shape)
        else:
            ydim = self.ydim
            yf = self.yf
        
        if self.implicit:
            fit_type = 'implicit-ODR'
        else:
            fit_type = 'explicit-ODR'
        if 'fit_type' in kw:
            fit_type = kw.pop('fit_type')
        if 'task' in kw:
            fit_type = kw.pop('task')
                    
        self.fit_output = odr_fit(self._f,self.xf,yf,self.p0,weight_x=weight_x,
                                  weight_y=weight_y,task=fit_type,
                                  jac_beta=self._jacp,jac_x=self._jacx)
        
        if self.fit_output.info > 5:
            raise RuntimeError('Error encountered during fit:\n' + str(self.fit_output.stopreason) + '\nODR info  = ' + str(self.fit_output.info))
        if self.fit_output.info == 4:
            warn('the iteration limit was reached',FitWarning)

        pf = self.fit_output.beta

        self.res = np.array([self.fit_output.delta,self.fit_output.eps])
        self.xres = self.fit_output.delta
        self.yres = self.fit_output.eps
        
        self.var = self.fit_output.res_var
        self.xvar = np.sum(self.xres**2)/(self.count - self.nparam)
        self.yvar = np.sum(self.yres**2)/(self.count - self.nparam)
        
        self.s = np.sqrt(self.var)
       
        self.dof = self.count - len(pf)
        
        self.rcov = self.fit_output.cov_beta
        if self.known_var is None or not self.variance_is_known:
            self.cov = self.var*self.fit_output.cov_beta
        elif self.variance_is_known:
            self.cov = self.known_var*self.fit_output.cov_beta
            
        if self.implicit or weight_y is None or _isscalar(weight_y):
            nwe = 1
        elif np.ndim(weight_y):
            nwe = np.size(weight_y)
        self.njacp,self.njacx = odr_jac(self.fit_output.rwork,self.count,
                                        self.xdim,self.nparam,ydim,nwe,
                                        fit_type)
        
        self.pf = pf
        self.p = self._create_p()
        
    def _shiftx(self,xav):
        # PolyFit implements this method which returns a matrix that transforms 
        # the fit parameters under: x -> x + xav
        
        # if this is impemented, the jacp method must also be implemented
        
        raise NotImplementedError()
                    
    def _ols(self,**kw):
        #ordinary least squares solver
        from scipy.linalg import pinvh
            
        try:
            # PolyFit implements _shiftx which returns a matrix that transforms 
            # the fit parameters under: x -> x + xav
            if self.xdim == 1:
                xav = np.mean(self.xf)
                xd = [self.xf - xav] # Bad data! Bad! (de-mean the data)
            else:
                xav = np.mean(self.xf,axis=1)
                xd = (self.xf.T - xav).T
            trp = self._shiftx(xav)

            x = self.jacp(*xd,*self.p0).T
            if self._jacp is None:
                raise NotImplementedError()
            self.njacp = self._jacp(self.xf,self.p0).T
        except NotImplementedError:
            trp = None
            if self._jacp is not None:
                x = self._jacp(self.xf,self.p0).T
            else:
                # numerically calculate the jacobian with repect to the fit 
                # parameters
                g0 = [gummy(p,1) for p in self.p0]
                if self.xdim == 1:
                    x = np.array([njacobian(lambda *p:self.f(i,*p),*g0) for i in self.xf])
                else:
                    x = np.array([njacobian(lambda *p:self.f(*i,*p),*g0) for i in self.xf.T]).T
            self.njacp = np.array(x)
        
        yf = self.yf
        w = self.sqrt_weights
        if np.ndim(w) == 2:
            x = w@x
            self.njacp = w@self.njacp
            yf = w@yf
        elif np.ndim(w) == 1:
            x = (x.T*w).T
            self.njacp = (w*self.njacp.T).T
            yf = yf*w
            
        sc = np.sqrt((x*x).sum(axis=0))
        x /= sc
        
        kw['return_rank'] = True
        
        cov,rank = pinvh(x.T@x,**kw)
        if rank != self.nparam:
            warn('the rank of the covariance matrix is not equal to the number of parameters',FitWarning)
        
        pf = cov@(yf@x)
        
        self.res = self.yres = self.yf - x@pf # residuals
        self.xres = None
        
        var = np.sum((yf - x@pf)**2)/(self.count - self.nparam)
        
        self.dof = self.count - self.nparam

        pf = (pf.T/sc).T
        cov /= np.outer(sc,sc)
        
        if trp is not None:
            cov = trp@cov@trp.T
            pf = trp@pf
        
        self.rcov = cov
            
        if self.variance_is_known and self.known_var is not None:
            self.cov = cov*self.known_var
        else:
            self.cov = cov*var
        
        self.var = self.yvar = var
        self.xweights = None
        self.sqrt_xweights = 1
        self.s = np.sqrt(var)
        self.pf = pf
        self.p = self._create_p()
        self.fit_output = None
            
    def ypredf(self,*x):
        """
        returns a float representing the value predicted by the fit at x
        """
        if len(x) != self.xdim:
            raise TypeError('the number of arguments to ypredf must be equal to xdim')
            
        if self.xdim == 1:
            x = x[0]
            scl = _isscalar(x)
            if scl:
                x = np.array([x])
        else:
            scl = np.all(_isscalar(i) for i in x)
            if scl:
                x = [np.array([i]) for i in x]

        ret = self._f(x,self.pf)
        
        if scl:
            if self.ydim is None or self.ydim == 1:
                return ret[0]
            else:
                return np.array([i[0] for i in ret])
        else:
            return np.asarray(ret)
        
    def ypred(self,*x):
        """
        returns a gummy representing the value predicted by the fit at x
        """
        if len(x) != self.xdim:
            raise TypeError('the number of arguments to ypred must be equal to xdim')
        
        if self.xdim == 1:
            scl = _isscalar(x[0])
            if scl:
                x = np.array([x])
        else:
            scl = np.all(_isscalar(i) for i in x)
            if scl:
                x = [np.array([i]) for i in x]
            
        p = [i/i.unit for i in self.p]
        
        if self.xdim == 1:
            xu = [self._xunit]
        else:
            xu = self._xunit
        x = [np.array([j.convert(xu[i])/xu[i] if isinstance(j,Quantity) else j for j in x[i]])
             for i in range(self.xdim)]
        
        try:
            r = gummy.apply(self.f,lambda *z:self.jac(*z),*x,*p)
        except NotImplementedError:
            r = gummy.napply(self.f,*x,*p)
            
        if self.ydim == 1:
            r = np.array([i*self._yunit for i in r])
        else:
            r = np.array([[j*self._yunit[i] for j in r[i]] for i in self._ydim])
            
        if scl:
            if self.ydim is None or self.ydim == 1:
                return r[0]
            else:
                return np.array([i[0] for i in r])
        else:
            return np.asarray(r)
                
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
                txt += self.flatex().strip('$') + ' \\\\ [10pt]'
            elif fmt == 'html':
                txt += self.fhtml() + '<br><br>'
            else:
                raise ValueError('format ' + str(fmt) + ' is not understood')
        except NotImplementedError:
            pass

        if fmt == 'unicode' or fmt == 'ascii':
            txt += 'best fit parameters:\n'
        elif fmt == 'latex':
            txt += type(self).latex_norm('best fit parameters:') + ' \\\\ '
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
        The Jacobian of the fit function.
            
        [df/dx1,df/dx2,...,df/dp1,df/dp2,...] 
        
        if f returns a scalar or:
        
        [[df1/dx1,df1/dx2,...,df1/dp1,df1/dp2,...],
         [df2/dx1,df2/dx2,...,df2/dp1,df2/dp2,...],...]
        
        if f returns a 1-d array [f1,f2,...].
        
        A NotImplementedError will be raised unless both the `jacx` and `jacp`
        methods are implemented.
        """
        return np.concatenate([self.jacx(*a),self.jacp(*a)])
        
    def jacp(self,*a):
        """
        The Jacobian of the fit function with respect to the fit parameters.
        
        This may be implemented in a derived class.  Note that this function
        can also be passed to the Fit initializer using the jacp keyword or
        can be omitted entirely in which case the Jacobian will be calculated
        numerically.
        
        It must either have signature:
        
        jacp(self,x,p1,p2,...,pn) 
        
        where there are p1 to pn are the n fit parameters and the independent
        variable x has one dimension, or:
        
        jacp(self,x1,x2,...,xm,p1,p2,...,pn)
        
        where the independent variable x has m dimensions at each observation.
        
        It must return an array (or array like object) with values:
            
        [df/dx1,df/dx2,...]
        
        if f returns a scalar or:
        
        [[df1/dx1,df1/dx2,...],
         [df2/dx1,df2/dx2,...],...]
        
        if f returns a 1-d array [f1,f2,...].
        """
        
        raise NotImplementedError()
    
    def jacx(self,*a):
        """
        The derivative fit function with respect to the x-values:
            
        This may be implemented in a derived class.  Note that this function
        can also be passed to the Fit initializer using the jacp keyword or
        can be omitted entirely in which case this Jacobian will be calculated
        numerically if it is needed.
        
        It must either have signature:
        
        jacp(self,x,p1,p2,...,pn) 
        
        where there are p1 to pn are the n fit parameters and the independent
        variable x has one dimension, or:
        
        jacp(self,x1,x2,...,xm,p1,p2,...,pn)
        
        where the independent variable x has m dimensions at each observation.
        
        It must return an array (or array like object) with values:
            
        [df/dp1,df/dp2,...] 
        
        if f returns a scalar or:
        
        [[df1/dp1,df1/dp2,...],
         [df2/dp1,df2/dp2,...],...]
        
        if f returns a 1-d array [f1,f2,...].
        """
        
        raise NotImplementedError()
        
    def get_p0(self):
        """
        Not implemented, may optionally be implemented by a derived class.
        
        Returns an initial guess for the the fit parameters [p1,p2,...] based
        on the input x and y data.  
        
        If this method is not implemented then the inital values must be passed 
        in the p0 parameter to the Fit initializer.
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
    def __init__(self,x,y,deg=1,ux=None,uy=None,variance_is_known=True,p0=None,
                 xunit=None,yunit=None,solver=None,**kw):
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
            the get_puints method is implemented in a subclass.  The maximum
            dimension for x is 3.
        y:  array_like
            The y-coordinates of the data (shape and type requirements are the
            same as for the x-coordinates).
        deg:  `int` or (for multi-dimensional x-data) array_like of `int`
            the degree of the polynomial
        p0: array_like of `float`, optional The initial values for the fit 
            parameters.  These are only used if solver is 'odr' or 'nls'.  If
            p0 is ommitted and the selected solver is 'odr' or 'nls', the 'ols' 
            solver will be used to find initial values for the selected solver.
        ux: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `x` values. This should not be specified if the `x`
            argument contains gummys.  If this is specified then only the odr
            solver may be used.  The default is `None`.
        uy: `float`, array_like of `float`  or `None`, optional
            Uncertainty in the `y` values. This should not be specified if the y
            argument contains gummys.  The default is `None`.
        variance_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`.
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
        if x.ndim > 3:
            raise TypeError('the maximum dimension for the x values is 3')
        if _isscalar(deg):
            if deg < 0:
                raise ValueError('deg must be >= 0')
            if x.ndim > 1:
                deg = np.array([int(deg)]*x.ndim)
            else:
                self.deg = int(deg)
            self._order = deg+1
        else:
            if x.ndim == 1:
                raise TypeError('if x is 1-d, then deg must be a scalar')
            try:
                self.deg = np.array([int(d) for d in deg])
                if np.any(self.deg < 0):
                    raise ValueError('all elements of deg must be >= 0')
                self._order = self.deg + 1
            except:
                raise TypeError('deg must be a scalar value or a 1-d array')

        super().__init__(x,y=y,p0=p0,ux=ux,uy=uy,variance_is_known=variance_is_known,xunit=xunit,
                 yunit=yunit,solver=solver,**kw)

    def _get_solver(self,solver,**kw):
        if self.xdim > 1:
            if _isscalar(self.deg) or len(self.deg) != self.xdim:
                raise TypeError('len(deg) != xdim')
                
        if solver is None:
            solver = 'ols'
                
        return super()._get_solver(solver,**kw)
            
    def _shiftx(self,xav):
        # returns a matrix that transforms the fit parameters under: x -> x + xav

        from scipy.special import comb
        
        if self.xdim == 1:
            ret = np.array([[0 if j < i else (-xav)**(j-i)*comb(j,i,exact=True) 
                            for j in range(self._order)] for i in range(self._order)])
        elif self.xdim == 2:
            ret = np.array([[0 if (j < i or l < k) else (-xav[0])**(j-i)*comb(j,i,exact=True)*(-xav[1])**(l-k)*comb(l,k,exact=True)
                            for l in range(self._order[1]) for j in range(self._order[0])] 
                            for k in range(self._order[1]) for i in range(self._order[0])])
        else:
            ret = np.array([[0 if (j < i or l < k) else (-xav[0])**(j-i)*comb(j,i,exact=True)*(-xav[1])**(l-k)*comb(l,k,exact=True)*(-xav[2])**(n-m)*comb(n,m,exact=True)
                            for n in range(self._order[2]) for l in range(self._order[1]) for j in range(self._order[0])] 
                            for m in range(self._order[2]) for k in range(self._order[1]) for i in range(self._order[0])])
        return ret
                                         
    def get_p0(self):
        if self.solver in ('ols','gls'):   
            # Linear fit will be used, no p0 needed.
            if self.xdim == 1:
                return [0]*self._order
            elif self.xdim == 2:
                return [0]*(self._order[0]*self._order[1])
            else:
                return [0]*(self._order[0]*self._order[1]*self._order[2])
        # If a nonlinear fit will be used, get p0 from a linear fit:
        try:
            return PolyFit(self.xf,self.yf,deg=self.deg,solver='ols').pf
        except:
            pass
            
        return [np.mean(self.yf)] + ([0]*np.array(self.deg).sum())
        
    def get_punits(self):
        if self.xdim == 1:
            return [self._yunit*self._xunit**-i for i in range(self._order)]
            
        ret = []
        if self.xdim == 2:
            for j in range(self._order[1]):
                for i in range(self._order[0]):
                    ret.append(self._yunit*self._xunit[0]**-i*self._xunit[1]**-j)
        else:
            for k in range(self._order[2]):
                for j in range(self._order[1]):
                    for i in range(self._order[0]):
                        ret.append(self._yunit*self._xunit[0]**-i*self._xunit[1]**-j*self._xunit[2]**-k)
        
        return ret
        
        
    def f(self,*a):
        if self.xdim == 1:
            x = np.polynomial.polynomial.polyvander(a[0],self.deg)
            y = x@a[1:]
        elif self.xdim == 2:
            x = np.polynomial.polynomial.polyvander2d(a[1],a[0],[self.deg[1],self.deg[0]])
            y = x@a[2:]
        else:
            x = np.polynomial.polynomial.polyvander3d(a[2],a[1],a[0],[self.deg[2],self.deg[1],self.deg[0]])
            y = x@a[3:]
            
        if _isscalar(a[0]):
            return y[0]
        return y
        
    def jacp(self,*a):
        if self.xdim == 1:
            return np.polynomial.polynomial.polyvander(a[0],self.deg).T
        elif self.xdim == 2:
            return np.polynomial.polynomial.polyvander2d(a[1],a[0],[self.deg[1],self.deg[0]]).T
        else:
            return np.polynomial.polynomial.polyvander3d(a[2],a[1],a[0],[self.deg[2],self.deg[1],self.deg[0]]).T
        
    def jacx(self,*a):
        scl = _isscalar(a[0])
        a = [i if _isscalar(i) else np.asarray(i) for i in a]
        
        if self.xdim == 1:
            if scl:
                d0 = 0
            else:
                d0 = np.zeros(a[0].shape)
            for i in range(1,self._order):
                d0 += i*a[i+1]*a[0]**(i-1)
            d0 = [d0]
        elif self.xdim == 2:
            k = 2
            if scl:
                d0 = [0,0]
            else:
                d0 = [np.zeros(a[0].shape),np.zeros(a[0].shape)]
            for i in range(1,self._order[0]):
                for j in range(1,self._order[1]):
                    d0[0] += i*a[k]*a[0]**(i-1)*a[1]**j
                    d0[1] += j*a[k]*a[0]**i*a[1]**(j-1)
                    k += 1
        else:
            k = 3
            if scl:
                d0 = [0,0,0]
            else:
                d0 = [np.zeros(a[0].shape),np.zeros(a[0].shape),np.zeros(a[0].shape)]
            for i in range(1,self._order[0]):
                for j in range(1,self._order[1]):
                    for m in range(1,self._order[2]):
                        d0[0] += i*a[k]*a[0]**(i-1)*a[1]**j*a[2]**m
                        d0[1] += j*a[k]*a[0]**i*a[1]**(j-1)*a[2]**m
                        d0[2] += m*a[k]*a[0]**i*a[1]**j*a[2]**(m-1)
                        k += 1
        
        if scl:
            d0 = [[d] for d in d0]
        return np.array(d0)
    
    def iiijac(self,*a):
        return np.concatenate([self.jacx(*a),self.jacp(*a)])
    
        scl = _isscalar(a[0])
        a = [i if _isscalar(i) else np.asarray(i) for i in a]
        
        if self.xdim == 1:
            if scl:
                d0 = 0
            else:
                d0 = np.zeros(a[0].shape)
            for i in range(1,self._order):
                d0 += i*a[i+1]*a[0]**(i-1)
            d0 = [d0]
            d = np.polynomial.polynomial.polyvander(a[0],self.deg).T
        elif self.xdim == 2:
            k = 2
            if scl:
                d0 = [0,0]
            else:
                d0 = [np.zeros(a[0].shape),np.zeros(a[0].shape)]
            for i in range(1,self._order[0]):
                for j in range(1,self._order[1]):
                    d0[0] += i*a[k]*a[0]**(i-1)*a[1]**j
                    d0[1] += j*a[k]*a[0]**i*a[1]**(j-1)
                    k += 1
            d = np.polynomial.polynomial.polyvander2d(a[1],a[0],[self.deg[1],self.deg[0]]).T
        else:
            k = 3
            if scl:
                d0 = [0,0,0]
            else:
                d0 = [np.zeros(a[0].shape),np.zeros(a[0].shape),np.zeros(a[0].shape)]
            for i in range(1,self._order[0]):
                for j in range(1,self._order[1]):
                    for m in range(1,self._order[2]):
                        d0[0] += i*a[k]*a[0]**(i-1)*a[1]**j*a[2]**m
                        d0[1] += j*a[k]*a[0]**i*a[1]**(j-1)*a[2]**m
                        d0[2] += m*a[k]*a[0]**i*a[1]**j*a[2]**(m-1)
                        k += 1
            d = np.polynomial.polynomial.polyvander3d(a[2],a[1],a[0],[self.deg[2],self.deg[1],self.deg[0]]).T
        
        if scl:
            d0 = [[d] for d in d0]
            return np.concatenate([d0,d]).T[0]
        
        return np.concatenate([d0,d])
        
    def funicode(self):
        if self.xdim == 1:
            txt = 'y = p(1)'
            if self._order >= 2:
                txt += ' + p(2)*x'
            for i in range(2,self._order):
                txt += ' + p('+str(i+1)+')*x**'+str(i)
            return txt
        
        txt = ''
        if self.xdim == 2:
            k = 1
            for j in range(self._order[1]):
                for i in range(self._order[0]):
                    txt += ' + p('+str(k)+')'
                    if i == 1:
                        txt += '*x(1)'
                    elif i > 1:
                        txt += '*x(1)**'+str(i)
                    if j == 1:
                        txt += '*x(2)'
                    elif j > 1:
                        txt += '*x(2)**'+str(j)
                    k += 1
        else:
            k = 1
            for m in range(self._order[2]):
                for j in range(self._order[1]):
                    for i in range(self._order[0]):
                        txt += ' + p('+str(k)+')'
                        if i == 1:
                            txt += '*x(1)'
                        elif i > 1:
                            txt += '*x(1)**'+str(i)
                        if j == 1:
                            txt += '*x(2)'
                        elif j > 1:
                            txt += '*x(2)**'+str(j)
                        if m == 1:
                            txt += '*x(3)'
                        elif m > 1:
                            txt += '*x(3)**'+str(m)
                        k += 1
                
        return 'y = ' + txt[3:]
    
    def flatex(self):
        if self.xdim == 1:
            txt = '$ y = p_{1}'
            if self._order >= 2:
                txt += ' + p_{2} x'
            for i in range(2,self._order):
                txt += ' + p_{'+str(i+1)+'} x^{'+str(i)+'}'
            return txt + ' $'
    
        txt = ''
        if self.xdim == 2:
            k = 1
            for j in range(self._order[1]):
                for i in range(self._order[0]):
                    txt += ' + p_{'+str(k)+'}'
                    if i == 1:
                        txt += ' x_{1}'
                    elif i > 1:
                        txt += ' x_{1}^{'+str(i)+'}'
                    if j == 1:
                        txt += ' x_{2}'
                    elif j> 1:
                        txt += ' x_{2}^{'+str(j)+'}'
                    k += 1
        else:
            for m in range(self._order[2]):
                for j in range(self._order[1]):
                    for i in range(self._order[0]):
                        txt += ' + p_{'+str(k)+'}'
                        if i == 1:
                            txt += ' x_{1}'
                        elif i > 1:
                            txt += ' x_{1}^{'+str(i)+'}'
                        if j == 1:
                            txt += ' x_{2}'
                        elif j > 1:
                            txt += ' x_{2}^{'+str(j)+'}'
                        if m == 1:
                            txt += ' x_{3}'
                        elif m > 1:
                            txt += ' x_{3}^{'+str(m)+'}'
                        k += 1
                            
        return '$ y = ' + txt[3:] + ' $'
    
    def fhtml(self):
        if self.xdim == 1:
            txt = '<i>y</i> = <i>p</i><sub>1</sub>'
            if self._order >= 2:
                txt += ' + <i>p</i><sub>2</sub> <i>x</i>'
            for i in range(2,self._order):
                txt += ' + <i>p</i><sub>'+str(i+1)+'</sub> <i>x</i><sup>'+str(i)+'</sup>'
            return txt
        
        txt = ''
        if self.xdim == 2:
            k = 1
            for j in range(self._order[1]):
                for i in range(self._order[0]):
                    txt += ' + <i>p</i><sub>'+str(k)+'</sub>'
                    if i == 1:
                        txt += ' <i>x</i><sub>1</sub>'
                    elif i > 1:
                        txt += ' <i>x</i><sub>1</sub><sup>'+str(i)+'</sup>'
                    if j == 1:
                        txt += ' <i>x</i><sub>2</sub>'
                    elif j> 1:
                        txt += ' <i>x</i><sub>2</sub><sup>'+str(j)+'</sup>'
                    k += 1
        else:
            k = 1
            for m in range(self._order[2]):
                for j in range(self._order[1]):
                    for i in range(self._order[0]):
                        txt += ' + <i>p</i><sub>'+str(k)+'</sub>'
                        if i == 1:
                            txt += ' <i>x</i><sub>1</sub>'
                        elif i > 1:
                            txt += ' <i>x</i><sub>1</sub><sup>'+str(i)+'</sup>'
                        if j == 1:
                            txt += ' <i>x</i><sub>2</sub>'
                        elif j> 1:
                            txt += ' <i>x</i><sub>2</sub><sup>'+str(j)+'</sup>'
                        if m == 1:
                            txt += ' <i>x</i><sub>3</sub>'
                        elif m > 1:
                            txt += ' <i>x</i><sub>3</sub><sup>'+str(m)+'</sup>'
                        k += 1
                
        return '<i>y</i> = ' + txt[3:]


class DoubleExpFit(Fit):
    """
    DoubleExpFit(x,y,p0=None,ux=None,uy=None,variance_is_known=True,xunit=None, yunit=None,
             solver=None,**keywords)
   
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
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
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
        return '$ y = p_{1}\\exp(-x/p_{2}) + p_{3}\\exp(-x/p_{4}) + p_{5} $'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub> exp(-<i>x</i>/<i>p</i><sub>2</sub>) + <i>p</i><sub>3</sub> exp(-<i>x</i>/<i>p</i><sub>4</sub>) + <i>p</i><sub>5</sub>'
        
        
class ExpFit(Fit):
    """
    ExpFit(x,y,p0=None,ux=None,uy=None,variance_is_known=True,xunit=None, yunit=None,
             solver=None,**keywords)
   
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
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
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
        return '$ y = p_{1}\\exp(x/p_{2}) + p_{3} $'
    
    def fhtml(self):
        return '<i>y</i> = <i>p</i><sub>1</sub> exp(<i>x</i>/<i>p</i><sub>2</sub>) + <i>p</i><sub>3</sub>'
        
        
        
class OneOverTFit(Fit):
    """
    DoubleExpFit(x,y,p0=None,ux=None,uy=None,variance_is_known=True,xunit=None, yunit=None,
             solver=None,**keywords)
   
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
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
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
    SinFit(x,y,p0=None,ux=None,uy=None,variance_is_known=True,xunit=None, yunit=None,
             solver=None,**keywords)
   
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
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
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
        
        