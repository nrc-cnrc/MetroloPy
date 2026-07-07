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

from ._ummy import njacobian
from .dof import DoF,_DoF_inf
from .util import _isscalar
from ._gummy import gummy
from .exceptions import FitWarning
from .printing import PrettyPrinter
from .abc import UncertainValue,AbcQuantity

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
            jx = np.reshape(jx,(count,xdim),order='F').T
        else:
            jx = np.reshape(jx,(count,xdim,ydim),order='F').T
    
    return (jp,jx)

def getvalues(txt,x,u,cov,w,units,ignore_corr):
    # unpack Fit x or y values
    
    from scipy.linalg import pinvh
    if x is None:
        return None,None,None,None,None,None,False,1,0
    elif _isscalar(x):
        if isinstance(x,UncertainValue):
            pcorr = True
            xf = float(x.x)
            u = x.u
        else:
            pcorr = False
            xf = float(x)
            u = None
        if units is None:
            if isinstance(x,AbcQuantity) and not x.unit_is_one:
                units = x.unit
            else:
                units = 1
        else:
            if units != 1:
                from ._unit import Unit
                units = Unit.unit(units)
        return x,xf,u,None,w,None,pcorr,units,0
    
    x = np.asarray(x)
    if u is not None and not _isscalar(u):
        u = np.asarray(u)
    if cov is not None:
        cov = np.asarray(cov)
    if w is not None and not _isscalar(w):
        w = np.asarray(w)
        
    shape = np.shape(x)
    if np.ndim(x) == 1:
        dim = 1
    else:
        dim = shape[0]
        
    
    if units is not None:
        if _isscalar(units):
            if units == 1:
                units = [1]*dim
            else:
                from ._unit import Unit
                units = [Unit.unit(units)]*dim
        else:
            if np.shape(units) != (dim,):
                raise TypeError(txt + 'units must be either None, a single unit or a list of units with length equal to xdim')
            if any([i != 1 for i in units]):
                from ._unit import Unit
                units = [Unit.unit(u) for u in units]
        if dim == 1:
            units = units[0]
        
    pcorr = False
    if x.dtype is np.dtype('O'):
        xx = np.empty(shape,dtype=np.dtype('O'))
        xf = np.empty(shape)
        uf = np.empty(shape)

        if dim == 1:
            units = 1
            if isinstance(x[0],AbcQuantity) and not x[0].unit_is_one:
                units = x[0].unit
        if units is None:
            units = [i[0].unit if (isinstance(i[0],AbcQuantity) and not i[0].unit_is_one) else 1 for i in x]
            
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
                if iunit == 1:
                    ii = ii.convert(1) if isinstance(i,AbcQuantity) else ii
                else:
                    if isinstance(ii,AbcQuantity):
                        if ii.unit_is_one:
                            ii *= iunit
                        else:
                            ii = ii.convert(iunit)
                    else:
                        ii = gummy(ii,unit=iunit)
                xx[it.multi_index] = ii
                ii = ii.value if isinstance(ii,AbcQuantity) else ii
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
                units = 1
            else:
                units = [1]*dim
            
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
            #if not np.any(np.any(cov[i] - np.diag((u.T[i])**2) for i in range(len(cov)))):
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
                    w = np.array([pinvh(i) for i in cov]).T
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
    if np.ndim(w) < 2 or np.shape(w)[0] != np.shape(w)[1]:
        return np.sqrt(w)
    elif np.ndim(w) == 2:
        return cholesky(w)
    else:
        return cholesky(w.T).T

def mmul(m,v):
    if np.ndim(m) > 1 and np.shape(m)[0] == 1:
        return mmul(m[0],v)
    if np.ndim(m) >= 2 and np.shape(m)[-1] == np.shape(m)[-2] == np.shape(v)[-1]:
        return m@v
    return m*v
    
def repl(a,p,fix):
    p = np.array(p)
    p[~fix] = a
    return p

class _Fit:
    # Base class for fitting.  This implements plotting and unpacks any
    # uncertainty or unit information in the data.
    
<<<<<<< HEAD
    plot_points = 500
=======
    plot_points = 100
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    over_plot = 0.05
    
    xlabel = None
    ylabel = None
    
    def __init__(self,x,y,ux=None,uy=None,xweights=None,weights=None,xcov=None,
                 ycov=None,variance_is_known=None,xunit=None,yunit=None,ignore_correlations=False):
        self.variance_is_known = variance_is_known
            
        if isinstance(x,np.ma.MaskedArray) or isinstance(y,np.ma.MaskedArray):
            if not isinstance(x,np.ma.MaskedArray):
                x = np.ma.array(x)
                
            if x.ndim > 1:
                mask = x[0].mask
                for i in  range(1,x.shape[0]):
                    mask = np.ma.mask_or(mask,x[i].mask)
            else:
                mask = x.mask
                
            if y is not None:
                if not isinstance(y,np.ma.MaskedArray):
                    y = np.ma.array(y)
                    
                if y.ndim > 1:
                    for i in  range(y.shape[0]):
                        mask = np.ma.mask_or(mask,y[i].mask)
                else:
                    mask = np.ma.mask_or(mask,y.mask)
            
            x.mask = mask
            x = x[~x.mask]
            
            if y is not None:
                y.mask = mask
                y = y[~y.mask]
                
        self.x,self.xf,self.ux,self.xcov,self.xweights,self.known_xvar,self.x_is_gummies,self._xunit,self.xdim = getvalues('x',x,ux,xcov,xweights,xunit,ignore_correlations)

        self.y,self.yf,self.uy,self.ycov,self.weights,self.known_yvar,self.y_is_gummies,self._yunit,self.ydim = getvalues('y',y,uy,ycov,weights,yunit,ignore_correlations)
            
        self.sqrt_weights = wsqrt(self.weights)
        self.sqrt_xweights = wsqrt(self.xweights)
            
        if self.known_yvar is None and self.known_xvar is None:
            if self.variance_is_known:
                warn('variance_is_known is set to True, but no uncertainties have been defined for the x or y values')
            self.variance_is_known = False
        elif self.variance_is_known is None:
            self.variance_is_known = True
                       
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
<<<<<<< HEAD
             xmin=None,xmax=None,xlabel=None,ylabel=None,title=None,hold=False,
             plot_points=None,fig_options={},subplot_options={}):
=======
             xmin=None,xmax=None,xlabel=None,ylabel=None,hold=False,
             plot_points=None):
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        """
        Plots the data points, fitted curve, as well as confidence limits and 
        control limits around the fitted curve.
        
        
        Parameters
        ----------
<<<<<<< HEAD
        show_data: `bool`, optional
            Whether or not to plot the data points. The default is `True`.
            
        show_fit: `bool`, optional
            Whether or not to plot the fitted curve. The default is `True`.
            
        error_bars: `bool`, optional
            Whether or not to plot error bars on the data points (if uncertainty
            values were defined for the data).  The default is `True`
            
        cik and cip: `Number`, optional
            Specifying a value for cip or cik will add uncertainty bands for the
            fitted curve to the plot.  These values give the coverage factor 
            (cik) or confidence level (cip) for the uncertainty bands in the 
            plot.  cik should be a float or int > 0 or cip should be a float 
            between 0 and 1.  Do not specify both `cik` and `cip`.
            
        clk and clp: `Number`, optional
            Specifying a value for clp or clk will add control limit bands for 
            the fitted curve to the plot.  The control limit band is the region
            where a new data point is expected to lie.  These values give the 
            coverage factor (cik) or confidence level (cip) for the conrol limit 
            in the plot.  clk should be a float or int > 0 or clp should be a 
            float between 0 and 1.  Do not specify both `clk` and `clp`.
            
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        xmin and xmax: `float`, optional
            The lower and upper limits of the fitted, confidence interval
            and control limit curves.  If this is None, the limits are equal
            to x1 +/- (x2 - x1)*`Fit.over_plot` where x1 is the x value of the
            first data point, x2 is the x value of the last data  point and
            `Fit.over_plot` is an attribute of the `Fit` object with  default
            value 0.05.
<<<<<<< HEAD
            
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        xlabel and ylabel:  `str`, optional
            Labels for the x and y axes. If units are defined for the x or y axes,
            the unit symbol will be added to the end of the labels defined here.
            If these are set to `None`, then the values of the `Fit.xlabel` and
            `Fit.ylabel` attributes will be used.  The default is `None`.
<<<<<<< HEAD
            
        hold: `bool`, optional
            If hold is `False` then ``pyplot.show()`` is executed just before this
            function returns.  The detault is `False`
            
        error_bar_k: `int` or `float`, optional
            The length of the error bars are determined by multiplying the
            uncertainty for each data point by this quantity. The default value
            is 1.
                
        data_format: `str`, optional
            The format string passed to `pyplot.plot` or `pyplot.errorbar` when
            plotting the data points.  The default is 'ko'.
            
        data_options: `dict`, optional
            A dictionary containing key words that are passed to `pyplot.plot` or
            `pyplot.errorbar` when plotting the data points.
            
        fit_format: `str`, optional
            The format string passed to `pyplot.plot` or` pyplot.errorbar` when
            plotting the fitted curve.  The default is 'k-'.
            
        fit_options: `dict`, optional
            A dictionary containing key words that are passed to `pyplot.plot`
            or `pyplot.errorbar` when plotting the fitted curve.
            
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        plot_points:  `int`, optional
            The number of points to use in each curve when plotting the fit,
            confidence interval, and control limit curves.  If this is set to
            `None`, then the value of the `Fit.plot_points` attribute will be used,
            which has a default value of 100.
<<<<<<< HEAD
            
        ciformat: `str`, optional
            Format string passes to the `pyplot.plot` command that plots the
            uncertainty bands. The default is 'g-'.
            
        cioptions:  `dict`, optional
            Keywork options passed to the `pyplot.plot` command that plots the
            uncertainty bands.
            
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        clk,clp,clformat, cloptions:  optional
            Control limit options, same as above for the uncertainty bands.  The
            control limit band if the control limit k factor multiplied by the
            RSS of the fit uncertainty and the standard deviation of the residuals.
<<<<<<< HEAD
            
        fig_options: `dict`, optional
            keywords passed to `pyplot.figure` when creating the figure
            
        subplot_options: `dict`, options
            keywords passed to `pyplot.figure.add_subplot` when creating the
            subplot
            
        Returns
        -------
        Figure, Axes
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
        
<<<<<<< HEAD
        if self.xunit == 1:
            xs = None
        else:
            xs = self.xunit.tostring(fmt='latex')
        xlabel = gummy._plotlabel(xlabel,symbol=xs)
            
        if self.yunit == 1:
            ys = None
        else:
            ys = self.yunit.tostring(fmt='latex')
        ylabel = gummy._plotlabel(ylabel,symbol=ys)
        
        if 'ylabel' not in subplot_options and ylabel is not None:
            subplot_options['ylabel'] = ylabel
            
        if 'xlabel' not in subplot_options and xlabel is not None:
            subplot_options['xlabel'] = xlabel
            
        if 'title' not in subplot_options and title is not None:
            subplot_options['title'] = title
            
        fig = plt.figure(**fig_options)
        ax = fig.add_subplot(**subplot_options)
        
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
<<<<<<< HEAD
                ax.errorbar(xf,yf,xerr=ux,yerr=uy,fmt=data_format,
                               **data_options)
            else:
                if data_format is None:
                    ax.plot(xf,yf,**data_options)
                else:
                    ax.plot(xf,yf,data_format,**data_options)            
=======
                plt.errorbar(xf,yf,xerr=ux,yerr=uy,fmt=data_format,
                               **data_options)
            else:
                if data_format is None:
                    plt.plot(xf,yf,**data_options)
                else:
                    plt.plot(xf,yf,data_format,**data_options)            
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            
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
<<<<<<< HEAD
                ax.plot(fx,fy,**fit_options)
            else:
                ax.plot(fx,fy,fit_format,**fit_options)
=======
                plt.plot(fx,fy,**fit_options)
            else:
                plt.plot(fx,fy,fit_format,**fit_options)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        
        if cik is not None:
            u = np.array([self.ypred(x).u for x in fx])
            up = fy + u*cik
            un = fy - u*cik
            if ciformat is None:
<<<<<<< HEAD
                ax.plot(fx,up,**cioptions)
                ax.plot(fx,un,**cioptions)
            else:
                ax.plot(fx,up, ciformat, **cioptions)
                ax.plot(fx,un, ciformat, **cioptions)
=======
                plt.plot(fx,up,**cioptions)
                plt.plot(fx,un,**cioptions)
            else:
                plt.plot(fx,up, ciformat, **cioptions)
                plt.plot(fx,un, ciformat, **cioptions)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        
        if clk is not None:
            u = np.array([self.control_limit(x, clk) for x in fx])
            upl = fy + u
            unl = fy - u
            if clformat is None:
<<<<<<< HEAD
                ax.plot(fx,upl,**cloptions)
                ax.plot(fx,unl,**cloptions)
            else:
                ax.plot(fx,upl,clformat,**cloptions)
                ax.plot(fx,unl,clformat,**cloptions)
=======
                plt.plot(fx,upl,**cloptions)
                plt.plot(fx,unl,**cloptions)
            else:
                plt.plot(fx,upl,clformat,**cloptions)
                plt.plot(fx,unl,clformat,**cloptions)
            
        if self.xunit == 1:
            xs = None
        else:
            xs = self.xunit.tostring(fmt='latex')
        xlabel = gummy._plotlabel(xlabel,symbol=xs)
        if xlabel is not None:
            plt.xlabel(xlabel)
            
        if self.yunit == 1:
            ys = None
        else:
            ys = self.yunit.tostring(fmt='latex')
        ylabel = gummy._plotlabel(ylabel,symbol=ys)
        if ylabel is not None:
            plt.ylabel(ylabel)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            
        if not hold:
            plt.show()
            
<<<<<<< HEAD
        return fig,ax
            
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    def control_limit(self,x,k=2):
        if self.known_var is None:
            s = self.s
        else:
            s = np.sqrt(self.known_var)
        return k*(np.sqrt(self.ypred(x).u**2 + s**2))
        
    
class Fit(_Fit,PrettyPrinter):
    
    latex_math = None
    
    def __init__(self,x,y=None,f=None,p0=None,ux=None,uy=None,
                 variance_is_known=None,xunit=None, yunit=None,solver=None,
                 xweights=None,weights=None,xcov=None,ycov=None,
                 ignore_correlations=False,fix=None, fargs=[],fkwds={},**kw):
        """
<<<<<<< HEAD
        Performs a least squares fit.  The function may be passed in the arguments
=======
        Performs a non-linear fit.  The function may be passed in the arguments
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        or may be specified by overriding the Fit.f(...) method in a subclass.
        The fitting is performed as soon as the instance is created.
        
        Parameters
        ----------
<<<<<<< HEAD
        x: array
            The indepenant variables x.  For n data points this should be an
            array of shape (n,) for a one-dimensional fit or shape (m,n) for an
            m-dimensional fit.
            
            The x-values may be all `float` or all `gummy`,  If uncertainties
            u are defined for the x-values, the xweights will be set to 1/u**2
            for each point.  To override this behavior set the keyword parameter
            "xweights = 1".  If uncertainties are defined for the x-values,
            the default solver is odr.  If the nls or ols solver is manually
            selected then the uncertainties in the x-values will be ignored.
            For an m-dimential fit, the odr solver will take into account 
            correlations between the m elements of the x-values for each data 
            point, but ignores correlations between different data points and 
            bewteen the x- and y-values.
            
        y:  array or `Number`, optional
            Response variables y. Usually this is an array of shape (n,) but 
            this can be omitted or a scalar `Number` if the nls or odr solvers 
            are used or an array of shape (k,n) if the odr solver is used.  The
            shape must match the shape of the array returned by f.
            
            The y-values may be all `float` or all `gummy`,  If uncertainties
            u are defined for the y-values, the weights will be set to 1/u**2
            for each point.  To override this behavior set the keyword parameter
            "weights = 1".
            
            If the nls or odr solver is used, the weighting will take into
            account correlations between the different points.  However, the odr
            solver ignores correlations between different data points.  If the
            y-values are multi-dimensional, the odr solver will take into
            account correlations between the different elements at each point.
            
        f:  function or method
            The fit function.  This function can be passed as an argument to
            the fit initializer or it can be provided by implementing the
            f method in a derived class.
            
            For a one-dimensional fit witn N parameters,
            this function will be called as f(x,p1,p2,...,pN,*fargs,**fkwds)
            where x is usually a numpy float array with shape (n,) The fit 
            parameters p1,p2,...,pN will always be scalar values.  For an M 
            dimensional fit the function will be called as 
            f(x1,x2,...,xM,p1,p2,...,pN,*fargs,**fkwds) and the x1,x2,..,xN all
            have shape (n,).
            
            Usually f will return an array of shape (n,), but if the nls sovler
            is used, the array length need not be n. 
            
            (The one exception where the function will be called with scalar 
            x-values is when the ols solver is used with no jacp defined; the
            function will be called with scalar values for x when numerically 
            calculating the Jacobian.) 
            
        p0: array of `float`
            The initial values for the fit parameters.  This parameter is 
            required unless the get_p0 method is overridden in a subclass to 
            provide an initial esitmate for the fit parameters.  If get_p0
            is implemented, then some elements of p0 can be set to `None` and
            they will be replaced with the estimated value from get_p0.
            
        fix: array of `bool`
            A mask for the fit parameters.  For any element in `fix` that is
            `True`, the corresponding fit parameter will be held constant at
            its initial value.
            
        solver:  {'ols','nls','odr'}, optional
            A non-linear least squares (nls), ordinary least squares (ols) or
            othogonal distance regression solver (odr) can be used to perform
            the fit.
            
            nls is the default solver unless uncertainties are specified for
            the x-values, in which case the default solver is odr.
            
            The nls solver uses `scipy.optimize.least_squares` while the odr
            solver used `odrpack.odr_fit`.  The ols solver uses 
            `scipy.linag.pinvh` to invert the normal equation matrix.
            
        ux: array or `float`, optional
            Uncertainties for the x-values.  Specifying ux is an alternative to 
            passing gummys with uncerainties defined as the the x-values.  ux 
            can be a number that applies to all the x-values or an array giving
            the uncertainty for each x-value.
            
            If uncertainties u are defined individually for the x-values, the 
            xweights will be set to 1/u**2 for each point.  To override this 
            behavior set the keyword parameter "xweights = 1".  If 
            uncertainties are defined for the x-values, the default solver is 
            odr.  If the nls or ols solver is manually selected then the 
            uncertainties in the x-values will be ignored.
            
        uy: array or `float`, optional
            Uncertainties for the y-values.  Specifying uy is an alternative to 
            passing gummys with uncerainties defined as the the y-values.  uy 
            can be a number that applies to all the x-values or an array giving
            the uncertainty for each y-value.
            
            If uncertainties u are defined individually for the y-values, the 
            weights will be set to 1/u**2 for each point.  To override this 
            behavior set the keyword parameter "weights = 1".
            
        jacp: function, optional
            The Jacobian of the fit function with respect to the fit parameters.
            If this is not provided the Jacobian will be calculated numerically.
            
            This function can be passed as an argument to the fit initializer 
            or it can be provided by implementing the jacp method in a derived 
            class.
            
            If provided jacp will be called as jacp(x,p1,p2,...,pN) if the
            x-values are a one-dimensional array and 
            jacp(x1,x2,...,xM,p1,p2,...,pN) if the x-values are M dimensional.
            For n data points and N fit parameters, jacp should return an
            array of shape (N,n).
            
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        variance_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`.
<<<<<<< HEAD
            
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        xunits, yunits: `str` or `None`, optional
            units for the x and y coordinates. These should not be specified
            if the `x` and `y` parameters contain gummys. These may only be
            specified if the `get_punits` method is overridden in a subclass.
<<<<<<< HEAD
            
        fargs:  `list`, optional
            Additional (fixed) arguments passed to the fit function after the 
            fit parameters.
            
        fkwds:  `list`, optional
            Fixed keyword arguments passed to the fit function.
        
=======
        solver:  {'ols','nls','odr'}, optional
            If this is 'nls' then `scipy.optimize.least_squares` is used to perform
            the fit.  If it is 'odr' then `odrpack.odr_fit` is used.  'nls' may 
            not be used if the y-coordinate is `None` or multi-dimensional or if
            there is uncertainty in the x-coordinates.  If this is `None`,
            then 'nls' will be used when possible.
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        other keywords:  optional
            Any additional keyword parameters will be passed to the solver.
        
        Attributes
        ----------
<<<<<<< HEAD
        p:  `numpy.array` of `gummy`
            The fitted values for the fit function parameters as gummys
            including uncertainties and units.
            
        pf:  `numpy.array` of `float`
            The fitted values for the fit function parameters as floats
            
        res:  `numpy.ndarray` of `float`
            the weighted fit residuals
            
        var:  `float`
            the variance of the weighted fit residuals
            
        s:  `float`
            the standard deviation of weighted the residuals
            
        cov:  `numpy.ndarray` of `float`
            the covariance matrix of the (non-fixed) parameters
            
        known_var:  `float`
            If uncertainties are defined for the y-values (and possible the
            x-values), this is the predicted variance of the weighted fit
            residuals.  If uncertainites are not defined for the y- or x-values,
            this is `None`.
            
        yvar:  `float`
            the variance of the residual differences between the fitted values
            and the weighted y-values. For the nls and ols solver this is the 
            same as var.
        
        xvar:  `float`
            the variance of the residual differences between the fitted values
            and the weighted x-values. For the nls and ols solver this is `None`.
            
        rcov:  `numpy.ndarray` of `float`
            the covariance matrix of the (non-fixed) parameters divided by the
            variance
        
        fit_output:
            the return value of `scipy.optimize.least_squares` for the nls
            solver of the return value of `odr_fit.odrpack` for the odr
            solver.  This is `None` if the ols solver is used.
            
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
            the number of (x) data points
            
        p0:  `list` of `float`
            the initial values for the fit function parameters
            
        solver:  `str`
            the solver used
            
        punits:  `list` of `Unit`
            the units of the fit parameters
            
        nparam:  `int`
            the number of (non-fixed) fit parameters
            
        dof: `float`
            degrees of freedom for the fit
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5

        Methods
        -------
        ypred(x1,x2,...):
<<<<<<< HEAD
            Takes `xdim` floats (or arrays of float) and returns a gummy 
            representing the predicted value(s) at that x-coordinate.
            
        ypredf(x1,x2,...):
            Takes `xdim` floats (or arrays of float) and returns a float giving 
            the  predicted value(s) at that x-coordinate.
            
        plot(...):
            plots the data (only available if x and y are one-dimensional)
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
        self.p_names = None
        self.p_names_html = None
        self.p_names_latex = None
        self.p_names_ascii = None
        
        if f is not None:
            self.f = f
            
        if 'jacp' in kw:
            self.jacp = kw.pop('jacp')
            
        if 'jacx' in kw:
            self.jacx = kw.pop('jacx')
            
        if 'p_names' in kw:
            self.p_names = kw.pop('p_names')
        if 'p_names_html' in kw:
            self.p_names_html = kw.pop('p_names_html')
        if 'p_names_latex' in kw:
            self.p_names_latex = kw.pop('p_names_latex')
        if 'p_names_ascii' in kw:
            self.p_names_html = kw.pop('p_names_ascii')
        
        if p0 is None:
            try:
                p0 = self.get_p0()
            except NotImplementedError:
                raise ValueError('initial values p0 must be specified for this fit')
<<<<<<< HEAD
        else:
            if np.any([i is None for i in p0]):
                try:
                    gp0 = self.get_p0()
                except NotImplementedError:
                    raise ValueError('all the initial values p0 must be specified for this fit (None is not allowed in p0)')
                    
                if _isscalar(p0) or len(gp0) < len(p0):
                    raise TypeError('p0 must be a list with length equal to the number of fit parameters')
                    
                p0 = [gp0[i] if p0[i] is None else p0[i] for i in range(len(p0))]
                
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        self.p0 = np.asarray(p0,dtype=float)

        getpu = False
        if self.xdim == 1:
            getpu = getpu or (self.xunit != 1)
        else:
            for u in self.xunit:
                getpu = getpu or (u != 1)
        if self.ydim <= 1:
            getpu = getpu or (self.yunit is not None and self.yunit != 1)
        else:
            for u in self.yunit:
                getpu = getpu or (u is not None and u != 1)
        if getpu:
            try:
                self.punits = self.get_punits()
            except NotImplementedError:
                raise ValueError('only dimensionless quanities are accepted for this fit')
        else:
            self.punits = [1]*len(self.p0)
            
            
        self.count = None
        yshape = np.shape(self.yf)
        if len(yshape) > 0:
            self.count = yshape[-1]
            
        try:
            # See if f will broadcast properly across the xf array...
            if self.xdim == 1:
                self.f0 = self.f(self.xf,*self.p0,*fargs,**fkwds)
            else:
                self.f0 = self.f(*self.xf,*self.p0,*fargs,**fkwds)
            
            if self.ydim > 0 and np.shape(self.f0)[-1] != self.count:
                raise TypeError
                
            #if not isinstance(self.f0,np.ndarray):
                #self.f = lambda *x: np.array(self.f(*x))
                #self.f0 = np.array(self.f0)
                
        except NotImplementedError:
            raise TypeError('the fit function f has not been sepecified')
<<<<<<< HEAD
        except Exception as e:
=======
        except:
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            try:
                if self.count is None:
                    raise
                self.f0 = None
                # ...if not vectorize it so that it does.
                if self.ydim <= 1:
                    ot = [np.float64]
                else:
                    ot = [np.float64]*self.ydim
                self.f = np.vectorize(self.f,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))
            except:
<<<<<<< HEAD
                raise e
=======
                raise TypeError('the calling the fit function with the initial parameters raises an error')
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            
        jp = True
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
            jp = False
        except:
            raise
            # ...if not vectorize it so that it does.
            ot = [np.float64]*self.nparam
            self.jacp = np.vectorize(self.jacp,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))
        
        jx = True
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
            jx = False
        except:
            # ...if not vectorize it so that it does.
            ot = [np.float64]*self.xdim
            self.jacx = np.vectorize(self.jacx,ot)#,excluded=list(range(self.xdim,self.xdim+self.nparam)))

        self._jacp = None
        self._jacx = None
        
        if fix is None:
            if self.xdim == 1:
<<<<<<< HEAD
                self._f = lambda x,p:np.asarray(self.f(x,*p,*fargs,**fkwds))
                if jp:
                    self._jacp = lambda x,p:np.asarray(self.jacp(x,*p))
                if jx:
                    self._jacx = lambda x,p:np.asarray(self.jacx(x,*p))
            else:
                self._f = lambda x,p:np.asarray(self.f(*x,*p,*fargs,**fkwds))
                if jp:
                    self._jacp = lambda x,p:np.asarray(self.jacp(*x,*p))
                if jx:
                    self._jacx = lambda x,p:np.asarray(self.jacx(*x,*p))
=======
                self._f = lambda x,p:self.f(x,*p,*fargs,**fkwds)
                if jp:
                    self._jacp = lambda x,p:self.jacp(x,*p)
                if jx:
                    self._jacx = lambda x,p:self.jacx(x,*p)
            else:
                self._f = lambda x,p:self.f(*x,*p,*fargs,**fkwds)
                if jp:
                    self._jacp = lambda x,p:self.jacp(*x,*p)
                if jx:
                    self._jacx = lambda x,p:self.jacx(*x,*p)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            self._p0 = self.p0
        else:
            fix = np.asarray(fix,dtype=bool)
            self._p0 = self.p0[~fix]
            if self.xdim == 1:
<<<<<<< HEAD
                self._f = lambda x,p:np.asarray(self.f(x,*repl(p,self.p0,fix),*fargs,**fkwds))
                if jp:
                    self._jacp = lambda x,p:np.asarray(self.jacp(x,*repl(p,self.p0,fix))[~fix])
                if jx:
                    self._jacx = lambda x,p:np.asarray(self.jacx(x,*repl(p,self.p0,fix)))
            else:
                self._f = lambda x,p:np.asarray(self.f(*x,*repl(p,self.p0,fix),*fargs,**fkwds))
                if jp:
                    self._jacp = lambda x,p:np.asarray(self.jacp(*x,*repl(p,self.p0,fix))[~fix])
                if jx:
                    self._jacx = lambda x,p:np.asarray(self.jacx(*x,*repl(p,self.p0,fix)))
=======
                self._f = lambda x,p:self.f(x,*repl(p,self.p0,fix),*fargs,**fkwds)
                if jp:
                    self._jacp = lambda x,p:self.jacp(x,*repl(p,self.p0,fix))[~fix]
                if jx:
                    self._jacx = lambda x,p:self.jacx(x,*repl(p,self.p0,fix))
            else:
                self._f = lambda x,p:self.f(*x,*repl(p,self.p0,fix),*fargs,**fkwds)
                if jp:
                    self._jacp = lambda x,p:self.jacp(*x,*repl(p,self.p0,fix))[~fix]
                if jx:
                    self._jacx = lambda x,p:self.jacx(*x,*repl(p,self.p0,fix))
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        self.fix = fix
                    
        self.nparam = len(self._p0)
        
        if self.f0 is None:
            self.f0 = self._f(self.xf,self._p0)
        if self.count is None:
            self.count = np.shape(self.f0)[-1]
        else:
            if np.shape(self.f0) != yshape:
                raise TypeError('the fit function does not return an array the same shape as y')
            
        self._solver(**kw)
        
        for p in self.p:
            if not p.unit_is_one and not p.unit.linear:
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
            
    @property
    def njacx(self):
        if self._njacx is None:
            try:
                if self.xdim <= 1:
                    self._njacx = mmul(self.sqrt_weights,self.jacx(self.xf,*self._pf))
                else:
                    self._njacx = mmul(self.sqrt_weights,self.jacx(*self.xf,*self._pf))
    
            except NotImplementedError:
                # numerically calculate the jacobian of the fit function with respect to x
                f = self.f
                x = self.x
                xf = self.xf
                pf = self._pf
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
                        self._njacx = (w*r)
                    else:
                        self._njacx = (w@r)
        return self._njacx
        
    @njacx.setter
    def njacx(self,j):
        self._njacx = j
            
    def _create_p(self):
        # Given fit parameters pf (float array) and the covariance matrix cov for the
        # parameters, return correlated gummys representing the fit parameters. jac,
        # if provided, should be the jacobian mutiplied by any weighting of the data
        # points.
    
        self.s = np.sqrt(self.var)
        self.dof = self.count - len(self._pf)
        
        if self.variance_is_known and self.known_var is None:
            self.variance_is_known = False
            
        if self.variance_is_known:
            dof = _DoF_inf
        else:
            dof = DoF(self.dof)
            
        if self.variance_is_known:
            self.cov = self.rcov*self.known_var
        elif self.var is not None:
            self.cov = self.rcov*self.var
        else:
            self.cov = None

        if self.y_is_gummies or self.x_is_gummies or self._ch_create_p:
            odr = self.solver.startswith('odr') and self.xvar > 0
            
            uy = self.uy
<<<<<<< HEAD
            if self.uy is None:
                if self.yvar is None:
                    yvar = self.var
                else:
                    yvar = self.yvar
                if yvar is not None:
                    uy = np.sqrt(yvar)
                
            if not self.y_is_gummies and uy is not None:
=======
            if self.yvar is None:
                yvar = self.var
            else:
                yvar = self.yvar
            if yvar is not None:
                uy = np.sqrt(yvar)
                
            
            y = self.y
            if not self.y_is_gummies and uy is not None:
                if self.ydim == 0:
                    y = [self.y]
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
                y = np.empty(np.shape(self.yf),dtype=np.dtype('O'))
                scl = _isscalar(uy)
                rd = not scl and np.shape(uy) == (self.ydim,)
                with np.nditer(y,flags=['multi_index','refs_ok'],op_flags=['readwrite']) as it:
                    for i in it:
                        if scl:
                            iuy = uy
                        elif rd:
                            iuy = uy[it.multi_index[0]]
                        else:
                            iuy = uy[it.multi_index]
                        i[...] = gummy(0,iuy,dof=dof)
<<<<<<< HEAD
            elif self.y_is_gummies:
                y = np.array([i/i.unit for i in self.y])
            else:
                y = self.y
=======
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
                        
            if self.weights is not None:
                y = mmul(self.sqrt_weights,y)
                            
<<<<<<< HEAD
            if odr:
                if not self.variance_is_known:
                    ux = np.sqrt(self.xvar)
                else:
                    ux = self.ux
            else:
                ux = None
=======
            if odr and not self.variance_is_known:
                ux = np.sqrt(self.xvar)
            else:
                ux = self.ux
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    
            if odr:
                if self.xdim == 1:
                    cjx = np.cos(np.arctan((self.njacx)[0]/self.sqrt_xweights))
                    y = cjx*y
                else:
                    cjx = np.cos(np.arctan((self.njacx)/self.sqrt_xweights))
                    y = np.sqrt(np.sum(cjx)**2,axis=0)*y
                
            p = self.rcov@(y@self.njacp).T
    
<<<<<<< HEAD
            if ux is not None:
                if not self.x_is_gummies:
=======
            if ux is not None and self.x is not None:
                x = self.x
                if not self.x_is_gummies and ux is not None and x is not None:
                    if self.xdim == 0:
                        x = [x]
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
                    x = np.empty(np.shape(self.xf),dtype=np.dtype('O'))
                    scl = _isscalar(ux)
                    rd = not scl and np.shape(ux) == (self.xdim,)
                    with np.nditer(x,flags=['multi_index','refs_ok'],op_flags=['readwrite']) as it:
                        for i in it:
                            if scl:
                                iux = ux
                            elif rd:
                                iux = ux[it.multi_index[0]]
                            else:
                                iux = ux[it.multi_index]
                            i[...] = gummy(0,iux,dof=dof)
<<<<<<< HEAD
                else:
                    x = np.array([i/i.unit for i in self.x])
                            
                #if self.xweights is not None:
                    #x = mmul(self.sqrt_xweights,x)

                if odr:
                    x = cjx*x
                    if self.xweights is not None:
                        x = mmul(self.sqrt_xweights,x)
=======
                            
                if self.xweights is not None:
                    x = mmul(self.sqrt_xweights,x)

                if odr:
                    x = cjx*x
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
                else:
                    x = mmul(self.njacx,x)
    
                px = self.rcov@(x@self.njacp).T
                
                if self.xdim > 1:
                    px = np.sum(px,axis=1)
                
                p = p + px
                    
            for g,gf in zip(p,self._pf):
                g.value._x = gf
            
        else:
            try:
                p = gummy.create(self._pf,covariance_matrix=self.cov,dof=dof)
            except np.linalg.LinAlgError:
                warn('the covariance matrix returned by the solver is not positive semidefinate; uncertainties cannot be calculated',FitWarning)
                p = gummy.create(self._pf)
<<<<<<< HEAD
=======
            # If we have a jacobian, calculate effiective degrees of freedom for each
            # parameter as sum(weights)**2/sum(weights**2) where the weights are the
            # square of the elements of the jacobian.
            #try:
                # sqrtm below if used to transform the jacobian from a basis with 
                # correlated parameters to one where the parameters are uncorrelated
                #cov = self.cov
                #m = np.array([[cov[i][j]/np.sqrt(cov[i][i]*cov[j][j]) 
                               #if cov[i][i] != 0 and cov[j][j] != 0 else 0 
                               #for i in range(self.nparam)] for j in range(self.nparam)])
                #u = np.sqrt(np.diag(cov))
                #val,vec = np.linalg.eig(m)
                #val = np.real(val)
                #val = np.clip(val,0,None)
                #vec = np.real(vec)
                #sqrtm = ((vec*np.sqrt(val)@np.linalg.inv(vec)).T*u).T
                #jact = sqrtm.T@self.njacp.T
                
                # calculate the effective degrees of freedom for the uncorrelated 
                # parameters then transform the resulting gummys back to the
                # correlated basis
                #df = (np.sum(jact**2,axis=1)**2/np.sum(jact**4,axis=1))*(self.count - 1)/self.count
                #df = np.array([i if i >= 1 else 1 for i in df])
                #g = np.array([gummy(0,1,dof=df[i]) for i in range(self.nparam)])
                #p = sqrtm@g + self._pf
            #except np.linalg.LinAlgError:
                #warn('unable to calculate the effective degrees of freedom for the fit parameters')
                #p = gummy.create(self._pf,cov=cov)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            
        self._p = p
        if self.fix is None:
            self.pd = p
            self.pf = self._pf
        else:
            self.pd = np.empty(np.shape(self.p0),dtype=np.dtype('O'))
            self.pd[self.fix] = [gummy(i) for i in self.p0[self.fix]]
            self.pd[~self.fix] = p
            self.pf = repl(self._pf,self.p0,self.fix)
        self.p = np.array([g*u for g,u in zip(self.pd,self.punits)])
    
    
    def _leastsq(self,**kw):
        # non-linear least square solver
        
        from scipy.optimize import least_squares
        from scipy.linalg import pinvh
        
        w = self.sqrt_weights
        if np.ndim(w) == 2:
            if self.y is None:
                def func(params,x,f):
                    return w@f(x,params)
            else:
                def func(params,x,y,f):
                    return w@(f(x,params) - y)
                
            if self._jacp is not None:
                def dfun(*a):
                    return (self._jacp(a[1],a[0])@w.T).T
            else:
                dfun = None
        else:
            if self.y is None:
                def func(params,x,f):
                    return w*f(x,params)
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
            
        if self.y is None:
            args = (self.xf,self._f)
        else:
            args = (self.xf,self.yf,self._f)
        self.fit_output = least_squares(func,self._p0,args=args,**kw)
        self._pf = self.fit_output.x
        
        if self.fit_output.status not in {1,2,3,4}:
            warn('the fit was not sucessful, status ' + str(self.fit_output.status) + ', ' + self.fit_output.message,FitWarning)

        if self.y is None:
            if self.xdim == 1:
                self.res = self.yres = mmul(w,self.ypredf(self.xf))
            else:
                self.res = self.yres = mmul(w,self.ypredf(*self.xf))
        else:
            if self.xdim == 1:
                self.res = self.yres = mmul(w,self.yf - self.ypredf(self.xf))
            else:
                self.res = self.yres = mmul(w,self.yf - self.ypredf(*self.xf))
        self.xres = None
        
        var = np.sum(self.res**2)/(self.count - self.nparam)
                    
        self.var = self.yvar = var
        
        self.rcov = pinvh(self.fit_output.jac.T@self.fit_output.jac)
        
        self.njacp = self.fit_output.jac

        self.xweights = None
        self.sqrt_xweights = 1
        self.known_var = self.known_yvar
        if self.known_var is None:
            if self.variance_is_known:
                warn('variance_is_known is set to True, but no uncertainties have been defined for the y values')
            self.variance_is_known = False
        self.knowm_xvar = None
        self.ux = None
        self.xcov = None
        self.x_is_gummies = False
        
        self._create_p()
        
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
                
        if self.y is None:
            ydim = 1
            yf = np.zeros(self.f0.shape)
            fit_type = 'implicit-ODR'
        else:
            ydim = self.ydim
            yf = self.yf
            fit_type = 'explicit-ODR'
        

        if 'fit_type' in kw:
            fit_type = kw.pop('fit_type')
        if 'task' in kw:
            fit_type = kw.pop('task')
            
        self.fit_output = odr_fit(self._f,self.xf,yf,self._p0,
                                  weight_x=weight_x,weight_y=weight_y,
                                  task=fit_type,jac_beta=self._jacp,
                                  jac_x=self._jacx)
        
        if self.fit_output.info > 5:
            raise RuntimeError('Error encountered during fit:\n' + str(self.fit_output.stopreason) + '\nODR info  = ' + str(self.fit_output.info))
        if self.fit_output.info == 4:
            warn('the iteration limit was reached',FitWarning)

        self._pf = self.fit_output.beta

<<<<<<< HEAD
        self.res = [self.fit_output.delta,self.fit_output.eps]
=======
        self.res = np.array([self.fit_output.delta,self.fit_output.eps])
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
        self.xres = self.fit_output.delta
        self.yres = self.fit_output.eps
        
        self.var = self.fit_output.res_var
        self.xvar = np.sum(self.xres**2)/(self.count - self.nparam)
        self.yvar = np.sum(self.yres**2)/(self.count - self.nparam)
        
        self.known_var = None
        if self.known_yvar is not None:
            self.known_var = self.known_yvar
            if self.known_xvar is not None:
                self.known_var +=self.known_xvar
        
        if self.y is None or weight_y is None or _isscalar(weight_y):
            nwe = 1
        elif np.ndim(weight_y):
            nwe = np.size(weight_y)
        self.njacp,self.njacx = odr_jac(self.fit_output.rwork,self.count,self.xdim,
                                 self.nparam,ydim,nwe,fit_type)

        self.rcov = self.fit_output.cov_beta

        self._create_p()
        
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
            if self.fix is not None or self._jacp is None:
                raise NotImplementedError()
                
            if self.xdim == 1:
                xav = np.mean(self.xf)
                xd = self.xf - xav # Bad data! Bad! (de-mean the data)
            else:
                xav = np.mean(self.xf,axis=1)
                xd = (self.xf.T - xav).T
            trp = self._shiftx(xav)
            
            x = self._jacp(xd,self._p0).T
            self.njacp = self._jacp(self.xf,self._p0).T
        except NotImplementedError:
            trp = None
            if self._jacp is not None:
                x = self._jacp(self.xf,self._p0).T
            else:
                # numerically calculate the jacobian with repect to the fit 
                # parameters
                g0 = [gummy(p,1) for p in self._p0]
                if self.xdim == 1:
                    x = np.array([njacobian(lambda *p:self._f(i,p),*g0) for i in self.xf])
                else:
                    x = np.array([njacobian(lambda *p:self._f(i,p),*g0) for i in self.xf.T])
                if len(self._p0) == 1:
                    x = np.reshape(x,np.shape(x) + (1,))
            self.njacp = np.array(x)
        
        yf = self.yf
        if self.fix is not None:
            yf = yf - self._f(self.xf,np.zeros(self.nparam))
            
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
        self.var = self.yvar = var

        pf = (pf.T/sc).T
        cov /= np.outer(sc,sc)
        
        if trp is not None:
            cov = trp@cov@trp.T
            pf = trp@pf
        
        self.rcov = cov
        self._pf = pf
        
        self.xweights = None
        self.sqrt_xweights = 1
        self.known_var = self.known_yvar
        self.knowm_xvar = None
        if self.known_var is None:
            if self.variance_is_known:
                warn('variance_is_known is set to True, but no uncertainties have been defined for the y values')
            self.variance_is_known = False
        self.ux = None
        self.xcov = None
        self.x_is_gummies = False
        
        self._create_p()
        self.fit_output = None
            
    def ypredf(self,*x):
        """
        returns a float representing the value predicted by the fit at x
        """
        if len(x) != self.xdim:
            raise TypeError('the number of arguments to ypredf must be equal to xdim')
            
        if self.xdim == 1:
            x = x[0]

        return self._f(x,self._pf)
        
    def ypred(self,*x):
        """
        returns a gummy representing the value predicted by the fit at x
        """
        if len(x) != self.xdim:
            raise TypeError('the number of arguments to ypred must be equal to xdim')
            
        ux = False
        x = list(x)
        for i,ix in enumerate(x):
            if self.xdim == 1:
                xun = self._xunit
            else:
                xun = self._xunit[i]
            if not _isscalar(ix):
                ix = np.array(ix)
                if ix.dtype is np.dtype('O'): 
                    with np.nditer(ix,flags=['multi_index','refs_ok'],op_flags=['readwrite']) as it:
                        for i in it:
                            if isinstance(i,AbcQuantity):
                                i[...] = i.convert(xun)/xun
                            if isinstance(i,UncertainValue):
                                ux = True
            else:
                if isinstance(ix,AbcQuantity):
                    x[i] = ix.convert(xun)/xun
                if isinstance(ix,UncertainValue):
                    ux = True
        if ux:
            try:
                r = gummy.apply(self.f,lambda *z:self.jac(*z),*x,*self.pd)
            except NotImplementedError:
                r = gummy.napply(self.f,*x,*self.p)
        else:
            try:
                r = gummy.apply(lambda *z:self.f(*x,*z),lambda *z:self.jacp(*x,*z),*self.pd)
            except NotImplementedError:
<<<<<<< HEAD
                r = gummy.napply(lambda *z:self.f(*x,*z),*self.pd)
=======
                r = gummy.napply(lambda *z:self.f(*x,*z),*self.p)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
            
        if self.ydim <= 1 and self._yunit != 1:
            if _isscalar(r):
                r = r*self._yunit
            else:
                r = np.array([i*self._yunit for i in r])
        elif self.ydim > 1 and any([i != 1 for i in self._yunit]):
            if np.ndim(r) == 1:
                r = np.array([r[i] for i in self._ydim])
            else:
                r = np.array([[j*self._yunit[i] for j in r[i]] for i in self._ydim])
            
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
                nm = 'p('+str(i+1)+')'
                try:
                    if self.p_names is not None:
                        nm = str(self.p_names[i])
                except:
                    pass
                txt += nm + ' = ' + p.tostring(fmt='unicode') + '\n'
            elif fmt == 'ascii' or fmt == 'utf-8':
                nm = 'p('+str(i+1)+')'
                try:
                    if self.p_names_ascii is not None:
                        nm = str(self.p_names_ascii[i])
                except:
                    pass
                txt += nm + ' = ' + p.tostring(fmt='ascii') + '\n'
            elif fmt == 'latex':
                if self.p_names_latex is not None:
                    nm = str(self.p_names_latex[i])
                elif self.p_names is not None:
                    if len(self.p_names[i]) == 1:
                        nm = str(self.p_names[i])
                    else:
                       nm = PrettyPrinter.latex_norm(str(self.p_names[i]))
                txt += 'p_{'+str(i+1)+'} &= ' + p.tostring(fmt='latex').strip('$') + '\\\\'
            elif fmt == 'html':
                nm = '<i>p</i><sub>'+str(i+1)+'</sub>'
                try:
                    if self.p_names_html is not None:
                        nm = str(self.p_names_html[i])
                    elif self.p_names is not None:
                        if len(self.p_names[i]) == 1:
                            nm = '<i>' + str(self.p_names[i]) + '</i>'
                        else:
                            nm = str(self.p_names[i])
                except:
                    pass
                txt += nm + ' = ' + p.tostring(fmt='html') + '<br>'
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
