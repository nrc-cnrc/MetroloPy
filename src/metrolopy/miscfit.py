# -*- coding: utf-8 -*-
<<<<<<< HEAD

# module miscfit

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
=======
"""
Created on Fri Jun 12 09:58:43 2026

@author: Parksh
"""
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5

import numpy as np
from .fit import Fit
from .polyfit import PolyFit

class DoubleExpFit(Fit):
    """
    DoubleExpFit(x,y,p0=None,ux=None,uy=None,variance_is_known=True,xunit=None, yunit=None,
             solver=None,**keywords)
   
    Fits the x,y data to a function of the form:
       
    p[0]*np.exp(x/p[1])+p[2]*np.exp(x/p[3])+p[4]
   
    Parameters
    ----------
<<<<<<< HEAD
    x: array
        The indepenant variables x.  For n data points this should be an
        array of shape (n,) for a one-dimensional fit or shape (m,n) for an
        m-dimensional fit.  The maximum dimension m is 3.
        
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
        
    y:  array
        Response variables y. 
        
        The y-values may be all `float` or all `gummy`,  If uncertainties
        u are defined for the y-values, the weights will be set to 1/u**2
        for each point.  To override this behavior set the keyword parameter
        "weights = 1".
        
        If the nls or odr solver is used, the weighting will take into
        account correlations between the different points.  However, the odr
        solver ignores correlations between different data points.  If the
        y-values are multi-dimensional, the odr solver will take into
        account correlations between the different elements at each point.
        
    fix: array of `bool`
        A mask for the fit parameters.  For any element in `fix` that is
        `True`, the corresponding fit parameter will be held constant at
        its initial value.
        
    solver:  {'ols','nls','odr'}, optional
        A non-linear least squares (nls), ordinary least squares (ols) or
        othogonal distance regression solver (odr) can be used to perform
        the fit.
        
        ols is the default solver unless uncertainties are specified for
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
        behavior set the keyword parameter "weights = 1"
        
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
<<<<<<< HEAD
        
    p0: array of `float`, optional
        initial values for the fit parameters.

    Attributes
    ----------
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
        
=======
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    plot(...):
        plots the data (only available if x and y are one-dimensional)
    """
    
    def get_p0(self):
        ft = ExpFit(self.xf,self.yf)
        r = np.abs(ft.pf[1])/(np.max(self.xf) - np.min(self.xf))
<<<<<<< HEAD
        count = len(self.yf)
        if r < 0.33:
            a = int(count*2*r)
            b = count
        elif r < 1:
            a = int(2*count/3)
            b = count
        else:
            a = 0
            b = int(count/3)
=======
        if r < 0.33:
            a = int(self.count*2*r)
            b = self.count
        elif r < 1:
            a = int(2*self.count/3)
            b = self.count
        else:
            a = 0
            b = int(self.count/3)
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
<<<<<<< HEAD
    x: array
        The indepenant variables x.  For n data points this should be an
        array of shape (n,) for a one-dimensional fit or shape (m,n) for an
        m-dimensional fit.  The maximum dimension m is 3.
        
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
        
    y:  array
        Response variables y. 
        
        The y-values may be all `float` or all `gummy`,  If uncertainties
        u are defined for the y-values, the weights will be set to 1/u**2
        for each point.  To override this behavior set the keyword parameter
        "weights = 1".
        
        If the nls or odr solver is used, the weighting will take into
        account correlations between the different points.  However, the odr
        solver ignores correlations between different data points.  If the
        y-values are multi-dimensional, the odr solver will take into
        account correlations between the different elements at each point.
        
    fix: array of `bool`
        A mask for the fit parameters.  For any element in `fix` that is
        `True`, the corresponding fit parameter will be held constant at
        its initial value.
        
    solver:  {'ols','nls','odr'}, optional
        A non-linear least squares (nls), ordinary least squares (ols) or
        othogonal distance regression solver (odr) can be used to perform
        the fit.
        
        ols is the default solver unless uncertainties are specified for
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
        behavior set the keyword parameter "weights = 1"
        
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
<<<<<<< HEAD
        
    p0: array of `float`, optional
        initial values for the fit parameters.

    Attributes
    ----------
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
        
=======
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
<<<<<<< HEAD
    x: array
        The indepenant variables x.  For n data points this should be an
        array of shape (n,) for a one-dimensional fit or shape (m,n) for an
        m-dimensional fit.  The maximum dimension m is 3.
        
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
        
    y:  array
        Response variables y. 
        
        The y-values may be all `float` or all `gummy`,  If uncertainties
        u are defined for the y-values, the weights will be set to 1/u**2
        for each point.  To override this behavior set the keyword parameter
        "weights = 1".
        
        If the nls or odr solver is used, the weighting will take into
        account correlations between the different points.  However, the odr
        solver ignores correlations between different data points.  If the
        y-values are multi-dimensional, the odr solver will take into
        account correlations between the different elements at each point.
        
    fix: array of `bool`
        A mask for the fit parameters.  For any element in `fix` that is
        `True`, the corresponding fit parameter will be held constant at
        its initial value.
        
    solver:  {'ols','nls','odr'}, optional
        A non-linear least squares (nls), ordinary least squares (ols) or
        othogonal distance regression solver (odr) can be used to perform
        the fit.
        
        ols is the default solver unless uncertainties are specified for
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
        behavior set the keyword parameter "weights = 1"
        
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
<<<<<<< HEAD
        
    p0: array of `float`, optional
        initial values for the fit parameters.

    Attributes
    ----------
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
    """
    
    def get_p0(self):
        count = len(self.yf)
        b = np.mean(self.yf[(count-2-int(count/10)):(count-1)])
=======
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
    plot(...):
        plots the data (only available if x and y are one-dimensional)    """
    
    def get_p0(self):
        b = np.mean(self.yf[(self.count-2-int(self.count/10)):(self.count-1)])
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
<<<<<<< HEAD
    x: array
        The indepenant variables x.  For n data points this should be an
        array of shape (n,) for a one-dimensional fit or shape (m,n) for an
        m-dimensional fit.  The maximum dimension m is 3.
        
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
        
    y:  array
        Response variables y. 
        
        The y-values may be all `float` or all `gummy`,  If uncertainties
        u are defined for the y-values, the weights will be set to 1/u**2
        for each point.  To override this behavior set the keyword parameter
        "weights = 1".
        
        If the nls or odr solver is used, the weighting will take into
        account correlations between the different points.  However, the odr
        solver ignores correlations between different data points.  If the
        y-values are multi-dimensional, the odr solver will take into
        account correlations between the different elements at each point.
        
    fix: array of `bool`
        A mask for the fit parameters.  For any element in `fix` that is
        `True`, the corresponding fit parameter will be held constant at
        its initial value.
        
    solver:  {'ols','nls','odr'}, optional
        A non-linear least squares (nls), ordinary least squares (ols) or
        othogonal distance regression solver (odr) can be used to perform
        the fit.
        
        ols is the default solver unless uncertainties are specified for
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
        behavior set the keyword parameter "weights = 1"
        
=======
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
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
    variance_is_known: `bool`, optional
        If this is `True` then any uncertainties in the data  (either as
        gummys in the `x` or `y` values or in the ux or uy parameters)
        are used to calculate the uncertainties in the fit.  Otherwise,
        the uncertainties are based on the standard deviation of the
        residuals and the uncertainties in the data are used only for
        weighting the data points.  The default value is `True`.
<<<<<<< HEAD
        
    p0: array of `float`, optional
        initial values for the fit parameters.

    Attributes
    ----------
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
        
=======
        Takes `xdim` floats and returns a gummy representing the predicted
        value at that x-coordinate.
    ypredf(x1,x2,...):
        Takes `xdim` floats and returns a float giving the  predicted value
        at that x-coordinate.
>>>>>>> 521c361ba2fc57e9677804d95b4bb16b2095dfa5
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
        return [self._yunit, 1/self._xunit, 1, self._yunit]
        
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

        