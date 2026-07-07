# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:54:45 2026

@author: Parksh
"""

import numpy as np
from ._ummy import _isscalar

from .fit import Fit

class PolyFit(Fit):
    def __init__(self,x,y,deg=1,**kw):
        """
        Fits the x,y data to a polynomial

        Parameters
        ----------
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
       
        deg:  `int` or list of `int`
            The degree of the polynomial.  For one-dimensional x-values
            this is scalar value, or a list of length two or three for a two
            or three dimensional fit.
            
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
            
        variance_is_known: `bool`, optional
            If this is `True` then any uncertainties in the data  (either as
            gummys in the `x` or `y` values or in the ux or uy parameters)
            are used to calculate the uncertainties in the fit.  Otherwise,
            the uncertainties are based on the standard deviation of the
            residuals and the uncertainties in the data are used only for
            weighting the data points.  The default value is `True`.

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

        Methods
        -------
        ypred(x1,x2,...):
            Takes `xdim` floats (or arrays of float) and returns a gummy 
            representing the predicted value(s) at that x-coordinate.
            
        ypredf(x1,x2,...):
            Takes `xdim` floats (or arrays of float) and returns a float giving 
            the  predicted value(s) at that x-coordinate.
            
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

        super().__init__(x,y=y,**kw)

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