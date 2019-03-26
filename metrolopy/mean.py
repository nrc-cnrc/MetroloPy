# -*- coding: utf-8 -*-

# module mean

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
Some miscellaneous functions that are useful for uncertainty analysis.
"""

import numpy as np
from .gummy import gummy
import datetime

def autocorrelation(x):
    """
    Find the de-meaned and normalized autocorrelation of the array x
    """
    z = x - np.mean(x)
    rho = np.correlate(z, z, mode='full')
    rho = rho[int(rho.size/2):]
    return rho/np.sum(z**2)

def n_eff(x, y = None):
    """Takes a series of data points which may be 
    correlated and uses the auto-correlation of the set to find an effective 
    number of data points such that an estimate of the standard uncertainty is
    (standard deviation)/sqrt(n_eff).
    
    See [N.F. Zhang, Metrologia, 43, S276 (2006)] for details.

    If only one argument is given it is assumed that the argument represents
    an array of evenly spaced values.  This list may be a masked array, in
    which case the masked values are filled in with a linear interpolation.

    If two arguments are passed to n_eff, then it is assumed that the two
    arrays represent an array of x values followed by an array representing
    the corresponding y values of a set (x,y) points.  A cubic spline
    interpolation is used to get a set of evently spaced points.

    Returns
    -------
    `float`:
        A float that represents the effective number of data points.
    """
    n = len(x)
    if n == 0:
        return 0
    if n == 1:
        return 1
        
    if y is None:
        # If y is None assume the data points x are evenly spaced.
        
        # If x is a masked array then fill in the missing values with a
        # linear interpolation.
        if isinstance(x,np.ma.masked_array):
            if np.ma.getmask(x) is np.ma.nomask:
                # Nothing is masked so we don't need to do anything.
                xx = x.data
                ct = n
            else:
                ct = x.count()
                if ct == 0: 
                    return 0
                
                # Trim off any masked values on the ends of the array.
                ep = np.ma.flatnotmasked_edges(x)
                a = x[ep[0]:(ep[1] + 1)]
                a.unshare_mask()
                n = len(a)
                if n == 1:
                    return 1
                for i in range(n):
                    if a[i] is np.ma.masked:
                        j = i - 1
                        while a[j] is np.ma.masked:
                            j -= 1
                        k = i + 1
                        while a[k] is np.ma.masked:
                            k += 1
                        a[i] = a[j] + (i-j)*(a[k]-a[j])/(k-j)
                xx = a.data
        else:
            xx = x
            ct = n
            
        # Now follow the recipe given in Zhang.
        nq = int(n/4)
        rho = autocorrelation(xx)
        if np.isnan(rho[0]):
            return 1.0
        
        a = 0
        nc = 0
        for i in range(1, nq):
            if np.abs(rho[i]) > 1.96*np.sqrt((1 + 2*a)/n):
                nc = i
            a += rho[i]**2
            
        r = 0
        for i in range(1, nc):
            r += (n - i)*rho[i]
        ne = n/(1 + 2*r/n)
        if ne > ct:
            ne = ct
        if ne < 1:
            return 1.0
        return ne
    
    #If y is not none then use a cubic spline interpolation to get evenly
    #spaced points and call this function again.
    from scipy import interpolate
    tck = interpolate.splrep(x, y, s=0)
    xx = np.arange(x[0], x[-1], (x[-1] - x[0])/n)
    yy = interpolate.splev(xx, tck, der=0)
    return n_eff(yy)

def wmean(x,chi_correct=False):
    """
    Takes an array of gummys and returns the weighted mean with weights that 
    minimize the uncertainty of the returned value (taking into account correlations)
    
    Parameters
    ----------
    x:  array_like of `gummy`
        the values to be averaged

    chi_correct:  `bool`, optional
        If this is `True`, the uncertainty of the returned value is multiplied by the
        square root of the reduced chi-squared of the residuals to take into account
        any under estimation of the uncertainties in `x`.  The default is `False`
            
    Returns
    -------
    `gummy`:
        a gummy representing the weighted mean
    """
    from scipy.linalg import inv
    
    unit = x[0].unit
    if not unit.linear:
        iunit = unit.conversion.unit
    else:
        iunit = unit
        
    x = [g.convert(iunit) for g in x]
    
    d = np.ones(len(x))
    v = gummy.covariance_matrix(x)
    vsc = np.diagonal(v).sum()/len(x)
    v /= vsc
    iv = inv(v)
    ret = np.sum((d.T@iv)*x)/(d.T@iv@d)
    
    if chi_correct:
        res = np.array([g.x - ret.x for g in x])
        chisq = res.T@iv@res/(vsc*(len(x)-1))
        ret._u *= np.sqrt(chisq)
        ret._set_U()
        ret._dof = len(x) - 1
    
    if iunit is not unit:
        ret.unit = unit
        
    return ret
    
def mean(x,n_sigma_trim=3,unit=1,ignore_nan=True,use_n_eff=None,bayesian=None):
    """
    Returns a gummy representing the mean of a float array.
    
    Parameters
    ----------
    x:  array_like of `float` or `int`
        the value to be averaged

    n_sigma_trim:  `int`, optional
        If this is not `None`, then ``sigma_trim(x, n_sigma_trim)`` is applied
        to the data before taking the mean.  Set this argument to `None` if you
        don't want `sigma_trim` to be applied.  The default value is 3.

    unit: `str`, `Unit` or 1
        The unit of returned gummy.  The default is 1.

    ignore_nan: `bool`, optional
        If this is `True`, elements with a ``float('nan')`` or `None` value will be ignored.
        The default value is `True`

    use_n_eff: `bool`, optional
        Whether to use the `n_eff` function to calculate an effective number of degrees of freedom.
        If use_n_eff is `None`, then `n_eff` will be used if the length of x is greater or equal
        to 20.  The default value is `None`.

    bayesian:  `bool`, optional
        If bayesian is `False` the standard uncertainty of the returned
        gummy is s/sqrt(n) where s is the standard deviation of x and n is
        the the number of samples (or n_eff).  If bayesian is `True` then
        the standard uncertainty is ((n-1)/(n-3))*s/sqrt(n).  If bayesian
        is `None` then the value of `gummy.bayesian` will be used.  The
        default value is `None`.
    """
    x = np.asanyarray(x)
    if ignore_nan:
        x = np.ma.masked_where([y is None or np.isnan(y) for y in x],x)
    if sigma_trim is not None:
        x = sigma_trim(x, n_sigma_trim)
    
    n = x.count()
    if use_n_eff or (use_n_eff is None and n >= 20):
        n = n_eff(x)

    m = float(np.mean(x))
    
    dof = n - 1
    if dof < 1:
        dof = 1
        
    if bayesian is None:
        bayesian = gummy.bayesian
    if bayesian:
        if dof <= 2:
            raise ValueError('dof is ' + str(dof) + '; it must be > 2')
        u = (dof/(dof - 2))*float(np.std(x,ddof=1))/np.sqrt(n)
        return gummy(m,u,unit=unit)
    else:
        u = float(np.std(x,ddof=1))/np.sqrt(n)
        return gummy(m,u,dof=dof,unit=unit)
    
def sigma_trim(x, n_sigma = 3):
    """
    Returns a masked array with data attribute equal to `x` and any elements
    more the `n_sigma` standard deviations from the mean masked.  (The
    standard deviation is calculated excluding the masked outliers.)
    """
    ma = np.ma.masked_where(np.abs(x - np.mean(x)) > n_sigma*np.std(x,ddof=1),x)
    n = 0
    
    # Iteratievly mask elements until no new elements are masked.
    while ma.count() != n:
        n = ma.count()
        ma = np.ma.masked_where(np.abs(ma - np.mean(ma)) > n_sigma*np.std(ma,ddof=1),ma)
    return ma
    
def delta_diff(x):
    """
    Returns a list containing differences of the type:
    
    ``x[1] - (x[0]+x[2])/2 , (x[1]+x[3])/2 - x[2], x[3] - (x[2]+x[4])/2, ...``
    
    Differences of this type are useful when we want the difference between 
    alternate points in a data set (e.g. we are alterately switching between
    a signal and a background).  This delta type difference removes the effect
    of a slow drift from the data.
    
    Parameters
    ----------
    x:  array_like of `float`
        the data
            
    Returns
    -------
    `numpy.ndarray`
    
    Notes
    -----
    We assume the signal is in the second element (array index 1 and
    all elements with odd indices) and the background is in the first
    element (array index 0 and all elements with even indices).
    """
    z = np.asanyarray(x)
    z = np.diff(z,2)/2
    z[0::2] = -z[0::2]
    return z

def delta_diff_mean(x,n_sigma=None,unit=1,bayesian=None):
    """
    Returns a gummy representing the mean value and uncertainty of a delta 
    type difference taken on the data.  A delta type difference removes a 
    linear drift from the data by taking differences:
    
    ``x[1] - (x[0]+x[2])/2 , (x[1]+x[3])/2 - x[2], x[3] - (x[2]+x[4])/2, ...``
    
    Parameters
    ----------
    x:  array_like of `float`
        the data

    n_sigma:  `int`, optional
        If this is not `None`, then ``sigma_trim(x, n_sigma_trim)`` is applied
        to the differences before taking the mean.  Set this argument to `None`
        if you don't want `sigma_trim` to be applied.  The default value is `None`.

    unit: `str`, `Unit` or 1
        The unit of returned gummy.  The default is 1.

    bayesian:  `bool`, optional
        If bayesian is `False` the standard uncertainty of the returned
        gummy is s/sqrt(n) where s is the standard deviation of x and n is
        the the number of samples (or n_eff).  If bayesian is `True` then
        the standard uncertainty is ((n-1)/(n-3))*s/sqrt(n).  If bayesian
        is `None` then the value of `gummy.bayesian` will be used.  The
        default value is `None`.
            
    Returns
    -------
    `gummy`:
        A gummy representing the mean of the differences along with the
        uncertainty and effective degrees of freedom.
    """
    diff = delta_diff(x)
    if n_sigma is not None:
        diff = sigma_trim(diff, n_sigma)
        n = diff.count() + 2
    else:
         n = len(diff) + 2
    mn = np.mean(diff)
    if n <= 3:
        return gummy(mn)
    sd = np.std(diff,ddof=1)
    u = (np.sqrt(4*n-11)/(n-2))*sd
    dof = (4*n-11)**2/(41/2+16*(n-4))
    
    if bayesian is None:
        bayesian = gummy.bayesian
    if bayesian:
        u = u*dof/(dof-2)
        return gummy(mn,u,unit=unit)
    
    return gummy(mn,u,dof=dof,unit=unit)
    
def delta_sum(x):
    """Returns a list containing differences of the type:
    
    ``x[1] + (x[0]+x[2])/2 , (x[1]+x[3])/2 + x[2], x[3] + (x[2]+x[4])/2, ...``
    
    Parameters
    ----------
    x:  array_like of `float`
        the data
            
    Returns
    -------
    `numpy.ndarray`
    """
    z = np.asanyarray(x)
    z[0::2] = -z[0::2]
    return delta_diff(z)/2.0
    
def delta_sum_mean(x,n_sigma=None,unit=1):
    """
    Returns a gummy representing the mean value
    and uncertainty of a delta type difference taken on the data.  A delta type 
    difference removes a linear drift from the data by taking differences:
    
    ``x[1] + (x[0]+x[2])/2 , (x[1]+x[3])/2 + x[2], x[3] + (x[2]+x[4])/2, ...``

    Parameters
    ----------
    x:  array_like of `float`
        the data

    n_sigma:  `int`, optional
        If this is not `None`, then ``sigma_trim(x, n_sigma_trim)`` is applied
        to the differences before taking the mean.  Set this argument to `None`
        if you don't want `sigma_trim` to be applied.  The default value is `None`.

    unit: `str`, `Unit` or 1, optional
        The unit of returned gummy.  The default is 1.

    Returns
    -------
    `gummy`:
        A gummy representing the mean of the differences along with the uncertainty and
        effective degrees of freedom.
    """
    z = np.asanyarray(x)
    z[0::2] = -z[0::2]
    return delta_diff_mean(z,n_sigma,unit)/2.0
    
def mean_datetime(*params):
    if len(params) == 1:
        p = params[0]
    else:
        p = params
    z = datetime.timedelta(0)
    for d in p[1:]:
        z += d - p[0]
    z = z/len(p)
    return p[0] + z