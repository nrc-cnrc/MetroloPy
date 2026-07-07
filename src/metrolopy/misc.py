# -*- coding: utf-8 -*-
# module misc

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
Misc functions.
"""
from ._gummy import gummy

def correlation_matrix(gummys):
    """
    Returns the correlation matrix of a list or array of gummys.
    
    This is an alias for the `gummy.correlation_matrix` static method.
    """
    return gummy.correlation_matrix(gummys)

def covariance_matrix(gummys):
    """
    Returns the variance-covariance matrix of a list or array of gummys.
    
    This is an alias for the `gummy.covariance_matrix` static method.
    """
    return gummy.covariance_matrix(gummys)

def correlation_matrix_sim(gummys):
    """
    The staticmethod takes a list of gummys an returns the correlation
    matrix calculated from Monte-Carlo data.  The return value is numpy 
    ndarray.
    
    See the method `gummy.correlation_matrix(gummys)` for the corresponding
    result based on first order error propagation.
    
    This is an alias for the `gummy.correlation_matrix_sim` static method.
    """
    return gummy.correlation_matrix_sim(gummys)

@staticmethod
def covariance_matrix_sim(gummys):
    """
    The staticmethod takes a list of gummys an returns the variance-covariance
    matrix calculated from Monte-Carlo data.  The return value is numpy
    ndarray.
    
    See the method gummy.covariance_matrix(gummys) for the corresponding
    result based on first order error propagation.
    
    This is an alias for the `gummy.covariance_matrix_sim` static method.
    """
    return gummy.covariance_matrix_sim(gummys)

def simulate(gummys,n=100000,ufrom=None):
    """
    Generates Monte-Carlo data for one or more gummys.  Calling this method
    erases previously generated Monte-Carlo data for all gummys.  See also
    the `gummy.sim()` method to generate data for one gummy only.
    
    Parameters
    ----------
    n:  `int` > 0, optional
        The number of samples to generate.  The default value is 100000.
        
    gummys: A list or array of `gummy` for which to generate the Monte-Carlo
        data.

    ufrom: `None`, `gummy`, `str` or array_like
        If this is not `None`, then only the gummys referenced here will be
        allowed to vary, and all other gummys will be held fixed at their
        mean values.  This can be a gummy, a string referencing a utype or
        a list containing  gummys and strings.  The default value is `None`.
        
    This is an alias for the `gummy.simulate` static method.
    """
    return gummy.simulate(gummys,n=n,ufrom=ufrom)

def clear_all_sim():
    """
    Clears Monte-Carlo data from all existing gummys.
    
    This is an alias for the `gummy.clear_all` static method.
    """
    gummy.clear_all()
    
def covplot(x,y,title=None,xlabel=None,ylabel=None,mean_marker=False,
            mean_marker_options={},hold=False,math=None,**plot_options):
    """
    Creates scatter plot showing the covariance between two gummys.
    
    Parameters
    ----------
    x:  `gummy`
        The gummy to plot on the horizontal axis.
    
    y:  `gummy`
        The gummy to plot on the vertical axis.
    
    title:  `str` or `None`, optional
        A title for the plot.  If this is omitted or set
        to None then the correlation will be displayed as the title.
       
    xlabel:  `str` or `None`, optional
        A label for the horizontal axis.  If this os
        omitted or None then that axis will be labeled either "x" or with
        the `x` gummy's unit.
       
    ylabel:  `str` or `None`, optional
        A label for the vertical axis.  If this os
        omitted or None then that axis will be labeled either "y" or with
        the `y` gummy's unit.
       
    mean_marker:  `bool`, optional
        Whether or not to display line markers at the mean
        values of `x` and `y`.  The default is `False`.
       
    mean_marker_options:  `dict`, optional
        A dictionary of options to be passed to the
        `pyplot.axvline` and `pyplot.axhline` methods that draw the `mean_marker`.
       
    hold:  `bool`, optional
        If this is `False` then ``pyplot.show()`` is called before this method
        exits.  If it is `True` ``pyplot.show()`` is not called.  The default is
        `False`.
       
    plot_options:  These are optional keyword arguments that are passed to
         the `pyplot.plot` method.  For example ``ms=0.1`` decreases the size of the
         dots in the plot.
         
    This is an alias for the `gummy.covplot` static method.
    """
    gummy.covplot(x,y,title=title,xlabel=xlabel,ylabel=ylabel,
                  mean_marker=mean_marker,
                  mean_marker_options=mean_marker_options,hold=hold,math=math,
                  **plot_options)
