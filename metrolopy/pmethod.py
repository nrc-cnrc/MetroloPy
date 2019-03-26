# -*- coding: utf-8 -*-

# module pmethod

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
Functions and classes for converting between coverage factor k and level of 
confidence p.
"""

import numpy as np

def coverage_factor(p,dof=float('inf'),bayesian=False):
    """Returns the k factor for a given level of confidence and degrees of freedom
    calculated from a normal or Student's t distribution
    
    Parameters:
        
    p:  `float`
        The level of confidence (as a fraction of 1).
    
    dof:  `float` or `int`,
        optional (default is float('inf')), the degrees of freedom.
        
    Returns
    -------
    `float`
    """
    from scipy.stats import t as students_t
    from scipy.stats import norm
    
    if np.isinf(dof):
        ret = norm.ppf(0.5 + p/2)
    else:
        ret = students_t.ppf(0.5 + p/2, dof)
        if bayesian:
            ret *= np.sqrt((dof-2)/dof)
               
    return float(ret)
    
def loc_from_k(k,dof=float('inf'),bayesian=False):
    """Returns the level of confidence given a coverage factor k and degrees of
    freedom for a Student's t distribution.
    
    Parameters:
        k:  `float`
            coverage factor.
        dof:  `float` or `int`, optional
            (default is float('inf')), the degrees of freedom.
        
    Returns:  float
    """
    from scipy.stats import t as students_t
    from scipy.stats import norm
    
    if np.isinf(dof):
        return float(1-2*(1-norm.cdf(k)))
    
    if bayesian:
        if dof <= 2:
            return 1
        k *= np.sqrt(dof/(dof-2))
    return float(1-2*(1-students_t.cdf(k,dof)))

def coverage_probability(p,dof=float('inf'),bayesian=False):
    k = 2/(3*(np.sqrt(1 - p)))
    if not bayesian and np.isfinite(dof):
        k *= np.sqrt(dof/(dof-2))
    return float(k)

def cp_from_k(k,dof=float('inf'),bayesian=False):
    if not bayesian and np.isfinite(dof):
        k *= np.sqrt((dof-2)/dof)
    cp = 1 - 4/(9*k**2)
    return float(cp)

def conservative_coverage_probability(p,dof=float('inf'),bayesian=False):
    k = 1/np.sqrt(1 - p)
    if not bayesian and np.isfinite(dof):
        k *= np.sqrt(dof/(dof-2))
    return float(k)

def ccp_from_k(k,dof=float('inf'),bayesian=False):
    if not bayesian and np.isfinite(dof):
        k *= np.sqrt((dof-2)/dof)
    ccp = 1 - 1/k**2
    return float(ccp)

class _Pmthd:
    def __init__(self,text):
        text = text.lower().strip()
        
        if text in {'loc','level of confidence'}:
            self.method = 'loc'
            self.fptok = coverage_factor
            self.fktop = loc_from_k
            self.text = 'level of confidence'
            return
        
        if text in {'cp','coverage probability','gauss'}:
            self.method = 'cp'
            self.fptok = coverage_probability
            self.fktop = cp_from_k
            self.text = 'coverage probability'
            return
        
        if text in {'ccp','conservative coverage probability','chebyshev'}:
            self.method = 'ccp'
            self.fptok = conservative_coverage_probability
            self.fktop = ccp_from_k
            self.text = 'conservative coverage probability'
            return