# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:08:11 2026

@author: Parksh
"""

import numpy as np
from .Fit import Fit

class DistFit(Fit):
    def __init__(self,x,cdf,p0,**kwds):
        self.cdf = cdf
        x = np.sort(x)
        super().__init__(x,p0=p0,**kwds)
        
    def f(self,x,*p,**kwds):
        z = len(x)*np.diff(self.cdf(x,*p,**kwds),prepend=0)
        return np.log(z*np.exp(z)) - 1
        