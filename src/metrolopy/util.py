# -*- coding: utf-8 -*-

# module util

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

from .abc import AbcQuantity,UncertainValue,UncertainComplexValue

def _isscalar(x):
    #try:
        #return np.ndim(x) == 0
    #except:
        #return True
    if isinstance(x,str):
        return True
    try:
        len(x)
        return False
    except:
        return True
    
def _replq(x):
    if _isscalar(x):
        if isinstance(x,np.ndarray):
            return _replq(x.item())
        if isinstance(x,AbcQuantity):
            return _replq(x.convert(1).value)
        return x
    
    x = np.array(x)
    with np.nditer(x,flags=['refs_ok'],op_flags=['readwrite']) as it:
        for r in it:
            r[...] = _replq(r)
    return x

def _replu(x):
    if _isscalar(x):
        if isinstance(x,np.ndarray):
            return _replu(x.item())
        if isinstance(x,(UncertainValue,UncertainComplexValue)):
            return x.x
        return x
    
    x = np.array(x)
    with np.nditer(x,flags=['refs_ok'],op_flags=['readwrite']) as it:
        for r in it:
            r[...] = _replu(r)
