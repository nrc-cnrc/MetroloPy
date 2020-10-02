# -*- coding: utf-8 -*-

# module codata2018

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
constants
"""

from .constant import GummyConstant
from .siunits import _const_c,_const_G

with GummyConstant._builtin():
    GummyConstant(_const_c,unit='m s**-1',name='speed of light in vacuum',
                  symbol='c',add_symbol=True,
                  description='SI defining constant, exact, CODATA 2018')
    
    _const_G = _const_G.copy().graft('m**3 kg**-1 s**-2')
    GummyConstant(_const_G,name='Newtonian constant of gravitation',symbol='G',
                  add_symbol=True,description='CODATA 2018')