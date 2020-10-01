# -*- coding: utf-8 -*-

# module constant

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
constant
"""

from .indexed import Indexed
from .gummy import gummy,jummy

class GummyConstant(gummy,Indexed):
    
    def tostring(self,fmt=None,style=None,k=None,p=None,show_k=None,
                 show_p=None,show_dof=None,show_name=None,name=None,
                 norm=None,raw=False,nsig=None,solidus=None,
                 mulsep=None,**kwds):
                     
        if name is None:
            name = super(Indexed,self).tostring(fmt)
        
        return super(gummy,self).tostring(fmt,style,k,p,show_k,show_p,show_dof,
                                          show_name,name,norm,raw,nsig,solidus,
                                          mulsep,**kwds)
