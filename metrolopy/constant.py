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
    
    def __new__(cls,x,u=0,unit=1,dof=float('inf'),k=1,p=None,uunit=None,
                 utype=None,name=None,symbol=None,short_name=None,
                 add_symbol=False, html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None):
        if symbol is None:
            ret = gummy.__new__(gummy,x,u=u,unit=unit,dof=dof,k=k,p=p,
                                uunit=uunit,utype=utype,name=name)
            ret.__init__(x,u=u,unit=unit,dof=dof,k=k,p=p,uunit=uunit,
                         utype=utype,name=name)
            return ret
        
        return super().__new__(cls)
    
    def __init__(self,x,u=0,unit=1,dof=float('inf'),k=1,p=None,uunit=None,
                 utype=None,name=None,symbol=None,short_name=None,
                 add_symbol=False, html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None):
        if name is None:
            name = symbol
        gummy.__init__(self,x,u=u,unit=unit,dof=dof,k=k,p=p,uunit=uunit,
                      utype=utype,name=name)
        Indexed.__init__(self,name,symbol=symbol,short_name=short_name,
                         add_symbol=add_symbol,html_symbol=html_symbol,
                         latex_symbol=latex_symbol,ascii_symbol=ascii_symbol,
                         description=description)
    
    def tostring(self,fmt=None,style=None,k=None,p=None,show_k=None,
                 show_p=None,show_dof=None,show_name=None,name=None,
                 norm=None,raw=False,nsig=None,solidus=None,
                 mulsep=None,**kwds):
                     
        if name is None:
            name = Indexed.tostring(self,fmt)
        
        return super().tostring(fmt=fmt,style=style,k=k,p=p,show_k=show_k,
                                show_p=show_p,show_dof=show_dof,
                                show_name=show_name,name=name,norm=norm,
                                raw=raw,nsig=nsig,solidus=solidus,
                                mulsep=mulsep)
    
class JummyConstant(jummy,Indexed):
    pass
