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
from .unit import unit
from .siunits import _const_c,_const_G,_const_h,_const_e,_const_k
import numpy as np

with GummyConstant._builtin():
    
    GummyConstant(9192631770,unit='Hz', 
                  name='hyperfine transition frequency of Cs',
                  symbol='\u0394\u03bd(Cs)',
                  html_symbol='&Delta;<i>&nu;</i><sub>Cs</sub>',
                  latex_symbol='\\Delta\\nu_{Cs}',
                  ascii_symbol='Delta-nu(Cs)',
                  add_symbol='True',
                  description='SI defining constant, SI Brochure 9th ed.,Caesium, Cesium')
    
    _h = GummyConstant(_const_h,unit='J s',symbol='h',name='Plank constant',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.')
    
    GummyConstant(_const_c,unit='m/s',name='speed of light in vacuum',
                  symbol='c',add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed.')
    
    GummyConstant(_const_e,unit='C',symbol='e',name='elementary charge',
                  add_symbol='True',
                  description='SI defining constant, SI Brochure 9th ed., electron charge')
    
    GummyConstant(6.02214076e23,unit='mol**-1',name='Avogadro constant',
                  symbol='N(A)',
                  html_symbol='<i>N</i><sub>A</sub>',
                  latex_symbol='N_{A}',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., mole')
    
    GummyConstant(_const_k,unit='J/K',symbol='k',name='Boltzmann constant',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., Kelvin')
    
    GummyConstant(683,unit='lm/W',
                  name='luminous efficacy of monochromatic radiation of frequency 540e12 Hz',
                  symbol='K(cd)',html_symbol='<i>K</i><sub>cd</sub>',
                  latex_symbol='K_{cd}',add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., candela')


    GummyConstant(_h/(2*np.pi),unit='J s',name='reduced Plank constant',
                  symbol='h-bar',html_symbol='&hbar;',latex_symbol='\\bar',
                  description='Plank constant over 2 pi')
    
    
    _G =  GummyConstant(_const_G*unit('m**3 kg**-1 s**-2'),
                        name='Newtonian constant of gravitation',symbol='G',
                        add_symbol=True,short_name='constant of gravitation',
                        description='CODATA 2018, gravity')
    