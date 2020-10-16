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
constant definitions
"""

from .constant import GummyConstant
from .unit import unit
from .ummy import MFraction
from .gummy import gummy
from .siunits import (pi,_const_c,_const_G,_const_h,_const_e,_const_k,
                      _const_dalton,_const_me,_const_alpha,
                      _const_rydberg,_const_a0,_const_proton_mass)

norm = GummyConstant.latex_norm

with GummyConstant._builtin():
    
    GummyConstant(9192631770,unit='Hz', 
                  name='hyperfine transition frequency of Cs',
                  symbol='\u0394\u03bd(Cs)',
                  html_symbol='&Delta;<i>&nu;</i><sub>Cs</sub>',
                  latex_symbol='\\Delta\\nu_{' + norm('Cs') + '}',
                  ascii_symbol='hf(Cs)',
                  add_symbol='True',
                  description='SI defining constant, SI Brochure 9th ed.,Caesium, Cesium')
    
    _h = GummyConstant(_const_h,unit='J s',symbol='h',name='Plank constant',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.')
    
    _c = GummyConstant(_const_c,unit='m/s',name='speed of light in vacuum',
                  symbol='c',add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed.')
    
    _e = GummyConstant(_const_e,unit='C',symbol='e',name='elementary charge',
                  add_symbol='True',
                  description='SI defining constant, SI Brochure 9th ed., electron charge')
    
    GummyConstant(MFraction('6.02214076e23'),unit='mol**-1',
                  name='Avogadro constant',
                  symbol='N(A)',
                  html_symbol='<i>N</i><sub>A</sub>',
                  latex_symbol='N_{' + norm('A') + '}',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., mole')
    
    GummyConstant(_const_k,unit='J/K',symbol='k',name='Boltzmann constant',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., Kelvin')
    
    GummyConstant(683,unit='lm/W',
                  name='luminous efficacy of monochromatic radiation of frequency 540e12 Hz',
                  symbol='K(cd)',html_symbol='<i>K</i><sub>cd</sub>',
                  latex_symbol='K_{' + norm('cd') + '}',add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., candela')


    _hbar = GummyConstant(_h/(2*pi),unit='J s',name='reduced Plank constant',
                          symbol='h-bar',html_symbol='&hbar;',latex_symbol='\\hbar',
                          add_symbol=True,
                          description='Plank constant over 2 pi')
    GummyConstant.alias('hbar',_hbar)
    
    
    _mu = GummyConstant(_const_dalton*unit('kg'),name='atomic mass constant',
                  symbol='m(u)',html_symbol='<i>m<sub>u</sub></i>',
                  latex_symbol='m_u',add_symbol=True,
                  description='CODATA 2018, Dalton')
    
    _me = GummyConstant(_const_me*unit('kg'),name='electron mass',symbol='m(e)',
                  html_symbol='<i>m<sub>e</sub></i>',
                  latex_symbol='m_{' + norm('e')+ '}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    _alpha = GummyConstant(_const_alpha,name='fine structure constant',
                           symbol='alpha',html_symbol='<i>&alpha;</i>',
                           latex_symbol='\\alpha',
                           add_symbol=True,
                           description='CODATA 2018')
    
    GummyConstant(_const_rydberg*unit('m**-1'),
                  name='Rydberg constant',symbol='R(\u221e)',
                  html_symbol='<i>R</i><sub>&infin;</sub>(e)',
                  latex_symbol='R_{\\infty}',
                  ascii_symbol='R(inf)',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(_const_a0*unit('m'),
                  name='Bohr radius',symbol='a(0)',
                  html_symbol='<i>a</i><sub>0</sub>',
                  latex_symbol='a_{0}',add_symbol=True,
                  description='CODATA 2018')
    
    _pm = GummyConstant(_const_proton_mass*unit('kg'),
                  name='proton mass',symbol='m(p)',
                  html_symbol='<i>m</i><sub>p</sub>',
                  latex_symbol='m_{' + norm('p') + '}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    _e0 = GummyConstant(_e**2/(2*_h*_c*_alpha),
                  name='vacuum electric permittivity',symbol='\u03b5(0)',
                  html_symbol='<i>&epsilon;</i><sub>0</sub>',
                  latex_symbol='\\epsilon_{0}',add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(1/(4*pi*_e0),name='Coulomb constant',symbol='k(e)',
                  html_symbol='<i>k</i><sub>e</sub>',
                  latex_symbol='k_{' + norm('e') + '}',add_symbol=True,
                  description='CODATA 2018').unit='N m**2/C**2'
    
    GummyConstant(2*_h*_alpha/(_c*_e**2),
                  name='vacuum magnetic permeability',
                  symbol='mu(0)',html_symbol='<i>&mu;</i><sub>0</sub',
                  latex_symbol='\\mu_{0}',add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(_h/_e**2,name='von Klitzing constant',symbol='R(K)',
                  html_symbol='<i>R</i><sub>K</sub>',
                  latex_symbol='R_{' + norm('K') + '}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(2*_e/_h,name='Josephson constant',symbol='K(J)',
                  html_symbol='<i>K</i><sub>J</sub>',
                  latex_symbol='K_{' + norm('J') + '}',
                  add_symbol=True,
                  description='CODATA 2018').unit = 'Hz/V'
    
    GummyConstant(MFraction('483597.9e9'),unit='Hz/V',
                  name='1990 conventional value of Josephson constant',
                  symbol='K(J-90)',
                  html_symbol='<i>K</i><sub>J-90</sub>',
                  latex_symbol='K_{' + norm('J-90') + '}',
                  add_symbol=True)
    
    GummyConstant(MFraction('25812.807'),unit='ohm',
                  name='1990 conventional value of von Klitzing constant',
                  symbol='R(K-90)',
                  html_symbol='<i>R</i><sub>k90</sub>',
                  latex_symbol='R_{' + norm('K-90') + '}',
                  add_symbol=True)
    
    GummyConstant(_h/(2*_e),name='magnetic flux quantum',symbol='Phi(0)',
                  html_symbol='<i>&Phi;</i><sub>0</sub>',
                  latex_symbol='\\Phi_{0}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(_h/(_me*_c),name='Compton wavelength',symbol='lambda(C)',
                  html_symbol='<i>&lambda;</i><sub>C</sub>',
                  latex_symbol='\\lambda_{' + norm('C') + '}',
                  add_symbol=True,
                  description='Compton wavelength of the electron, CODATA 2018'
                  ).unit = 'm'
    
    GummyConstant(_e**2/(4*pi*_e0*_me*_c**2),name='classical electron radius',
                  symbol='r(e)',
                  html_symbol='<i>r</i><sub>r</sub>',
                  latex_symbol='r_{' + norm('e') + '}',
                  add_symbol=True,
                  description='CODATA 2018').unit='m'
    
    GummyConstant(2*_h*_alpha/_e**2,
                  name='vacuum impedance',
                  symbol='Z(0)',
                  html_symbol='<i>Z</i><sub>0</sub>',
                  latex_symbol='Z_{0}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    _aralpha = gummy(4.001506179127,0.000000000063) 
    GummyConstant(_aralpha*_mu,
                  name='alpha particle mass',
                  symbol='m(alpha)',
                  html_symbol='<i>m</i><sub>&alpha;</sub>',
                  latex_symbol='m_{' + norm('\\alpha') + '}',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(1.00001495e-10,0.00000090e-10,unit='m',
                  name='Angstrom star',
                  symbol='\u212B*',
                  html_symbol='&#8491;<sup>*</sup>',
                  latex_symbol='\u212B^*',
                  ascii_symbol='A*',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(_e*_hbar/(2*_me),
                  name='Bohr magneton',
                  symbol='\u03bc(B)',
                  html_symbol='<i>&mu;</i><sub>B</sub>',
                  latex_symbol='\\mu_{' + norm('B') + '}',
                  ascii_symbol='mu(B)',
                  add_symbol=True,
                  description='CODATA 2018').unit = 'J/T'
    
    
    _G =  GummyConstant(_const_G*unit('m**3 kg**-1 s**-2'),
                        name='Newtonian constant of gravitation',symbol='G',
                        add_symbol=True,short_name='constant of gravitation',
                        description='CODATA 2018, gravity')
    