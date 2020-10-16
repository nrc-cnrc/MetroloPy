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
from .ummy import MFraction,_getfinfo
from .gummy import gummy,pi
from .siunits import (_const_c,_const_G,_const_h,_const_e,_const_k,
                      _const_dalton,_const_me,_const_a0,_const_proton_mass)
from .alpha_cor import alpha,rinf,ae,gh,Arh

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
    
    _NA = GummyConstant(MFraction('6.02214076e23'),unit='mol**-1',
                        name='Avogadro constant',
                        symbol='N(A)',
                        html_symbol='<i>N</i><sub>A</sub>',
                        latex_symbol='N_{' + norm('A') + '}',
                        add_symbol=True,
                        description='SI defining constant, SI Brochure 9th ed., mole')
    
    _k = GummyConstant(_const_k,unit='J/K',symbol='k',name='Boltzmann constant',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed., Kelvin')
    
    GummyConstant(683,unit='lm/W',
                  name='luminous efficacy of monochromatic radiation of frequency 540e12 Hz',
                  symbol='K(cd)',html_symbol='<i>K</i><sub>cd</sub>',
                  latex_symbol='K_{' + norm('cd') + '}',add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., candela')

    _fhbar = _h.x/(2*pi)
    if not gummy.rounding_u:
        try:
            fi,_ = _getfinfo(_fhbar)
            _fhbar = gummy(_fhbar,_fhbar*fi.rel_u,unit='J s')
        except:
            _fhbar = gummy(_fhbar,unit='J s')
    else:
        _fhbar = gummy(_fhbar,unit='J s')
    _hbar = GummyConstant(_fhbar,name='reduced Plank constant',
                          symbol='hbar',html_symbol='&hbar;',latex_symbol='\\hbar',
                          add_symbol=True,
                          description='Plank constant over 2 pi')
    GummyConstant.alias('h-bar',_hbar)
    
    
    _mu = GummyConstant(_const_dalton*unit('kg'),name='atomic mass constant',
                  symbol='m(u)',html_symbol='<i>m<sub>u</sub></i>',
                  latex_symbol='m_u',add_symbol=True,
                 description='calculated from CODATA 2018 values, Dalton')
    
    _me = GummyConstant(_const_me*unit('kg'),name='electron mass',symbol='m(e)',
                  html_symbol='<i>m<sub>e</sub></i>',
                  latex_symbol='m_{' + norm('e')+ '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values')
    
    _alpha = GummyConstant(alpha,name='fine structure constant',
                           symbol='alpha',html_symbol='<i>&alpha;</i>',
                           latex_symbol='\\alpha',
                           add_symbol=True,
                           description='calculated from CODATA 2018 values')
    
    GummyConstant(rinf*unit('m**-1'),
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
                  description='calculated from CODATA 2018 values')
    
    _pm = GummyConstant(_const_proton_mass*unit('kg'),
                  name='proton mass',symbol='m(p)',
                  html_symbol='<i>m</i><sub>p</sub>',
                  latex_symbol='m_{' + norm('p') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values')
    
    _e0 = GummyConstant(_e**2/(2*_h*_c*_alpha),
                  name='vacuum electric permittivity',symbol='\u03b5(0)',
                  html_symbol='<i>&epsilon;</i><sub>0</sub>',
                  latex_symbol='\\epsilon_{0}',
                  ascii_symbol = 'epsilon(0)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values')
    _e0.unit = 'F/m'
    
    GummyConstant(1/(4*pi*_e0),name='Coulomb constant',symbol='k(e)',
                  html_symbol='<i>k</i><sub>e</sub>',
                  latex_symbol='k_{' + norm('e') + '}',add_symbol=True,
                  description='calculated from CODATA 2018 values'
                  ).unit='N m**2/C**2'
    
    GummyConstant(2*_h*_alpha/(_c*_e**2),
                  name='vacuum magnetic permeability',
                  symbol='\u03bc(0)',
                  html_symbol='<i>&mu;</i><sub>0</sub>',
                  latex_symbol='\\mu_{0}',
                  ascii_symbol='mu(0)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values'
                  ).unit = 'N/A**2'
    
    GummyConstant(_h/_e**2,name='von Klitzing constant',symbol='R(K)',
                  html_symbol='<i>R</i><sub>K</sub>',
                  latex_symbol='R_{' + norm('K') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values')
    
    GummyConstant(2*_e/_h,name='Josephson constant',symbol='K(J)',
                  html_symbol='<i>K</i><sub>J</sub>',
                  latex_symbol='K_{' + norm('J') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values'
                  ).unit = 'Hz/V'
    
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
                  description='calculated from SI Brochure 9th ed. values'
                  ).unit = 'Wb'
    
    GummyConstant(_h/(_me*_c),name='Compton wavelength',symbol='lambda(C)',
                  html_symbol='<i>&lambda;</i><sub>C</sub>',
                  latex_symbol='\\lambda_{' + norm('C') + '}',
                  add_symbol=True,
                  description='Compton wavelength of the electron, calculated from CODATA 2018 values'
                  ).unit = 'm'
    
    GummyConstant(_e**2/(4*pi*_e0*_me*_c**2),name='classical electron radius',
                  symbol='r(e)',
                  html_symbol='<i>r</i><sub>r</sub>',
                  latex_symbol='r_{' + norm('e') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values'
                  ).unit='m'
    
    GummyConstant(2*_h*_alpha/_e**2,
                  name='vacuum impedance',
                  symbol='Z(0)',
                  html_symbol='<i>Z</i><sub>0</sub>',
                  latex_symbol='Z_{0}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 2*h*alpha/e**2'
                  ).unit = 'ohm'
    
    _aralpha = gummy(4.001506179127,0.000000000063) 
    GummyConstant(_aralpha*_mu,
                  name='alpha particle mass',
                  symbol='m(alpha)',
                  html_symbol='<i>m</i><sub>&alpha;</sub>',
                  latex_symbol='m_{' + norm('\\alpha') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as Ar(alpha)*m(u)'
                  )
    
    GummyConstant(1.00001495e-10,0.00000090e-10,unit='m',
                  name='Angstrom star',
                  symbol='\u212B*',
                  html_symbol='&#8491;<sup>*</sup>',
                  latex_symbol='\u212B^*',
                  ascii_symbol='A*',
                  add_symbol=True,
                  description='CODATA 2018')
    
    _muB = GummyConstant(_e*_hbar/(2*_me),
                         name='Bohr magneton',
                         symbol='\u03bc(B)',
                         html_symbol='<i>&mu;</i><sub>B</sub>',
                         latex_symbol='\\mu_{' + norm('B') + '}',
                         ascii_symbol='mu(B)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as e*hbar/(2*m(e))')
    _muB.unit = 'J/T'
                  
    _emma = GummyConstant(ae,
                          name='electron magnetic moment anomaly',
                          symbol='a(e)',
                          html_symbol='<i>a</i><sub>e</sub>',
                          latex_symbol='a_{' + norm('e') + '}',
                          add_symbol=True,
                          description='CODATA 2018')
    
    _ge = GummyConstant(-2*(1+_emma),
                        name='electron g factor',
                        symbol='g(e-)',
                        html_symbol='<i>g</i><sub>e<sup>-</sup></sub>',
                        latex_symbol='g_{' + norm('e') + '^{-}}',
                        add_symbol=True,
                        description='calculated from CODATA 2018 values as -2*(1-a(e))')
    
    GummyConstant(_ge*_muB/_hbar,
                  name='electron gyromagnetic ratio',
                  symbol='\u03B3(e)',
                  html_symbol='<i>&gamma;</i><sub>e</sub>',
                  latex_symbol='\\gamma_{' + norm('e') + '}',
                  ascii_symbol='gamma(e)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as g(e)*mu(B)/hbar where g(e) = -2*(1-a(e)) and mu(B) = e*hbar/(2*m(e))')
    
    GummyConstant(_ge*_muB/2,
                  name='electron magnetic moment',
                  symbol='\u03bc(e)',
                  html_symbol='<i>&mu;</i><sub>e</sub>',
                  latex_symbol='\\mu_{' + norm('e') + '}',
                  ascii_symbol = 'mu(e)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as g(e)*mu(B)/2 where g(e) = -2*(1-a(e)) and mu(B) = e*hbar/(2*m(e))'
                  ).unit = 'J/T'
                  
    GummyConstant(2*_e**2/_h,
                  name='conductance quantum',
                  symbol='G(0)',
                  html_symbol='<i>G</i><sub>0</sub>',
                  latex_symbol='G{0}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*e**2/h'
                  ).unit = 'S'
                  
                  
    (
     _angstar,
     _copperx
     ) = gummy.create([1.00001495e-10,1.00207697e-13],
                      [0.00000090e-10,0.00000028e-13],
                      correlation_matrix = [[1,0.00039],
                                            [0.00039,1]])
    GummyConstant(_angstar*unit('m'),unit='m',
                  name='Angstrom star',
                  symbol='\u212B*',
                  html_symbol='&#8491;<sup>*</sup>',
                  latex_symbol='\u212B^*',
                  ascii_symbol='A*',
                  add_symbol=True,
                  description='CODATA 2018')
    
    GummyConstant(_copperx*unit('m'),unit='m',
                  name='Copper x unit',
                  symbol='xu(CuK\u03B11)',
                  html_symbol='xu(CuK&alpha;<sub>1</sub>)',
                  latex_symbol='xu(CuK\\alpha_{1})',
                  ascii_symbol='xu(CuKalpha1)',
                  short_name='xu(Cu)',
                  add_symbol=True,
                  description='CODATA 2018')
    
    _ard = gummy(2.013553212745,0.000000000040)
    GummyConstant(_mu*_ard,
                  name='deuteron mass',
                  symbol='m(d)',
                  html_symbol='<i>m</i><sub>d</sub>',
                  latex_symbol='m_{' + norm('d') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(u)*Ar(d)')
    
    GummyConstant(0.8574382338,0.0000000022,
                  name='deuteron g factor',
                  symbol='g(d)',
                  html_symbol='<i>g</i><sub>d</sub>',
                  latex_symbol='g_{' + norm('d') + '}',
                  add_symbol=True,
                  description='CODATA 2018')    
    
    _G =  GummyConstant(_const_G*unit('m**3 kg**-1 s**-2'),
                        name='Newtonian constant of gravitation',symbol='G',
                        add_symbol=True,short_name='constant of gravitation',
                        description='CODATA 2018, gravity')
    
    GummyConstant(_e*_NA,
                  name='Faraday constant',
                  symbol='F',
                  html_symbol='<i>F</i>',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as e*N(A)')
    
    GummyConstant(2*_h*_c**2,
                  name='first radiation constant for spectral radiance',
                  symbol='c(1L)',
                  html_symbol='<i>c</i><sub>1L</sub>',
                  latex_symbol='c_{1' + norm('L') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*h*c**2'
                  ).unit = 'W m**2/sr'
    
    GummyConstant(2*pi*_h*_c**2,
                  name='first radiation constant',
                  symbol='c(1)',
                  html_symbol='<i>c</i><sub>1</sub>',
                  latex_symbol='c_{1}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*pi*h*c**2'
                  ).unit = 'W m**2'
    
    GummyConstant(_h*_c/_k,
                  name='second radiation constant',
                  symbol='c(2)',
                  html_symbol='<i>c</i><sub>2</sub>',
                  latex_symbol='c_{2}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as h*c/k'
                  ).unit = 'm K'
                  
    GummyConstant(_me*_c**2*_alpha**2,
                  name='Hartree energy',
                  symbol='E(h)',
                  html_symbol='<i>E</i><sub>h</sub>',
                  latex_symbol='E_{' + norm('h') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(e)*c**2*alpha**2 where m(e) = m(u)*Ar(e) and m(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)'
                  ).unit = 'J'
    
    GummyConstant(gh,
                  name='helion g factor',
                  symbol='g(h)',
                  html_symbol='<i>g</i><sub>h</sub>',
                  latex_symbol='g_{' + norm('h') + '}',
                  add_symbol=True,
                  description='CODATA 2018'
                  )
    
    GummyConstant(_mu*Arh,
                  name='helion mass',
                  symbol='m(h)',
                  html_symbol='<i>m</i><sub>h</sub>',
                  latex_symbol='m_{' + norm('h') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(u)*Ar(h) where m(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)'
                  )
    
    GummyConstant(5.996743e-5,0.000010e-5,
                  name='helion shielding shift',
                  symbol='\u03c3(h)',
                  html_symbol='<i>&sigma;</i><sub>h</sub>',
                  latex_symbol='\\sigma_{' + norm('h') + '}',
                  ascii_symbol='sigma(h)',
                  add_symbol=True,
                  description='CODATA 2018'
                  )