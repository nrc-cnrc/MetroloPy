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

from numpy import sqrt
from .constant import GummyConstant
from .unit import unit
from .ummy import MFraction
from .gummy import gummy
from .constcom import h,c,e,hbar,k,dalton,me,mp,a0,G,pi,euler,sqrt2
from .constcom import (alph,aral,ryd,ae,rd,ghn,arh,sigmah,mmu,amu,
                       arn,gnn,ard,gdn,gp,sigmapp,mtu,gtn)
                       

def norm(x):
    return GummyConstant.latex_norm(x)

with GummyConstant._builtin():
    
    GummyConstant(pi, 
                  name='float64 representation of pi',
                  symbol='\u03c0',
                  html_symbol='<i>&pi;</i>',
                  latex_symbol='\\pi',
                  ascii_symbol='pi',
                  add_symbol=True,
                  description='numpy.pi',
                  )
    
    _euler = GummyConstant(euler, 
                           name='float64 representation of Euler\'s number e',
                           symbol='e',
                           html_symbol='<i>e</i>',
                           add_symbol=False,
                           short_name='Euler',
                           description='numpy.e',
                           )
    GummyConstant.alias('euler',_euler)
    GummyConstant.alias('math e',_euler)
    GummyConstant.alias('math_e',_euler)
    
    GummyConstant(sqrt2, 
                  name='float64 representation of sqrt(2)',
                  symbol='sqrt(2)',
                  latex_symbol='\\sqrt{2}',
                  add_symbol=True,
                  description='numpy.sqrt(2)',
                  )
    
    GummyConstant(9192631770,unit='Hz', 
                  name='hyperfine transition frequency of Cs-133',
                  symbol='\u0394\u03bd(Cs)',
                  html_symbol='&Delta;<i>&nu;</i><sub>Cs</sub>',
                  latex_symbol='\\Delta\\nu_{' + norm('Cs') + '}',
                  ascii_symbol='hf(Cs)',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed.,Caesium, Cesium',
                  _key='nucs')
    
    _h = GummyConstant(h,unit='J s',
                       name='Planck constant',
                       symbol='h',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.',
                       _key='h')
    
    _c = GummyConstant(c,unit='m/s',
                       name='speed of light in vacuum',
                       symbol='c',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.',
                       _key='c')
    
    _e = GummyConstant(e,unit='C',
                       name='elementary charge',
                       symbol='e',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.',
                       _key='e')
    
    _NA = GummyConstant(MFraction('6.02214076e23'),unit='mol**-1',
                        name='Avogadro constant',
                        symbol='N(A)',
                        html_symbol='<i>N</i><sub>A</sub>',
                        latex_symbol='N_{' + norm('A') + '}',
                        add_symbol=True,
                        description='SI defining constant, SI Brochure 9th ed.',
                        _key='na')
    
    _k = GummyConstant(k,unit='J/K',
                       name='Boltzmann constant',
                       symbol='k',
                       add_symbol=True,
                       description='SI defining constant, SI Brochure 9th ed.',
                       _key='k')
    
    GummyConstant(683,unit='lm/W',
                  name='luminous efficacy',
                  symbol='K(cd)',
                  html_symbol='<i>K</i><sub>cd</sub>',
                  latex_symbol='K_{' + norm('cd') + '}',
                  add_symbol=True,
                  description='SI defining constant, SI Brochure 9th ed., candela',
                  _key='kcd')

    _hbar = GummyConstant(hbar*unit('J s'),
                          name='reduced Planck constant',
                          symbol='hbar',
                          html_symbol='&hbar;',
                          latex_symbol='\\hbar',
                          add_symbol=True,
                          description='calculated from SI Brochure 9th ed. value as h/(2*pi)',
                          _key='hbar')
    GummyConstant.alias('h-bar',_hbar)
    
    
    _mu = GummyConstant(dalton*unit('kg'),
                        name='atomic mass constant',
                        symbol='m(u)',
                        html_symbol='<i>m<sub>u</sub></i>',
                        latex_symbol='m_u',
                        add_symbol=True,
                        description='calculated from CODATA 2018 values as 2*R(inf)*h/(Ar(e)*c*alpha**2\nsee also the Unit dalton',
                        _key='u')
    
    _me = GummyConstant(me*unit('kg'),
                        name='electron mass',
                        symbol='m(e)',
                        html_symbol='<i>m<sub>e</sub></i>',
                        latex_symbol='m_{' + norm('e')+ '}',
                        add_symbol=True,
                        description='calculated from CODATA 2018 values as 2*R(inf)*h/c*alpha**2',
                        _key='me')
    
    _alpha = GummyConstant(alph,
                           name='fine-structure constant',
                           symbol='alpha',
                           html_symbol='<i>&alpha;</i>',
                           latex_symbol='\\alpha',
                           add_symbol=True,
                           description='CODATA 2018',
                           _key='alph')
    
    GummyConstant(ryd*unit('m**-1'),
                   name='Rydberg constant',
                   symbol='R(\u221e)',
                   html_symbol='<i>R</i><sub>&infin;</sub>',
                   latex_symbol='R_{\\infty}',
                   ascii_symbol='R(inf)',
                   add_symbol=True,
                   description='CODATA 2018',
                   _key='ryd')
    
    _a0= GummyConstant(a0*unit('m'),
                       name='Bohr radius',
                       symbol='a(0)',
                       html_symbol='<i>a</i><sub>0</sub>',
                       latex_symbol='a_{0}',
                       add_symbol=True,
                       description='calculated from CODATA 2018 values as alpha/(4*pi*R(inf))',
                       _key='bohrrada0')
        
    _mp = GummyConstant(mp*unit('kg'),
                  name='proton mass',
                  symbol='m(p)',
                  html_symbol='<i>m</i><sub>p</sub>',
                  latex_symbol='m_{' + norm('p') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as p_e_mass_ratio*m(e) where\nm(e) = 2*R(inf)*h/c*alpha**2',
                  _key='mp')
    
    _e0 = GummyConstant(_e**2/(2*_h*_c*_alpha),
                  name='vacuum electric permittivity',
                  symbol='\u03b5(0)',
                  html_symbol='<i>&epsilon;</i><sub>0</sub>',
                  latex_symbol='\\epsilon_{0}',
                  ascii_symbol = 'epsilon(0)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as e**2/(s*h*c*alpha)',
                  _key='ep0')
    _e0.unit = 'F/m'
    
    GummyConstant(1/(4*pi*_e0),name='Coulomb constant',
                  symbol='k(e)',
                  html_symbol='<i>k</i><sub>e</sub>',
                  latex_symbol='k_{' + norm('e') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 1/(4*pi*e(0)) where\ne(0) = e**2/(s*h*c*alpha)'
                  ).unit='N m**2/C**2'
    
    GummyConstant(2*_h*_alpha/(_c*_e**2),
                  name='vacuum magnetic permeability',
                  symbol='\u03bc(0)',
                  html_symbol='<i>&mu;</i><sub>0</sub>',
                  latex_symbol='\\mu_{0}',
                  ascii_symbol='mu(0)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 2*h*alpha/(c*e**2)',
                  _key='mu0'
                  ).unit = 'N/A**2'
    
    GummyConstant(_h/_e**2,name='von Klitzing constant',
                  symbol='R(K)',
                  html_symbol='<i>R</i><sub>K</sub>',
                  latex_symbol='R_{' + norm('K') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as h/e**2',
                  _key='rk')
    
    GummyConstant(2*_e/_h,name='Josephson constant',
                  symbol='K(J)',
                  html_symbol='<i>K</i><sub>J</sub>',
                  latex_symbol='K_{' + norm('J') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*e/h',
                  _key='kjos'
                  ).unit = 'Hz/V'
    
    GummyConstant(MFraction('483597.9e9'),unit='Hz/V',
                  name='1990 conventional value of Josephson constant',
                  symbol='K(J-90)',
                  html_symbol='<i>K</i><sub>J-90</sub>',
                  latex_symbol='K_{' + norm('J-90') + '}',
                  add_symbol=True,
                  _key='kj90')
    
    GummyConstant(MFraction('25812.807'),unit='ohm',
                  name='1990 conventional value of von Klitzing constant',
                  symbol='R(K-90)',
                  html_symbol='<i>R</i><sub>K-90</sub>',
                  latex_symbol='R_{' + norm('K-90') + '}',
                  add_symbol=True,
                  _key='rk90')
    
    GummyConstant(_h/(2*_e),
                  name='magnetic flux quantum',
                  symbol='Phi(0)',
                  html_symbol='<i>&Phi;</i><sub>0</sub>',
                  latex_symbol='\\Phi_{0}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as h/(2*e)',
                  _key='flxquhs2e'
                  ).unit = 'Wb'
    
    GummyConstant(_h/(_me*_c),
                  name='Compton wavelength',
                  symbol='\u03bb(C)',
                  html_symbol='<i>&lambda;</i><sub>C</sub>',
                  latex_symbol='\\lambda_{' + norm('C') + '}',
                  ascii_symbol='lambda(C)',
                  add_symbol=True,
                  description='Compton wavelength of the electron, calculated from CODATA 2018 values as\nh/(m(e)*c) where m(e) = 2*R(inf)*h/c*alpha**2',
                  _key='ecomwl'
                  ).unit = 'm'
    
    GummyConstant(_e**2/(4*pi*_e0*_me*_c**2),
                  name='classical electron radius',
                  symbol='r(e)',
                  html_symbol='<i>r</i><sub>r</sub>',
                  latex_symbol='r_{' + norm('e') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as e**2/(4*pi*e(0)*m(e)*c**2) where\nm(e) = 2*R(inf)*h/c*alpha**2 and e(0) = e**2/(s*h*c*alpha)',
                  _key='re'
                  ).unit='m'
    
    GummyConstant(2*_h*_alpha/_e**2,
                  name='vacuum impedance',
                  symbol='Z(0)',
                  html_symbol='<i>Z</i><sub>0</sub>',
                  latex_symbol='Z_{0}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 2*h*alpha/e**2',
                  _key='z0'
                  ).unit = 'ohm'
    
    GummyConstant(aral*_mu,
                  name='alpha particle mass',
                  symbol='m(\u03b1)',
                  html_symbol='<i>m</i><sub>&alpha;</sub>',
                  latex_symbol='m_{' + norm('\u03b1') + '}',
                  ascii_symbol='m(alpha)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as Ar(alpha)*m(u)',
                  _key='mal'
                  )
    
    _muB = GummyConstant(_e*_hbar/(2*_me),
                         name='Bohr magneton',
                         symbol='\u03bc(B)',
                         html_symbol='<i>&mu;</i><sub>B</sub>',
                         latex_symbol='\\mu_{' + norm('B') + '}',
                         ascii_symbol='mu(B)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as e*hbar/(2*m(e))',
                         _key='mub')
    _muB.unit = 'J/T'
                  
    _emma = GummyConstant(ae,
                          name='electron magnetic moment anomaly',
                          symbol='a(e)',
                          html_symbol='<i>a</i><sub>e</sub>',
                          latex_symbol='a_{' + norm('e') + '}',
                          add_symbol=True,
                          description='CODATA 2018',
                          _key='ae')
    
    _ge = GummyConstant(-2*(1+_emma),
                        name='electron g factor',
                        symbol='g(e-)',
                        html_symbol='<i>g</i><sub>e<sup>-</sup></sub>',
                        latex_symbol='g_{' + norm('e') + '^{-}}',
                        add_symbol=True,
                        description='calculated from CODATA 2018 values as -2*(1-a(e))',
                        _key='gem')
    
    GummyConstant(-_ge*_muB/_hbar,
                  name='electron gyromagnetic ratio',
                  symbol='\u03B3(e)',
                  html_symbol='<i>&gamma;</i><sub>e</sub>',
                  latex_symbol='\\gamma_{' + norm('e') + '}',
                  ascii_symbol='gamma(e)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as -g(e)*mu(B)/hbar where\ng(e) = -2*(1-a(e)) and mu(B) = e*hbar/(2*m(e))',
                  _key='gammae')
    
    GummyConstant(_ge*_muB/2,
                  name='electron magnetic moment',
                  symbol='\u03bc(e)',
                  html_symbol='<i>&mu;</i><sub>e</sub>',
                  latex_symbol='\\mu_{' + norm('e') + '}',
                  ascii_symbol = 'mu(e)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as g(e)*mu(B)/2 where\ng(e) = -2*(1-a(e)) and mu(B) = e*hbar/(2*m(e))',
                  _key='muem'
                  ).unit = 'J/T'
                  
    GummyConstant(2*_e**2/_h,
                  name='conductance quantum',
                  symbol='G(0)',
                  html_symbol='<i>G</i><sub>0</sub>',
                  latex_symbol='G_{0}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*e**2/h',
                  _key='conqu2e2sh'
                  ).unit = 'S'
                  
                  
    (
     _angstar,
     _copperx,
     _molyx,
     _a
     ) = gummy.create([1.00001495e-10,
                       1.00207697e-13,
                       1.00209952e-13,
                       5.431020511e-10],
                      [0.00000090e-10,
                       0.00000028e-13,
                       0.00000053e-13,
                       0.000000089e-10],
                      correlation_matrix = [[1,0.00039,0.00100,0.00818],
                                            [0.00039,1,0.00067,0.01272],
                                            [0.00100,0.00067,1,0.01398],
                                            [0.00818,0.01272,0.01398,1]])
    GummyConstant(_angstar*unit('m'),
                  name='Angstrom star',
                  symbol='\u212B*',
                  html_symbol='&#8491;<sup>*</sup>',
                  latex_symbol='\u212B^{*}',
                  ascii_symbol='A*',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='angstar')
    
    GummyConstant(_copperx*unit('m'),
                  name='copper x unit',
                  symbol='xu(CuK\u03B11)',
                  html_symbol='xu(CuK&alpha;<sub>1</sub>)',
                  latex_symbol=norm('xu') + '(' + norm('CuK\u03b1') + '_{1})',
                  ascii_symbol='xu(CuKalpha1)',
                  short_name='xu(Cu)',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='xucukalph1'
                  )
    
    GummyConstant(_molyx*unit('m'),
                  name='molybdenum x unit',
                  symbol='xu(MoK\u03B11)',
                  html_symbol='xu(MoK&alpha;<sub>1</sub>)',
                  latex_symbol=norm('xu') + '(' + norm('MoK\u03b1') + '_{1})',
                  ascii_symbol='xu(MoKalpha1)',
                  short_name='xu(Mo)',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='xumokalph1')
    
    GummyConstant(_a*unit('m'),
                  name='lattice parameter of silicon',
                  symbol='a',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='asil')
    
    GummyConstant((_a/sqrt(8))*unit('m'),
                  name='lattice spacing of ideal Si (220)',
                  symbol='d(220)',
                  html_symbol='d<sub>220</sub>)',
                  latex_symbol='d_{220}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as a/sqrt(8)',
                  _key='d220sil')
    
    GummyConstant(_mu*ard,
                  name='deuteron mass',
                  symbol='m(d)',
                  html_symbol='<i>m</i><sub>d</sub>',
                  latex_symbol='m_{' + norm('d') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(u)*Ar(d)',
                  _key='md')
    
    GummyConstant(gdn,
                  name='deuteron g factor',
                  symbol='g(d)',
                  html_symbol='<i>g</i><sub>d</sub>',
                  latex_symbol='g_{' + norm('d') + '}',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='gdn')   
    
    _muN = GummyConstant(_e*_hbar/(2*_mp),
                         name='nuclear magneton',
                         symbol='\u03bc(N)',
                         html_symbol='<i>&mu;</i><sub>N</sub>',
                         latex_symbol='\\mu_{' + norm('N') + '}',
                         ascii_symbol='mu(N)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as e*hbar/(2*m(p))',
                         _key='mun')
    _muN.unit = 'J/T'

    GummyConstant(gdn*_muN,
                  name='deuteron magnetic moment',
                  symbol='\u03bc(d)',
                  html_symbol='<i>&mu;</i><sub>d</sub>',
                  latex_symbol='\\mu_{' + norm('d') + '}',
                  ascii_symbol='mu(d)',
                  add_symbol=True,
                  description='calculate from CODATA 2018 values as g(d)*mu(N)',
                  _key='mud') 
    
    GummyConstant(rd,
                  name='deuteron rms charge radius',
                  symbol='r(d)',
                  html_symbol='<i>r</i><sub>d</sub>',
                  latex_symbol='r_{' + norm('d') + '}',
                  add_symbol=True,
                  _key='rd'
                  )
    
    _G =  GummyConstant(G*unit('m**3 kg**-1 s**-2'),
                        name='Newtonian constant of gravitation',symbol='G',
                        add_symbol=True,short_name='constant of gravitation',
                        description='CODATA 2018, gravity',
                        _key='bg')
    
    GummyConstant(_e*_NA,
                  name='Faraday constant',
                  symbol='F',
                  html_symbol='<i>F</i>',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as e*N(A)',
                  _key='f')
    
    GummyConstant(2*_h*_c**2,
                  name='first radiation constant for spectral radiance',
                  symbol='c(1L)',
                  html_symbol='<i>c</i><sub>1L</sub>',
                  latex_symbol='c_{1' + norm('L') + '}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*h*c**2',
                  _key='c1l'
                  ).unit = 'W m**2/sr'
    
    GummyConstant(2*pi*_h*_c**2,
                  name='first radiation constant',
                  symbol='c(1)',
                  html_symbol='<i>c</i><sub>1</sub>',
                  latex_symbol='c_{1}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as 2*pi*h*c**2',
                  _key='eqc11strc'
                  ).unit = 'W m**2'
    
    GummyConstant(_h*_c/_k,
                  name='second radiation constant',
                  symbol='c(2)',
                  html_symbol='<i>c</i><sub>2</sub>',
                  latex_symbol='c_{2}',
                  add_symbol=True,
                  description='calculated from SI Brochure 9th ed. values as h*c/k',
                  _key='c22ndrc'
                  ).unit = 'm K'
                  
    GummyConstant(_me*_c**2*_alpha**2,
                  name='Hartree energy',
                  symbol='E(h)',
                  html_symbol='<i>E</i><sub>h</sub>',
                  latex_symbol='E_{' + norm('h') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(e)*c**2*alpha**2\nwhere m(e) = m(u)*Ar(e) and m(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)',
                  _key='hr'
                  ).unit = 'J'
    
    _gh = GummyConstant(ghn,
                        name='helion g factor',
                        symbol='g(h)',
                        html_symbol='<i>g</i><sub>h</sub>',
                        latex_symbol='g_{' + norm('h') + '}',
                        add_symbol=True,
                        description='CODATA 2018',
                        _key='ghn'
                        )
    
    GummyConstant(_gh*_muN/2,
                  name='helion magnetic moment',
                  symbol='\u03bc(h)',
                  html_symbol='<i>&mu;</i><sub>h</sub>',
                  latex_symbol='\\mu_{' + norm('h') + '}',
                  ascii_symbol = 'mu(h)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as g(h)*mu(B)/2 where\nmu(B) = e*hbar/(2*m(e))',
                  _key='muh'
                  ).unit = 'J/T'
    
    GummyConstant(_mu*arh,
                  name='helion mass',
                  symbol='m(h)',
                  html_symbol='<i>m</i><sub>h</sub>',
                  latex_symbol='m_{' + norm('h') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(u)*Ar(h) where\nm(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)',
                  _key='mh'
                  )
    
    GummyConstant(sigmah,
                  name='helion shielding shift',
                  symbol='\u03c3(h)',
                  html_symbol='<i>&sigma;</i><sub>h</sub>',
                  latex_symbol='\\sigma_{' + norm('h') + '}',
                  ascii_symbol='sigma(h)',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='sigmah'
                  )
    
    GummyConstant(_mu*_NA,
                  name='molar mass constant',
                  symbol='M(u)',
                  html_symbol='<i>M</i><sub>u</sub>',
                  latex_symbol='M_{' + norm('u') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as m(u)*N(A)',
                  _key='mu'
                  )
    
    GummyConstant(_NA*_k,
                  name='molar gas constant',
                  symbol='R',
                  html_symbol='<i>R</i>',
                  latex_symbol='R',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as N(A)*k',
                  _key='r'
                  )
    
    GummyConstant(12*_NA*_mu,
                  name='molar mass of carbon-12',
                  symbol='M(12C)',
                  html_symbol='<i>M</i>(<sup>12</sup>C)',
                  latex_symbol='M(^{12}' + norm('C') + ')',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 12*m(u)',
                  _key='mm12c'
                  )
    
    GummyConstant(_NA*_a**3/8,
                  name='molar volume of silicon',
                  symbol='Vm(Si)',
                  html_symbol='<i>V</i><sub>m</sub>(Si)',
                  latex_symbol='V_{' + norm('m') + '}(' + norm('Si') + ')',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as N(A)*a**3/8',
                  _key='mvolsil'
                  )
    
    GummyConstant(100000,unit='Pa',
                  name='standard-state pressure',
                  symbol='ssp',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='stdspr'
                  )
    
    GummyConstant(101325,unit='Pa',
                  name='standard atmosphere',
                  symbol='atm',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='stdatm'
                  )
        
    _mmu = GummyConstant(mmu*unit('kg'),
                         name='muon mass',
                         symbol='m(\u03bc)',
                         html_symbol='<i>m</i><sub>&mu;</sub>',
                         latex_symbol='m_{' + norm('\u03bc') + '}',
                         ascii_symbol='m(mu)',
                         add_symbol=True,
                         description='CODATA 2018',
                         _key='mmu'
                         )
    
    GummyConstant(amu,
                  name='muon magnetic moment anomaly',
                  symbol='a(\u03bc)',
                  html_symbol='<i>a</i><sub>&mu;</sub>',
                  latex_symbol='a_{' + norm('\u03bc') + '}',
                  ascii_symbol='a(mu)',
                  add_symbol=True,
                  description='CODATA 2018',
                  _key='amu'
                  )
    
    _gmu = GummyConstant(-2*(1+amu),
                  name='muon g factor',
                  symbol='g(\u03bc-)',
                  html_symbol='<i>g</i><sub>&mu;<sup>-</sup></sub>',
                  latex_symbol='g_{' + norm('\u03bc') + '^-}',
                  ascii_symbol='g(mu-)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as -2*(1 + a(mu))',
                  _key='gmum'
                  )
    
    GummyConstant(_gmu*_e*_hbar/(4*_mmu),
                  name='muon magnetic moment',
                  symbol='\u03bc(\u03bc)',
                  html_symbol='<i>&mu;</i><sub>&mu;</sub>',
                  latex_symbol='\\mu_{' + norm('\u03bc') + '}',
                  ascii_symbol = 'mu(mu)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as g(mu-)*e*hbar/(4*m(mu)))',
                  _key='mumum'
                  ).unit = 'J/T'
    
    GummyConstant(_h/(_mmu*_c),
                  name='muon Compton wavelength',
                  symbol='\u03bb(\u03bc)',
                  html_symbol='<i>&lambda;</i><sub>&mu;</sub>',
                  latex_symbol='\\lambda_{' + norm('\u03bc') + '}',
                  ascii_symbol = 'lambda(mu)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as\nh/(m(mu)*c)',
                  _key='mcomwl'
                  ).unit = 'm'
    
    _mn = GummyConstant(arn*_mu,
                        name='neutron mass',
                        symbol='m(n)',
                        html_symbol='<i>m</i><sub>n</sub>',
                        latex_symbol='m_{' + norm('n') + '}',
                        add_symbol=True,
                        description='CODATA 2018 calculated from CODATA 2018 values as Ar(n)*m(u) where\nm(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)',
                        _key='mn'
                        )
    
    _gn = GummyConstant(gnn,
                        name='neutron g factor',
                        symbol='g(n)',
                        html_symbol='<i>g</i><sub>n</sub>',
                        latex_symbol='g_{' + norm('g') + '}',
                        add_symbol=True,
                        description='CODATA 2018',
                        _key='gnn'
                        )
        
    _mun = GummyConstant(_gn*_muN/2,
                         name='neutron magnetic moment',
                         symbol='\u03bc(n)',
                         html_symbol='<i>&mu;</i><sub>n</sub>',
                         latex_symbol='\\mu_{' + norm('n') + '}',
                         ascii_symbol = 'mu(n)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as g(n)*mu(N)/2 where\nmu(N) = e*hbar/(2*m(p))',
                         _key='munn'
                         )
    _mun.unit = 'J/T'
    
    GummyConstant(-2*_mun/_hbar,
                  name='neutron gyromagnetic ratio',
                  symbol='\u03b3(n)',
                  html_symbol='<i>&gamma;</i><sub>n</sub>',
                  latex_symbol='\\gamma_{' + norm('n') + '}',
                  ascii_symbol = 'gamma(n)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as -2*mu(n)/hbar where\nmu(n) =  g(n)*mu(N)/2 and mu(N) = e*hbar/(2*m(p))',
                  _key='gamman'
                  )
    
    GummyConstant(_h/(_mn*_c),
                  name='neutron Compton wavelength',
                  symbol='\u03bb(n)',
                  html_symbol='<i>&lambda;</i><sub>n</sub>',
                  latex_symbol='\\lambda_{' + norm('n') + '}',
                  ascii_symbol = 'lambda(n)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as\nh/(m(n)*c) where m(e) = Ar(n)*m(u)',
                  _key='ncomwl'
                  ).unit = 'm'
    
    _gp = GummyConstant(gp,
                        name='proton g factor',
                        symbol='g(p)',
                        html_symbol='<i>g</i><sub>p</sub>',
                        latex_symbol='g_{' + norm('p') + '}',
                        add_symbol=True,
                        description='CODATA 2018',
                        _key='gp'
                        )
    
    _mup = GummyConstant(_gp*_muN/2,
                         name='proton magnetic moment',
                         symbol='\u03bc(p)',
                         html_symbol='<i>&mu;</i><sub>p</sub>',
                         latex_symbol='\\mu_{' + norm('p') + '}',
                         ascii_symbol = 'mu(p)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as g(p)*mu(N)/2 where\nmu(N) = e*hbar/(2*m(p))',
                         _key='mup'
                         )
    _mup.unit = 'J/T'
    
    GummyConstant(2*_mup/_hbar,
                  name='proton gyromagnetic ratio',
                  symbol='\u03b3(p)',
                  html_symbol='<i>&gamma;</i><sub>p</sub>',
                  latex_symbol='\\gamma_{' + norm('p') + '}',
                  ascii_symbol = 'gamma(p)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 2*mu(p)/hbar where\nmu(p) =  g(p)*mu(N)/2 and mu(N) = e*hbar/(2*m(p))',
                  _key='gammap'
                  )
    
    GummyConstant(sigmapp,
                  name='proton magnetic shielding correction',
                  symbol='\u03c3\u2032(p)',
                  html_symbol='<i>&sigma;</i>&prime;<sub>p</sub>',
                  latex_symbol='\\sigma\'_{' + norm('p') + '}',
                  ascii_symbol='sigma\'(p)',
                  add_symbol=True,
                  _key='sigmapp'
                  )
    
    GummyConstant(mtu*_mu,
                  name='triton mass',
                  symbol='m(t)',
                  html_symbol='<i>m</i><sub>t</sub>',
                  latex_symbol='m_{' + norm('t') + '}',
                  add_symbol=True,
                  description='CODATA 2018 calculated from CODATA 2018 values as Ar(t)*m(u) where\nm(u) = 2*R(inf)*h/(Ar(e)*c*alpha**2)',
                  _key='mt'
                  )
    
    _gtn = GummyConstant(gtn,
                         name='triton g factor',
                         symbol='g(t)',
                         html_symbol='<i>g</i><sub>t</sub>',
                         latex_symbol='g_{' + norm('t') + '}',
                         add_symbol=True,
                         _key='gtn'
                         )
    
    _mut = GummyConstant(_gtn*_muN/2,
                         name='triton magnetic moment',
                         symbol='\u03bc(t)',
                         html_symbol='<i>&mu;</i><sub>t</sub>',
                         latex_symbol='\\mu_{' + norm('t') + '}',
                         ascii_symbol = 'mu(t)',
                         add_symbol=True,
                         description='calculated from CODATA 2018 values as g(t)*mu(N)/2 where\nmu(N) = e*hbar/(2*m(t))',
                         _key='mut'
                         )
    _mut.unit = 'J/T'
    
    GummyConstant((8*pi/3)*_alpha**4*_a0**2,
                  name='Thomson cross section',
                  symbol='\u03c3(e)',
                  html_symbol='<i>&sigma;</i><sub>e</sub>',
                  latex_symbol='\\sigma_{' + norm('e') + '}',
                  ascii_symbol = 'sigma(e)',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as (8*pi/3)*alpha**4*a0**2',
                  _key='eqsigmae'
                  ).unit = 'm**2'
    
    _wm =  GummyConstant(0.22290,0.00030,
                         name='weak mixing angle',
                         symbol='sin2(\u03f4(W))',
                         html_symbol='sin<sup>2</sup> <i>&Theta;</i><sub>W</sub>',
                         latex_symbol='sin^{2}\\Theta_{' + norm('W') + '}',
                         ascii_symbol='sin2(Theta(W))',
                         add_symbol=True,
                         _key='sin2th'
                         )
    
    GummyConstant(sqrt(1 - _wm),
                  name='W to Z mass ratio',
                  symbol='m(W)/m(Z)',
                  html_symbol='<i>m</i><sub>W</sub>/<i>m</i><sub>Z</sub>',
                  latex_symbol='m_{' + norm('W') + '}/m_{' + norm('Z') + '}',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as sqrt(1 - sin2(Theta(W)))',
                  _key='rmwmz'
                  )
    
    GummyConstant(2*pi**5*_k**4/(15*_h**3*_c**2),
                  name='Stefan-Boltzmann constant',
                  symbol='\u03c3',
                  html_symbol='<i>&sigma;</i>',
                  latex_symbol='\\sigma',
                  ascii_symbol = 'sigma',
                  add_symbol=True,
                  description='calculated from CODATA 2018 values as 2*pi*5*k**4/(15*h**3*c**2)',
                  _key='sigma'
                  ).unit = 'W m**-2 K**-4'

    GummyConstant(1.1663787e-05,6e-12,unit='GeV**-2',
                  name='Fermi coupling constant',
                  symbol='G(F)/(hbar*c)**3',
                  html_symbol='<i>G</i><sub>F</sub>/(&hbar;c)<sup>3</sup>',
                  latex_symbol='G_{' + norm('F') + '}/(\\hbar c)^{3}',
                  add_symbol=True,
                  _key='gf'
                  )