# -*- coding: utf-8 -*-

# module test_misc

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
Simple unit tests for the gummy module
"""

import metrolopy as uc
import numpy as np

def test_ku():
    from metrolopy.gummy import _ku
    from decimal import Decimal
    
    assert _ku(2.1,3.3) == 2.1*3.3
    assert _ku(2.1,Decimal(3.3)) == Decimal(float(2.1))*Decimal(3.3)
    
def test_lg10():
    from metrolopy.gummy import _lg10
    from mpmath import mpf,log10
    from decimal import Decimal
    from fractions import Fraction
    import math
    
    assert _lg10(mpf(2.2)) == log10(mpf(2.2))
    assert _lg10(Decimal(2.2)) == Decimal(2.2).log10()
    assert _lg10(2.2) == np.log10(2.2)
    assert _lg10(Fraction(2,3)) == math.log10(float(Fraction(2,3)))
    
def test_meta_properties():
    a = uc.gummy.cimethod
    assert isinstance(a,str)
    uc.gummy.cimethod = 'symmetric'
    assert uc.gummy.cimethod == 'symmetric'
    g = uc.gummy(1,1)
    assert g.cimethod == 'symmetric'
    
    try:
        uc.gummy.cimethod = 'asdfkjhasd'
        assert False
    except ValueError:
        pass
    
    uc.gummy.cimethod = a

    b = uc.gummy.bayesian
    assert isinstance(b,bool)
    uc.gummy.bayesian = not b
    assert uc.gummy.bayesian is not b
    g = uc.gummy(1,1)
    assert g.bayesian is not b
    
    uc.gummy.bayesian = b
    
    b = uc.gummy.style
    assert isinstance(b,str)
    uc.gummy.style = 'xf'
    assert uc.gummy.style == 'xf'
    g = uc.gummy(1,1)
    assert g.style == 'xf'
    
    try:
        uc.gummy.style = 'asdfjkhasdf'
        assert False
    except ValueError:
        pass
    
    uc.gummy.style = b
    
    k = uc.gummy.cmp_k
    p = uc.gummy.cmp_p
    if p is None:
        assert k >= 0
    else:
        assert k is None
        assert p >= 0
        assert p <= 1
        
    
    uc.gummy.cmp_k = 2.1
    assert uc.gummy.cmp_p is None
    assert uc.gummy.cmp_k == 2.1
    
    if k is None:
        uc.gummy.cmp_p = p
    else:
        uc.gummy.cmp_k = k
        
    pm = uc.gummy.p_method
    assert isinstance(pm,str)
    uc.gummy.p_method = 'ccp'
    assert uc.gummy.p_method == 'ccp'
    uc.gummy.p_method = None
    assert uc.gummy.p_method == 'loc'
    uc.gummy.p_method = pm
    
    d = uc.gummy.max_dof
    uc.gummy.max_dof = 1021
    assert uc.gummy.max_dof == 1021
    uc.gummy.max_dof = d
    
    n = uc.gummy.nsig
    uc.gummy.nsig = n + 1
    assert uc.gummy.nsig == n + 1
    g = uc.gummy(1,1)
    assert g.nsig == n + 1
    uc.gummy.nsig = n
    
    ts = uc.gummy.thousand_spaces
    assert isinstance(ts,bool)
    uc.gummy.thousand_spaces = not ts
    assert uc.gummy.thousand_spaces is not ts
    g = uc.gummy(1,1)
    assert g.thousand_spaces is not ts
    uc.gummy.thousand_spaces = ts
    
def test_init():
    g = uc.gummy(1.2,3.4)
    assert g.x == 1.2
    assert g.u == 3.4
    assert g.U == g.u
    assert g.unit is uc.one
    
    gg = uc.gummy(g)
    assert gg is not g
    assert gg.x == 1.2
    assert gg.u == 3.4
    assert g.correlation(gg) == 1
    
    q = uc.Quantity(5.6,unit='m')
    gq = uc.gummy(q)
    assert gq.x == 5.6
    assert gq.u == 0
    assert gq.unit is uc.unit('m')
    
    u = uc.ummy(7.8,9.1)
    gu = uc.gummy(u)
    assert gu.x == 7.8
    assert gu.u == 9.1
    assert gu.unit is uc.one
    
    ung = uc.gummy(2.1,unit='mm')
    gun = uc.gummy(3.4,u=ung,unit='m')
    assert gun.unit is uc.unit('m')
    assert gun.x == 3.4
    assert gun.U == 2.1
    assert gun.u == 2.1/1000
    assert gun.uunit is uc.unit('mm')
    
    ung = uc.Quantity(2.1,unit='mm')
    gun = uc.gummy(3.4,u=ung,unit='m')
    assert gun.unit is uc.unit('m')
    assert gun.x == 3.4
    assert gun.U == 2.1
    assert gun.u == 2.1/1000
    assert gun.uunit is uc.unit('mm')
    
    guu = uc.gummy(1,1,unit='m',uunit='m')
    assert guu.unit is uc.unit('m')
    assert guu.uunit is None
    
    guuz = uc.gummy(1,unit='m',uunit='mm')
    assert guuz.unit is uc.unit('m')
    assert guuz.uunit is None
    
    gp = uc.gummy(2,1,p=0.95)
    assert gp.x == 2
    assert gp.p == 0.95
    assert abs(gp.k - 1.96) < 0.01
    assert abs(gp.u - 1/1.96) < 0.01
    
    gk = uc.gummy(3,1,k=2)
    assert gk.x == 3
    assert gk.k == 2
    assert abs(gk.p - 0.954) < 0.01
    
    gk = uc.gummy(3,1,k=1.96)
    assert gk.x == 3
    assert gk.k == 1.96
    assert abs(gk.p - 0.95) < 0.01
    
    gd = uc.gummy(uc.TriangularDist(2.2,0.73))
    assert gd.x == 2.2
    assert abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12
    assert abs(gd.U - np.sqrt(0.73**2/6)) < 1e-12
    
    gd = uc.gummy(uc.TriangularDist(2.2,0.73),k=2.6)
    assert gd.x == 2.2
    assert abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12
    assert abs(gd.U - 2.6*np.sqrt(0.73**2/6)) < 1e-12
    
    gd = uc.gummy(uc.TriangularDist(2.2,0.73),k=2.6,unit='m',uunit='cm')
    assert gd.x == 2.2
    assert abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12
    assert abs(gd.U - 2.6*100*np.sqrt(0.73**2/6)) < 1e-12
    
    try:
        uc.gummy(3.3,1.3,unit='%',uunit='m')
        assert False
    except uc.NoUnitConversionFoundError:
        pass
    
    gg = uc.gummy(3.3,1.4,unit='%',uunit=uc.one)
    assert abs(gg.u - 100*1.4) < 1e-12
    
    gg = uc.gummy(3.3,1.4,unit='%',uunit='ppm')
    assert abs(gg.u - 3.3*1.4e-6) < 1e-12
    
    gg = uc.gummy(3.3,1.4,unit='m',uunit='%')
    assert abs(gg.u - 3.3*0.014) < 1e-12
    
    gg = uc.gummy(1,1,dof=3.3,utype='abc',name='def')
    assert gg.dof == 3.3
    assert gg.utype == 'abc'
    assert gg.name == 'def'
    
def test_U():
    a = uc.gummy(1.1,2.2)
    assert a.U == 2.2
    b = uc.gummy(3,4.4)
    c = a + b
    c.ubreakdown = [a,b]
    assert c.U[0] == 2.2
    assert c.U[1] == 4.4
    
    a = uc.gummy(1.1,2.2,'m')
    assert a.U == 2.2
    b = uc.gummy(3,4.4,'s')
    c = a/b
    c.ubreakdown = [a,b]
    assert abs(c.U[0] - 2.2/3) < 1e-12
    assert abs(c.U[1] - 4.4*1.1/3**2) < 1e-12
    
def test_set_U():
    g = uc.gummy(1,1,k=2,unit='cm',uunit='mm')
    assert g.U == 1
    assert g.uunit is uc.unit('mm')
    g.unit = 'm'
    assert g.U == 0.001
    
    a = uc.gummy(1.1,2.2)
    b = uc.gummy(3,4.4)
    c = a + b
    c.ubreakdown = [a,b]
    assert c.U[0] == 2.2
    assert c.U[1] == 4.4
    
def test_iset_U():
    g = uc.gummy(1)
    assert g._iset_U() == 0
    assert g._iset_U(unit='m').unit is uc.unit('m')
    assert g._iset_U(unit='m').value == 0
    assert g._iset_U(u=float('inf')) == float('inf')
    assert g._iset_U(u=float('inf'),unit='m').value == float('inf')
    assert g._iset_U(u=float('inf'),unit='m').unit is uc.unit('m')
    
    g = uc.gummy(1,1,k=2,unit='m')
    assert g._iset_U(unit='m').value == 1
    assert g._iset_U(k=4) == 2
    assert g._iset_U(unit='%').value == 100
    assert g._iset_U(unit='%').unit is uc.unit('%')
    assert g._iset_U(unit='cm').value == 100
    assert g._iset_U(unit='cm').unit is uc.unit('cm')
    
    try:
        assert g._iset_U(unit='kg')
        assert False
    except uc.NoUnitConversionFoundError:
        pass
    
    g = uc.gummy(0,1,unit='m')
    g._iset_U(unit='%').value == float('inf')
    g._iset_U(unit='%').unit is uc.unit('%')
    
    g = uc.gummy(1,1,unit='degC')
    assert g._iset_U(k=2) == 2
    
def test_Usim():
    g = uc.gummy(-8.5,1,unit='m')
    g.sim()
    assert abs(g.Usim[0] - 1) < 0.1
    assert abs(g.Usim[1] - 1) < 0.1
    g.uunit = '%'
    assert abs(g.Usim[0] - 100) < 10
    assert abs(g.Usim[1] - 100) < 10
    g.uunit = 'mm'
    assert abs(g.Usim[0] - 1000) < 100
    assert abs(g.Usim[1] - 1000) < 100
    
def test_xusim():
    """
    tests gummy.xsim and gummy.usim
    """
    g = uc.gummy(-3.3,0.11)
    g.sim()
    assert abs(g.xsim + 3.3) < 0.1
    assert abs(g.usim - 0.11) < 0.01
    
def test_cisim():
    """
    tests gummy.cisim and gummy.cimethod
    """
    g = uc.gummy(-3.3,0.11)
    g.sim()
    
    g.cimethod = 'symmetric'
    ci = g.cisim
    assert abs(ci[0] + 3.41) < 0.01
    assert abs(ci[1] + 3.19) < 0.01
    
    g.cimethod = 'shortest'
    ci = g.cisim
    assert abs(ci[0] + 3.41) < 0.01
    assert abs(ci[1] + 3.19) < 0.01
    
    try:
        g.cimethod = 'asdfhasd'
        assert False
    except ValueError:
        pass
    
    a = uc.gummy(uc.UniformDist(center=0,half_width=1))
    b = a**2
    b.p = 0.8
    b.cimethod = 'shortest'
    b.sim()
    ci = b.cisim
    assert abs(ci[0]) < 0.001
    assert abs(ci[1] - 0.64) < 0.01
    b.cimethod = 'symmetric'
    ci = b.cisim
    assert abs(ci[0]- 0.01) < 0.001
    assert abs(ci[1] - 0.81) < 0.01
    
def test_simdata():
    """
    tests gummy.simdata and gummy.simsorted
    """
    g = uc.gummy(1,1)
    g.sim(n=10)
    assert len(g.simdata) == 10
    assert len(g.simsorted) == 10
    for i in range(9):
        assert g.simsorted[i+1] >= g.simsorted[i]
        
def test_distribution():
    g = uc.gummy(1,1)
    assert isinstance(g.distribution,uc.Distribution)
   
def test_ksim():
    g = uc.gummy(1,1)
    g.p = 0.9545
    assert abs(g.ksim - 2) < 0.1
    
def test_independant():
    a = uc.gummy(1,1)
    b = uc.gummy(-2,2)
    c = a + b
    assert a.independant
    assert not c.independant
    d = uc.gummy(1)
    assert not d.independant
    
def test_name():
    g = uc.gummy(1,1)
    assert g.name is None
    g.name = 'abc'
    assert g.name == 'abc'
    g.name = ['def','g','hi','jkl']
    assert g.name == 'def'
    assert g.get_name() == 'def'
    assert g.get_name(fmt='html') == 'g'
    assert g.get_name(fmt='latex') == 'hi'
    assert g.get_name(fmt='ascii') == 'jkl'
    g.name = 'f'
    assert g.get_name(fmt='html') == '<i>f</i>'
    g.name = 'ff'
    assert g.get_name(fmt='latex') == '\\text{ff}'
    assert g.get_name(fmt='latex',norm=lambda x:'*'+x+'*') == '*ff*'
    
    try:
        g.get_name(fmt='asdgfads')
        assert False
    except ValueError:
        pass
    
    g.name = None
    assert g.name is None
    
    try:
        g.name = 3
        assert False
    except ValueError:
        pass
    
    try:
        g.name = ['a','b']
        assert False
    except ValueError:
        pass
    
    try:
        g.name = [1,2]
        assert False
    except ValueError:
        pass
    
def test_unit():
    g = uc.gummy(1,1,unit='m')
    assert g.unit is uc.unit('m')
    assert str(g.unit) == 'm'
    g.unit = 'cm'
    assert str(g.unit) = 'cm'
    assert g.x = 100
    assert g.u = 100
    g.unit = uc.unit('mm')
    assert str(g.unit) = 'mm'
    assert g.x = 1000
    assert g.u = 1000
    
def test_uunit():
    g = uc.gummy(1,1)
    assert g.uunit is None
    g = uc.gummy(1,1,unit='m',uunit='cm')
    assert str(g.uunit) == 'cm'
    
    g = uc.gummy(1,unit='m',uunit='cm')
    assert uunit is None
    g.uunit = 'mm'
    assert g.uunit is None
    
    g = uc.gummy(1,1,unit='m')
    assert g.uunit is None
    g.uunit = 'cm'
    assert str(g.uunit) == 'mm'
    assert g.U = 1000
    g.uunit = None
    assert g.uunit is None
    