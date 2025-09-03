# -*- coding: utf-8 -*-

# module test_misc

# Copyright (C) 2025 National Research Council Canada
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
import unittest

class TestGummy(unittest.TestCase):
    """
    Various unit tests of gummy properties and methods.
    """
    def test_ku(self):
        from metrolopy.gummy import _ku
        from decimal import Decimal
        
        self.assertTrue(_ku(2.1,3.3) == 2.1*3.3)
        self.assertTrue(_ku(2.1,Decimal(3.3)) == Decimal(float(2.1))*Decimal(3.3))
        
    def test_meta_properties(self):
        """
        Test class properties.

        """
        a = uc.gummy.cimethod
        self.assertTrue(isinstance(a,str))
        uc.gummy.cimethod = 'symmetric'
        self.assertTrue(uc.gummy.cimethod == 'symmetric')
        g = uc.gummy(1,1)
        self.assertTrue(g.cimethod == 'symmetric')
        
        try:
            uc.gummy.cimethod = 'asdfkjhasd'
            self.assertTrue(False)
        except ValueError:
            pass
        
        uc.gummy.cimethod = a
    
        b = uc.gummy.bayesian
        self.assertTrue(isinstance(b,bool))
        uc.gummy.bayesian = not b
        self.assertTrue(uc.gummy.bayesian is not b)
        g = uc.gummy(1,1)
        self.assertTrue(g.bayesian is not b)
        
        uc.gummy.bayesian = b
        
        b = uc.gummy.style
        self.assertTrue(isinstance(b,str))
        uc.gummy.style = 'xf'
        self.assertTrue(uc.gummy.style == 'xf')
        g = uc.gummy(1,1)
        self.assertTrue(g.style == 'xf')
        
        try:
            uc.gummy.style = 'asdfjkhasdf'
            self.assertTrue(False)
        except ValueError:
            pass
        
        uc.gummy.style = b
        
        k = uc.gummy.cmp_k
        p = uc.gummy.cmp_p
        if p is None:
            self.assertTrue(k >= 0)
        else:
            self.assertTrue(k is None)
            self.assertTrue(p >= 0)
            self.assertTrue(p <= 1)
            
        
        uc.gummy.cmp_k = 2.1
        self.assertTrue(uc.gummy.cmp_p is None)
        self.assertTrue(uc.gummy.cmp_k == 2.1)
        
        if k is None:
            uc.gummy.cmp_p = p
        else:
            uc.gummy.cmp_k = k
            
        try:
            uc.gummy.p_method = 'ccp'
            self.assertTrue(uc.gummy.p_method == 'ccp')
            uc.gummy.p_method = None
            self.assertTrue(uc.gummy.p_method == 'loc')
        finally:
            uc.gummy.p_method = 'loc'
        
        d = uc.gummy.max_dof
        uc.gummy.max_dof = 1021
        self.assertTrue(uc.gummy.max_dof == 1021)
        uc.gummy.max_dof = d
        
        n = uc.gummy.nsig
        uc.gummy.nsig = n + 1
        self.assertTrue(uc.gummy.nsig == n + 1)
        g = uc.gummy(1,1)
        self.assertTrue(g.nsig == n + 1)
        uc.gummy.nsig = n
        
        ts = uc.gummy.thousand_spaces
        self.assertTrue(isinstance(ts,bool))
        uc.gummy.thousand_spaces = not ts
        self.assertTrue(uc.gummy.thousand_spaces is not ts)
        g = uc.gummy(1,1)
        self.assertTrue(g.thousand_spaces is not ts)
        uc.gummy.thousand_spaces = ts
        
    def test_init(self):
        """
        Test gummy creation
        """
        g = uc.gummy(1.2,3.4)
        self.assertTrue(g.x == 1.2)
        self.assertTrue(g.u == 3.4)
        self.assertTrue(g.U == g.u)
        self.assertTrue(g.unit is uc.one)
        
        gg = uc.gummy(g)
        self.assertTrue(gg is not g)
        self.assertTrue(gg.x == 1.2)
        self.assertTrue(gg.u == 3.4)
        self.assertTrue(g.correlation(gg) == 1)
        
        q = uc.Quantity(5.6,unit='m')
        gq = uc.gummy(q)
        self.assertTrue(gq.x == 5.6)
        self.assertTrue(gq.u == 0)
        self.assertTrue(gq.unit is uc.unit('m'))
        
        u = uc.ummy(7.8,9.1)
        gu = uc.gummy(u)
        self.assertTrue(gu.x == 7.8)
        self.assertTrue(gu.u == 9.1)
        self.assertTrue(gu.unit is uc.one)
        
        ung = uc.gummy(2.1,unit='mm')
        gun = uc.gummy(3.4,u=ung,unit='m')
        self.assertTrue(gun.unit is uc.unit('m'))
        self.assertTrue(gun.x == 3.4)
        self.assertTrue(gun.U == 2.1)
        self.assertTrue(gun.u == 2.1/1000)
        self.assertTrue(gun.uunit is uc.unit('mm'))
        
        ung = uc.Quantity(2.1,unit='mm')
        gun = uc.gummy(3.4,u=ung,unit='m')
        self.assertTrue(gun.unit is uc.unit('m'))
        self.assertTrue(gun.x == 3.4)
        self.assertTrue(gun.U == 2.1)
        self.assertTrue(gun.u == 2.1/1000)
        self.assertTrue(gun.uunit is uc.unit('mm'))
        
        guu = uc.gummy(1,1,unit='m',uunit='m')
        self.assertTrue(guu.unit is uc.unit('m'))
        self.assertTrue(guu.uunit is None)
        
        guuz = uc.gummy(1,unit='m',uunit='mm')
        self.assertTrue(guuz.unit is uc.unit('m'))
        self.assertTrue(guuz.uunit is None)
        
        gp = uc.gummy(2,1,p=0.95)
        self.assertTrue(gp.x == 2)
        self.assertTrue(gp.p == 0.95)
        self.assertTrue(abs(gp.k - 1.96) < 0.01)
        self.assertTrue(abs(gp.u - 1/1.96) < 0.01)
    
        gk = uc.gummy(3,1,k=2)
        self.assertTrue(gk.x == 3)
        self.assertTrue(gk.k == 2)
        self.assertTrue(abs(gk.p - 0.954) < 0.1)
        
        gk = uc.gummy(3,1,k=1.96)
        self.assertTrue(gk.x == 3)
        self.assertTrue(gk.k == 1.96)
        self.assertTrue(abs(gk.p - 0.95) < 0.1)
        
        gd = uc.gummy(uc.TriangularDist(2.2,0.73))
        self.assertTrue(gd.x == 2.2)
        self.assertTrue(abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12)
        self.assertTrue(abs(gd.U - np.sqrt(0.73**2/6)) < 1e-12)
        
        gd = uc.gummy(uc.TriangularDist(2.2,0.73),k=2.6)
        self.assertTrue(gd.x == 2.2)
        self.assertTrue(abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12)
        self.assertTrue(abs(gd.U - 2.6*np.sqrt(0.73**2/6)) < 1e-12)
        
        gd = uc.gummy(uc.TriangularDist(2.2,0.73),k=2.6,unit='m',uunit='cm')
        self.assertTrue(gd.x == 2.2)
        self.assertTrue(abs(gd.u - np.sqrt(0.73**2/6)) < 1e-12)
        self.assertTrue(abs(gd.U - 2.6*100*np.sqrt(0.73**2/6)) < 1e-12)
        
        try:
            uc.gummy(3.3,1.3,unit='%',uunit='m')
            self.assertTrue(False)
        except uc.IncompatibleUnitsError:
            pass
        
        gg = uc.gummy(3.3,1.4,unit='%',uunit=uc.one)
        self.assertTrue(abs(gg.u - 100*1.4) < 1e-12)
        
        gg = uc.gummy(3.3,1.4,unit='%',uunit='ppm')
        self.assertTrue(abs(gg.u - 3.3*1.4e-6) < 1e-12)
        
        gg = uc.gummy(3.3,1.4,unit='m',uunit='%')
        self.assertTrue(abs(gg.u - 3.3*0.014) < 1e-12)
        
        gg = uc.gummy(1,1,dof=3.3,utype='abc',name='def')
        self.assertTrue(gg.dof == 3.3)
        self.assertTrue(gg.utype == 'abc')
        self.assertTrue(gg.name == 'def')
        
    def test_U(self):
        """
        Test expanded uncertainty.
        """
        a = uc.gummy(1.1,2.2)
        self.assertTrue(a.U == 2.2)
        b = uc.gummy(3,4.4)
        c = a + b
        c.ubreakdown = [a,b]
        self.assertTrue(c.U[0] == 2.2)
        self.assertTrue(c.U[1] == 4.4)
        
        a = uc.gummy(1.1,2.2,'m')
        self.assertTrue(a.U == 2.2)
        b = uc.gummy(3,4.4,'s')
        c = a/b
        c.ubreakdown = [a,b]
        self.assertTrue(abs(c.U[0] - 2.2/3) < 1e-12)
        self.assertTrue(abs(c.U[1] - 4.4*1.1/3**2) < 1e-12)
        
    def test_set_U(self):
        """
        Additional tests for expanded uncertainty.
        """
        g = uc.gummy(1.0,1.0,k=2,unit='cm',uunit='mm')
        self.assertTrue(g.U == 1)
        self.assertTrue(g.uunit is uc.unit('mm'))
        g.unit = 'm'
        self.assertTrue(g.U == 0.001)
        
        a = uc.gummy(1.1,2.2)
        b = uc.gummy(3,4.4)
        c = a + b
        c.ubreakdown = [a,b]
        self.assertTrue(c.U[0] == 2.2)
        self.assertTrue(c.U[1] == 4.4)
        
    def test_iset_U(self):
        """
        Additional tests for expanded uncertainty.
        """
        g = uc.gummy(1)
        self.assertTrue(g._iset_U() == 0)
        self.assertTrue(g._iset_U(unit='m').unit is uc.unit('m'))
        self.assertTrue(g._iset_U(unit='m').value == 0)
        self.assertTrue(g._iset_U(u=float('inf')) == float('inf'))
        self.assertTrue(g._iset_U(u=float('inf'),unit='m').value == float('inf'))
        self.assertTrue(g._iset_U(u=float('inf'),unit='m').unit is uc.unit('m'))
        
        g = uc.gummy(1,1,k=2,unit='m')
        self.assertTrue(g._iset_U(unit='m').value == 1)
        self.assertTrue(g._iset_U(k=4) == 2)
        self.assertTrue(g._iset_U(unit='%').value == 100)
        self.assertTrue(g._iset_U(unit='%').unit is uc.unit('%'))
        self.assertTrue(g._iset_U(unit='cm').value == 100)
        self.assertTrue(g._iset_U(unit='cm').unit is uc.unit('cm'))
        
        try:
            self.assertTrue(g._iset_U(unit='kg'))
            self.assertTrue(False)
        except uc.IncompatibleUnitsError:
            pass
        
        g = uc.gummy(0,1,unit='m')
        g._iset_U(unit='%').value == float('inf')
        g._iset_U(unit='%').unit is uc.unit('%')
        
        g = uc.gummy(1,1,unit='degC')
        self.assertTrue(g._iset_U(k=2) == 2)
        
    def test_Usim(self):
        """
        Test the gummy.Usim property
        """
        g = uc.gummy(-8.5,1,unit='m')
        g.sim()
        self.assertTrue(abs(g.Usim[0] - 1) < 0.5)
        self.assertTrue(abs(g.Usim[1] - 1) < 0.5)
        g.uunit = '%'
        self.assertTrue(abs(g.Usim[0] - 100*1/8.5) < 10)
        self.assertTrue(abs(g.Usim[1] - 100*1/8.5) < 10)
        g.uunit = 'mm'
        self.assertTrue(abs(g.Usim[0] - 1000) < 500)
        self.assertTrue(abs(g.Usim[1] - 1000) < 500)
        
    def test_xusim(self):
        """
        tests gummy.xsim and gummy.usim
        """
        g = uc.gummy(-3.3,0.11)
        g.sim()
        self.assertTrue(abs(g.xsim + 3.3) < 0.1)
        self.assertTrue(abs(g.usim - 0.11) < 0.1)
        
    def test_cisim(self):
        """
        tests gummy.cisim and gummy.cimethod
        """
        g = uc.gummy(-3.3,0.11)
        g.sim()
        
        g.cimethod = 'symmetric'
        ci = g.cisim
        self.assertTrue(abs(ci[0] + 3.41) < 0.1)
        self.assertTrue(abs(ci[1] + 3.19) < 0.1)
        
        g.cimethod = 'shortest'
        ci = g.cisim
        self.assertTrue(abs(ci[0] + 3.41) < 0.1)
        self.assertTrue(abs(ci[1] + 3.19) < 0.1)
        
        try:
            g.cimethod = 'asdfhasd'
            self.assertTrue(False)
        except ValueError:
            pass
        
        a = uc.gummy(uc.UniformDist(center=0,half_width=1))
        b = a**2
        b.p = 0.8
        b.cimethod = 'shortest'
        b.sim()
        ci = b.cisim
        self.assertTrue(abs(ci[0]) < 0.001)
        self.assertTrue(abs(ci[1] - 0.64) < 0.01)
        b.cimethod = 'symmetric'
        ci = b.cisim
        self.assertTrue(abs(ci[0]- 0.01) < 0.001)
        self.assertTrue(abs(ci[1] - 0.81) < 0.01)
        
    def test_simdata(self):
        """
        tests gummy.simdata and gummy.simsorted
        """
        g = uc.gummy(1,1)
        g.sim(n=10)
        self.assertTrue(len(g.simdata) == 10)
        self.assertTrue(len(g.simsorted) == 10)
        for i in range(9):
            self.assertTrue(g.simsorted[i+1] >= g.simsorted[i])
            
    def test_distribution(self):
        """
        Test gummy.distribution
        """
        g = uc.gummy(1,1)
        self.assertTrue(isinstance(g.distribution,uc.Distribution))
       
    def test_ksim(self):
        """
        Test gummy.ksim
        """
        g = uc.gummy(1,1)
        g.p = 0.9545
        g.sim()
        self.assertTrue(abs(g.ksim - 2) < 0.1)
        
    def test_independent(self):
        """
        Test gummy.independant
        """
        a = uc.gummy(1,1)
        b = uc.gummy(-2,2)
        c = a + b
        self.assertTrue(a.independent)
        self.assertTrue(not c.independent)
        d = uc.gummy(1)
        self.assertTrue(not d.independent)
        
    def test_name(self):
        """
        Test setting and getting gummy.name
        """
        g = uc.gummy(1,1)
        self.assertTrue(g.name is None)
        g.name = 'abc'
        self.assertTrue(g.name == 'abc')
        g.name = ['def','g','hi','jkl']
        self.assertTrue(g.name == 'def')
        self.assertTrue(g.get_name() == 'def')
        self.assertTrue(g.get_name(fmt='html') == 'g')
        self.assertTrue(g.get_name(fmt='latex') == 'hi')
        self.assertTrue(g.get_name(fmt='ascii') == 'jkl')
        g.name = 'f'
        self.assertTrue(g.get_name(fmt='html') == '<i>f</i>')
        g.name = 'ff'
        self.assertTrue(g.get_name(fmt='latex') == '\\text{ff}')
        self.assertTrue(g.get_name(fmt='latex',norm=lambda x:'*'+x+'*') == '*ff*')
        
        try:
            g.get_name(fmt='asdgfads')
            self.assertTrue(False)
        except ValueError:
            pass
        
        g.name = None
        self.assertTrue(g.name is None)
        
        try:
            g.name = 3
            self.assertTrue(False)
        except ValueError:
            pass
        
        try:
            g.name = ['a','b']
            self.assertTrue(False)
        except ValueError:
            pass
        
        try:
            g.name = [1,2]
            self.assertTrue(False)
        except ValueError:
            pass
        
    def test_unit(self):
        """
        Test setting and getting gummy.unit
        """
        g = uc.gummy(1,1,unit='m')
        self.assertTrue(g.unit is uc.unit('m'))
        self.assertTrue(str(g.unit) == 'm')
        g.unit = 'cm'
        self.assertTrue(str(g.unit) == 'cm')
        self.assertTrue(g.x == 100)
        self.assertTrue(g.u == 100)
        g.unit = uc.unit('mm')
        self.assertTrue(str(g.unit) == 'mm')
        self.assertTrue(g.x == 1000)
        self.assertTrue(g.u == 1000)
        
    def test_uunit(self):
        """
        Test setting and getting gummy.uunit
        """
        g = uc.gummy(1,1)
        self.assertTrue(g.uunit is None)
        g = uc.gummy(1,1,unit='m',uunit='cm')
        self.assertTrue(str(g.uunit) == 'cm')
        
        g = uc.gummy(1,unit='m',uunit='cm')
        g.uunit = 'mm'
        self.assertTrue(g.uunit is None)
        
        g = uc.gummy(1,1,unit='m')
        self.assertTrue(g.uunit is None)
        g.uunit = 'cm'
        self.assertTrue(str(g.uunit) == 'cm')
        self.assertTrue(g.U == 100)
        g.uunit = None
        self.assertTrue(g.uunit is None)
        
    def test_uunit_is_rel(self):
        """
        Test gummy.uunit_is_rel
        """
        g = uc.gummy(1,1)
        self.assertTrue(not g.uunit_is_rel)
        g.uunit='%'
        self.assertTrue(g.uunit_is_rel)
        
        g = uc.gummy(1,1,unit='m')
        self.assertTrue(not g.uunit_is_rel)
        g.uunit = 'mm'
        self.assertTrue(not g.uunit_is_rel)
        g.uunit = '%'
        self.assertTrue(g.uunit_is_rel)
        
    def test_k(self):
        """
        Test setting and getting gummy.k
        """
        g = uc.gummy(1)
        self.assertTrue(g.k is None)
        
        g = uc.gummy(1,1,k=2)
        self.assertTrue(g.k == 2)
        g.k = 3
        self.assertTrue(g.k == 3)
        try:
            g.k = -1
            self.assertTrue(False)
        except ValueError:
            pass
        
        uc.gummy.p_method = None
        g.p = 0.68
        self.assertTrue(abs(g.k - 0.99) < 0.01)
        g.p = 0.95
        self.assertTrue(abs(g.k - 1.96) < 0.01)
        
    def test_p(self):
        """
        Test setting and getting gummy.p
        """
        uc.gummy.p_method = None
        g = uc.gummy(1,1)
        g.p = 0.95
        self.assertTrue(g.p == 0.95)
        g.k = 2
        self.assertTrue(abs(g.p - 0.9545) < 0.0001)
        
        try:
            g.p = -0.1
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            g.p = 1.1
            self.assertTrue(False)
        except ValueError:
            pass
        
    def test_correlation(self):
        """
        Test gummy.correlation
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        self.assertTrue(a.correlation(b) == 0)
        self.assertTrue(a.correlation(1) == 0)
        self.assertTrue(abs(c.correlation(a) - 1/c.u) < 1e-4)
        self.assertTrue(abs(a.correlation(c) - 1/c.u) < 1e-4)
        self.assertTrue(abs(c.correlation(b) - 2/c.u) < 1e-4)
        self.assertTrue(abs(b.correlation(c) - 2/c.u) < 1e-4)
        self.assertTrue(a.correlation(a) == 1)
        
    def test_covariance(self):
        """
        Test gummy.covariance
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        self.assertTrue(a.covariance(b) == 0)
        self.assertTrue(a.covariance(1) == 0)
        self.assertTrue(abs(c.covariance(a) - 1) < 1e-4)
        self.assertTrue(abs(a.covariance(c) - 1) < 1e-4)
        self.assertTrue(abs(c.covariance(b) - 4) < 1e-4)
        self.assertTrue(abs(b.covariance(c) - 4) < 1e-4)
        self.assertTrue(a.covariance(a) == 1)
        self.assertTrue(b.covariance(b) == 4)
        
    def test_correlation_sim(self):
        """
        Test gummy.correlation_sim
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        uc.gummy.simulate([a,b,c])
        self.assertTrue(abs(a.correlation_sim(b)) < 1e-1)
        self.assertTrue(abs(c.correlation_sim(a) - 1/c.u) < 1e-1)
        self.assertTrue(abs(a.correlation_sim(c) - 1/c.u) < 1e-1)
        self.assertTrue(abs(c.correlation_sim(b) - 2/c.u) < 1e-1)
        self.assertTrue(abs(b.correlation_sim(c) - 2/c.u) < 1e-1)
        self.assertTrue(abs(a.correlation_sim(a) - 1) < 1e-6)
        
    def test_covariance_sim(self):
        """
        Test gummy.covariance_sim
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        uc.gummy.simulate([a,b,c])
        self.assertTrue(abs(a.covariance_sim(b)) < 1e-1)
        self.assertTrue(abs(c.covariance_sim(a) - 1) < 1e-1)
        self.assertTrue(abs(a.covariance_sim(c) - 1) < 1e-1)
        self.assertTrue(abs(c.covariance_sim(b) - 4) < 1e-1)
        self.assertTrue(abs(b.covariance_sim(c) - 4) < 1e-1)
        self.assertTrue(abs(a.covariance_sim(a) - 1) < 1e-1)
        self.assertTrue(abs(b.covariance_sim(b) - 4) < 1e-1)
        
    def test_correlation_matrix(self):
        """
        Test gummy.correlation_matrix
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        uc.gummy.simulate([a,b,c])
        m = uc.gummy.correlation_matrix_sim([a,b,c])
        self.assertTrue(abs(m[0][0] - 1) < 1e-1)
        self.assertTrue(abs(m[1][1] - 1) < 1e-1)
        self.assertTrue(abs(m[2][2] - 1) < 1e-1)
        self.assertTrue(abs(m[0][1]) < 1e-1)
        self.assertTrue(abs(m[1][0]) < 1e-1)
        self.assertTrue((m[2][0] - 1/c.u) < 1e-1)
        self.assertTrue((m[0][2] - 1/c.u) < 1e-1)
        self.assertTrue((m[2][1] - 2/c.u) < 1e-1)
        self.assertTrue((m[1][2] - 2/c.u) < 1e-1)
        
    def test_covariance_matrix(self):
        """
        Test gummy.covariance_matrix
        """
        a = uc.gummy(1,1)
        b = uc.gummy(1,2)
        c = a + b
        uc.gummy.simulate([a,b,c])
        m = uc.gummy.covariance_matrix_sim([a,b,c])
        self.assertTrue(abs(m[0][0] - 1) < 1e-1)
        self.assertTrue(abs(m[1][1] - 4) < 1e-1)
        self.assertTrue(abs(m[2][2] - 5) < 1e-1)
        self.assertTrue(abs(m[0][1]) < 1e-1)
        self.assertTrue(abs(m[1][0]) < 1e-1)
        self.assertTrue((m[2][0] - 1) < 1e-1)
        self.assertTrue((m[0][2] - 1) < 1e-1)
        self.assertTrue((m[2][1] - 4) < 1e-1)
        self.assertTrue((m[1][2] - 4) < 1e-1)
        
    def test_finfo(self):
        """
        Test gummy.finfo
        """
        a = uc.gummy(1,1)
        self.assertTrue(a.finfo.rel_u == 0)
        
    def test_real(self):
        """
        Test gummy.real
        """
        a = uc.gummy(1.1,0.3)
        b = a.real
        self.assertTrue(not a is b)
        self.assertTrue(a.x == b.x)
        self.assertTrue(a.u == b.u)
        self.assertTrue(a.correlation(b) == 1)
        
    def test_conjugate(self):
        """
        Test gummy.conjugate
        """
        a = uc.gummy(1.1,0.3)
        b = a.conjugate()
        self.assertTrue(not a is b)
        self.assertTrue(a.x == b.x)
        self.assertTrue(a.u == b.u)
        self.assertTrue(a.correlation(b) == 1)
        
    def test_angle(self):
        """
        Test gummy.angle
        """
        self.assertTrue(uc.gummy(1.1,0.3).angle() == 0)
        self.assertTrue(abs(abs(uc.gummy(-1.1,0.3).angle()) - np.pi) < 1e-10)
        
    def test_utype(self):
        """
        Test gummy.utype
        """
        a = uc.gummy(1.2,1,1,utype='A')
        self.assertTrue(a.utype == 'A')
        b = uc.gummy(1.2,1,1,utype='xyz')
        self.assertTrue(b.utype == 'xyz')
        
    def test_ufrom(self):
        """
        Test gummy.ufrom
        """
        a = uc.gummy(1.2,0.2,utype='A')
        b = uc.gummy(3.2,0.5,utype='A')
        c = uc.gummy(0.9,0.2,utype='B')
        d = a + b + c
        self.assertTrue(abs(d.ufrom('A') - np.sqrt(0.2**2+0.5**2)) < 1e-8)
        uc.gummy.simulate([d,b])
        self.assertTrue(abs(d.ufromsim('A') - np.sqrt(0.2**2+0.5**2)) < 1e-2)
        self.assertTrue(d.ufrom(c) == 0.2)
        
    def test_dof_from(self):
        """
        Test gummy.doffrom
        """
        a = uc.gummy(1.2,0.2,dof=5,utype='A')
        b = uc.gummy(3.2,0.5,dof=7,utype='A')
        c = uc.gummy(0.9,0.2,utype='B')
        d = a + b + c
        self.assertTrue(abs(d.doffrom('A') - (0.2**2+0.5**2)**2/(0.2**4/5 + 0.5**4/7)) < 1e-10)
        self.assertTrue(d.doffrom(a) == 5)
        
            
if __name__ == '__main__':
    unittest.main()