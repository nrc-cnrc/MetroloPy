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
Check that known bugs are fixed.
"""

import metrolopy as uc
import numpy as np
import unittest


class TestMisc(unittest.TestCase):

    def test_reduce(self):
        a = uc.gummy(12.367,0.22,unit='mm/m')
        a.reduce_unit()
        self.assertTrue(abs(a.x - 0.012367) < 1e-15)
        self.assertTrue(abs(a.u - 0.00022) < 1e-15)
        self.assertTrue(a.unit is uc.one)
        
    def test_zero_derivative(self):
        u = 0.00223
        a = uc.gummy(0,u)
        b = a**2
        
        self.assertTrue(b.u == 0)
        b.sim()
        self.assertTrue(abs((np.sqrt(2)*u**2-b.usim)/b.usim) < 0.1)
        b.style = 'pmsim'
        self.assertTrue('?' not in b.tostring())
        
    def test_triangular(self,plot=False):
        a = uc.gummy(uc.TriangularDist(1,1))
        a.sim()
        self.assertTrue(abs((a.usim - 1/np.sqrt(6))/a.usim) < 0.1)
        if plot:
            print('This should be a triangle with mode 1, lower limit 0 and upper limit 1:')
            a.hist()
            
        a = uc.gummy(uc.TriangularDist(-2.5,left_width=0.5,right_width=1.5))
        a.sim()
        if plot:
            print('this should be a triangle with mode -2.5, lower limit -3 and upper limit -1:')
            a.hist()
            
        a = uc.gummy(uc.TriangularDist(1,lower_limit=0.5,upper_limit=1.5))
        a.sim()
        if plot:
            print('this should be a triangle with mode 1, lower limit 0.5 and upper limit 1.5:')
            a.hist()
            
    def test_zero_degc(self):
        a = uc.gummy(0, 0.15,unit='degC')
        a.style='xf'
        self.assertTrue(a.tostring() == '0.00')
        
    def test_ufrom_multiletter(self):
        a = uc.gummy(10.0, 1, utype='A')
        b = uc.gummy(20.0, 2, utype='DUT')
        y = a - b
        self.assertTrue(abs(y.ufrom('A') - 1) < 1e-15)
        self.assertTrue(abs(y.ufrom('DUT') - 2) < 1e-15)
        
    def test_uniform_params(self):
        # issue 33
        g = uc.gummy(uc.UniformDist(lower_limit=6,upper_limit=9),unit='pA')
        g.sim()
        self.assertTrue(abs(min(g.simdata) - 6) < 0.01)
        self.assertTrue(abs(max(g.simdata) - 9) < 0.01)
        
        g = uc.gummy(uc.UniformDist(center=7.5,half_width=1.5),unit='pA')
        g.sim()
        self.assertTrue(abs(min(g.simdata) - 6) < 0.01)
        self.assertTrue(abs(max(g.simdata) - 9) < 0.01)
        
    def test_mean_b(self):
        # issue 32
        import numpy as np
        from scipy.stats import t
        from scipy import special
    
        con = special.erf(1/np.sqrt(2))
        def u(x):
            return t.ppf(1-(1-con)/2, len(x)-1)*x.std(ddof=1)/np.sqrt(len(x))
    
        x = np.random.rand(5)
        uc.gummy.bayesian = False
        g = uc.mean(x)
        g.p = 'ssd'
        self.assertTrue(abs(con-g.p) < 1e-6)
        self.assertTrue(abs(g.U - u(x)) < 1e-6)
        try:
            uc.gummy.bayesian = True
            gg = uc.mean(x)
            self.assertTrue(abs(gg.u - np.sqrt(2)*g.u) < 1e-6)
            gg.p = 'ssd'
            self.assertTrue(abs(gg.U - u(x)) < 1e-6)
        finally:
            uc.gummy.bayesian = False
            
    def test_uniform_hw(self):
        # issue 37
        a = uc.gummy(uc.UniformDist(upper_limit=1,lower_limit=-1))
        b = uc.gummy(uc.UniformDist(center=0,half_width=1))
        self.assertTrue(abs(a.u - b.u) < 1e-15)
        
        
if __name__ == '__main__':
    unittest.main()