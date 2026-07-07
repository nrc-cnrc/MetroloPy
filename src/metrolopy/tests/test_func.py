# -*- coding: utf-8 -*-

# module test_func

# Copyright (C) 2026 National Research Council Canada
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

import metrolopy as uc
import numpy as np
import unittest

class TestFunc(unittest.TestCase):
    def tstfunc(self,fa,fb,*p):
        if isinstance(fa(*p),(tuple,np.ndarray)):
            for i in range(len(fa(*p))):
                ifa = lambda *x:fa(*x)[i]
                ifb = lambda *x:fb(*x)[i]
                self.tstfunc(ifa,ifb,*p)
        else:
            if isinstance(p[0],list):
                pf = [[i.x for i in p[0]]]
                c = None
            else:
                pf = [i.x if isinstance(i,uc.gummy) else i for i in p]
                c = uc.gummy.napply(fb,*p)

            if fa(*pf) != 0:
                self.assertTrue(np.abs((fa(*p).x - fa(*pf))/fa(*pf)) < 1e-10)
                self.assertTrue(np.abs((fb(*p).x - fa(*pf))/fa(*pf)) < 1e-10)
                if c is not None:
                    self.assertTrue(np.abs((fa(*p).x - c.x)/fa(*p).x) < 1e-10)
            else:
                self.assertTrue(np.abs((fa(*p).x - fa(*pf))) < 1e-10)
                self.assertTrue(np.abs((fb(*p).x - fa(*pf))) < 1e-10)
                if c is not None:
                    self.assertTrue(np.abs((fa(*p).x - c.x)) < 1e-10)
    
            if fa(*p).u != 0:
                self.assertTrue(np.abs((fa(*p).u - fb(*p).u)/fa(*p).u) < 1e-6)
                if c is not None:
                    self.assertTrue(np.abs((fa(*p).u - c.u)/fa(*p).u) < 1e-6)
            else:
                self.assertTrue(np.abs((fa(*p).u - fb(*p).u)) < 1e-6)
                if c is not None:
                    self.assertTrue(np.abs((fa(*p).u - c.u)) < 1e-6)

    def test_sin(self):
        self.tstfunc(uc.sin,np.sin,uc.gummy(1.1,0.1))

    def test_cos(self):
        self.tstfunc(uc.cos,np.cos,uc.gummy(1.1,0.1))

    def test_tan(self):
        self.tstfunc(uc.tan,np.tan,uc.gummy(1.1,0.1))

    def test_arcsin(self):
        self.tstfunc(uc.arcsin,np.arcsin,uc.gummy(0.55,0.1))

    def test_arccos(self):
        self.tstfunc(uc.arccos,np.arccos,uc.gummy(0.55,0.1))

    def test_arctan(self):
        self.tstfunc(uc.arctan,np.arctan,uc.gummy(0.55,0.1))

    def test_arctan2(self):
        self.tstfunc(uc.arctan2,np.arctan2,uc.gummy(0.55,0.1),uc.gummy(0.54,0.2))

    def test_sinh(self):
        self.tstfunc(uc.sinh,np.sinh,uc.gummy(0.55,0.1))

    def test_cosh(self):
        self.tstfunc(uc.cosh,np.cosh,uc.gummy(3.2,0.1))

    def test_tanh(self):
        self.tstfunc(uc.tanh,np.tanh,uc.gummy(0.55,0.1))

    def test_arcsinh(self):
        self.tstfunc(uc.arcsinh,np.arcsinh,uc.gummy(0.55,0.1))

    def test_arccosh(self):
        self.tstfunc(uc.arccosh,np.arccosh,uc.gummy(3.2,0.1))

    def test_arctanh(self):
        self.tstfunc(uc.arctanh,np.arctanh,uc.gummy(0.55,0.1))

    def test_exp(self):
        self.tstfunc(uc.exp,np.exp,uc.gummy(1.1,0.1))

    def test_exp2(self):
        self.tstfunc(uc.exp2,np.exp2,uc.gummy(1.1,0.1))

    def test_expm1(self):
        self.tstfunc(uc.expm1,np.expm1,uc.gummy(1.1,0.1))

    def test_log(self):
        self.tstfunc(uc.log,np.log,uc.gummy(1.1,0.1))

    def test_log2(self):
        self.tstfunc(uc.log2,np.log2,uc.gummy(1.1,0.1))

    def test_log10(self):
        self.tstfunc(uc.log10,np.log10,uc.gummy(1.1,0.1))

    def test_log1p(self):
        self.tstfunc(uc.log1p,np.log1p,uc.gummy(1.1,0.1))

    def test_logaddexp(self):
        self.tstfunc(uc.logaddexp,np.logaddexp,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_logaddexp2(self):
        self.tstfunc(uc.logaddexp2,np.logaddexp2,uc.gummy(1.1,0.1),uc.gummy(1.1,0.1))

    def test_sqrt(self):
        self.tstfunc(uc.sqrt,np.sqrt,uc.gummy(1.1,0.1))

    def test_cbrt(self):
        self.tstfunc(uc.cbrt,np.cbrt,uc.gummy(1.1,0.1))

    def test_square(self):
        self.tstfunc(uc.square,np.square,uc.gummy(1.1,0.1))

    def test_add(self):
        self.tstfunc(uc.add,np.add,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_subtract(self):
        self.tstfunc(uc.subtract,np.subtract,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_negative(self):
        self.tstfunc(uc.negative,np.negative,uc.gummy(1.1,0.1))

    def test_multiply(self):
        self.tstfunc(uc.multiply,np.multiply,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_divide(self):
        self.tstfunc(uc.divide,np.divide,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_true_divide(self):
        self.tstfunc(uc.true_divide,np.true_divide,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_floor_divide(self):
        self.tstfunc(uc.floor_divide,np.floor_divide,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_reciprocal(self):
        self.tstfunc(uc.reciprocal,np.reciprocal,uc.gummy(1.1,0.1))

    def test_power(self):
        self.tstfunc(uc.power,np.power,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_absolute(self):
        self.tstfunc(uc.absolute,np.absolute,uc.gummy(1.1,0.1))

    def test_mod(self):
        self.tstfunc(uc.mod,np.mod,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_remainder(self):
        self.tstfunc(uc.remainder,np.remainder,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_divmod(self):
        self.tstfunc(uc.divmod,np.divmod,uc.gummy(1.1,0.1),uc.gummy(1.2,0.2))

    def test_modf(self):
        self.tstfunc(uc.modf,np.modf,uc.gummy(1.1,0.1))

    def test_angle(self):
        self.tstfunc(uc.angle,np.angle,uc.gummy(1.1,0.1))

    def test_real(self):
        self.tstfunc(uc.real,np.real,uc.gummy(1.1,0.1))

    def test_imag(self):
        self.tstfunc(uc.imag,np.imag,uc.gummy(1.1,0.1))

    def test_conj(self):
        self.tstfunc(uc.conj,np.conj,uc.gummy(1.1,0.1))

    def test_around(self):
        self.tstfunc(uc.around,np.around,uc.gummy(1.1,0.1))

    #def test_fix(self):
        #self.tstfunc(uc.fix,np.fix,uc.gummy(1.1,0.1))

    def test_floor(self):
        self.tstfunc(uc.floor,np.floor,uc.gummy(1.1,0.1))

    def test_ceil(self):
        self.tstfunc(uc.ceil,np.ceil,uc.gummy(1.1,0.1))

    def test_trunc(self):
        self.tstfunc(uc.trunc,np.trunc,uc.gummy(1.1,0.1))

    def test_heaviside(self):
        self.tstfunc(uc.heaviside,np.heaviside,uc.gummy(1.1,0.1),1)

    def test_sign(self):
        self.tstfunc(uc.sign,np.sign,uc.gummy(1.1,0.1))

    def test_sum(self):
        self.tstfunc(uc.sum,np.sum,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_prod(self):
        self.tstfunc(uc.prod,np.prod,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_cumsum(self):
        self.tstfunc(uc.cumsum,np.cumsum,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_cumprod(self):
        self.tstfunc(uc.cumprod,np.cumprod,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_diff(self):
        self.tstfunc(uc.diff,np.diff,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_ediff1d(self):
        self.tstfunc(uc.ediff1d,np.ediff1d,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    #def test_cross(self):
        #self.tstfunc(uc.cross,np.cross,[uc.gummy(1.1,0.1),uc.gummy(1.2,0.2),uc.gummy(1.3,0.3)])

    def test_misc(self):
        g = uc.gummy(-3,1)
        self.assertTrue(np.fmod(g,2).correlation(np.remainder(g,2)) == -1)
        self.assertTrue(np.fmod(g,2).x == -np.remainder(g,2).x)

        g = np.array([[uc.gummy(1.1,0.1) for i in range(3)] for j in range(2)])
        r = uc.arctan2(g,uc.gummy(2,2))
        self.assertTrue(r.shape == (2,3))
        r = uc.gummy.napply(np.arctan2,g,uc.gummy(2,2))
        self.assertTrue(r.shape == (2,3))
        r = uc.gummy.apply(np.sin,np.cos,g)
        self.assertTrue(r.shape == (2,3))