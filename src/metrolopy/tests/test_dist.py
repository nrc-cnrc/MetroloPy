# -*- coding: utf-8 -*-

# module test_dist

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
import scipy.stats as stats

class TestDist(unittest.TestCase):
    def tstdist(self,mdist,sdist,params,perm=None,fperm=None,n=10000,prnt=False,ch_u=True,uabs=False,sfix=None,fix=None,p0=None,ch_x=True):
        params = np.asarray(params)
        sp0=None
        if fperm is None:
            sparams = params[np.asarray(perm)]
            if p0 is not None:
                sp0 = p0[np.asarray(perm)]
        else:
            sparams = fperm(*params)
            if p0 is not None:
                sp0 = fperm(*p0)
        a = uc.gummy(mdist(*params))
        b = uc.gummy(sdist(*sparams))

        if sp0 is None and sfix is not None:
            sp0 = [p if fx else None for p,fx in zip(sparams,sfix)]

        if p0 is None and fix is not None:
            p0 = [p if fx else None for p,fx in zip(params,fix)]

        if ch_x:
            self.assertTrue(np.abs(a.x - b.x)/a.u < 10**-8)
        if ch_u:
            self.assertTrue(np.abs(a.u - b.u)/a.u < 0.1)

        uc.simulate([a,b],n=n)

        fa = uc.DistFit(a.simdata,sdist,fix=sfix,p0=sp0)
        faa = uc.DistFit(a.simdata,mdist(*params),fix=fix,p0=p0)
        fb = uc.DistFit(b.simdata,mdist,fix=fix,p0=p0)
        fbb = uc.DistFit(b.simdata,sdist(*sparams),fix=sfix,p0=sp0)

        if fperm is None:
            fx = faa.p[perm]
        else:
            fx = fperm(*faa.p)
        for p1,p2 in zip(fa.p,fx):
            if isinstance(p1,uc.gummy) and p1.u != 0:
                self.assertTrue(np.abs(p1.x - p2.x)/(4*p1.u) < 1)

        if fperm is None:
            fx = fb.p[perm]
        else:
            fx = fperm(*fb.p)
        for p1,p2 in zip(fx,fbb.p):
            if isinstance(p2,uc.gummy) and p2.u != 0:
                self.assertTrue(np.abs(p1.x - p2.x)/(4*p2.u) < 1)

        if uabs:
            for p,p0 in zip(faa.p,params):
                self.assertTrue(np.abs(p.x - p0) < 0.1)
    
            for p,p0 in zip(fbb.p,sparams):
                self.assertTrue(np.abs(p.x - p0) < 0.1)
        else:
            for p,p0 in zip(faa.p,params):
                if isinstance(p,uc.gummy) and p.u != 0:
                    self.assertTrue(np.abs(p.x - p0)/(4*p.u) < 1)
    
            for p,p0 in zip(fbb.p,sparams):
                if isinstance(p,uc.gummy) and p.u != 0:
                    self.assertTrue(np.abs(p.x - p0)/(4*p.u) < 1)

        if prnt:
            print('Distrbution params: ',params)
            print('Distrbution value: ',a)
            print('stats params: ',sparams)
            print('stats value: ',b)
            print('---')
            print('stats fit:')
            print(fa)
            print('var:  ',fa.var)
            fa.plot()
            fa.plot_cdf()
            print('---')
            print('Distribution fit:')
            print(fb)
            print('var:  ',fb.var)
            fb.plot()
            fb.plot_cdf()

    def test_t(self,prnt=False):
        self.tstdist(uc.TDist,stats.t,[3.2,1.2,5],perm=[2,0,1],n=10000,prnt=prnt,ch_u=False)

    def test_norm(self,prnt=False):
        self.tstdist(uc.NormalDist,stats.norm,[3.2,1.2],perm=[0,1],n=10000,prnt=prnt)

    def test_uniform(self,prnt=False):
        self.tstdist(uc.UniformDist,stats.uniform,[3.2,1.2],fperm=lambda *a:[a[0]-a[1],2*a[1]],n=10000,prnt=prnt,uabs=True)

    def test_gamma(self,prnt=False):
        self.tstdist(uc.GammaDist,stats.gamma,[3.2,1.2],fperm=lambda *a:[a[0],0,a[1]],n=10000,prnt=prnt,sfix=[False,True,False],ch_u=False,p0=[1,1])

    def test_laplace(self,prnt=False):
        self.tstdist(uc.LaplaceDist,stats.laplace,[3.2,1.2],perm=[0,1],n=10000,prnt=prnt)

    #def test_triangular(self,prnt=False):
        #self.tstdist(uc.TriangularDist,stats.triang,[3.2,1.2,2.2],fperm=lambda *a:[a[1]/(a[1]+a[2]),a[0]-a[1],a[1]+a[2]],n=10000,prnt=prnt,ch_x=False)

    def test_exponential(self,prnt=False):
        self.tstdist(uc.ExponentialDist,stats.expon,[3.2],fperm=lambda *a:[0,a[0]],sfix=[True,False],n=10000,prnt=prnt)

    def test_trapezoidal(self,prnt=False):
        self.tstdist(uc.TrapezoidalDist,stats.trapezoid,[0.4,2.2,0.6],fperm=lambda *a:[(1-a[2])/2,(1+a[2])/2,a[0],a[1]-a[0]],n=10000,prnt=prnt,ch_x=False)

    def test_lognormal(self,prnt=False):
        self.tstdist(uc.LogNormalDist,stats.lognorm,[0,1.2],perm=[1,0],n=10000,prnt=prnt,sfix=[False,True],fix=[True,False])

    def test_weibull(self,prnt=False):
        self.tstdist(uc.WeibullDist,stats.weibull_min,[3,1.2],fperm=lambda *a:[a[0],0,a[1]],sfix=[False,True,False],n=10000,prnt=prnt,ch_x=False,ch_u=False)

    def test_arcsin(self,prnt=False):
        p = [1.2,0.3]
        a = uc.gummy(uc.ArcSinDist(*p))
        a.sim()
        self.assertTrue(np.abs((a.xsim - p[0])/p[0]) < 1e-2)
        self.assertTrue(np.min(a.simdata) >= p[0] - p[1])
        self.assertTrue(np.max(a.simdata) <= p[0] + p[1])
        if prnt:
            print('params: ',p)
            a.hist()

    def test_ctrap(self,prnt=False):
        p = [1.2,0.3,0.1]
        a = uc.gummy(uc.CurvlinearTrapDist(*p))
        a.sim()
        self.assertTrue(np.abs((a.xsim - p[0])/p[0]) < 1e-2)
        self.assertTrue(np.min(a.simdata) >= p[0] - p[1] - p[2])
        self.assertTrue(np.max(a.simdata) <= p[0] + p[1] + p[2])
        if prnt:
            print('params: ',p)
            a.hist()

    def test_poisson(self,prnt=False):
        p = 8.2
        a = uc.gummy(uc.PoissonDist(p))
        b = uc.gummy(stats.poisson(p))
        uc.simulate([a,b])

        self.assertTrue(np.abs((a.xsim - a.x)/a.x <= 0.01))
        self.assertTrue(np.abs((a.usim - a.u)/a.u <= 0.01))
        self.assertTrue(np.abs((b.xsim - b.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.usim - b.u)/b.u <= 0.01))
        self.assertTrue(np.abs((b.x - a.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.u - a.u)/b.u <= 0.01))
        
        if prnt:
            print('param: ',p)
            a.hist()
            b.hist()

    def test_binom(self,prnt=False):
        p = [8,0.4]
        a = uc.gummy(uc.BinomialDist(*p))
        b = uc.gummy(stats.binom(*p))
        
        uc.simulate([a,b])

        self.assertTrue(np.abs((a.xsim - a.x)/a.x <= 0.01))
        self.assertTrue(np.abs((a.usim - a.u)/a.u <= 0.01))
        self.assertTrue(np.abs((b.xsim - b.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.usim - b.u)/b.u <= 0.01))
        self.assertTrue(np.abs((b.x - a.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.u - a.u)/b.u <= 0.01))
                        
        if prnt:
            print('params: ',p)
            a.hist()
            b.hist()

    def test_averaged(self,prnt=False):
        a = uc.gummy(uc.AveragedDist(uc.NormalDist(1,0.3),6))
        b = uc.gummy(uc.AveragedDist(stats.norm(1,0.3),6))

        uc.simulate([a,b])

        self.assertTrue(np.abs((a.x - 1)/1 <= 1e-6))
        u = 0.3/np.sqrt(6)
        self.assertTrue(np.abs((a.u - u)/u <= 1e-6))
        self.assertTrue(np.abs((a.xsim - a.x)/a.x <= 0.01))
        self.assertTrue(np.abs((a.usim - a.u)/a.u <= 0.01))
        self.assertTrue(np.abs((b.xsim - b.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.usim - b.u)/b.u <= 0.01))
        self.assertTrue(np.abs((b.x - a.x)/b.x <= 0.01))
        self.assertTrue(np.abs((b.u - a.u)/b.u <= 0.01))

        if prnt:
            print('params: ','normal',[1,0.3])
            a.hist()
            b.hist()