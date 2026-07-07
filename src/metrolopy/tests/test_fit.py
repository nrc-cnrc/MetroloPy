# -*- coding: utf-8 -*-

# module test_fit

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

def check_cov(cov1,cov2,dabs=None,drel = 0.0001):
    d = np.abs(cov1 - cov2)
    if dabs is not None:
        assert np.max(d) < dabs
        return
    d = np.array([[d[i][j]/np.sqrt(cov1[i][i]*cov1[j][j]) for i in range(len(cov1[0]))] for j in range(len(cov1))])
    assert np.max(d) < drel

class TestFit(unittest.TestCase):
        
    def test_poly(self,plot=False):
        rng = np.random.default_rng()
        x = np.linspace(0,1,20) - 0.2
        p = np.array([3.0,0.5,2.0])
        sd = 0.01
        y = p[0] + p[1]*x + p[2]*x**2 + rng.normal(0,sd,len(x))
        
        fit0 = uc.PolyFit(x,y,0)
        mn = uc.mean(y,use_n_eff=False)
        self.assertTrue(np.abs((fit0.pf[0] - mn.x)/mn.x) < 1e-6)
        self.assertTrue(np.abs((fit0.p[0].u - mn.u)/mn.u) < 1e-6)
        self.assertTrue(np.abs((fit0.p[0].dof - mn.dof)/mn.dof) < 0.01)
        
        fit1 = uc.PolyFit(x,y,2,xunit='s',yunit='m')
        fit2 = uc.PolyFit(x,y,2,solver='nls')
        fit2b = uc.PolyFit(x,y,2,_ch_create_p=True)
        fit3 = uc.PolyFit(x,y,2,solver='odr')
        fit3b = uc.PolyFit(x,y,2,solver='odr',_ch_create_p=True)

        p0 = p + rng.normal(0,0.1,3)
        def f(x,p0,p1,p2):
                return p0 + p1*x + p2*x**2
        class pfit(uc.Fit):
            def f(self,*p):
                return f(*p)
        fit4 = pfit(x,y,p0=p0)
        fit5 = uc.Fit(x,y,p0=p0,f=f)
        fit5a = uc.Fit(x,y,p0=p0,f=f,solver='ols')
        fit6 = uc.Fit(x,y,p0=p0,f=f,solver='odr')
        fit7 = uc.Fit(x,y,p0=p0,f=f,solver='odr',fit_type='OLS')
        fit8 = uc.Fit(x,y,p0=p0,f=f,solver='odr',xweights=1,weights=0.0001)

        self.assertTrue(np.shape(fit1.ypredf(1)) == ())
        self.assertTrue(np.shape(fit1.ypredf([1,2,3])) == (3,))
        self.assertTrue(isinstance(fit1.ypred(1),uc.gummy))
        self.assertTrue(np.shape(fit1.ypred([1,2,3])) == (3,))
        self.assertTrue(np.max(np.abs((fit1.pf - p)/np.array([x.u for x in fit1.p])) < 4))
        self.assertTrue(np.max(np.abs((fit1.pf - p)/p) < 0.1))
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit2.pf)/(fit1.pf + fit2.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit3.pf)/(fit1.pf + fit3.pf)) < 0.001))
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit4.pf)/(fit1.pf + fit4.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit5.pf)/(fit1.pf + fit5.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit5a.pf)/(fit1.pf + fit5.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit3.pf - fit6.pf)/(fit1.pf + fit6.pf)) < 0.0001))
        check_cov(fit1.cov,fit2.cov)
        check_cov(fit1.cov,uc.gummy.covariance_matrix(fit1.p))
        check_cov(fit3b.cov,uc.gummy.covariance_matrix(fit3b.p))
        check_cov(fit5.cov,uc.gummy.covariance_matrix(fit5.p))
        check_cov(fit6.cov,uc.gummy.covariance_matrix(fit6.p),drel=0.01)
        check_cov(fit1.cov,fit4.cov)
        check_cov(fit1.cov,fit5.cov)
        check_cov(fit1.cov,fit5a.cov)
        check_cov(fit3.cov,fit6.cov,drel=0.01)
        check_cov(fit1.cov,fit7.cov)
        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit7.pf)/(fit1.pf + fit7.pf)) < 0.0001))
        check_cov(fit1.cov,fit8.cov,drel=0.01)
        check_cov(fit8.cov,uc.gummy.covariance_matrix(fit8.p),drel=0.01)
        
        if plot:
            odof = uc.gummy.show_dof
            try:
                uc.gummy.show_dof = True
                fit1.plot(clp=0.999,cik=4)
                print('PolyFit(OLS):')
                print(fit1)
                print('PolyFit(NLS):')
                print(fit2)
                print('PolyFit(_ch_create_p=True):')
                print(fit2b)
                print('PolyFit(ODR):')
                print(fit3)
                print('Fit subclassed(NLS):')
                print(fit4)
                print('Fit(NLS):')
                print(fit5)
                print('Fit(OLS):')
                print(fit5a)
                print('Fit(ODR):')
                print(fit6)
                print('Fit(ODR,fit_type=OLS):')
                print(fit7)
                print('Fit(ODR,xweights=1,weights=0.0001):')
                print(fit8)
            finally:
                uc.gummy.show_dof = odof
                
    def test_poly_2d(self,prnt=False):
        rng = np.random.default_rng()
        x = np.mgrid[-1:1:0.05, -1:1:0.05].reshape(2,-1)
        x = (x.T + np.array([2,-3])).T
        p = np.array([3.0,0.5,2.0,1.1,2.2,1.5])
        y = p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[1] + p[4]*x[1]**2 + p[5]*x[0]*x[1] + rng.normal(0,0.01,len(x[0]))
        p0 = p + rng.normal(0,0.1,6)
        def f(*x):
                return x[2] + x[3]*x[0] + x[4]*x[0]**2 + x[5]*x[1] + x[6]*x[1]**2 + x[7]*x[0]*x[1]
            
        def f2(*x):
                return x[2] + x[3]*x[0] + x[4]*x[0]**2 + x[5]*x[1] + x[6]*x[0]*x[1] + x[7]*x[0]**2*x[1] + x[8]*x[1]**2 + x[9]*x[0]*x[1]**2 + x[10]*x[0]**2*x[1]**2
        
        fit1 = uc.PolyFit(x,y,deg=[2,2])
        fit2 = uc.PolyFit(x,y,deg=[2,2],solver='nls')
        fit3 = uc.Fit(x,y,p0=p0,f=f)
        p0_2 = fit2.pf + rng.normal(0,0.1,9)
        fit4 = uc.Fit(x,y,p0=p0_2,f=f2)

        p2 = np.array([3.0,0.5,2.0,1.1,1.5,0.7])
        y2 = p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[1] + p[4]*x[0]*x[1] + p[5]*x[0]**2*x[1] + rng.normal(0,0.01,len(x[0]))
        fit5 = uc.PolyFit(x,y2,deg=[2,1])
        fit6 = uc.PolyFit(x,y2,deg=[2,1],xunit=['m','s'],yunit='kg')

        xg = [[uc.gummy(i,unit='m') for i in x[0]],[uc.gummy(i,unit='s') for i in x[1]]]
        yg = [uc.gummy(i,unit='kg') for i in y2]
        fit7 =  uc.PolyFit(xg,yg,deg=[2,1])

        self.assertTrue(np.shape(fit1.ypredf(1,2)) == ())
        self.assertTrue(np.shape(fit1.ypredf([1,2,3],[1,2,3])) == (3,))
        
        check_cov(fit1.cov,fit2.cov)
        check_cov(fit1.cov,fit4.cov)

        self.assertTrue(np.max(np.abs(2*(fit1.pf - fit2.pf)/(fit1.pf + fit2.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit4.pf - fit2.pf)/(fit4.pf + fit2.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs((fit3.pf - p)/p) < 0.01))
        self.assertTrue(np.max(np.abs((fit5.pf - p2)/p2) < 0.01))

        self.assertTrue(fit6.p[0].unit is uc.unit('kg'))
        self.assertTrue(fit7.p[0].unit is uc.unit('kg'))
        self.assertTrue(fit6.p[1].unit is uc.unit('kg/m'))
        self.assertTrue(fit7.p[1].unit is uc.unit('kg/m'))
        self.assertTrue(fit6.p[3].unit is uc.unit('kg/s'))
        self.assertTrue(fit7.p[3].unit is uc.unit('kg/s'))

        if prnt:
            print('PolyFit(OLS,deg=[2,2]):')
            print(fit1)
            print('PolyFit(NLS):')
            print(fit2)
            print('Fit(non-zero p only):')
            print(fit3)
            print('Fit:')
            print(fit4)
            print('PolyFit(deg=[2,1]):')
            print(fit5)
            print('PolyFit(xunit=("m","s"),yunit="kg")')
            print(fit6)
            print('PolyFit with x and y constant gummies')
            print(fit7)

    def test_poly_3d(self,prnt=False):
        rng = np.random.default_rng()
        x = np.mgrid[-1:1:0.05, -1:1:0.05,-1:1:0.05].reshape(3,-1)
        x = (x.T + np.array([2,3,-2])).T
        p = np.array([3.0,0.5,2.0,1.1,2.2,1.5,-1.7,1.3])
        y = p[0] + p[1]*x[0] + p[2]*x[1] + p[3]*x[0]*x[1] + p[4]*x[2] + p[5]*x[0]*x[2] + p[6]*x[1]*x[2] + p[7]*x[0]*x[1]*x[2] + rng.normal(0,0.01,len(x[0]))

        fit1 = uc.PolyFit(x,y,deg=[1,1,1])

        def f(*x):
            return x[3] + x[4]*x[0] + x[5]*x[1] + x[6]*x[0]*x[1] + x[7]*x[2] + x[8]*x[0]*x[2] + x[9]*x[1]*x[2] + x[10]*x[0]*x[1]*x[2]
        p0 = p + rng.normal(0,0.1,8)
        fit2 = uc.Fit(x,y,f=f,p0=p0)

        self.assertTrue(np.max(np.abs((fit1.pf - p)/p) < 0.01))
        check_cov(fit1.cov,fit2.cov)
        
        if prnt:
            print('PolyFit(OLS):')
            print(fit1)
            print('Fit:')
            print(fit2)

    def test_ygummies(self,prnt=False):
        rng = np.random.default_rng()
        x = np.linspace(0.1,2,8)
        p = np.array([3.0,2.0])
        uy = 0.3*x + 0.3
        y = p[0] + p[1]*x + rng.normal(0,1,len(x))*uy
        y = [y[i] + uc.gummy(0,uy[i]) for i in range(len(x))]

        fit1 = uc.PolyFit(x,y,0)
        fit2 = uc.PolyFit(x,y,0,solver='nls')
        def f(x,*p):
            return(p[0])
        fit3 = uc.Fit(x,y,f=f,p0=[0])
        
        wm = uc.wmean(y)

        fit4 = uc.PolyFit(x,y,1)
        fit4b = uc.PolyFit(x,y,1,ux=0.2)
        fit5 = uc.PolyFit(x,y,1,solver='nls')

        def f(x,*p):
            return(p[0] + x*p[1])
        fit6 = uc.Fit(x,y,f=f,p0=[0,0])

        xo = [uc.gummy(i,u=0.0001) for i in x]
        fit6b = uc.Fit(xo,y,f=f,p0=[0,0],solver='odr')

        y2 = [y[i].x + uc.gummy(0,5*uy[i]) for i in range(len(x))]

        fit7 = uc.PolyFit(x,y2,1)

        x = np.mgrid[-1:1:0.1, -1:1:0.1].reshape(2,-1)
        x = (x.T + np.array([2,-3])).T
        p2d = np.array([3.0,0.5,2.0,1.1])
        uy = 0.03*np.sqrt(x[0]**2 + x[1]**2) + 0.03
        yi = p2d[0] + p2d[1]*x[0] + p2d[2]*x[1] + p2d[3]*x[0]*x[1] #+ rng.normal(0,1,len(x[0]))*uy
        y = [yi[i] + uc.gummy(0,uy[i]) for i in range(len(x[0]))]

        fit8 = uc.PolyFit(x,yi,deg=[1,1])
        fit9 = uc.PolyFit(x,y,deg=[1,1])

        ux = [0.01,0.02,0.03]
        xg = [[x[j][i] + uc.gummy(0,ux[j]) for i in range(len(x[0]))] for j in range(len(x))]
        fit9b = uc.PolyFit(xg,y,deg=[1,1])
        fit9c = uc.PolyFit(x,yi,deg=[1,1],ux=0.01,uy=0.03)
        fit9d = uc.PolyFit(x,yi,deg=[1,1],ux=0.01,uy=0.03,_ch_create_p=True)

        x = np.linspace(0.1,2,8)
        p = np.array([3.0,2.0])
        uy = 0.5
        ux = 0.2
        y = p[0] + p[1]*x + rng.normal(0,1,len(x))*uy
        x = x + rng.normal(0,1,len(x))*ux
        yg = [uc.gummy(i,uy) for i in y]
        xg = [uc.gummy(i,ux) for i in x]

        fit10 = uc.PolyFit(x,y,1,ux=ux,uy=uy,solver='odr',_ch_create_p=True)
        fit11 = uc.PolyFit(xg,yg,1,solver='odr')
        fit12 = uc.PolyFit(y,x,1,ux=uy,uy=ux,solver='odr',_ch_create_p=True)
        fit13 = uc.PolyFit(x,y,1,ux=ux,uy=uy)
        fit14 = uc.PolyFit(x,y,1,ux=ux,uy=uy,solver='nls')
        fit15 = uc.PolyFit(x,y,1)

        self.assertTrue(np.abs(fit1.p[0].x - wm.x) < 1e-6)
        self.assertTrue(np.abs((fit1.p[0].u - wm.u)/wm.u) < 1e-6)
        self.assertTrue(np.abs(fit2.p[0].x - wm.x) < 1e-6)
        self.assertTrue(np.abs((fit2.p[0].u - wm.u)/wm.u) < 1e-6)
        self.assertTrue(np.abs(fit3.p[0].x - wm.x) < 1e-6)
        self.assertTrue(np.abs((fit3.p[0].u - wm.u)/wm.u) < 1e-6)
        check_cov(fit1.cov,fit2.cov)
        check_cov(fit4.cov,fit5.cov)
        check_cov(fit4.cov,fit6.cov)
        check_cov(fit6b.cov,uc.gummy.covariance_matrix(fit6b.p))
        self.assertTrue(np.max(np.abs(2*(fit4.pf - fit5.pf)/(fit4.pf + fit5.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit4.pf - fit6.pf)/(fit4.pf + fit6.pf)) < 0.0001))
        self.assertTrue(np.max(np.abs(2*(fit4.pf - fit7.pf)/(fit4.pf + fit7.pf)) < 0.0001))
        check_cov(fit1.cov,uc.gummy.covariance_matrix(fit1.p))
        check_cov(fit2.cov,uc.gummy.covariance_matrix(fit2.p))
        check_cov(25*fit6.cov,fit7.cov)
        self.assertTrue(np.max(np.abs((fit4.pf - p)/np.array([x.u for x in fit4.p])) < 4))
        self.assertTrue(np.max(np.abs((fit5.pf - p)/np.array([x.u for x in fit5.p])) < 4))
        self.assertTrue(np.max(np.abs((fit6.pf - p)/np.array([x.u for x in fit6.p])) < 4))
        self.assertTrue(np.max(np.abs((fit8.pf - p2d)/np.array([x.u for x in fit8.p])) < 4))

        a = fit10.ypred(0)
        b = -fit12.p[0]/fit12.p[1]
        self.assertTrue(np.abs(a.x - b.x) < 1e-2)
        self.assertTrue(np.abs(a.u - b.u) < 1e-2)

        check_cov(fit4b.cov,uc.gummy.covariance_matrix(fit4b.p))
        check_cov(fit9c.cov,uc.gummy.covariance_matrix(fit9c.p))
        check_cov(fit9c.cov,uc.gummy.covariance_matrix(fit9d.p))
        
        if prnt:
            print('PolyFit(deg = 0)')
            fit1.plot(cik=2)
            print(fit1)
            print('PolyFit(deg = 1)')
            fit4.plot(cik=2,clk=2)
            print(fit4)
            print(fit6b)
            print(fit10)
            print(fit11)
            print(fit12)
            print(fit13)
            print(fit14)
            print(fit15)
            print(fit9b)

    def test_odr2dy(self,prnt=False):
        rng = np.random.default_rng()
        x = np.linspace(0.1,2,8)
        p1 = [3.0,2.0]
        p2 = [4.1,1.2]
        y1 = p1[0] + p1[1]*x + rng.normal(0,0.1,len(x))
        y2 = p2[0] + p2[1]*x + rng.normal(0,0.1,len(x))

        p0 = [3.1,2.2,4.0,1.0]
        y = np.array([y1,y2])

        def f(x,*p):
            return [p[0] + p[1]*x,p[2] + p[3]*x]

        fit1 = uc.Fit(x,y,f=f,p0 = p0,solver='odr')
        if prnt:
            print(fit1)