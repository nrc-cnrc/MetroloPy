# -*- coding: utf-8 -*-

# module test_operations

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

import metrolopy as uc
import numpy as np

from metrolopy.tests.common import rand,make_gummy,make_number

def binary_func(f,df,sim=False,exp=None,fionly=False,uexp=-6,allowazero=True,
                allowbzero=True,allowlargeu=True):
    if sim:
        dof = dof=float('inf')
    else:
        dof = None
    a,ax,au,adof,_,_ = make_gummy(unit=1,dof=dof,exp=exp,mpf=(not fionly),
                                  uexp=uexp,allowzero=allowazero,
                                  allowlargeu=allowlargeu)
    if rand.randint(2):
        b = make_number(exp=exp,fionly=fionly,allowzero=allowbzero)
        if isinstance(b,uc.gummy):
            bx = b.x
        else:
            bx = b
        bu = 0
        bdof = float('inf')
    else:
        b,bx,bu,bdof,_,_ = make_gummy(unit=1,dof=dof,exp=exp,mpf=(not fionly),
                                      uexp=uexp,allowzero=allowbzero,
                                      allowlargeu=allowlargeu)
        
    ax = float(ax)
    bx = float(bx)
    au = float(au)
    bu = float(bu)
    x = f(ax,bx)
    d = df(ax,bx)
    da = float(d[0])
    db = float(d[1])
    u = np.sqrt((da*au)**2 + (db*bu)**2)
    if dof is None:
        dof = ((da*au)**4/adof + (db*bu)**4/bdof)
        if dof == 0:
            dof = float('inf')
        else:
            dof = u**4/dof
            
        if dof > uc.gummy.max_dof:
            dof = float('inf')
        if dof < 1:
            dof = 1
            
        dof = float(dof)
            
    g = f(a,b)
    
    if x == 0:
        assert abs(g.x) < 1e-15
    else:
        assert abs((g.x - x)/x) < 1e-14
        
    if u == 0:
        assert abs(u) < 1e-6
    else:
        assert abs((g.u - u)/u) < 1e-3
    
    if dof is not None and not np.isinf(dof):
        assert abs((g.dof - dof)/dof) < 0.01
    
    assert a.correlation(b) == 0
    
    if sim and x != 0 and abs(u/x) > 1e-10:
        uc.gummy.simulate([a,b,g])
        
        assert abs((g.xsim - x)/(g.u/np.sqrt(1e5))) < 5
        assert abs((g.u - g.usim)/g.u) < 0.1
        
        if g.correlation(a) > 0.1:
            assert abs((g.correlation_sim(a)-g.correlation(a))/g.correlation(a)) < 0.1
            
        if isinstance(b,uc.gummy) and g.correlation(b) > 0.1:
            assert abs((g.correlation_sim(b)-g.correlation(b))/g.correlation(b)) < 0.1
    
def test_add(n=1000,sim=False):
    def f(a,b):
        return a + b
    def df(ax,bx):
        return (1,1)
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_sub(n=1000,sim=False):
    def f(a,b):
        return a - b
    def df(ax,bx):
        return (1,-1)
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_mul(n=1000,sim=False):
    def f(a,b):
        return a*b
    def df(ax,bx):
        return (bx,ax)
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_div(n=1000,sim=False):
    def f(a,b):
        return a/b
    def df(ax,bx):
        return (1/bx,-ax/bx**2)
    if sim:
        lu = False
    else:
        lu = True
    for i in range(n):
        binary_func(f,df,sim=sim,allowbzero=False,allowlargeu=lu)
        
def test_pow(n=1000,sim=False):
    def f(a,b):
        return abs(a)**b
    def df(ax,bx):
        ax = abs(ax)
        return (bx*ax**(bx-1),np.log(ax)*ax**bx)
    if sim:
        lu = False
    else:
        lu = True
    for i in range(n):
        binary_func(f,df,sim=sim,exp=0,allowazero=False,allowlargeu=lu)
        
def test_mod(n=1000,sim=False):
    def f(a,b):
        return a%b
    def df(ax,bx):
        return (1,np.sign(bx)*abs(ax//bx))
    for i in range(n):
        binary_func(f,df,sim=sim,fionly=True)
    
def test_abs(n=1000,sim=False):
    def f(a,b):
        return abs(a)*abs(b)
    def df(ax,bx):
        return (abs(bx),abs(ax))
    if sim:
        lu = False
    else:
        lu = True
    for i in range(n):
        binary_func(f,df,sim=sim,allowlargeu=lu,allowazero=lu,
                    allowbzero=lu)
        
def test_neg(n=1000,sim=False):
    def f(a,b):
        return (-a)*b
    def df(ax,bx):
        return (-bx,-ax)
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_pos(n=1000,sim=False):
    def f(a,b):
        return (+a)*b
    def df(ax,bx):
        return (bx,ax)
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_sincos(n=1000,sim=False):
    def f(a,b):
        return uc.sin(a) + uc.cos(b)
    def df(ax,bx):
        return (uc.cos(ax),-uc.sin(bx))
    for i in range(n):
        binary_func(f,df,sim=sim)
        
def test_npsincos(n=1000,sim=False):
    def f(a,b):
        return np.sin(a) + np.cos(b)
    def df(ax,bx):
        return (np.cos(ax),-np.sin(bx))
    for i in range(n):
        binary_func(f,df,sim=sim,fionly=True)
        
def test_ap1sincos(n=1000,sim=False):
    def f(a,b):
        return uc.gummy.apply(np.sin,np.cos,a) + uc.gummy.napply(np.cos,b)
    def df(ax,bx):
        return (np.cos(ax),-np.sin(bx))
    for i in range(n):
        binary_func(f,df,sim=sim,fionly=True,exp=0,allowlargeu=False)
        
def test_ap2sincos(n=1000,sim=False):
    def ff(a,b):
        return np.sin(a) + np.cos(b)
    def df(ax,bx):
        return (np.cos(ax),-np.sin(bx))
    def f(a,b):
        return uc.gummy.apply(ff,df,a,b)

    for i in range(n):
        binary_func(f,df,sim=sim,fionly=True)
        
def test_apnsincos(n=1000,sim=False):
    def ff(a,b):
        return np.sin(a) + np.cos(b)
    def df(ax,bx):
        return (np.cos(ax),-np.sin(bx))
    def f(a,b):
        return uc.gummy.napply(ff,a,b)

    for i in range(n):
        binary_func(f,df,sim=sim,fionly=True,exp=0,allowlargeu=False)
        
def test_addxmul(n=1000):
    def f(*x):
        r = x[0]
        for i in range(len(x)-1):
            r += x[i+1]
        return r
        
    for j in range(n):
        k = rand.randint(10) + 2
        x = make_gummy(unit=1)[0]
        a = f(*(k*[x]))
        b = k*x
        c = uc.gummy.apply(f,lambda *x: len(x)*[1],*(k*[x]))
        d = uc.gummy.napply(f,*(k*[x]))
        
        if b.x == 0:
            assert abs(a.x - b.x) < 1e-14
        else:
            assert abs((a.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert a.u == 0
            assert b.u == 0
        else:
            assert (a.u - b.u)/b.u < 1e-3
            assert a.correlation(b) == 1
            assert a.correlation(x) == 1
        
        if b.x == 0:
            assert abs(c.x - b.x) < 1e-14
        else:
            assert abs((c.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert c.u == 0
        else:
            assert (c.u - b.u)/b.u < 1e-3
            assert c.correlation(b) == 1
        
        if b.x == 0:
            assert abs(d.x - b.x) < 1e-14
        else:
            assert abs((d.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert d.u == 0
        else:
            assert abs((d.u - b.u)/b.u) < 1e-3
            assert d.correlation(b) == 1
        
def test_mulxpow(n=1000):
    def f(*x):
        r = x[0]
        for i in range(len(x)-1):
            r *= x[i+1]
        return r
        
    for j in range(n):
        k = rand.randint(10) + 2
        x = make_gummy(unit=1)[0]
        a = f(*(k*[x]))
        b = x**k
        xx = x.x
        c = uc.gummy.apply(f,lambda *x: len(x)*[xx**(k-1)],*(k*[x]))
        d = uc.gummy.napply(f,*(k*[x]))
        
        if b.x == 0:
            assert abs(a.x - b.x) < 1e-14
        else:
            assert abs((a.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert a.u == 0
            assert b.u == 0
        else:
            assert abs((a.u - b.u)/b.u) < 1e-3
            assert a.correlation(b) == 1
        
        if k%2 == 0 and x.x < 0:
            assert a.correlation(x) == -1
        else:
            assert a.correlation(x) == 1
            
        if b.x == 0:
            assert abs(c.x - b.x) < 1e-14
        else:
            assert abs((c.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert c.u == 0
        else:
            assert abs((c.u - b.u)/b.u) < 1e-3
            assert c.correlation(b) == 1
        
        if b.x == 0:
            assert abs(d.x - b.x) < 1e-14
        else:
            assert abs((d.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert d.u == 0
        else:
            assert abs((d.u - b.u)/b.u) < 1e-3
            assert d.correlation(b) == 1

        
def test_mulxnpow(n=1000):
    def f(*x):
        r = 1/x[0]
        for i in range(len(x)-1):
            r *= 1/x[i+1]
        return r
        
    for j in range(n):
        k = rand.randint(10) + 2
        x = make_gummy(unit=1,allowzero=False,allowlargeu=False)[0]
        a = f(*(k*[x]))
        b = x**-k
        d = uc.gummy.napply(f,*(k*[x]))
        
        if b.x == 0:
            assert abs(a.x - b.x) < 1e-14
        else:
            assert abs((a.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert a.u == 0
            assert b.u == 0
        else:
            assert abs((a.u - b.u)/b.u) < 1e-3
            assert a.correlation(b) == 1
        
        if k%2 == 0 and x.x < 0:
            assert a.correlation(x) == 1
        else:
            assert a.correlation(x) == -1
        
        if b.x == 0:
            assert abs(d.x - b.x) < 1e-14
        else:
            assert abs((d.x - b.x)/b.x) < 1e-14
        if x.u == 0:
            assert d.u == 0
        else:
            assert abs((d.u - b.u)/b.u) < 1e-3
            assert d.correlation(b) == 1