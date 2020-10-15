# -*- coding: utf-8 -*-

# module test_complex

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

import numpy as np
from metrolopy import immy,ummy,jummy,gummy
from metrolopy.tests.common import display

rand = np.random.RandomState()

def make_immy(prnt=False):
    if rand.randint(2):
        imake = immy
        make = ummy
    else:
        imake = jummy
        make = gummy
    r = 4*rand.rand() + 1
    if rand.randint(2):
        r = -r
    ru = 0.1*rand.rand()*(abs(r))
    i = 4*rand.rand() + 1
    if rand.randint(2):
        i = -i
    iu = 0.1*rand.rand()*(abs(i))
    c = (2*rand.rand() - 1)*ru*iu

    cov = [[ru**2,c],[c,iu**2]]
    real,imag = make.create([r,i],covariance_matrix=cov)
    
    if rand.randint(2):
        ret = imake(real=real,imag=imag)
    else:
        if rand.randint(2):
            ret = imake(real=r,imag=i,cov=cov)
            assert abs(ret.real.covariance(ret.imag) - cov[0][1]) < 1e-6
            if imake is immy:
                real._ref = ret.real._ref
                imag._ref = ret.imag._ref
            else:
                real.value._ref = ret.real.value._ref
                imag.value._ref = ret.imag.value._ref
        else:
            r = np.sqrt(real**2 + imag**2)
            phi = np.arctan2(imag,real)
            ret = imake(r=r.x,phi=phi.x,cov=make.covariance_matrix([r,phi]))
            assert abs(ret.real.covariance(ret.imag) - cov[0][1]) < 1e-6
            if imake is immy:
                real._ref = ret.real._ref
                imag._ref = ret.imag._ref
            else:
                real.value._ref = ret.real.value._ref
                imag.value._ref = ret.imag.value._ref
                
    if prnt:
        display(ret)
    
    return (ret,real,imag)

def assert_ummy_close(u1,u2):
        assert abs(u1.correlation(u2)) > 1 - 1e-4
        assert abs((u1.x - u2.x)/(u1.x)) < 1e-10
        assert abs((u1.u - u2.u)/(u1.u)) < 1e-4
        
        if u1.dof == float('inf'):
            assert u2.dof == float('inf')
        else:
            assert abs((u2.dof - u1.dof)/u1.dof) < 1e-2
        
def assert_immy_close(i1,i2):
        assert_ummy_close(i1.real,i2.real)
        assert_ummy_close(i1.imag,i2.imag)

def test_immy_init(n=1000,prnt=False):
    for m in range(n):
        x,xr,xi = make_immy(prnt=prnt)
        
        if rand.randint(2):
            if rand.randint(2):
                x = type(x)(x)
            else:
                x = x.copy(formatting=False)
        
        assert_ummy_close(x.real,xr)
        assert_ummy_close(x.imag,xi)
        assert_ummy_close(x.angle(),np.arctan2(xi,xr))
        assert_ummy_close(abs(x),np.sqrt(xr**2 + xi**2))
        assert_immy_close(x.conjugate(),immy(real=xr,imag=-xi))
        assert_immy_close(-x,immy(real=-xr,imag=-xi))
        assert_immy_close(+x,immy(real=xr,imag=xi))
        
        if prnt:
            if rand.randint(2):
                y =1e12*x
            else:
                y = x/1e12
            display(y)

def _test_immy_bop(f,nf,n=1000,prnt=False):
    for m in range(n):
        a,ar,ai = make_immy()
        if True:#rand.randint(2):
            b,br,bi = make_immy()
        else:
            if rand.randint(2):
                b = 4*rand.rand() + 1
                if rand.randint(2):
                    b = -b
                if rand.randint(2):
                    bu = 0.1*rand.rand()*(abs(b))
                    b = ummy(b,u=bu)
            else:
                br = 4*rand.rand() + 1
                if rand.randint(2):
                    br = -br
                bi = 4*rand.rand() + 1
                if rand.randint(2):
                    bi = -bi
                b = complex(br,bi)
            if rand.randint(2):
                bb = b
                b = a
                a = bb
                
        if isinstance(a,(immy,ummy)):
            ax = a.x
        else:
            ax = a
        if isinstance(b,(immy,ummy)):
            bx = b.x
        else:
            bx = b
        
            
        cx = f(ax,bx)
        c = f(a,b)
        cn = type(c).napply(nf,a,b)
        
        if prnt:
            display(a)
            display(b)
            display(c)
            display(cn)
            print('---')

        assert abs((c.real.x - cx.real)/cx.real) < 1e-10
        assert abs((c.imag.x - cx.imag)/cx.imag) < 1e-10
        assert_immy_close(c,cn)
    
def test_immy_add(n=1000,prnt=False):
    _test_immy_bop(lambda a,b: a + b,np.add,n,prnt)
    
    for m in range(10):
        i = make_immy()[0]
        assert i + 0 == i
        assert 0 + i == i
    
def test_immy_sub(n=1000,prnt=False):
    _test_immy_bop(lambda a,b: a - b,np.subtract,n,prnt)
    
    for m in range(10):
        i = make_immy()[0]
        assert i - 0 == i
        assert 0 - i == -i
    
def test_immy_mul(n=1000,prnt=False):
    _test_immy_bop(lambda a,b: a*b,np.multiply,n,prnt)
    
    for m in range(10):
        i = make_immy()[0]
        assert i*0 == 0
        assert i*immy(0) == 0
        assert i*1 == i
        assert i*immy(1) == i
        assert 1*i == i
        assert immy(1)*i == i
    
def test_immy_div(n=1000,prnt=False):
    _test_immy_bop(lambda a,b: a/b,np.divide,n,prnt)
    
    for m in range(10):
        i = make_immy()[0]
        assert i/1 == i
    
def test_immy_pow(n=1000,prnt=False):
    _test_immy_bop(lambda a,b: a**b,np.power,n,prnt)
    
    for m in range(10):
        i = make_immy()[0]
        assert 0**i == immy(0)
        assert 0**i == ummy(0)
        assert 0**i == 0
        assert i**0 == 1
        
        assert_immy_close(i**-1,1/i)
        assert_immy_close(i**2,i*i)
        assert_immy_close(i**3,i*i*i)
    
