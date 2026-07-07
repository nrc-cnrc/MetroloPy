# -*- coding: utf-8 -*-

import metrolopy as uc
import numpy as np
from fractions import Fraction

try:
    import mpmath as mp
except:
    mp = None

rand = np.random.RandomState()

def display(g):
    w = rand.randint(5)
    
    if isinstance(g,uc.immy):
        if rand.randint(2):
            g.style = 'polar'
            
    assert '?' not in g.tostring()
            
    if w == 0:
        print(g)
    elif w == 1:
        assert '?' not in g.tolatex()
        g.latex()
    elif w == 2:
        assert '?' not in g.tohtml()
        g.html()
    elif w == 3:
        assert '?' not in g.tounicode()
        g.unicode()
    else:
        assert '?' not in g.toascii()
        g.ascii()

def make_gummy(sign=None,exp=None,uexp=-6,sometimes_small=True,dof=None,
               unit=None,units=None,mpf=False,allowzero=True,allowlargeu=True):
    
    utypes = ['A','B','D',None]
    if units is None:
        units = ['m','lb','m**2 s**3/kg**4','m**3','s**-2','1']
        
    if rand.randint(10) == 0:
        if exp is None:
            x = rand.randint(6000)
        else:
            x = rand.randint(int(10**exp))
    else:
        x = rand.rand()*0.9 + 0.1
        
        if exp is None:
            exp = rand.randint(1,6)
            if rand.randint(2):
                exp = -exp
                
        x = x*10**exp
    
    if sign is None:
        if rand.randint(2):
            x = -x
    else:
        x *= sign
        
    if not allowzero and x == 0:
        x = x + 1
        
    if allowlargeu and rand.randint(10) == 0:
        u = rand.randint(20) + 1
    else:
        u = rand.rand()*0.5 + 0.5
        if sometimes_small and not rand.randint(10):
            uexp -= 4
        u = abs(x)*u*10**uexp
    
    if mpf and mp is not None and not rand.randint(4):
        x = mp.mpf(x)
        if rand.randint(2):
            u = mp.mpf(u)
        
    if dof is None:
        if not rand.randint(4):
            dof = float('inf')
        else:
            dof = rand.randint(5,15)
            
    if unit is None:
        unit = units[rand.randint(len(units))]

    utype = utypes[rand.randint(len(utypes))]
    
    g = uc.gummy(x,u=u,dof=dof,unit=unit,utype=utype)
    
    return (g,x,u,dof,unit,utype)

def make_number(sign=None,gconst=True,exp=None,fionly=True,allowzero=True):
    q = rand.randint(4)
    if fionly and q > 1:
        q = 1
    if q == 0:
        if exp is not None:
            x = rand.randint(10**(exp+1))
        else:
            x = rand.randint(100000)
    if q == 1 or q == 2:
        x = rand.rand()*0.9 + 0.1
        if exp is None:
            exp = rand.randint(1,6)
            if rand.randint(2):
                exp = -exp

        x = x*10**exp
    if q == 2 and mp is not None:
        x = mp.mpf(x)
    elif q == 3:
        if exp is not None:
            a = rand.randint(10**(exp+1))
            b = rand.randint(10**(exp+1))
        else:
            a = rand.randint(100000)
            b = rand.randint(100000)
        if b == 0:
            x = Fraction(0)
        else:
            x = Fraction(a,b)
        
    
    if sign is None:
        if rand.randint(2):
            x = -x
    else:
        x *= sign
    
    if gconst and rand.randint(2):
        x = uc.gummy(x)
    
    if not allowzero and x == 0:
        x = x + 1
        
    return x