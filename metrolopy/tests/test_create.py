# -*- coding: utf-8 -*-

import metrolopy as uc
import numpy as np
import warnings
from fractions import Fraction

try:
    import mpmath as mp
except:
    mp = None

rand = np.random.RandomState()

def test_gummy_init(n=None,exception_on_warning=True,prnt=False,plot=False):
    print('uc.gummy.p_method may be modified\nset uc.gummy.p_method = \'loc\' to recover the default value')
    from scipy.stats import norm, t
    
    if n is None:
        if prnt or plot:
            n = 100
        else:
            n=1000000
    
    units=['m','lb','m**2 s**3/kg**4','degF','dB(SPL)','%','Np']
    uunits = {'%':100, 'ppm':1e6, 'ppb':1e9, 'ms/s':1000}
    
    with warnings.catch_warnings():
        if exception_on_warning:
            warnings.simplefilter('error')
        else:
            warnings.simplefilter('ignore')
    
        for i in range(n):
            
            if rand.randint(2):
                unit = units[rand.randint(len(units))]
            else:
                unit = uc.one
                
            if rand.randint(2):
                if unit == 'm':
                    if rand.randint(2):
                        uunit = 'mm'
                    elif rand.randint(2):
                        uunit = 'm'
                    else:
                        uunit = 'in'
                else:
                    if unit in ['degF','dB(SPL)','Np']:
                        uunit = None
                    else:
                        uunit = list(uunits.keys())[rand.randint(len(uunits))]
            else:
                uunit = None
            
            bayesian = bool(rand.randint(2))
            
            if rand.randint(2):
                dof = rand.randint(3,20)
            else:
                dof = float('inf')
                
            if rand.randint(2):
                if rand.randint(2):
                    k = 4*rand.rand() + 0.1
                    p = None
                else:
                    p = rand.rand()/2 + 0.49
                    k = 1
            else:
                k = 1
                p = None
            
            if rand.randint(2):
                pmethods = ['loc','level of confidence','cp','coverage probability','gauss','ccp','ccp','conservative coverage probability','chebyshev']
                uc.gummy.p_method = pmethods[rand.randint(len(pmethods))]
            else:
                uc.gummy.p_method = None
            
            if unit == 'dB(SPL)':
                if rand.randint(10) == 0:
                    x = rand.randint(-10,10)
                else:
                    x = 100.0*(2*rand.rand() - 1)
            elif unit == 'Np':
                if rand.randint(10) == 0:
                    x = rand.randint(-10,10)
                else:
                    x = 20.0*(2*rand.rand() - 1)
            else:
                if rand.randint(10) == 0:
                    x = rand.randint(-60000,60000)
                else:
                    x = (2*rand.rand() - 1)*10.0**rand.randint(-10,10)
            if rand.randint(20) == 0:
                x = 0
                    
            if uunit in ['%','ppm','ppb','ms/s'] and not (unit == '%' and uunit == '%'):
                U = 1.0e3*rand.rand()
                u = abs(x)*U/uunits[uunit]
            elif unit == 'm' and uunit == 'mm':
                if x == 0:
                    U = 1.0e3*rand.rand()
                else:
                    U = abs(x)*1.0e3*rand.rand()
                u = U/1000
            elif unit == 'm' and uunit == 'in':
                if x == 0:
                    U = 100*rand.rand()
                else:
                    U = abs(x)*100*rand.rand()
                u = U*0.0254 
            else:
                if x == 0:
                    U = 1.0e12*rand.rand()
                else:
                    if rand.randint(20) == 0:
                        U = 1e6*abs(x)*rand.rand()
                    elif rand.randint(20) == 0:
                        U = rand.randint(20)+1
                    else:
                        U = 1.1*abs(x)*rand.rand()
                u = U
                    
            if rand.randint(2):
                name = chr(rand.randint(20) + 65) + chr(rand.randint(20) + 65) + chr(rand.randint(20) + 65)
            else:
                name = None
                
            if rand.randint(2):
                utype = chr(rand.randint(20) + 65) + chr(rand.randint(20) + 65) + chr(rand.randint(20) + 65)
            else:
                utype = None
                
            if rand.randint(20) == 0:
                assert uc.gummy(x) == x
                continue
            
            if rand.randint(20) == 0:
                const = True
                U = 0
                u = 0
            else:
                if u == 0:
                    const = True
                else:
                    const = False
                
            uc.gummy.bayesian = bayesian
            if uunit is not None and rand.randint(2):
                uu = uc.gummy(U,unit=uunit)
                g = uc.gummy(x,u=uu,unit=unit,dof=dof,k=k,p=p,name=name,
                             utype=utype)
            else:
                g = uc.gummy(x,u=U,unit=unit,dof=dof,k=k,p=p,uunit=uunit,
                             name=name,utype=utype)
                
            if not const:
                u = u/g.k
                
            orig_g = g

            if rand.randint(10) == 0:
                n = rand.randint(12)
                if bayesian and n > 0:
                    dof = float('inf')
                if n == 0:
                    g = uc.gummy(g)
                    assert g.name is None
                elif n == 1:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        g = 1*g
                elif n == 2:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        g = g*1
                elif n == 3:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        g = g/1
                elif  n == 4:
                    if unit == 'degF':
                        v = uc.gummy(0,unit='degF-i')
                    elif unit is not uc.one or rand.randint(2):
                        v = uc.gummy(0,unit=unit)
                    else:
                        v = 0
                    g = g + v
                elif n == 5:
                    if unit == 'degF':
                        v = uc.gummy(0,unit='degF-i')
                    elif unit is not uc.one or rand.randint(2):
                        v = uc.gummy(0,unit=unit)
                    else:
                        v = 0
                    g = v + g
                elif n == 6:
                    if unit == 'degF':
                        v = uc.gummy(0,unit='degF-i')
                    elif unit is not uc.one or rand.randint(2):
                        v = uc.gummy(0,unit=unit)
                    else:
                        v = 0
                    g = g - v
                elif n == 7:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    elif unit in ['dB(SPL)','Np']:
                        g =  g + uc.gummy(0,unit=unit)
                    else:
                        g = g**1
                elif n == 8:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        g = +g
                elif n == 9:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        g = 2.1*g
                        g = g/2.1
                elif n == 10:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        if unit is not uc.one or rand.randint(2):
                            v = uc.gummy(3.5*x,u=0.45*u,unit=unit)
                        else:
                            v = 3.3*x
                        g = v - g
                        g = g - v
                        g = -g
                else:
                    if unit == 'degF':
                        g =  g + uc.gummy(0,unit='degF-i')
                    else:
                        if x >= 0:
                            g = abs(g)
                        else:
                            g = -abs(g)
    
                k = 1
                p = None
                reinit = True
                uunit = None
                name = None
            else:
                reinit = False
    
            if x != 0 and g.x != 0:
                assert abs(g.x - x)/abs(x) < 1e-10
            else:
                assert abs(g.x - x) < 1e-10
                
            if x != 0 and g._x != 0:
                assert abs(g._x - x)/abs(x) < 1e-10
            else:
                assert abs(g._x - x) < 1e-10
                            
            if unit is 'degF':
                assert g.unit in [uc.Unit.unit('degF'),uc.Unit.unit('degF-i')]
            else:
                assert g.unit is uc.Unit.unit(unit)
            assert g.name == name
            
            if const:
                assert g.const
                assert g._u == 0
                if unit is uc.one:
                    if x == 0:
                        assert abs(g - x) < 1e-10
                    else:
                        assert abs(g - x)/x < 1e-10
                elif unit == '%':
                    if x == 0:
                         assert abs(g - x/100) < 1e-10
                    else:
                        assert abs((g - x/100)/(x/100)) < 1e-10
                elif unit == 'Np':
                    assert abs(g.convert(1) - np.exp(x))/np.exp(x) < 1e-10
                continue
            
            if uunit in ['%','ppm','ppb','dB','ms/s'] and unit != uunit:
                assert g.uunit_is_rel
            else:
                assert not g.uunit_is_rel
                
            assert (g._x - x)/u < 1e-6
            
            assert not g.const
            
            if dof == float('inf'):
                assert g.dof == float('inf')
            else:
                assert (g.dof - dof)/dof < 1e-6
            
            if uunit is not None and uunit != unit:
                assert g.uunit is uc.Unit.unit(uunit)
                
            if p is None:
                assert g.k == k
            else:
                assert g.p == p
            
            if uc.gummy.p_method in ['loc','level of confidence',None]:
                if dof == float('inf'):
                    assert abs(g.k - norm.ppf(0.5 + g.p/2)) < 1e-6
                else:
                    if bayesian:
                        if g.p is not None:
                            assert abs(g.k - np.sqrt((dof-2)/dof)*t.ppf(0.5 + g.p/2,dof)) < 1e-6
                    else:
                        assert abs(g.k - t.ppf(0.5 + g.p/2,dof)) < 1e-6
                        
            elif uc.gummy.p_method in ['cp','coverage probability','gauss']:
                if dof == float('inf') or bayesian:
                    if g.p is not None:
                        assert abs(g.k - 2/(3*(np.sqrt(1 - g.p)))) < 1e-6
                else:
                    if g.p is not None:
                        assert abs(g.k - np.sqrt(dof/(dof-2))*(2/(3*(np.sqrt(1 - g.p))))) < 1e-6                 
            else:
                if dof == float('inf') or bayesian:
                    if g.p is not None:
                        assert abs(g.k - 1/np.sqrt(1 - g.p)) < 10.0**-6
                else:
                    if g.p is not None:
                        assert abs(g.k - np.sqrt(dof/(dof-2))*(1/np.sqrt(1 - g.p))) < 1e-6
                    
            assert abs((g._u - u)/u) < 1e-10
            
            if not reinit:
                assert abs(g.U - U)/U < 1e-10
            
            assert abs((g._u - u)/u) < 1e-10
            
            assert '?' not in g.tostring()
            
            if prnt:
                print('------')
                print()
                if rand.randint(2):
                    styles = ['pm','pmi','concise','ueq','x','xf','u','uf','xunit','uunit']
                    g.style = styles[rand.randint(10)]
                    
                    g.show_k = bool(rand.randint(2))
                    g.show_dof = bool(rand.randint(2))
                    g.show_p = bool(rand.randint(2))
                    g.show_name = bool(rand.randint(2))
                    g.solidus = bool(rand.randint(2))
                    g.mulsep = bool(rand.randint(2))
                    if not rand.randint(5):
                        g.nsig = rand.randint(5) + 1
                    
                    assert '?' not in g.tostring()
                    print(g.style,', show_k=',g.show_k,', show_dof=',g.show_dof,', show_p=',g.show_p,', show_name=',g.show_name,', solidus=',g.solidus,', mulsep=',g.mulsep,', nsig=',g.nsig,':')
                    print()
                display(g)
                print()
                    
            if plot:
                print('------')
                print()
                print(g)
                g.sim()
                if rand.randint(2):
                    styles = ['pmsim','pmsimi','mcisim','cisim','usim','ufsim']
                    g.style = styles[rand.randint(6)]
                g.slashaxis = bool(rand.randint(2))
                print(g.style,', slashaxis=',g.slashaxis,':')
                print()
                display(g)
                print()
                g.hist()
                print()
                
            uc.gummy.bayesian = False
            
def display(g):
    w = rand.randint(5)
    assert '?' not in g.tostring()
    if w == 0:
        print(g)
    elif w == 1:
        g.latex()
    elif w == 2:
        g.html()
    elif w == 3:
        g.unicode()
    else:
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
        return (-bx,-ax)
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
        
        assert abs((a.x - b.x)/b.x) < 1e-14
        assert (a.u - b.u)/b.u < 1e-3
        assert a.correlation(b) == 1
        assert a.correlation(x) == 1
        
        assert abs((c.x - b.x)/b.x) < 1e-14
        assert (c.u - b.u)/b.u < 1e-3
        assert c.correlation(b) == 1
        
        assert abs((d.x - b.x)/b.x) < 1e-14
        assert (d.u - b.u)/b.u < 1e-3
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
        
        assert abs((a.x - b.x)/b.x) < 1e-14
        assert abs((a.u - b.u)/b.u) < 1e-3
        assert a.correlation(b) == 1
        
        if k%2 == 0 and x.x < 0:
            assert a.correlation(x) == -1
        else:
            assert a.correlation(x) == 1
            
        assert abs((c.x - b.x)/b.x) < 1e-14
        assert abs((c.u - b.u)/b.u) < 1e-3
        assert c.correlation(b) == 1
        
        assert abs((d.x - b.x)/b.x) < 1e-14
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
        
        assert abs((a.x - b.x)/b.x) < 1e-14
        assert abs((a.u - b.u)/b.u) < 1e-3
        assert a.correlation(b) == 1
        
        if k%2 == 0 and x.x < 0:
            assert a.correlation(x) == 1
        else:
            assert a.correlation(x) == -1
        
        assert abs((d.x - b.x)/b.x) < 1e-14
        assert abs((d.u - b.u)/b.u) < 1e-3
        assert d.correlation(b) == 1

def test_ubreakdown(n=None,prnt=False,budget=False):
    
    if n is None:
        if prnt:
            n = 100
        else:
            n = 10000
        
    for i in range(n):
        k = rand.randint(1,6)
        
        if rand.randint(2):
            unt = '1'
        else:
            unt = None
        g,x,u,dof,unit,utype = zip(*[make_gummy(unit=unt) for y in range(k)])
        
        for i,y in enumerate(g):
            assert y.utype is utype[i]
        
        c = rand.rand()*0.9 + 0.1
        if rand.randint(2):
            c = -c
        
        def f(*p):
            r = c
            for y in p:
                r *= y
            return r
        
        def d(k,*p):
            r = c
            for i,y in enumerate(p):
                if i != k:
                    r *= y
            return r
        
        gr = f(*g)
        
        ut = None
        for q in g:
            if gr._ref is q._ref:
                ut = q.utype
        assert gr.utype is ut
        
        x = list(x)
        u = list(u)
        lb = 'lb' in str(gr.unit)
        kg = 'kg' in str(gr.unit)
        for i,un in enumerate(unit):
            if 'lb' in un and not lb:
                x[i] *= 0.45359237
                u[i] *= 0.45359237
            if 'kg' in un and not kg and lb:
                x[i] *= 0.45359237**4
                u[i] *= 0.45359237**4
                
        xr = f(*x)
        
        assert abs(gr.x - xr)/abs(xr) < 1e-10
        
        ur = 0
        for i,y in enumerate(u):
            ur += y**2*d(i,*x)**2
        ur = np.sqrt(ur)
        assert abs(gr.u - ur)/ur < 1e-10
        
        dofr = 0
        for i,(v,idof) in enumerate(zip(u,dof)):
            dofr += v**4*d(i,*x)**4/idof
        if dofr == 0:
            dofr = float('inf')
        else:
            dofr = ur**4/dofr
        if dofr > uc.gummy.max_dof:
            dofr = float('inf')
            
        # Due to the fact that any dof larger than gummy.mox_dof is rounded to 
        # infinity in intermediate calculations, gr.dof may not be exactly dofr
        # particularly for large dofr.  
        if dofr > 100:
            assert gr.dof > 100
        else:
            assert abs(gr.dof - dofr)/dofr < 1e-2
            
        if unt == '1':
            if k == 1:
                dr = lambda *x: d(0,*x)
            else:
                dr = lambda *x: [d(i,*x) for i in range(k)]
            ga = uc.gummy.apply(f,dr,*g)
            assert abs(ga.x - xr)/abs(xr) < 1e-10
            assert abs(ga.u - ur)/ur < 1e-6
            if dofr > 100:
                assert ga.dof > 100
            else:
                assert abs(ga.dof - dofr)/dofr < 1e-2
                
            gn = uc.gummy.napply(f,*g)
            assert abs(gn.x - xr)/abs(xr) < 1e-10
            assert abs(gn.u - ur)/ur < 1e3
            if dofr > 100:
                assert gn.dof > 100
            else:
                assert abs(gn.dof - dofr)/dofr < 1e-2
                
            if rand.randint(2):
                if rand.randint(2):
                    gf = ga
                else:
                    gf = gn
            else:
                gf = gr
        else:
            gf = gr
            
        gl = list(g)
        if not rand.randint(4):
            gl.append(make_gummy(unit=unt)[0])
            
        if rand.randint(2):
            e = gl[rand.randint(len(gl))]
        else:
            e = [gl[rand.randint(len(gl))] for i in range(rand.randint(len(gl)))]
            
        ua = 0
        ub = 0
        ud = 0
        ue = 0
        dofa = 0
        dofb = 0
        dofd = 0
        dofe = 0
        for i,(ig,iu,idof) in enumerate(zip(g,u,dof)):
            if gf.ufrom(ig) == 0:
                assert iu*d(i,*x)/ur < 1e-3
            else:
                if ig.utype == 'A':
                    ua += iu**2*d(i,*x)**2
                    dofa += iu**4*d(i,*x)**4/idof
                elif ig.utype == 'B':
                    ub += iu**2*d(i,*x)**2
                    dofb += iu**4*d(i,*x)**4/idof
                elif ig.utype == 'D':
                    ud += iu**2*d(i,*x)**2
                    dofd += iu**4*d(i,*x)**4/idof
                
                try:
                    if e is ig or ig in e:
                        ue += iu**2*d(i,*x)**2
                        dofe += iu**4*d(i,*x)**4/idof
                except TypeError:
                    pass
            
        ua = np.sqrt(ua)
        ub = np.sqrt(ub)
        ud = np.sqrt(ud)
        ue = np.sqrt(ue)
        
        if dofa == 0:
            dofa = float('inf')
        else:
            dofa  = ua**4/dofa
        if dofb == 0:
            dofb = float('inf')
        else:
            dofb  = ub**4/dofb
        if dofd == 0:
            dofd = float('inf')
        else:
            dofd  = ud**4/dofd
        if dofe == 0:
            dofe = float('inf')
        else:
            dofe  = ue**4/dofe
                    
        if ua == 0:
            assert gf.ufrom('A') == 0
        else:
            assert abs(gf.ufrom('A') - ua)/ur < 1e-3
        if ua/ur > 0.01:
            if dofa > 100 or gf.ufrom('A') == 0:
                assert gf.doffrom('A') > 100
            else:
                assert abs(gf.doffrom('A') - dofa)/dofa < 1e-2
            
        if ub == 0:
            assert gf.ufrom('B') == 0
        else:
            assert abs(gf.ufrom('B') - ub)/ur < 1e-3
        if ub/ur > 0.01:
            if dofb > 100 or gf.ufrom('B') == 0:
                assert gf.doffrom('B') > 100
            else:
                assert abs(gf.doffrom('B') - dofb)/dofb < 1e-2
            
        if ud == 0:
            assert gf.ufrom('D') == 0
        else:
            assert abs(gf.ufrom('D') - ud)/ur < 1e-3
        if ud/ur > 0.01:
            if dofd > 100 or gf.ufrom('D') == 0:
                assert gf.doffrom('D') > 100
            else:
                assert abs(gf.doffrom('D') - dofd)/dofd < 1e-2
            
        if ue == 0:
            assert gf.ufrom(e) == 0
        else:
            assert abs(gf.ufrom(e) - ue)/ur < 1e-3
        #if ue/ur > 0.01:
            #if dofe > 100 or gf.ufrom(e) == 0:
                #assert gf.doffrom(e) > 100
            #else:
                #assert abs(gf.doffrom(e) - dofe)/dofe < 1e-2
                
        if not rand.randint(4):
            gf.uunit = '%'
                
        styles = ['pm','pmi','concise','ueq','u','uf']
        gf.style = styles[rand.randint(6)]
                
        assert '?' not in gf.tostring()
        
        ub = {i.utype for i in g if i.utype is not None}
        if not rand.randint(8) and len(ub) > 1:
            ub.pop()
        ub = list(ub)
        
        gfub = gf.copy(formatting=True)
        if len(ub) > 0:
            gfub.ubreakdown = ub
        
        assert '?' not in gfub.tostring()
                
        if prnt:
            w = rand.randint(5)
            if w == 0:
                print(gf)
                print(gfub)
            elif w == 1:
                gf.latex()
                gfub.latex()
            elif w == 2:
                gf.html()
                gfub.html()
            elif w == 3:
                gf.unicode()
                gfub.unicode()
            else:
                gf.ascii()
                gfub.ascii()
        
        if budget:
            show_expanded_u = bool(rand.randint(2))
            show_subtotals = bool(rand.randint(2))
            
            k = None
            p = None
            if bool(rand.randint(2)):
                if bool(rand.randint(2)):
                    if bool(rand.randint(2)):
                        gf.k = 2
                    else:
                        k = 2
                else:
                    if bool(rand.randint(2)):
                        gf.p = 0.95
                    else:
                        p = 0.95
                        
            if bool(rand.randint(2)):
                gl[0].name = 'gl0'
                
            if bool(rand.randint(2)) and len(gl) > 1:
                gl[1].name = 'gl1'
                
            if bool(rand.randint(2)):
                gf.name = 'gf'
            
            if not bool(rand.randint(4)):
                description = ['description ' + str(i) for i in range(len(gl)+1)]
            else:
                description = None
                
            if not bool(rand.randint(4)):
                custom = [str(i) for i in range(len(gl)+1)]
                custom_heading = 'n'
            else:
                custom = None
                custom_heading = None
                
            show_s = None
            show_c = None
            show_d = None
            if bool(rand.randint(2)):
                show_s = bool(rand.randint(2))
                show_c = bool(rand.randint(2))
                show_d = bool(rand.randint(2))
                
            b = gf.budget(gl,show_subtotals=show_subtotals,
                          show_expanded_u=show_expanded_u,k=k,p=p,
                          description=description,show_s=show_s,show_c=show_c,
                          show_d=show_d,custom=custom,
                          custom_heading=custom_heading)
            
            w = rand.randint(5)
            if w == 0:
                print(b)
            elif w == 1:
                b.latex()
            elif w == 2:
                b.html()
            else:
                b.unicode()
        
        if prnt or budget:
            print()
            print('---')
            print()