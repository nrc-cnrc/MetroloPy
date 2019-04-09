# -*- coding: utf-8 -*-

import metrolopy as uc
import numpy as np
import warnings

from metrolopy.tests.common import rand,display


def test_gummy_init(n=None,exception_on_warning=True,prnt=False,plot=False):
    print('uc.gummy.p_method may be modified\nset uc.gummy.p_method = \'loc\' to recover the default value')
    from scipy.stats import norm, t
    
    uc.gummy.p_method = None
    
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

