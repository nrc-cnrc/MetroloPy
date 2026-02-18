# -*- coding: utf-8 -*-

# module test_ubreakdown

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

import metrolopy as uc
import numpy as np
import unittest

from metrolopy.tests.common import rand,make_gummy

class TestUbreakdown(unittest.TestCase):
    def test_ubreakdown(self,n=None,prnt=False,budget=False):
        """
        Test gummy.ubreakdown and the budget class.
        """
        uc.gummy.p_method = None
        
        if n is None:
            if prnt or budget:
                n = 100
            else:
                n = 1000
            
        for i in range(n):
            k = rand.randint(1,6)
            
            if rand.randint(2):
                unt = '1'
            else:
                unt = None
            g,x,u,dof,unit,utype = zip(*[make_gummy(unit=unt) for y in range(k)])
            
            for i,y in enumerate(g):
                if y.u != 0:
                    self.assertTrue(y.utype is utype[i])
            
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
                if abs(gr.correlation(q)) == 1:
                    ut = q.utype
            if ut is not None:
                self.assertTrue(gr.utype is ut)
            
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
            
            if xr == 0:
                self.assertTrue(abs(gr.x - xr) < 1e-10)
            else:
                self.assertTrue(abs(gr.x - xr)/abs(xr) < 1e-10)
            
            ur = 0
            for i,y in enumerate(u):
                ur += y**2*d(i,*x)**2
            ur = np.sqrt(ur)
            if ur == 0:
                self.assertTrue(abs(gr.u) < 1e-10)
            else:
                self.assertTrue(abs(gr.u - ur)/ur < 1e-10)
            
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
            #if dofr > 100:
                #self.assertTrue(gr.dof > 100
            #else:
                #self.assertTrue(abs(gr.dof - dofr)/dofr < 1e-2
                
            if unt == '1':
                if k == 1:
                    dr = lambda *x: d(0,*x)
                else:
                    dr = lambda *x: [d(i,*x) for i in range(k)]
                ga = uc.gummy.apply(f,dr,*g)
                if xr == 0:
                    self.assertTrue(abs(ga.x) < 1e-10)
                else:
                    self.assertTrue(abs(ga.x - xr)/abs(xr) < 1e-10)
                if ur == 0:
                    self.assertTrue(abs(ga.u) < 1e-6)
                else:
                    self.assertTrue(abs(ga.u - ur)/ur < 1e-6)
                if dofr > 100:
                    self.assertTrue(ga.dof > 100)
                else:
                    self.assertTrue(abs(ga.dof - dofr)/dofr < 1e-2)
                    
                gn = uc.gummy.napply(f,*g)
                if xr == 0:
                    self.assertTrue(abs(gn.x) < 1e-10)
                else:
                    self.assertTrue(abs(gn.x - xr)/abs(xr) < 1e-10)
                if ur == 0:
                    self.assertTrue(abs(gn.u) < 1e-3)
                else:
                    self.assertTrue(abs(gn.u - ur)/ur < 1e-3)
                if dofr > 100:
                    self.assertTrue(gn.dof > 100)
                else:
                    self.assertTrue(abs(gn.dof - dofr)/dofr < 1e-2)
                    
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
                
            if ur > 0:
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
                        self.assertTrue(iu*d(i,*x)/ur < 1e-3)
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
                    self.assertTrue(gf.ufrom('A') == 0)
                else:
                    self.assertTrue(abs(gf.ufrom('A') - ua)/ur < 1e-3)
                if ua/ur > 0.01:
                    if dofa > 100 or gf.ufrom('A') == 0:
                        self.assertTrue(gf.doffrom('A') > 100)
                    else:
                        self.assertTrue(abs(gf.doffrom('A') - dofa)/dofa < 1e-2)
                    
                if ub == 0:
                    self.assertTrue(gf.ufrom('B') == 0)
                else:
                    self.assertTrue(abs(gf.ufrom('B') - ub)/ur < 1e-3)
                if ub/ur > 0.01:
                    if dofb > 100 or gf.ufrom('B') == 0:
                        self.assertTrue(gf.doffrom('B') > 100)
                    else:
                        self.assertTrue(abs(gf.doffrom('B') - dofb)/dofb < 1e-2)
                    
                if ud == 0:
                    self.assertTrue(gf.ufrom('D') == 0)
                else:
                    self.assertTrue(abs(gf.ufrom('D') - ud)/ur < 1e-3)
                if ud/ur > 0.01:
                    if dofd > 100 or gf.ufrom('D') == 0:
                        self.assertTrue(gf.doffrom('D') > 100)
                    else:
                        self.assertTrue(abs(gf.doffrom('D') - dofd)/dofd < 1e-2)
                    
                if ue == 0:
                    self.assertTrue(gf.ufrom(e) == 0)
                else:
                    self.assertTrue(abs(gf.ufrom(e) - ue)/ur < 1e-3)
                #if ue/ur > 0.01:
                    #if dofe > 100 or gf.ufrom(e) == 0:
                        #self.assertTrue(gf.doffrom(e) > 100
                    #else:
                        #self.assertTrue(abs(gf.doffrom(e) - dofe)/dofe < 1e-2
                        
            if not rand.randint(4):
                gf.uunit = '%'
                    
            styles = ['pm','pmi','concise','ueq','u','uf']
            gf.style = styles[rand.randint(6)]
                    
            self.assertTrue('?' not in gf.tostring())
            
            ub = {i.utype for i in g if i.utype is not None}
            if not rand.randint(8) and len(ub) > 1:
                ub.pop()
            ub = list(ub)
            
            gfub = gf.copy(formatting=True)
            if len(ub) > 0:
                gfub.ubreakdown = ub
            
            self.assertTrue('?' not in gfub.tostring())
                    
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
    

if __name__ == '__main__':
    unittest.main()