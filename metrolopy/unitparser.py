# -*- coding: utf-8 -*-

# module unitparser

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
unit parser
"""

from .exceptions import UnitLibError
import numpy as np

class _UnitParser:
    def __init__(self,txt):  
        self.txt= txt
        self.n = len(txt)
        self.i = 0
        
    def parse(self):
        return self.ppar(0,1)
                
    def ppar(self,lev,es):
        cu = ''
        ul = []
        sol = False
        while self.i < self.n:
            c = self.txt[self.i]
            if c == ')':
                cu = cu.strip()
                if cu != '':
                    ul += [[cu,es]]
                if len(ul) == 0:
                    raise UnitLibError('syntax error in: ' + self.txt)
                self.i += 1
                e = self.look_for_pow()
                if e != 1:
                    for u in ul:
                        u[1] *= e
                        if np.modf(u[1]) == 0:
                            u[1] = int(u[1])
                return ul
            elif c == '(':
                cu = cu.strip()
                p = self.pch(cu,')',']',self.i)
                if p is None:
                    if cu != '':
                        ul += [[cu,es]]
                        cu = ''
                    self.i += 1
                    ul += self.ppar(lev+1,es)
                    if sol:
                        es = -es
                        sol = False
                else:
                    cu += p
                    self.i += len(p)
            elif c == '[':
                cu = cu.strip()
                p = self.pch(cu,']',')',self.i)
                if p is None:
                    if cu != '':
                        ul += [[cu,es]]
                        cu = ''
                    self.i += 1
                    ul += [self.pbr(es,cu)]
                    cu = ''
                    if sol:
                        es = -es
                        sol = False
                else:
                    cu += p
                    self.i += len(p)
            elif c.isspace() or c == '*':
                cu = cu.strip()
                if cu != '':
                    ul += [[cu,es*self.look_for_pow()]]
                cu = ''
                if sol:
                    es = -es
                    sol = False
            elif c == '/':
                cu = cu.strip()
                if cu != '':
                    ul += [[cu,es]]
                    cu = ''
                 #if sol:
                    # raise UnitLibError('ambugus use of several soliudses (forward slashes) in: ' + self.txt)
                sol = not sol
                es = -es
                self.i += 1
            else:
                cu += c
                self.i += 1
        if lev != 0:
            raise UnitLibError('unmatched ")" found in: ' + self.txt)
        cu = cu.strip()
        if cu != '':
            ul += [[cu,es]]
        return ul
            
    def look_for_pow(self):
        while self.i < self.n:
            c = self.txt[self.i]
            if c == '*':
                if self.i == (self.n - 1):
                    raise UnitLibError(self.txt + ' ends with *')
                self.i += 1
                c = self.txt[self.i]
                if not c == '*':
                    return 1
                par = -1
                et = ''
                self.i += 1
                cl = ''
                stnum = False
                while self.i < self.n:
                    c = self.txt[self.i]
                    if c == '(':
                        if par == 0:
                            if cl.ispace() or cl == '*':
                                f =  float(et)
                                if np.modf(f)[0] == 0:
                                    return int(f)
                                else:
                                    return f
                            raise UnitLibError('syntax error in: ' + self.txt)
                        else:
                            if par == -1:
                                par = 1
                            else:
                                par += 1
                    elif c == '[':
                        if cl.ispace() or cl == '*':
                            f =  float(et)
                            if np.modf(f)[0] == 0:
                                return int(f)
                            else:
                                return f
                        raise UnitLibError('syntax error in: ' + self.txt )
                    elif c == ']':
                        raise UnitLibError('unmatched bracket: ' + self.txt)
                    elif c.isalpha():
                        if par > 0:
                            raise UnitLibError('syntax error in: ' + self.txt)
                        f =  float(et)
                        if np.modf(f)[0] == 0:
                            return int(f)
                        else:
                            return f
                    elif c == '/':
                        if par <= 0:
                            f =  float(et)
                            if np.modf(f)[0] == 0:
                                return int(f)
                            else:
                                return f
                    elif c == ')':
                        if par <= 0:
                            f =  float(et)
                            if np.modf(f)[0] == 0:
                                return int(f)
                            else:
                                return f
                        else:
                            par -= 1
                    elif c.isnumeric():
                        stnum = True
                    elif c == '+'  or c == '-':
                        if stnum and par <= 0:
                            raise UnitLibError('syntax error in: ' + self.txt)
                    elif c == '*':
                        if self.i == self.n - 1:
                            raise UnitLibError('syntax error in: ' + self.txt)
                        if par <= 0 and cl != '*' and self.txt[self.i+1] != '*':
                            self.i += 1
                            f =  float(et)
                            if np.modf(f)[0] == 0:
                                return int(f)
                            else:
                                return f
                    elif c.isspace():
                        if stnum and par <= 0:
                            self.i += 1
                            f =  float(et)
                            if np.modf(f)[0] == 0:
                                return int(f)
                            else:
                                return f
                    else:
                        if par <= 0:
                            raise UnitLibError('syntax error in: ' + self.txt)
                    self.i += 1
                    cl = c
                    et += c
                f = float(et)
                if np.modf(f)[0] == 0:
                    return int(f)
                else:
                    return f
    
            elif not c.isspace():
                return 1
            self.i += 1
        return 1
            
    def pbr(self,es,cu):
        br = 0
        while self.i < self.n:
            c = self.txt[self.i]
            if c == '[':
                br += 1
                cu += c
            if c == ']':
                if br == 0:
                    self.i += 1
                    return [cu,es*self.look_for_pow()]
                br -= 1
            cu += c
            self.i += 1
        raise UnitLibError('syntax error in: ' + self.txt)
        
    def pch(self,cu,ech,och,i):
        if cu == '':
            return None
        ret = self.txt[i]
        i += 1
        while i < self.n:
            c = self.txt[i]
            
            if c == ech:
                return ret + c
            if c.isspace():
                return None
            if c in {och,'*','/'}:
                return None
            
            if c == '(':
                r = self.pch('a',')',']',i)
                if r is None:
                    return None
                else:
                    ret += r
                    i += len(r) - 1
            elif c == '[':
                r = self.pch('a',']',')',i)
                if r is None:
                    return None
                else:
                    ret += r
                    i += len(r) - 1
            else:
                ret += c
                
            i += 1
            
        return None
        
        