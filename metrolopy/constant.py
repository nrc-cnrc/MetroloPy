# -*- coding: utf-8 -*-

# module constant

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
constant
"""

from .indexed import Indexed
from .gummy import gummy,jummy
from .unitutils import _mrtxt
from importlib import import_module
from .printing import PrettyPrinter,print_markdown,print_html,ipython_installed
from .exceptions import ConstantNotFoundError
from .printing import _latex_math

class GummyConstant(gummy,Indexed):
    
    return_toummy = False
    return_splonk = False
    
    _builtins_to_import = ['..codata2018']
    
    _builtin_lib = {}
    _lib = {}
    _open_lib = _lib
    
    @classmethod
    def _raise_not_found(cls,name):
        raise ConstantNotFoundError('constant "' + str(name) + '" was not found')
    
    def __new__(cls,x,u=0,unit=1,dof=float('inf'),k=1,p=None,uunit=None,
                 utype=None,name=None,symbol=None,short_name=None,
                 add_symbol=False, html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None,_key=None):
        if symbol is None or name is None:
            ret = gummy.__new__(gummy,x,u=u,unit=unit,dof=dof,k=k,p=p,
                                uunit=uunit,utype=utype,name=name)
            ret.__init__(x,u=u,unit=unit,dof=dof,k=k,p=p,uunit=uunit,
                         utype=utype,name=name)
            return ret
        
        return super().__new__(cls)
    
    def __init__(self,x,u=0,unit=1,dof=float('inf'),k=1,p=None,uunit=None,
                 utype=None,name=None,symbol=None,short_name=None,
                 add_symbol=False, html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None,_key=None):
        """
        A subclass of `gummy` that can be retrieved by name or alias from the
        constant library with the `constant` function.
        
        `GummyConstant` takes all the parameters that `gummy` takes, though 
        the `name` parameter is now required.  In addition it takes a required
        `symobol` parameter and several other optional parameters:
            
        Parameters
        ----------
        symbol:  `str`, required
            A unicode symbol.  If the `add_symbol` parameter is set to `True`, 
            then this symbol can also be used to access the unit with the 
            `constant` function.
            
        short_name: `str` or `None`
            a short name which can be used as an additional alias for the 
            constant in the constant library
        
        add_symbol: `bool`, optional
            If this is `True`, then the symbol can be used to look up the constant
            in the constant library.  The default is `False`
            
        html_symbol, latex_symbol, ascii_symbol:  `str` or `Mone`, optional
            html, latex, and ascii versions of the symbol if they are different
            from the unicode representation of the symbol.
            
        description:  `str` or `None`, optional
            a description of the unit
        """
        if name is None:
            name = symbol
            
        if html_symbol is None:
            html_symbol = '<i>' + symbol + '</i>'
            
        gummy.__init__(self,x,u=u,unit=unit,dof=dof,k=k,p=p,uunit=uunit,
                      utype=utype,name=name)
        Indexed.__init__(self,name,symbol=symbol,short_name=short_name,
                         add_symbol=add_symbol,html_symbol=html_symbol,
                         latex_symbol=latex_symbol,ascii_symbol=ascii_symbol,
                         description=description)
        
        self._key = _key
    
    def tostring(self,fmt=None,style=None,k=None,p=None,show_k=None,
                 show_p=None,show_dof=None,show_name=None,name=None,
                 norm=None,raw=False,nsig=None,solidus=None,
                 mulsep=None,**kwds):
                     
        if name is None:
            name = Indexed.tostring(self,fmt)
        
        return super().tostring(fmt=fmt,style=style,k=k,p=p,show_k=show_k,
                                show_p=show_p,show_dof=show_dof,
                                show_name=show_name,name=name,norm=norm,
                                raw=raw,nsig=nsig,solidus=solidus,
                                mulsep=mulsep)
    
    @property
    def unit(self):
        return self._unit
    @unit.setter
    def unit(self,u):
        name = self.name
        gummy.unit.fset(self,u)
        self.name = name
    
class JummyConstant(jummy,Indexed):
    _builtins_to_import = ['..codata2018']
    
    _builtin_lib = {}
    _lib = {}
    _open_lib = _lib
    
    @classmethod
    def _raise_not_found(cls,name):
        raise ConstantNotFoundError('constant "' + str(name) + '" was not found')
    
    def __new__(cls,real=None,imag=None,r=None,phi=None,cov=None,unit=1,
                 name=None,symbol=None,short_name=None,add_symbol=False,
                 html_symbol=None,latex_symbol=None,ascii_symbol=None,
                 description=None):
        if symbol is None:
            ret = jummy.__new__(jummy,real=real,imag=imag,r=r,phi=phi,cov=cov,
                                unit=unit,name=name)
            ret.__init__(real=real,imag=imag,r=r,phi=phi,cov=cov,unit=unit,
                         name=name)
            return ret
        
        return super().__new__(cls)
    
    def __init__(self,real=None,imag=None,r=None,phi=None,cov=None,unit=1,
                 name=None,symbol=None,short_name=None,add_symbol=False,
                 html_symbol=None,latex_symbol=None,ascii_symbol=None,
                 description=None):
        if name is None:
            name = symbol
            
        if html_symbol is None:
            html_symbol = '<i>' + symbol + '</i>'
            
        jummy.__init__(self,real=real,imag=imag,r=r,phi=phi,cov=cov,unit=1,
                       name=name)
        Indexed.__init__(self,name,symbol=symbol,short_name=short_name,
                         add_symbol=add_symbol,html_symbol=html_symbol,
                         latex_symbol=latex_symbol,ascii_symbol=ascii_symbol,
                         description=description)
    
    def tostring(self,fmt='unicode',norm=None,nsig=None,solidus=None,
                 mulsep=None,show_name=False,name=None,**kwds):
                     
        if name is None:
            name = Indexed.tostring(self,fmt)
        
        return super().tostring(fmt=fmt,norm=norm,nsig=nsig,solidus=solidus,
                                mulsep=mulsep,show_name=show_name,name=name)

def constant(name,toummy=None,splonk=None):
    """
    Finds an returns a constants from the constant library.
    
    Parameter
    ----------
        name: `str`, `Unit` or 1
            The name or alias of the constant.  
            
    Returns
    -------
    A `GummyConstant` or `JummyConstant` instance 
    
    ummy/immy Quantity instances can bre retrieved by setting `toummy` to True
    or GummyConstant.toummy to True.  splonk the returned value by setting
    `splonk` to True or GummyConstant.splonk to True
    """
    if isinstance(name,GummyConstant) or isinstance(name,JummyConstant):
        return name
    ret = JummyConstant._lib.get(name)
    if ret is None:
        ret = GummyConstant.get(name,exception=False)
    if ret is None:
        ret = JummyConstant.get(name,exception=False)
    if ret is None:
        raise ConstantNotFoundError('constant "' + str(name) + '" was not found')
        
    if splonk is None:
        splonk = GummyConstant.return_splonk
    if splonk:
        ret = ret.splonk()
    else:
        if toummy is None:
            toummy = GummyConstant.return_toummy
        if toummy:
            ret = ret.toummy()
        
    return ret

class _search_display(PrettyPrinter):
    def __init__(self,search,constants):
        self.search = search
        self.constants = constants
        
    def tostring(self,fmt='unicode',**kwds):
        return search_constants(self.search,fmt=fmt,constants=self.constants,
                                prnt=False)

def search_constants(search=None,fmt=None,constants=None,prnt=True):
    """
    Prints a list of all loaded constant or all constants that match the search 
    terms.
    
    Parameters
    ----------
    search: `str` or `None`, optional
        A space separated list of search terms to case insentively match.
        If this is omitted or set equal to `None` then a list of all loaded
        constants will be printed.  The default is `None`.

    fmt: {'html','latex','unicode','ascii',`None`},optional
        The output format.  If `None`, then the `gummy.printer` value is used.
        If latex output is selected, Markdown is actually used with the unit
        symbols and conversion displayed using inline LaTeX.

    constants: `list` of `str`,optional
        A list of constants to print.  If this parameter is specified the values
        of the search and `show_all` parameters are ignored.

    prnt: `bool`, optional
        If this is `True`, the results are printed.  If it is `False` the results
        are returned as a string.  The default is `True`.
    """
    
    if fmt is None and prnt:
        return _search_display(search,constants)

    fmt = fmt.lower().strip()
    if fmt == 'utf-8':
        fmt = 'unicode'
    
    if constants is None:
        while len(GummyConstant._builtins_to_import) > 0:
            import_module(GummyConstant._builtins_to_import.pop(),
                          GummyConstant.__module__)
            
        constants = {id(c):c for c in GummyConstant._builtin_lib.values()}
        constants.update({id(c):c for c in GummyConstant._lib.values()})
        constants.update({id(c):c for c in JummyConstant._builtin_lib.values()})
        constants.update({id(c):c for c in JummyConstant._lib.values()})
        constants = constants.values()
        
        if search is None:
            if len(constants) == 0:
                if prnt:
                    print('no constants are loaded')
                    return
                return ''
        else:
            uf = set()
            for u in constants:
                s = set()
                for a in u.aliases:
                    s = s.union(set(a.lower().split()))
                for a in u.shadowed_aliases:
                    s = s.union(set(a.lower().split()))
                s = s.union(u.name.lower().split())
                s = s.union(set(u.symbol.lower().split()))
                s = s.union(set(u.symbol.lower().split()))
                if u.ascii_symbol is not None:
                    s = s.union(set(u.ascii_symbol.lower().split()))
                if u.description is not None:
                    s = s.union(set(u.description.lower().split()))
                        
                s = {i.strip(',.;:()') for i in s}
                
                srch = search.lower().split()
                ad = True
                for a in srch:
                    if a.strip(',.;') not in s:
                        ad = False
                        break
                if ad:
                    uf.add(u)
                    
            constants = uf
            if len(constants) == 0:
                if prnt:
                    print('no constants found matching "' + search + '"')
                    return
                return ''
            
        uf = []
        for u in constants:
            uf.append((u.name,u))
        constants = uf
    else:
        if (isinstance(constants,str) or isinstance(constants,GummyConstant) 
            or isinstance(constants,JummyConstant)):
            constants  = [constants]
        constants = ([constant(u) for u in constants])
                    
        constants = [(u.name,u) for u in constants]

    constants = sorted(constants,key=lambda u:u[0].lower())
        
    if fmt == 'latex':
        txt = "<ul style=\"font-family: 'Times New Roman', Times, serif;font-size:1.2em\">\n"
    elif fmt == 'html':
        txt = "<ul>\n"
    else:
        txt = ''
    for u in constants:
        if fmt in ['latex','html']:
            txt += "<li>"
        try:
            u = u[1]
            txt += u.name + ' '
            
            if fmt == 'latex':
                txt += _latex_math(u.tostring(fmt=fmt,show_name=True))
            else:
                txt += u.tostring(fmt=fmt,show_name=True)
                
            aliases = u.aliases
            if len(aliases) == 1:
                txt += ', alias: ' + _mrtxt(aliases.pop(),fmt)
            if len(aliases) > 1:
                txt += ', aliases: '
                if u.short_name in aliases:
                    txt += _mrtxt(u.short_name,fmt) + ', '
                    aliases.remove(u.short_name)
                txt += _mrtxt(', '.join(sorted(aliases,key=str.lower)),fmt)
                
            saliases = u.shadowed_aliases
            if len(saliases) > 0:
                if len(aliases) > 0:
                    txt += '; '
                else:
                    txt += ', '
                if len(saliases) == 1:
                    txt += 'shadowed alias: ' + _mrtxt(saliases.pop(),fmt)
                else:
                    txt += 'shadowed aliases: ' + _mrtxt(', '.join(sorted(saliases,key=str.lower)),fmt)
        except:
            raise
            txt += '??'
                
        if fmt == 'html' or fmt == 'latex':
            txt += '</li>\n'
        else:
            txt += '\n'
            
    txt = txt[:-1]
    if fmt in ['latex','html']:
        txt += '</ul>'
    
    if not prnt:
        return txt
    
    if fmt == 'latex' and ipython_installed:
            print_markdown(txt)
    elif fmt == 'html' and ipython_installed:
            print_html(txt)
    else:
        print(txt)
        
        
def shadowed_constants(fmt=None,prnt=True):
    """
    Lists any constants which have a shadowed name or alias.  Constants may be 
    shadowed if the user has defined a new unit with the same name or alias as 
    an existing unit.
    
    Parameters
    ---------
    fmt: {'html','latex','unicode','ascii',`None`},optional
        The output format.  If `None`, then the `gummy.printer` value is used.
        If latex output is selected, Markdown is actually used with the unit
        symbols and conversion displayed using inline LaTeX.

    prnt: `bool`, optional
        If this is `True`, the results are printed.  If it is `False` the results
        are returned as a string.  The default is `True`.
    """
    constants = {id(c):c for c in GummyConstant._builtin_lib.values()}
    constants.update({id(c):c for c in GummyConstant._lib.values()})
    constants.update({id(c):c for c in JummyConstant._builtin_lib.values()})
    constants.update({id(c):c for c in JummyConstant._lib.values()})
    constants = constants.values()
    constants = [u for u in constants if len(u.shadowed_aliases) > 0]
    return search_constants(fmt=fmt,constants=constants,prnt=prnt)
