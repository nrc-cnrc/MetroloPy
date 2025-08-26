# -*- coding: utf-8 -*-

# module unitutils

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
the search_units, shadow_units, and convert functions are defined here
"""
from .gummy import gummy
from .ummy import ummy
from .unit import Unit,one,Quantity
from .printing import PrettyPrinter


def _mrtxt(txt,fmt):
    if fmt == 'latex':
        txt = txt.replace('*','\\*')
        txt = txt.replace('_','\\_')
        txt = txt.replace('#','\\#')
    return txt

class search_units_result(PrettyPrinter):
    """
    A `search_units_result` instance emulates a list of units returned 
    from a 'search_units` function call, and pretty-prints the results to 
    the output
    """
    
    def __init__(self,units,show_all):
        self.show_all = show_all
        self._units = units
        self.units = [u[1] for u in units]
        
    def tostring(self,fmt='unicode',**kwds):
        #return search_units(self.search,fmt=fmt,show_all=self.show_all,
                            #units=self.units,converts_to=self.converts_to,
                            #prnt=False)
                            
        from .prefixedunit import PrefixedUnit
        from .nonlinearunit import NonlinearUnit
        from .offsetunit import OffsetUnit
        from .logunit import LogUnit
        from .unit import _CompositeUnit
        
        if fmt == 'latex':
            txt = "<ul style=\"font-family: 'Times New Roman', Times, serif;font-size:1.2em\">\n"
        elif fmt == 'html':
            txt = "<ul>\n"
        else:
            txt = ''
        for u in self._units:
            if fmt in ['latex','html']:
                txt += "<li>"
            try:
                if u[0] != u[1].name:
                    txt += u[0] + ', alias for: '
                    if isinstance(u[1],_CompositeUnit):
                        stxt = u[1].tostring(fmt=fmt)
                        if fmt == 'latex':
                            txt += '$ ' + stxt + ' $'
                        else:
                            txt += stxt
                    else:
                        txt += _mrtxt(u[1].name,fmt)
                else:
                    u = u[1]
                    txt += u.name
                    ttxt = ''
                    if isinstance(u,PrefixedUnit) and not self.show_all:
                        if not self.show_all:
                            if len(u.prefixes) == 1:
                                ttxt = '1 prefix'
                            else:
                                ttxt = str(len(u.prefixes)) + ' prefixes'
                    elif isinstance(u,LogUnit):
                        bse = u.conversion.log_base
                        if isinstance(bse,Quantity):
                            bse = bse.value
                        if isinstance(bse,ummy):
                            bse = bse.x
                        if (round(bse,3) - 2.718) < 0.001:
                            if fmt == 'latex':
                                ttxt += ', log base ' + PrettyPrinter.latex_math('e')
                            if fmt == 'html':
                                ttxt += ', log base <i>e</i>'
                            else:
                                ttxt += ', log base e'
                        else:
                            ttxt += ', log base ' + str(u.conversion.log_base)
                        ttxt += ', multiplier = ' + str(u.conversion.multiplier)
                    elif isinstance(u,NonlinearUnit):
                        ttxt += 'non-linear unit'
                        
                    if ttxt.startswith(', '):
                        ttxt = ttxt[2:]
                    if ttxt != '':
                        txt += ' (' + ttxt + ')'
                    
                    if u.conversion is not None:
                        try:
                            if isinstance(u,OffsetUnit):
                                g = gummy(0,unit=u)
                            elif isinstance(u,LogUnit):
                                g = gummy(u.conversion.offset,unit=u)
                            else:
                                g = gummy(1,unit=u)
                            gc = g.convert(u.conversion.unit)
                            ctxt = g.tostring(fmt=fmt)
                            if fmt == 'html':
                                ctxt += '&nbsp;=&nbsp;' 
                            else:
                                ctxt += ' = '
                            ctxt += gc.tostring(fmt=fmt)
                            if fmt == 'latex':
                                txt += ', $ ' + ctxt + ' $'
                            else:
                                txt += ', ' + ctxt
                        except:
                            raise
                            txt += ', ?? = ??'
                    elif u is not one:
                        if fmt == 'latex':
                            txt += ', symbol: $ ' + u.tostring(fmt=fmt) + ' $'
                        else:
                            txt += ', symbol: ' + u.tostring(fmt=fmt)
                        
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
            
        return txt
    
    def __len__(self):
        return len(self.units)
    
    def __getitem__(self,i):
        return self.units[i]
    
    def __iter__(self):
        return iter(self.units)
    
    def __reversed__(self):
        return reversed(self.units)
    
    def __contains__(self,item):
        return item in self.units
    
def search_units(search=None,fmt=None,show_all=False,units=None,
                        converts_to=None):
    """
    Prints a list of all loaded units or all units that match the search terms.
    
    Parameters
    ----------
    search: `str` or `None`, optional
        A space separated list of search terms to case insentively match.
        If this is omitted or set equal to `None` then a list of all loaded
        units will be printed.  The default is `None`.

    fmt: {'html','latex','unicode','ascii',`None`},optional
        The output format.  If `None`, then the `gummy.printer` value is used.
        If latex output is selected, Markdown is actually used with the unit
        symbols and conversion displayed using inline LaTeX.

    show_all: `bool`, optional
        If `True` units are shown with each prefix listed on a separate line
        (e.g. the millisecond and the microsecond are listed in addition to
        the second) and interval units are shown. If `False` only the base
        unit is shown.  The default is `False`.

    units: `list` of `str` or `Unit`,optional
        A list of units to print.  If this parameter is specified the values
        of the search and `show_all` parameters are ignored.

    Returns
    -------
    A `search_units_result` instance which emulates a list of the returned
    constants and pretty-prints the results to the output or `None` if no
    units are found.
    """
    from importlib import import_module

    from .unit import _CompositeUnit
    
    if fmt is not None:
        fmt = fmt.lower().strip()
        if fmt == 'utf-8':
            fmt = 'unicode'
    
    if units is None:
        while len(Unit._builtins_to_import) > 0:
            import_module(Unit._builtins_to_import.pop(),Unit.__module__)
            
        units = set(Unit._builtin_lib.values()).union(set(Unit._lib.values()))
        
        if converts_to is not None:
            uf = set()
            bunit = Unit.unit(converts_to)
            if isinstance(bunit,_CompositeUnit):
                raise TypeError('converts_to may not be a composite unit')
            bunit = bunit.base
            for u in units:
                if not isinstance(u,_CompositeUnit) and u.base is bunit:
                    uf.add(u)
            units = uf
        elif search is None:
            if len(units) == 0:
                print('no units are loaded')
                return None
        else:
            uf = set()
            for u in units:
                s = set()
                for a in u.aliases:
                    s = s.union(set(a.lower().split()))
                for a in u.shadowed_aliases:
                    s = s.union(set(a.lower().split()))
                if not isinstance(u,_CompositeUnit):
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
                    
            uf = uf.union({u for u in units if u.parent in uf})
                    
            units = uf
            if len(units) == 0:
                print('no units found matching "' + search + '"')
                return None
            
        if not show_all:
            units = [u for u in units if u.parent is None or u.parent not in units]
            
        uf = []
        for u in units:
            if isinstance(u,_CompositeUnit):
                uf += [(a,u) for a in u._aliases]
            else:
                uf.append((u.name,u))
        units = uf
    else:
        if isinstance(units,str) or isinstance(units,Unit):
            units = [units]
        units = [Unit.unit(u) for u in units]
        
        if show_all:
            uf = set(units)
            allunits = set(Unit._builtin_lib.values()).union(set(Unit._lib.values()))
            for u in units:
                if u.parent is None:
                    p = u
                else:
                    p = u.parent
                for a in allunits:
                    if a is p or a.parent is p:
                        uf.add(a)
            units = uf
                    
        units = [(u.name,u) for u in units]

    units = sorted(units,key=lambda u:u[0].lower())
        

    
    return search_units_result(units,show_all)
        
        
def shadowed_units(fmt=None):
    """
    Lists any units which have a shadowed name or alias.  Units may be shadowed
    if the user has defined a new unit with the same name or alias as an
    existing unit.
    
    Parameters
    ---------
    fmt: {'html','latex','unicode','ascii',`None`},optional
        The output format.  If `None`, then the `gummy.printer` value is used.
        If latex output is selected, Markdown is actually used with the unit
        symbols and conversion displayed using inline LaTeX.
    """
    units = set(Unit._builtin_lib.values()).union(set(Unit._lib.values()))
    units = [u for u in units if len(u.shadowed_aliases) > 0]
    return search_units(fmt=fmt,units=units)


def convert(amount,from_unit,to_unit):
    """
    Performs a unit conversion of `amount` in units of `from_unit` to units
    of `to_unit`.
    
    equivalent to ``gummy(amount,from_unit).convert(to_unit)``
    """
    return gummy(amount,unit=from_unit).convert(to_unit)
    
            
