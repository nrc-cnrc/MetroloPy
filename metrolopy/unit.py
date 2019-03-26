# -*- coding: utf-8 -*-

# module Unit

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
This module contains the procedures for defining units and unit algebera used
by the gummy class.
"""

import numpy as np
import weakref
import warnings
from importlib import import_module
from .exceptions import (UnitLibError,CircularUnitConversionError,
                         NoUnitConversionFoundError,UnitNotFoundError,
                         IncompatibleUnitsError,UnitLibNotFoundError,
                         UnitWarning)
from .printing import PrettyPrinter
from .nummy import nummy

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
        
        
class Conversion:
    """
    Represents a unit conversion.  This class should only be used as arguments
    to `Unit` object initializers.  Each conversion should be associated with one
    and only one parent Unit.
    
    Parameters
    ----------
    unit:  `str` or `Unit`
        the Unit that the parent Unit will be converted to.
        
    factor:  `float`, optional
        The conversion factor between the parent Unit and the new Unit:
        [value with new Unit] = factor * [value with parent Unit]
        The default value is 1
    """
    linear = True
    
    def __init__(self,unit,factor=1):
        self._unit = unit
        self.factor = factor
        
    def _to(self,g):
        return g*self.factor
        
    def to(self,g):
        unit_in = None
        try:
            try:
                unit_in = g._unit
                g._unit = one
            except:
                pass
            ret = self._to(g)
            try:
                if unit_in is not None:
                    ret._unit = self.unit
            except:
                pass
            return ret
        finally:
            if unit_in is not None:
                try:
                    g._unit = unit_in
                except:
                    pass
    
    def _frm(self,g):
        return g/self.factor
        
    def frm(self,g):
        unit_in = None
        try:
            try:
                unit_in = g._unit
                g._unit = one
            except:
                pass
            ret = self._frm(g)
            try:
                if unit_in is not None:
                    ret._unit = self.parent
            except:
                pass
            return ret
        finally:
            if unit_in is not None:
                try:
                    g._unit = unit_in
                except:
                    pass
        
    def chain(self,c):
        if c.linear:
            ret = Conversion(c.unit,self.factor*c.factor)
            ret.parent = self.parent
            return ret
        else:
            return c.rchain(self)
        
    @property
    def unit(self):
        if not isinstance(self._unit,Unit):
            self._unit = Unit.unit(self._unit)
        return self._unit
        
    def copy(self):
        r = Conversion(self.unit,self.factor)
        r.parent = self.parent
        return r
        
    def pow(self,e):
        r = Conversion(self.unit**e,self.factor**e)
        r.parent = self.parent**e
        return r


class Unit(PrettyPrinter):
    """
    Creating an instance of this class creates a representation of a physical
    unit and adds it to the unit library.  Units already in the unit library or 
    derived units made up of other units in the unit library can be accessed 
    by passing a text string with the unit name or symbol to the static method
    `Unit.unit`.
    
    Parameters
    ----------
        
    name:  `str`
        The name of the unit.  The name can be used access the unit with the
        `Unit.unit` method, but note that if you define a Unit with an
        identical name to a previously defined unit then the older name will
        be shadowed.
    
    symbol:  `str`
        A unicode symbol used when displaying the unit.  If the `add_symbol`
        parameter is set to `True`, then this symbol can also be used
        to access the unit with the `Unit.unit` method.
    
    conversion: `Conversion` or `None`, optional
        A conversion to another unit.  When creating units be careful to avoid
        circular conversions, i.e. you can define:
            
        Unit('inch','in',conversion=None)
        Unit('foot','ft',conversion=Conversion('in',12))
        Unit('yard','yd',conversion=Conversion('ft',3))
        
        but not:
            
        Unit('inch','in',conversion=Conversion('yd',1/36))
        Unit('foot','ft',conversion=Conversion('in',12))
        Unit('yard','yd',conversion=Conversion('ft',3))
        
        Note that an equivent and allowable way of defining the first set of 
        units above is:
            
        Unit('inch','in',conversion=None)
        Unit('foot','ft',conversion=Conversion('in',12))
        Unit('yard','yd',conversion=Conversion('in',36))
        
        Either way will allow the free conversion between inches, feet, and
        yards.  You can also define the inch as:
            
        Unit('inch','in',conversion=('cm',2.54)
        
        (The inch is actually defined this way in the builtin unit library.)
    
    short_name: `str` or `None`
        a short name
    
    add_symbol: `bool`, optional
        If this is `True`, then the symbol can be used to look up the unit
        in the unit library.  The default is `False`
        
    html_symbol, latex_symbol, ascii_symbol:  `str` or `Mone`, optional
        html, latex, and ascii versions of the symbol if they are different
        from the unicode representation of the symbol.
        
    description:  `str` or `None`, optional
        a description of the unit
    
    order:  `int` , optional
        When displaying composite derived units, the symbols with a lower
        `order` value will be displayed.  The default if -1.
        
        
    See Also
    --------
    The following are sub-classes of `Unit`:

    PrefixedUnit:  Creates a set of units with SI prefixes (..., kilo, mega, giga, ...)
    
    BinaryPrefixedUnit:  Creates a set of unit with binary prefixes (..., kibi, mebi, gibi, ...)
    
    LogUnit:  Logrithmic units (e.g. decibel or neper)
    
    OffsetUnit:  Units with an offset origin (degree Celsius or degree Fahrenheit)
    """
    # to speed loading of the module, only wait to load units until they are
    # needed.  If a unit cannot be found the modules below are loaded in 
    # reverse order.
    _builtins_to_import = ['..usunits','..siunits','..relunits']
    
    # _lib and _builtin_lib are dictionarys with string keys that reference
    # Unit object instances.  _lib is searched first so definitions in _lib
    # may shadow definitions in _builtin_lib.
    
    _lib = {} # Contains Unit instances created or loaded by the user.
    _builtin_lib = {} # Contains Unit instances created in the gummy module.
    
    _open_lib = _lib
    _used_units = {}
        
    @staticmethod
    def unit(txt,exception=True):
        """
        Finds an returns a Unit object from the Unit library.
        
        Parameters:
            txt: `str`, `Unit` or 1
                This may be a string representing the unit.  The string can
                contain the name, short name or (if the unit was created
                with `add_symbol` set to `True`) the symbol of the unit or a
                combination of names and/or symbols of several different
                units.  Spaces or the character '*' represent multiplication,
                the character '/' represents division and the string '**'
                represents the power operator. For example txt can be:
                    
                'kg m**2/s'
                
                or equivalently:
                    
                'kilogram*metre*metre*second**-1' or '(kg/s)*m**2'.
                
                If a unit name contains a space, '*' or '/' character then the
                name must be enclosed in square brackets, e.g:
                    
                [light year]
                
                If txt is a Unit instance that instance is returned.
                
            exception: `bool`, optional
                If this is `True` then a `UnitNotFoundError` or UnitLibError
                is raised if a unit is not found that matches `txt`.  If it is
                `False` and a unit is not found, then `Unit.unit` returns
                `None` without raising an exception.  The default is `True`.
        """
        try:
            if isinstance(txt,Unit):
                return txt
            if txt is None or txt is 1:
                return one
            try:
                txt = txt.strip()
            except AttributeError:
                raise UnitNotFoundError('a unit may only be represented by a Unit object, the integer "1", or a string')
            if txt == '' or txt == '1':
                return one
            uu = Unit._used_units.get(txt)
            if uu is not None:
                return uu 
            ul = _UnitParser(txt).parse()
            if len(ul) == 1 and ul[0][1] == 1:
                ret = Unit._single_unit(ul[0][0])
                return ret
            ull = []
            for u in ul:
                un = Unit._single_unit(u[0])
                if hasattr(un,'_getme'):
                    un = un._getme(ul,u[1])
                if isinstance(un,_CompositeUnit):
                    ull += [(uu[0],uu[1]*u[1]) for uu in un._units]
                else:
                    ull.append((un,u[1]))
            ret = _CompositeUnit(ull)
            Unit._used_units[txt] = ret
            return ret
        except UnitNotFoundError or UnitLibError:
            if (not txt.startswith('[') and not txt.endswith(']')):
                return Unit.unit('[' + txt + ']')
            elif exception:
                raise
            return None
        
    @staticmethod
    def reorder(txt):
        """
        This changes the order in which the symbols of composite derived units
        are displayed.

        Examples
        --------
        
        >>> print(Unit.unit('ft lb'))
        ft lb
        >>> print(Unit.unit('lb ft'))  #This is the same unit as above and displays identically
        ft lb
        >>> Unit.reorder('lb ft')  #Now the order will be changed when the unit is displayed
        >>> print(Unit.unit('ft lb'))
        lb ft
            
        """
        unit = Unit.unit(txt)
        if not isinstance(unit,_CompositeUnit):
            return
        
        ul = _UnitParser(txt).parse()
        for u in ul:
            un = Unit._single_unit(u[0])
            if hasattr(un,'_getme'):
                u[0] = un._getme(ul,u[1])
            else:
                u[0] = un
            if isinstance(un,_CompositeUnit):
                raise ValueError('aliases for composite units may not be used in Unit.reorder')
        unit._make_symbols(ul)
                
    @staticmethod
    def _single_unit(txt):
        if txt == '1' or txt == '' or txt is 1:
            return one

        un = Unit._lib.get(txt)
        if un is None:
            un = Unit._builtin_lib.get(txt)
            if un is None:
                while len(Unit._builtins_to_import) > 0:
                    import_module(Unit._builtins_to_import.pop(),Unit.__module__)
                    un = Unit._builtin_lib.get(txt)
                    if un is not None:
                        break    

        if un is None:
            raise UnitNotFoundError('unit "' + txt + '" was not found')
        return un
        
    @staticmethod
    def alias(alias,unit):
        """
        Creates an alias that can be used to reference a Unit.
        
        Parameters
        ----------
        alias:  `str`
            a string containing the new alias

        unit:  `str` or `Unit`
            A string referencing the `Unit` that will be assigned
            the alias or the `Unit` instance its self.
        """
        unit = Unit.unit(unit) # I like this line.
        Unit._open_lib[alias] = unit
        unit._add_alias(alias)
        
    @staticmethod
    def load(library_name):
        import_module(library_name)
        
    @staticmethod
    def format_latex(text):
        if text == '':
            return ''
        text = text.replace(' ','\\,')
        if text.startswith('\t\t'):
            return text[2:]
        if text.startswith('\\'):
            return text
        if text.startswith('\t'):
            return '\t\\mathrm{' + text[1:] + '}'
        return '\\mathrm{' + text + '}'
              
    _gconv = None
    _conversion = None
    linear = True
    
    def __init__(self,name,symbol,conversion=None,short_name = None,
                 add_symbol=False,html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None,order = -1,**kwds):
        if conversion is not None and not isinstance(conversion,Conversion):
            raise ValueError('conversion is not an instance of the Conversion class')
        self.name = name
        self.conversion = conversion
        self.description=description
        self.order = order
        
        if symbol.startswith('\t\t'):
            raise ValueError('the symbol may not start with more than one tab character; only the latex_symbol may start with double tabs')
        self.symbol = symbol
        if html_symbol is None:
            self.html_symbol = symbol
        else:
            self.html_symbol=html_symbol
        if latex_symbol is None:
            self.latex_symbol = Unit.format_latex(symbol)
        else:
            self.latex_symbol = Unit.format_latex(latex_symbol)
        if ascii_symbol is None:
            self.ascii_symbol=symbol
        else:
            self.ascii_symbol=ascii_symbol
            
        self._aliases = set()
        
        self.parent = kwds.get('parent')
        
        #reset _used_units incase we are shadowing any unit definitions
        Unit._used_units = {} 
        
        Unit._open_lib[name] = self
        if add_symbol:
            if symbol.strip() != name:
                Unit._open_lib[symbol.strip()] = self
                self._aliases.add(symbol.strip())
            if (ascii_symbol is not None and ascii_symbol.strip() != name and 
                   ascii_symbol.strip() != symbol.strip()):
                Unit._open_lib[ascii_symbol.strip()] = self
                self._aliases.add(ascii_symbol.strip())
        if short_name is not None:
            Unit._open_lib[short_name] = self
            self.short_name = short_name
            self._aliases.add(short_name)
        else:
            if add_symbol:
                if ascii_symbol is not None:
                    self.short_name = ascii_symbol
                else:
                    self.short_name = symbol
            else:
                self.short_name = name

    @property
    def conversion(self):
        """
        Gets or sets the `Conversion` instance for the unit.
        
        ***This property is not intended to be used directly and setting this 
        property may cause problems***
        """
        return self._conversion
    @conversion.setter
    def conversion(self,c):
        self._conversion = c
        if c is not None:
            c.parent = self
            
    @property
    def aliases(self):
        """read-only
        
        Returns a set of the unshadowed aliases of this unit.  To add aliases
        to the unit use the `Unit.alias` static method.
        """
        return {a for a in self._aliases if (Unit.unit('[' + a + ']',exception=False) is self and a != self.name)}
    
    @property
    def shadowed_aliases(self):
        """read-only
        
        Returns a set of the shadowed aliases of this unit.
        """
        sha = {a for a in self._aliases if Unit.unit('[' + a + ']',exception=False) is not self}
        if Unit.unit('[' + self.name + ']') is not self:
            sha.add(self.name)
        return sha
    
    def _add_alias(self,alias):
        self._aliases.add(alias)
    
    def tostring(self,fmt=None,**kwds):
        """
        Returns a string containing the symbol for the unit the format given by
        the keyword `fmt` which may be set to a string the values 'html', 'latex',
        'ascii', or 'unicode'.  
        """
        if fmt is None or fmt == 'unicode':
            return self.symbol
        if fmt is None:
            return self.symbol
        if fmt == 'html':
            return self.html_symbol
        if fmt == 'latex':
            return self.latex_symbol
        if fmt == 'ascii':
            return self.ascii_symbol
        raise ValueError('format ' + str(fmt) + ' is not recognized')

    def _convs(self,d,c,h):          
        if c is None:
            c = self._conversion.copy()
        else:
            if self in h:
                if self is d:
                    raise CircularUnitConversionError('a circular conversion chain was found starting with unit ' + d.name + ' in library ' + d.library)
                else:
                    raise CircularUnitConversionError('a circular conversion chain was found;\nstarted search with unit ' + d.name + ' in library ' + d.library + ' and the loop started with unit ' + self.name + ' in library ' + self.library)
            h.add(self)
            d._gconv = c
            if self._conversion is None:
                return
            c = c.chain(self._conversion)
            
        try:
             u = self._conversion.unit
        except UnitNotFoundError:
            warnings.warn('while searching conversions:  unit ' + self.name + ' should have a conversion to ' + self._conversion.unit_name + ' in library ' + str(self._conversion.unit_library) + ' but the second unit was not found',UnitWarning,stacklevel=3)
            return
        except UnitLibNotFoundError:
            warnings.warn('while searching conversions:  unit ' + self.name + ' should have a conversion to ' + self._conversion.unit_name + ' in library ' + str(self._conversion.unit_library) + ' but the library containing the second unit has not been loaded',UnitWarning,stacklevel=3)
            return
        
        u._convs(d,c,h)
                
    def _convtree(self):
        if self._conversion is not None:
            self._convs(self,None,set([self]))
            
    def _convert(self,g,un):
        if self._conversion is not None and self._conversion.unit is un:
            return self._conversion.to(g)
            
        if self._gconv is not None and self._gconv._unit is un:
            return self._gconv.to(g)
            
        if un._conversion is not None and un._conversion.unit is self:
            return un._conversion.frm(g)
                
        if un._gconv is not None and un._gconv._unit is self:
            return un._gconv.frm(g)
        
        if self._gconv is not None and un._gconv is not None and un._gconv._unit is self._gconv._unit:
            g = self._gconv.to(g)
            g = un._gconv.frm(g)
            return g
            
        return None
        
    def convert(self,g,unit):
        """
        Converts a number or gummy from a quanitity with the units self to unit.
        
        This is not intended be used directly, instead use the `gummy.convert`
        method.
        """
        unit = Unit.unit(unit)
        if unit is self:
            return g
            
        ret = self._convert(g,unit)
        if ret is not None:
            return ret
        self._convtree()
        unit._convtree()
        ret = self._convert(g,unit)

        unm = unit.short_name
        if ret is None:
            raise NoUnitConversionFoundError('no conversion found from unit ' + self.short_name + ' to unit ' + unm)
        return ret
    
    @property
    def is_dimensionless(self):
        """
        `bool`, read-only
        
        Returns `True` if a conversion exists between `self` and `one`, and `False`
        if not.
        """
        if self.conversion is not None and self.conversion.unit is one:
            return True
        if self._gconv is not None and self._gconv._unit is one:
            return True
        return False
        
    @property
    def units(self):
        """read-only
        
        Returns a list of the constituent units and their exponents, e.g. for
        kg m**2/s units would be [(kg, 1), (m, 2), (s, -1)].
        """
        return [(self,1)]
    
    @staticmethod
    def _cancel(a,b=[],asign=1,bsign=1):
        # Takes .units lists a**asign and b**bsign and attempts to convert 
        # each unit in a into another unit either in a or b, and if successful 
        # cancels that unit.  Returns a tuple containing the  new .units 
        # list for the reduced composite unit a**asign * b**bsign as well as
        # the conversion factor between the unreduced and reduced composite unit.
        # In addition if both a an b contain dimensionless units then all
        # dimensionless units are converted to one.
        c = 1
        b = {i: bsign*j for i,j in b}
        for k,v in a:
            v *= asign
            if k in b:
                b[k] += v
            else:
                converted = False
                if not k.is_dimensionless:
                    for f in b.keys():
                        try:
                            c *= k.convert(1,f)**v
                            b[f] += v
                            converted = True
                            break
                        except NoUnitConversionFoundError:
                            pass
                if not converted:
                    b[k] = v
        un = []
        for k,v in b.items():
            if v != 0 and k is not one:
                uun,cc = k._cpow(v)
                un.append(uun[0])
                c *= cc
        return (un,c)
        
    def _mul_cancel(self,b):
        # self*b canceling units when possible
        c = Unit._cancel(b.units)
        d = Unit._cancel(self.units,c[0])
        return (_CompositeUnit(d[0]),c[1]*d[1])
        
    def _div_cancel(self,b):
        # self/b canceling units when possible        
        c = Unit._cancel(b.units)
        d = Unit._cancel(self.units,b=c[0],asign=1,bsign=-1)
        return (_CompositeUnit(d[0]),c[1]*d[1])
        
    def _rdiv_cancel(self,b):
        # b/self canceling units when possible
        c = Unit._cancel(b.units)
        d = Unit._cancel(self.units,b=c[0],asign=-1,bsign=1)
        return (_CompositeUnit(d[0]),c[1]*d[1])
    
    def _nummy_add(self,a,b):
        return nummy._add(a,b)

    def _add(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_add(a,b),self)
        if not bunit.linear:
            return bunit._radd(b,self,a,not aconv)
        try:
            if aconv:
                return (self._nummy_add(self.convert(a,bunit),b),bunit)
            return (self._nummy_add(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())

    def _nummy_radd(self,a,b):
        # ummy._radd assumes that b is not an ummy.  
        if isinstance(b,nummy):
            return nummy._add(b,a)
        return nummy._radd(a,b)
    
    def _radd(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_radd(a,b),self)
        try:
            if aconv:
                return (self._nummy_radd(self.convert(a,bunit),b),bunit)
            return (self._nummy_radd(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())
    
    def _nummy_sub(self,a,b):
        return nummy._sub(a,b)
    
    def _sub(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_sub(a,b),self)
        try:
            if not bunit.linear:
                return bunit._rsub(b,self,a,not aconv)
            if aconv:
                return (self._nummy_sub(self.convert(a,bunit),b),bunit)
            return (self._nummy_sub(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + bunit.tostring() + ' may not be subtracted from a quantity with unit ' + self.tostring())

    def _nummy_rsub(self,a,b):
        if isinstance(b,nummy):
            return nummy._sub(b,a)
        return nummy._rsub(a,b)
    
    def _rsub(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_rsub(a,b),self)
        try:
            if aconv:
                return (self._nummy_rsub(self.convert(a,bunit),b),bunit)
            return (self._nummy_rsub(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be subtracted from a quantity with unit ' + bunit.tostring())
        
    def _nummy_mod(self,a,b):
        return nummy._mod(a,b)
    
    def _mod(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_mod(a,b),self)
        if not bunit.linear:
            return bunit._rmod(b,self,a,not aconv)
        try:
            if aconv:
                return (self._nummy_mod(self.convert(a,bunit),b),bunit)
            return (self._nummy_mod(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())

    def _nummy_rmod(self,a,b):
        if isinstance(b,nummy):
            return nummy._mod(b,a)
        return nummy._rmod(a,b)
    
    def _rmod(self,a,bunit,b,aconv):
        if self is bunit:
            return (self._nummy_rmod(a,b),self)
        try:
            if aconv:
                return (self._nummy_rmod(self.convert(a,bunit),b),bunit)
            return (self._nummy_rmod(a,bunit.convert(b,self)),self)
        except NoUnitConversionFoundError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())

    def _nummy_mul(self,a,b):
        return nummy._mul(a,b)
    
    def _mul(self,a,bunit,b,aconv):
        if not bunit.linear:
            return bunit._rmul(b,self,a,not aconv)

        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return (self._nummy_mul(self._nummy_mul(a,b),c),un)
    
    def _nummy_rmul(self,a,b):
        if isinstance(b,nummy):
            return nummy._mul(b,a)
        return nummy._rmul(a,b)

    def _rmul(self,a,bunit,b,aconv):
        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return (self._nummy_mul(self._nummy_rmul(a,b),c),un)
    
    def _nummy_truediv(self,a,b):
        return nummy._truediv(a,b)
    
    def _truediv(self,a,bunit,b,aconv):
        if not bunit.linear:
            return bunit._rtruediv(b,self,a,not aconv)

        if aconv:
            un,c = self._div_cancel(bunit)
        else:
            un,c = bunit._rdiv_cancel(self)
        
        return (self._nummy_mul(self._nummy_truediv(a,b),c),un)
    
    def _nummy_rtruediv(self,a,b):
        if isinstance(b,nummy):
            return nummy._truediv(b,a)
        return nummy._rtruediv(a,b)

    def _rtruediv(self,a,bunit,b,aconv):
        if aconv:
            un,c = self._rdiv_cancel(bunit)
        else:
            un,c = bunit._div_cancel(self)
        
        return (self._nummy_mul(self._nummy_rtruediv(a,b),c),un)
    
    def _cpow(self,v):
        return ([[self,v]],1)
    
    def _nummy_pow(self,a,b):
        return nummy._pow(a,b)
    
    def _pow(self,a,bunit,b,aconv):
        if not bunit.linear:
            return bunit._rpow(b,self,a,not aconv)
        
        if bunit is not one:
            try:
                b = bunit.convert(b,one)
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
                
        bb = float(b)
        if bb != b:
            try:
                a = self.convert(a,one)
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('a gummy that is not dimensonless may not be raised to a power with an uncertainty')
            return (self._nummy_pow(a,b),one)
        b = float(b)
        if int(b) == b:
            b = int(b)
        
        if b == 0:
            return (self._nummy_pow(a,0),one)
        
        un,c = self._cpow(b)
        un = _CompositeUnit(un)
        
        return (self._nummy_mul(self._nummy_pow(a,b),c),un)

    def _nummy_rpow(self,a,b):
        if isinstance(b,nummy):
            return nummy._pow(b,a)
        return nummy._rpow(a,b)
    
    def _rpow(self,a,bunit,b,aconv):
        if self is not one:
            try:
                b = bunit.convert(b,one)
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
                
        aa = float(a)
        if aa != a:
            try:
                b = self.convert(b,one)
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('a gummy that is not dimensonless may not be raised to a power with an uncertainty')
            return (self._nummy_rpow(a,b),one)
        a = float(a)
        if int(a) == a:
            a = int(a)
        
        if a == 0:
            return (self._nummy_pow(b,0),one)
        
        un,c = _CompositeUnit(bunit._cpow(a))
        
        return (self._nummy_mul(self._nummy_rpow(a,b),c),un)
    
    def _nummy_neg(self,a):
        return nummy.__neg__(a)
    
    def _neg(self,a):
        return (self._nummy_neg(a),self)
    
    def _nummy_pos(self,a):
        return nummy.__pos__(a)
    
    def _pos(self,a):
        return (self._nummy_pos(a),self)
    
    def _nummy_abs(self,a):
        return nummy.__abs__(a)
    
    def _abs(self,a):
        return (self._nummy_abs(a),self)
    
    def __mul__(self,v):          
        if v is one:
            return self
        
        if not isinstance(v,Unit):
            try:
                if v == 1:
                    return self
            except:
                pass
            raise TypeError('unsupported operand type(s) for *: Unit and ' + str(v))
            
        if not v.linear:
            return v.__rmul__(self)
            
        return _CompositeUnit(self.units + v.units)
        
    def __truediv__(self,v):
        if v is one:
            return self
        
        if not isinstance(v,Unit):
            try:
                if v == 1:
                    return self
            except:
                pass
            raise TypeError('unsupported operand type(s) for /: Unit and ' + str(v))
            
        if not v.linear:
            return v.__rtruediv__(self)
            
        vi = [(e[0],-e[1]) for e in v.units]
        return _CompositeUnit(self.units + vi)

    def __rtruediv__(self,v):
        if v is one or v == 1:
            return _CompositeUnit([(e[0],-e[1]) for e in self.units])
        
    def __pow__(self,v):
        try:
            if v == 0:
                return one
            
            if v == 1:
                return self
        except:
            raise TypeError('unsupported operand type(s) for **: Unit and ' + str(v))
        
        if v != float(v):
            raise TypeError('unsupported operand type(s) for **: Unit and ' + str(v))

        return _CompositeUnit([(e[0],v*e[1]) for e in self.units])
        
    # Define comparison operators and a hash function so, for the 
    # _CompositeUnit subclass we can sort the Unit order and store them as keys 
    # in a unique OrderedDict.
    def __lt__(self, b):
        if self.order == b.order:
            return self.name < b.name
        return self.order < b.order

    def __le__(self, b):
        if self.__eq__(b):
            return True
        if self.order == b.order:
            return self.name < b.name
        return self.order < b.order

    def __eq__(self, b):
       return (self is b)

    def __ge__(self, b):
        if self.__eq__(b):
            return True
        if self.order == b.order:
            return self.name > b.name
        return self.order > b.order

    def __gt__(self, b):
        if self.order == b.order:
            return self.name > b.name
        return self.order > b.order

    def __ne__(self, b):
        return not (self is b)
        
    def __hash__(self):
        return id(self)
        
    
class _CompositeUnit(Unit):
    # Instances of this class represent composite derived unit.
    
    living_units = weakref.WeakValueDictionary()
    
    parent = None
    
    @staticmethod
    def _unicode_super(e):
        if np.modf(e)[0] == 0:
            se = str(int(e))
            so = ''
            for s in se:
                if s == '-':
                    so += '\u207B'
                if s == '0':
                    so += '\u2070'
                if s == '1':
                    so += '\u00B9'
                if s == '2':
                    so += '\u00B2'
                if s == '3':
                    so += '\u00B3'
                if s == '4':
                    so += '\u2074'
                if s == '5':
                    so += '\u2075'
                if s == '6':
                    so += '\u2076'
                if s == '7':
                    so += '\u2077'
                if s == '8':
                    so += '\u2078'
                if s == '9':
                    so += '\u2079'
            return so
        return '**' + '{:.2f}'.format(e).rstrip('0')
    
    @staticmethod
    def _format_sscript(e):
        mf = np.modf(e)
        if mf[0] == 0:
            return '{:.0f}'.format(e)
        if mf[0] == 0.5:
            return '{:.0f}'.format(2*e) + '/2'
        if mf[0] == 0.25:
            return '{:.0f}'.format(4*e) + '/4'
        if np.modf(3*e)[0] < 0.01:
            return '{:.0f}'.format(3*e) + '/3'
        return '{:.2f}'.format(e).rstrip('0')
        
    @staticmethod
    def _iaddsym(txt,s,e,fmt,f):
        if txt.startswith('\t'):
            txt = txt[1:]
        if e == 0:
            return txt
        if f == 1 and e < 0:
            return txt
        if f == 2:
            if e > 0:
                return txt
            e = -e
            
        if fmt == 'html':
            if txt != '':
                txt += '&middot;'
            txt += s
            if e != 1:
                txt += '<sup>' + _CompositeUnit._format_sscript(e) + '</sup>'
        elif fmt == 'latex':
            if txt != '':
                txt += r'\cdot'
            txt += s
            if e != 1:
                txt += '^{' + _CompositeUnit._format_sscript(e) + '}'
        elif fmt == 'ascii':
            if txt != '':
                txt += '*'
            txt += s
            if e != 1:
                txt += '**' + '{:.2f}'.format(e).rstrip('0').rstrip('.')
        else:
            if txt != '':
                txt += ' '
            txt += s
            if e != 1:
                txt += _CompositeUnit._unicode_super(e)
            
        return txt
        
    def _addsym(txt,u,fmt=None,f=0):
        return _CompositeUnit._iaddsym(txt,u[0].tostring(fmt=fmt),u[1],fmt,f)
    
    def __new__(cls,ul):
        ul = list(ul)
        
        if len(ul) == 0:
            return one
        
        if len(ul) == 1 and ul[0][1] == 1:
            return ul[0][0]
        
        ul = sorted(ul,key = lambda x: (x[0],-x[1]))
        
        uid = tuple([(id(u[0]),u[1]) for u in ul])
        ret = _CompositeUnit.living_units.get(uid)
        if ret is not None:
            return ret
        
        linear = True
        for u,e in ul:
            if not u.linear:
                linear = False
                try:
                    cls = u.get_composite(ul)
                except NotImplementedError:
                    raise IncompatibleUnitsError('unit ' + u.tostring() + ' may not be combined with other units')
                break
        
        ret = object.__new__(cls)
        
        ret.linear = linear
        ret.order = 100
        ret._units = ul
        _CompositeUnit.living_units[uid] = ret
        ret._oconv = None
        ret._have_conversion = False
        ret._aliases = set()
        
        nm = ''
        snm = ''
        for u in ul:          
            nm += ' '
            snm += ' '
            nmi = u[0].name
            if ' ' in nmi and (not nmi[0].startswith('[') or not nmi.endswith(']')):
                nmi = '[' + nmi + ']'
            snmi = u[0].short_name
            if ' ' in snmi and (not snmi[0].startswith('[') or not snmi.endswith(']')):
                snmi = '[' + snmi + ']'
            if u[1] == 1:
                    nm += nmi
                    snm += snmi
            else:
                    nm += nmi + '**' + str(u[1])
                    snm += snmi + '**' + str(u[1])
        
        ret.name = nm.strip()
        ret.short_name = snm.strip()
        
        ret._make_symbols(ul)
        
        return ret
    
    def __init__(self,ul):
        pass
        
    def _make_symbols(self,ul):
        sym = ''
        hsym = ''
        lsym = ''
        asym = ''
        symsn = ''
        symsd = ''
        hsymsn = ''
        hsymsd = ''
        lsymsn = ''
        lsymsd = ''
        asymsn = ''
        asymsd = ''
        
        nd = 0
        for u in ul:
            sym = _CompositeUnit._addsym(sym,u)
            hsym = _CompositeUnit._addsym(hsym,u,fmt = 'html')
            lsym = _CompositeUnit._addsym(lsym,u,fmt = 'latex')
            asym = _CompositeUnit._addsym(asym,u,fmt = 'ascii')
            
            symsn = _CompositeUnit._addsym(symsn,u,f=1)
            hsymsn = _CompositeUnit._addsym(hsymsn,u,fmt = 'html',f=1)
            lsymsn = _CompositeUnit._addsym(lsymsn,u,fmt = 'latex',f=1)
            asymsn = _CompositeUnit._addsym(asymsn,u,fmt = 'ascii',f=1)
            
            symsd = _CompositeUnit._addsym(symsd,u,f=2)
            hsymsd = _CompositeUnit._addsym(hsymsd,u,fmt = 'html',f=2)
            lsymsd = _CompositeUnit._addsym(lsymsd,u,fmt = 'latex',f=2)
            asymsd = _CompositeUnit._addsym(asymsd,u,fmt = 'ascii',f=2)
            
            if u[1] < 0:
                nd += 1
            
        self.sym = sym
        self.hsym = hsym
        self.lsym = lsym
        self.asym = asym
        
        if symsn == '' and nd > 0:
            symsn = '1'
            hsymsn = '1'
            lsymsn = '1'
            asymsn = '1'
        
        if nd > 0:
            if nd > 1:
                symsd = '(' + symsd + ')'
                hsymsd = '(' + hsymsd + ')'
                lsymsd = '(' + lsymsd + ')'
                asymsd = '(' + asymsd + ')'
            syms = symsn + '/' + symsd
            hsyms = hsymsn + '/' + hsymsd
            lsyms = lsymsn + '/' + lsymsd
            asyms = asymsn + '/' + asymsd
        else:
            syms = symsn
            hsyms = hsymsn
            lsyms = lsymsn
            asyms = asymsn
        self.syms = syms
        self.hsyms = hsyms
        self.lsyms = lsyms
        self.asyms = asyms
        
    @property
    def units(self):
        """
        read-only

        Returns a list of the constituent units and their exponents, e.g. for
        kg m**2/s units would be [(kg, 1), (m, 2), (s, -1)].
        """
        return self._units
        
    def tostring(self,fmt=None,solidus=True,mulsep=False,**kwds):
        """
        Returns a string containing the symbol for the unit the format given by
        the keyword fmt which may be set to a string the values 'html', 'latex', 
        'ascii', or 'unicode'.
        
        If solidus is `True` a slash will be used instead of negative exponents
        and if `mulsep` is `True` a mid-dot or asterix will be used instead of a space
        to separate units.
        """
        if fmt is None or fmt == 'utf-8' or fmt == 'unicode':
            if solidus:
                ret = self.syms
            else:
                ret = self.sym
            if not mulsep:
                ret = ret.replace('\u00B7',' ')
            return ret
        elif fmt == 'html':
            if solidus:
                ret = self.hsyms
            else:
                ret = self.hsym
            if not mulsep:
                ret = ret.replace('&middot;','&thinsp;')
        elif fmt == 'latex':
            if solidus:
                ret = self.lsyms
            else:
                ret = self.lsym
            if not mulsep:
                ret = ret.replace(r'\cdot',r'\,')
            return ret
        elif fmt == 'ascii':
            if solidus:
                ret = self.asyms
            else:
                ret = self.asym
            if not mulsep:
                ret = ret.replace('**','#@@#')
                ret = ret.replace('*',' ')
                ret = ret.replace('#@@#','**')
        else:
            raise ValueError('format ' + str(fmt) + ' is unrecognized')
        return ret
        
    def _find_conv(self):    
        c = None
        cun = {}
        cg = None
        cgun = {}
        for un,e in self._units:
            if un.conversion is None:
                if un in cun:
                    cun[un] += e
                else:
                    cun[un] = e
            else:
                p = un.conversion.pow(e)
                if c is None:
                    c = p
                else:
                    c = c.chain(p)
                for t in p.unit.units:
                    if t[0] in cun:
                        cun[t[0]] += t[1]
                    else:
                        cun[t[0]] = t[1]
                    
            un._convtree()
                    
            if un._gconv is None:
                if un in cgun:
                    cgun[un] += e
                else:
                    cgun[un] = e
            else:
                p = un._gconv.pow(e)
                if cg is None:
                    cg = p
                else:
                    cg = c.chain(p)
                for t in p.unit.units:
                    if t[0] in cgun:
                        cgun[t[0]] += t[1]
                    else:
                        cgun[t[0]] = t[1]
                        
        cun = [(k,v) for k,v in cun.items() if v != 0]
        cgun = [(k,v) for k,v in cgun.items() if v != 0]
        
        cun = _CompositeUnit(cun)
        if cun is not self:
            if c is None:
                c = Conversion(cun,1)
            else:
                c._unit = cun
            c.parent = self
            self._conversion = c
        else:
            self._conversion = None
        
        cgun = _CompositeUnit(cgun)
        if cgun is not self:
            if cg is None:
                cg = Conversion(cgun,1)
            else:
                cg._unit = cgun
            cg.parent = self
            self._gconv = cg
        else:
            self._gconv = None

        self._have_conversion = True

    
    def _cpow(self,v):
        c = 1
        un = []
        for u,e in self.units:
            uun,cc = u._cpow(v*e)
            un.append(uun[0])
            c *= cc
        return (un,c)

    @property
    def conversion(self):
        """
        Gets or sets the conversion object for the unit.  
        
        ***This property is not intended to be used directly and setting this 
        property may cause problems***
        """
        if not self._have_conversion:
            self._find_conv()
        return self._conversion
    
    @property
    def shadowed_aliases(self):
        """
        Returns a set of the shadowed aliases of this unit.
        """
        sha = {a for a in self._aliases if Unit.unit('[' + a + ']',exception=False) is not self}
        return sha
        
    def _convtree(self):
        if not self._have_conversion:
            self._find_conv()
        Unit._convtree(self)
        
    def _convs(self,d,c,h):
        if not self._have_conversion:
            self._find_conv()
        Unit._convs(self,d,c,h)

    
class _BuiltinLib:
    """
    Allows Unit objects to be added to Unit._builtin_lib rather than Unit._lib
    using a "with" statement.
    
    Example
    -------
    
    >>> with _BuiltinLib():
    ...     one = Unit('1','',add_symbol=False)
    
    """
    def  __enter__(self):
        Unit._open_lib = Unit._builtin_lib
    def __exit__(self,type, value, traceback):
        Unit._open_lib = Unit._lib
    
class _One(Unit):
    """
    The only instance of this class should be the unit `one`.
    """
    def __mul__(self,v):            
        return v
        
    def __truediv__(self,v):
        return 1/v

    def __rtruediv__(self,v):
        return v
        
    def __pow__(self,v):
        return self
    
    @property
    def is_dimensionless(self):
        return True
    
with _BuiltinLib():
    one = _One('1','',add_symbol=False)