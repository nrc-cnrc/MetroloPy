# -*- coding: utf-8 -*-

# module unit

# Copyright (C) 2026 National Research Council Canada
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
Classes that for handling units are defined here:  Conversion, Unit,
_CompositeUnit, one, Quantity and QuantityArray
"""

import numpy as np
import weakref
import warnings
from numbers import Rational,Integral
from fractions import Fraction
from abc import ABCMeta

from .exceptions import (UnitLibError,CircularUnitConversionError,
                         UnitNotFoundError,IncompatibleUnitsError,
                         UnitLibNotFoundError,UnitWarning)
from .printing import PrettyPrinter,MetaPrettyPrinter
from .indexed import Indexed
from .unitparser import _UnitParser
from .mfraction import MFraction
from .abc import AbcQuantity,AbcUnit


def unit(name,exception=True):
    """
    Finds an returns a `Unit` from the unit library.  This function is an 
    alias for the `Unit.unit` static method.
    
    Parameters
    ----------
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
            
    Returns
    -------
    A `Unit` instance or possibly None if the `exception` parameter is
    set to `True`.
    """

    return Unit.unit(name,exception=exception)

class Conversion:
    """
    Represents a unit conversion.  This class should only be used as arguments
    to `Unit` object initializers.  Each conversion should be associated with
    one and only one parent Unit; do not re-use conversion instances with more
    than one `Unit` instance.  
    
    The base `Conversion` class is defined by a single conversion factor.  
    More comlicated conversions overriding the `to` and `frm` methods should 
    inherit from the `NonlinearConversion` subclass and be used with units 
    that inherit from the `NonlinearUnit` subclass.
    
    Parameters
    ----------
    unit:  `str` or `Unit`
        the Unit that the parent Unit will be converted to.
        
    factor:  `float`, optional
        The conversion factor between the parent Unit and the new Unit:
        [value with new Unit] = factor * [value with parent Unit].  The
        conversion factor may not have a unit, but may be an ummy of gummy
        with not unit one.  The default value is 1.
    """
    linear = True
    
    def __init__(self,unit,factor=1):
        self._unit = unit
        if isinstance(factor,Fraction):
            factor = MFraction(factor)
        self.factor = factor
        
    def to(self,g):
        return g*self.factor
    
    def frm(self,g):
        return g/self.factor
        
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


def _f_bop(f,rf,x,y):
    if isinstance(x,AbcQuantity) and x.unit:
        if isinstance(y,AbcQuantity):
            yunit = y.unit
            y = y.value
        else:
            yunit = one
        return f(x.unit,x.value,yunit,y,x.autoconvert)
    else:
        if isinstance(x,AbcQuantity):
            xunit = x.unit
            x = x.value
        else:
            xunit = one
        return rf(y.unit,y.value,xunit,x,y.autoconvert)
    
def _f_ratio(f,x,y):
    if not isinstance(y,AbcQuantity):
        x = x.convert(one).value
    elif not isinstance(x,AbcQuantity):
        y = y.convert(one).value
    else:
        y = y.convert(x.unit).value
        x = x.value
    return (f(x,y),one)

class MetaUnit(MetaPrettyPrinter,ABCMeta):
    pass

class Unit(AbcUnit,PrettyPrinter,Indexed,metaclass = MetaUnit):
    """
    Creating an instance of this class creates a representation of a physical
    unit and adds it to the unit library.  Once created the unit intance may
    be retrived by passing a string with the unit name or alias to the `unit`
    or `Unit.unit` functions.  Units can be multiplied and divided by other 
    Units or Quantities and raised to numerical powers.  Multiplying or 
    dividing a numerical value by a Unit will create a Quantity instance.
    
    Parameters
    ----------
        
    name:  `str`
        The name of the unit.  The name can be used access the unit with the
        `unit` function, but note that if you define a Unit with an
        identical name to a previously defined unit then the older name will
        be shadowed.
    
    symbol:  `str`
        A unicode symbol used when displaying the unit.  If the `add_symbol`
        parameter is set to `True`, then this symbol can also be used
        to access the unit with the `unit` function.
    
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
        a short name which can be used as an additional alias for the unit in
        the unit library
    
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
    _builtins_to_import = ['..usunits','..miscunits','..siunits','..relunits']
    _builtin_lib = {}
    _lib = {}
    _open_lib = _lib
    
    _used_units = {}
    
    _ufunc_dict = {np.square: lambda x: (np.square(x.value),x.unit**2),
                   np.sqrt: lambda x: (np.sqrt(x.value),x.unit**0.5),
                   np.cbrt: lambda x: (np.cbrt(x.value),x.unit**(1/3)),
                   np.reciprocal: lambda x: (np.reciprocal(x.value),x.unit**-1),
                   np.absolute: lambda x: (np.absolute(x.value),x.unit),
                   np.negative: lambda x: (np.negative(x.value),x.unit),
                   np.positive: lambda x: (np.positive(x.value),x.unit),
                   np.conjugate: lambda x: (np.conjugate(x.value),x.unit),
                   np.rint: lambda x: (np.rint(x.value),x.unit),
                   np.floor: lambda x: (np.floor(x.value),x.unit),
                   np.ceil: lambda x: (np.ceil(x.value),x.unit),
                   np.trunc: lambda x: (np.trunc(x.value),x.unit),
                   np.add: lambda x,y: _f_bop(Unit._add,Unit._radd,x,y),
                   np.subtract: lambda x,y: _f_bop(Unit._sub,Unit._rsub,x,y),
                   np.multiply: lambda x,y: _f_bop(Unit._mul,Unit._rmul,x,y),
                   np.divide: lambda x,y: _f_bop(Unit._truediv,Unit._rtruediv,x,y),
                   np.floor_divide: lambda x,y: _f_bop(Unit._floordiv,Unit._rfloordiv,x,y),
                   np.mod: lambda x,y: _f_bop(Unit._mod,Unit._rmod,x,y),
                   np.power: lambda x,y: _f_bop(Unit._pow,Unit._rpow,x,y),
                   np.arctan2: lambda x,y: _f_ratio(np.arctan2,x,y)}
        
    solidus = True
    mulsep = False
    
    @staticmethod
    def unit(txt,exception=True):
        """
        Finds an returns a `Unit` from the unit library.
        
        Parameters
        ----------
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
                
        Returns
        -------
        A `Unit` instance or possibly None if the `exception` parameter is
        set to `True`.
        """
        try:
            if isinstance(txt,Unit):
                return txt
            if txt is None or txt == 1:
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
                if ul[0][0] == '1' or ul[0][0] == '' or ul[0][0] == 1:
                    return one
                else:
                    return Unit.get(ul[0][0])
            ull = []
            for u in ul:
                if u[0] == '1' or u[0] == '' or u[0] == 1:
                    un = one
                else:
                    un = Unit.get(u[0])
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
            if u[0] == '1' or u[0] == '' or u[0] == 1:
                un = one
            else:
                un = Unit.get(u[0])
            if hasattr(un,'_getme'):
                u[0] = un._getme(ul,u[1])
            else:
                u[0] = un
            if isinstance(un,_CompositeUnit):
                raise ValueError('aliases for composite units may not be used in Unit.reorder')
        unit._make_symbols(ul)
        
    @classmethod
    def _raise_not_found(cls,name):
        raise UnitNotFoundError('unit "' + str(name) + '" was not found')
        
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
    _linear = True
    
    def __init__(self,name,symbol,conversion=None,short_name = None,
                 add_symbol=False,html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,description=None,order = -1,**kwds):
        if conversion is not None and not isinstance(conversion,Conversion):
            raise ValueError('conversion is not an instance of the Conversion class')
        
        if symbol.startswith('\t\t'):
            raise ValueError('the symbol may not start with more than one tab character; only the latex_symbol may start with double tabs')
        
        super().__init__(name,symbol,short_name,add_symbol,html_symbol,
                         latex_symbol,ascii_symbol,description)
        
        self.conversion = conversion
        self.order = order
        
        self.parent = kwds.get('parent')
        
        #reset _used_units incase we are shadowing any unit definitions
        Unit._used_units = {} 

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
            
    def _convs(self,h):
        if self in h:
            raise CircularUnitConversionError('a circular conversion chain was found nwith unit ' + self.name)
        h.add(self)
        
        try:
            c = self._conversion
            u = c.unit
            if isinstance(u,_CompositeUnit):
                una = u._units
                unb = {}
                for k,v in una:
                    if k in unb:
                        unb[k] += v
                    else:
                        unb[k] = v
                unb = {k:v for k,v in unb.items() if v != 0 and k is not one}
                if len(unb) == 0:
                    u = one
                else:
                    u = _CompositeUnit([(k,v) for k,v in unb.items()])
                    u._find_conv()
                c = c.chain(Conversion(u,1))
                     
        except UnitNotFoundError:
            warnings.warn('while searching conversions:  unit ' + self.name + ' should have a conversion to ' + self._conversion.unit_name + ' in library ' + str(self._conversion.unit_library) + ' but the second unit was not found',UnitWarning,stacklevel=3)
            return
        except UnitLibNotFoundError:
            warnings.warn('while searching conversions:  unit ' + self.name + ' should have a conversion to ' + self._conversion.unit_name + ' in library ' + str(self._conversion.unit_library) + ' but the library containing the second unit has not been loaded',UnitWarning,stacklevel=3)
            return
        
        if u._conversion is None:
            self._gconv = c
        else:
            r = u._convs(h)
            self._gconv = c.chain(r)
            
        return self._gconv
        
    @property
    def linear(self):
        """
        Gets a `bool` value indicating whether the Unit is linear.  If the
        unit is linear then any associated values will follow the standard 
        rules of arithmatic for Quantaties and the unit's Conversion is 
        defined by multiplying or dividing by a single conversion factor.
        Nonlinear units may have a more complicated conversion and may
        override the unusal operator methods.
        """

        return self._linear
                
    def _convtree(self):
        if self.conversion is not None and self._gconv is None:
            self._convs(set())
            
    def _convert(self,g,un):
        if self.conversion is not None and self.conversion.unit is un:
            return self.conversion.to(g)
            
        if self._gconv is not None and self._gconv._unit is un:
            return self._gconv.to(g)
            
        if un.conversion is not None and un.conversion.unit is self:
            return un.conversion.frm(g)
                
        if un._gconv is not None and un._gconv._unit is self:
            return un._gconv.frm(g)
        
        if self._gconv is not None and un._gconv is not None and un._gconv._unit is self._gconv._unit:
            g = self._gconv.to(g)
            g = un._gconv.frm(g)
            return g
            
        return None
        
    def convert(self,g,unit):
        """
        Converts a number from a quantity with the units `self` to `unit`.
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
            raise  IncompatibleUnitsError('no conversion found from unit ' + self.short_name + ' to unit ' + unm)
        return ret
    
    @property
    def is_dimensionless(self):
        """
        `bool`, read-only
        
        Returns `True` if a conversion exists between `self` and `one`, and
        `False` if not.
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
    
    @property
    def base(self):
        """read-only
        
        Returns the base unit.  That is we follow the chain of conversions
        starting with this units `Conversion` instance until we get to a unit
        that has no `Conversion` instance.  If this unit has no `Conversion`
        instance `self` is returned.
        """
        if self.conversion is None:
            return self
        self._convtree()
        return self._gconv._unit
    
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
                        except  IncompatibleUnitsError:
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
    
    def _add(self,s,vunit,v,aconv):
        if self is vunit:
            return (s + v,self)
        
        if not vunit.linear:
            return vunit._radd(v,self,s,not aconv)
        
        if aconv:
            s = self.convert(s,vunit)
            return (s + v,vunit)
        else:
            v = vunit.convert(v,self)
            return (s + v,self)
        
    def _radd(self,s,vunit,v,aconv):
        if self is vunit:
            return (v + s,self)
        
        if not vunit.linear:
            return vunit._add(v,self,s,not aconv)
        
        if aconv:
            s = self.convert(s,vunit)
            return (v + s,vunit)
        else:
            v = vunit.convert(v,self)
            return (v + s,self)
        
    def _sub(self,s,vunit,v,aconv):
        if self is vunit:
            return (s - v,self)
        
        if not vunit.linear:
            return vunit._rsub(v,self,s,not aconv)
        
        if aconv:
            s = self.convert(s,vunit)
            return (s - v,vunit)
        else:
            v = vunit.convert(v,self)
            return (s - v,self)
        
    def _rsub(self,s,vunit,v,aconv):
        if self is vunit:
            return (v - s,self)
        
        if not vunit.linear:
            return vunit._sub(v,self,s,not aconv)
        
        if aconv:
            s = self.convert(s,vunit)
            return (v - s,vunit)
        else:
            v = vunit.convert(v,self)
            return (v - s,self)
        
    def _mul_cancel(self,b):
        # self*b canceling units when possible
        c = Unit._cancel(b.units)
        d = Unit._cancel(self.units,c[0])
        return (_CompositeUnit(d[0]),c[1]*d[1])
    
    def _mul(self,a,bunit,b,aconv):        
        if isinstance(b,Unit):
            bunit = b
            b = 1
        
        if not bunit.linear:
            return bunit._rmul(b,self,a,not aconv)

        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return ((a*b)*c,un)
    
    def _rmul(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if not bunit.linear:
            return bunit._mul(b,self,a,not aconv)

        if aconv:
            un,c = self._mul_cancel(bunit)
        else:
            un,c = bunit._mul_cancel(self)
        
        return ((b*a)*c,un)
        
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
        
    def _truediv(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if not bunit.linear:
            return bunit._rtruediv(b,self,a,not aconv)

        if aconv:
            un,c = self._div_cancel(bunit)
        else:
            un,c = bunit._rdiv_cancel(self)
        
        return ((a/b)*c,un)
    
    def _rtruediv(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if aconv:
            un,c = self._rdiv_cancel(bunit)
        else:
            un,c = bunit._div_cancel(self)
        
        return ((b/a)*c,un)
    
            
    def _floordiv(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if not bunit.linear:
            return bunit._rtruediv(b,self,a,not aconv)

        if aconv:
            un,c = self._div_cancel(bunit)
        else:
            un,c = bunit._rdiv_cancel(self)
        
        return ((a // b)*c,un)
    
    def _rfloordiv(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if aconv:
            un,c = self._rdiv_cancel(bunit)
        else:
            un,c = bunit._div_cancel(self)
        
        return ((b // a)*c,un)
        
    def _mod(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if self is bunit:
            return (a % b,self)
        if not bunit.linear:
            return bunit._rmod(b,self,a,not aconv)
        try:
            if aconv:
                return (self.convert(a,bunit) % b,bunit)
            return (a % bunit.convert(b,self),self)
        except IncompatibleUnitsError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())

    def _rmod(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if self is bunit:
            return (b % a,self)
        try:
            if aconv:
                return (b % self.convert(a,bunit),bunit)
            return (bunit.convert(b,self) % a,self)
        except IncompatibleUnitsError:
            raise IncompatibleUnitsError('a quantity with unit ' + self.tostring() + ' may not be added to a quantity with unit ' + bunit.tostring())

    def _cpow(self,v):
        return ([[self,v]],1)
    
    def _pow(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if not bunit.linear:
            return bunit._rpow(b,self,a,not aconv)
        
        if bunit is not one:
            try:
                b = bunit.convert(b,one)
            except IncompatibleUnitsError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
                
        un,c = self._cpow(b)
        un = _CompositeUnit(un)
        return ((a**b)*c,un)
    
    def _rpow(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if self is not one:
            try:
                a = self.convert(a,one)
            except IncompatibleUnitsError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
        un,c = bunit._cpow(a)
        un = _CompositeUnit(un)
        return ((b**a)*c,un)
    
    def _ufunc(self,func,*args,**kwds):
        if self.linear:
            try:
                nl = next(a for a in args 
                          if isinstance(a,AbcQuantity) and not a.unit.linear)
                return nl.unit._ufunc(func,*args,**kwds)
            except StopIteration:
                pass
            
        try:
            f = type(self)._ufunc_dict[func]
        except KeyError:
            raise NotImplementedError()
        return f(*args,**kwds)
    
    def _neg(self,a):
        return (-a,self)
    
    def _pos(self,a):
        return (+a,self)
    
    def _abs(self,a):
        return (abs(a),self)
    
    def __mul__(self,v):
        if v == 1:
            return self
        
        if not isinstance(v,Unit):
           raise TypeError('unsupported operand type(s) for *: Unit and ' + str(v))
        
        if not v.linear:
            return v.__rmul__(self)
        
        return _CompositeUnit(self.units + v.units)
    
    def __rmul__(self,v):
        if v == 1:
            return self
        
        if not isinstance(v,Unit):
            raise TypeError('unsupported operand type(s) for *: Unit and ' + str(v))
        
        if not v.linear:
            return v.__mul__(self)
        
        return _CompositeUnit(self.units + v.units)
        
    def __truediv__(self,v):
        if v == 1:
            return self
        
        if not isinstance(v,Unit):
            raise TypeError('unsupported operand type(s) for /: Unit and ' + str(v))
            
        if not v.linear:
            return v.__rtruediv__(self)
            
        vi = [(e[0],-e[1]) for e in v.units]
        return _CompositeUnit(self.units + vi)

    def __rtruediv__(self,v):
        if v == 1:
            return _CompositeUnit([(e[0],-e[1]) for e in self.units])
        
        raise TypeError('unsupported operand type(s) for /: Unit and ' + str(v))
        
        #return Quantity._make(v,unit=self**-1)
        
    def __pow__(self,v):
        if v == -1:
            return _CompositeUnit([(e[0],-e[1]) for e in self.units])
        
        if v == 0:
            return one
        
        if v == 1:
            return self
        
        if v != float(v):
            raise TypeError('unsupported operand type(s) for **: Unit and ' + str(v))

        return _CompositeUnit([(e[0],v*e[1]) for e in self.units])
        
    
class _CompositeUnit(Unit):
    # Instances of this class represent composite derived unit.
    
    living_units = weakref.WeakValueDictionary()
    
    parent = None
    
    @staticmethod
    def _unicode_super(e):
        if int(e) == e:
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
        return '**' + _CompositeUnit._format_sscript(e,parenth=True)
    
    @staticmethod
    def _format_sscript(e,parenth=False):
        if isinstance(e,Rational) and not isinstance(e,Integral):
            ret = str(e.numerator) + '/' + str(e.denominator)
            if parenth:
                ret = '(' + ret + ')'
            return ret
        mf = np.modf(e)
        if mf[0] == 0:
            return '{:.0f}'.format(e)
        if mf[0] == 0.5:
            return '{:.0f}'.format(2*e) + '/2'
        if mf[0] == 0.25:
            return '{:.0f}'.format(4*e) + '/4'
        return str(e).rstrip('0.')
        
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
                txt += '**' + _CompositeUnit._format_sscript(e,parenth=True)
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
        
        # ul is a list of tupels in each tuple a (non _CompositeUnit) Unit
        # instance followed by an exponent.
        ul = [u for u in ul if u[0] is not one and u[1] != 0]
        
        if len(ul) == 0:
            return one
        
        if len(ul) == 1 and ul[0][1] == 1:
            return ul[0][0]
        
        ul = sorted(ul,key = lambda x: (x[0].order,-x[1]))
        
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
        
        ret._linear = linear
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
        
    def tostring(self,fmt=None,solidus=None,mulsep=None,**kwds):
        """
        Returns a string containing the symbol for the unit the format given by
        the keyword fmt which may be set to a string the values 'html', 'latex', 
        'ascii', or 'unicode'.
        
        If solidus is `True` a slash will be used instead of negative exponents
        and if `mulsep` is `True` a mid-dot or asterix will be used instead of a space
        to separate units.
        """
        if solidus is None:
            solidus = self.solidus
            
        if mulsep is None:
            mulsep = self.mulsep
            
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
        cg = None
        cgun = {}
        for un,e in self._units:
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
                    cg = cg.chain(p)
                for t in p.unit.units:
                    if t[0] in cgun:
                        cgun[t[0]] += t[1]
                    else:
                        cgun[t[0]] = t[1]
                        
        cgun = [(k,v) for k,v in cgun.items() if v != 0]
        
        cgun = _CompositeUnit(cgun)
        if cgun is not self:
            if cg is None:
                cg = Conversion(cgun,1)
            else:
                cg._unit = cgun
            cg.parent = self
            self._gconv = self._conversion = cg
        else:
            self._gconv = self._conversion = None

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
    def aliases(self):
        """read-only
        
        Returns a set of the unshadowed aliases of this unit.  To add aliases
        to the unit use the `Unit.alias` static method.
        """
        return {a for a in self._aliases if (Unit.unit('[' + a + ']',exception=False) is self and a != self.name)}
    
    @property
    def shadowed_aliases(self):
        """
        Returns a set of the shadowed aliases of this unit.
        """
        sha = {a for a in self._aliases if Unit.unit('[' + a + ']',exception=False) is not self}
        return sha
    
class _One(Unit):
    """
    The only instance of this class should be the unit `one`.
    """
    def __mul__(self,v):
        if isinstance(v,Unit) or v == 1:            
            return v
        raise TypeError('a unit only be multipled by 1 or another unit')
        #return Quantity._make(v)
        
    def __truediv__(self,v):
        if isinstance(v,Unit):            
            return v**-1
        elif v == 1:
            return self
        raise TypeError('a unit only by divided by 1 or another unit')
        #return Quantity._make(1/v)

    def __rtruediv__(self,v):
        if isinstance(v,Unit):            
            return v
        elif v == 1:
            return self
        raise TypeError('a unit divide 1 or another unit')
        #return Quantity._make(v)
        
    def __pow__(self,v):
        return self
    
    def _cpow(self,v):
        return ([[self,1]],1)
    
    def _repr_html_(self):
        return 'one'
            
    def _repr_latex_(self):
        return 'one'
    
    def __str__(self):
        return 'one'
            
    def __repr__(self):
        return 'one'
    
    def __eq__(self,x):
        if x is self:
            return True
        return x == 1
    
    @property
    def is_dimensionless(self):
        return True
    
with Unit._builtin():
    one = _One('1','',add_symbol=False)

