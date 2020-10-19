# -*- coding: utf-8 -*-

# module unit

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
Classes that for handling units are defined here:  Conversion, Unit,
_CompositeUnit, one, Quantity and QuantityArray
"""

import numpy as np
import weakref
import warnings
from .exceptions import (UnitLibError,CircularUnitConversionError,
                         NoUnitConversionFoundError,UnitNotFoundError,
                         IncompatibleUnitsError,UnitLibNotFoundError,
                         UnitWarning)
from .printing import PrettyPrinter
from .indexed import Indexed
from .unitparser import _UnitParser

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
    if isinstance(x,Quantity) and x.unit:
        if isinstance(y,Quantity):
            yunit = y.unit
            y = y.value
        else:
            yunit = one
        return f(x.unit,x.value,yunit,y,x.autoconvert)
    else:
        if isinstance(x,Quantity):
            xunit = x.unit
            x = x.value
        else:
            xunit = one
        return rf(y.unit,y.value,xunit,x,y.autoconvert)
    
def _f_ratio(f,x,y):
    if not isinstance(y,Quantity):
        x = x.convert(one).value
    elif not isinstance(x,Quantity):
        y = y.convert(one).value
    else:
        y = y.convert(x.unit).value
        x = x.value
    return (f(x,y),one)

class Unit(PrettyPrinter,Indexed):
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
            raise NoUnitConversionFoundError('no conversion found from unit ' + self.short_name + ' to unit ' + unm)
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
        except NoUnitConversionFoundError:
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
        except NoUnitConversionFoundError:
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
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
                
        un,c = self._cpow(b)
        un = _CompositeUnit(un)
        return ((a**b)*c,un)
        #bb = float(b)
        #if bb != b:
            #try:
                #a = self.convert(a,one)
            #except NoUnitConversionFoundError:
                #raise IncompatibleUnitsError('a gummy that is not dimensonless may not be raised to a power with an uncertainty')
            #return (a**b,one)
        #b = bb
        #if int(b) == b:
            #b = int(b)
        
        #if b == 0:
            #return (a**0,one)
        
        #un,c = self._cpow(b)
        #un = _CompositeUnit(un)
        
        #return ((a**b)*c,un)
    
    def _rpow(self,a,bunit,b,aconv):
        if isinstance(b,Unit):
            bunit = b
            b = 1
            
        if self is not one:
            try:
                a = self.convert(a,one)
            except NoUnitConversionFoundError:
                raise IncompatibleUnitsError('the exponent is not dimensionless')
        un,c = bunit._cpow(a)
        un = _CompositeUnit(un)
        return ((b**a)*c,un)
        #aa = float(a)
        #if aa != a:
            #try:
                #b = self.convert(b,one)
            #except NoUnitConversionFoundError:
                #raise IncompatibleUnitsError('a gummy that is not dimensonless may not be raised to a power with an uncertainty')
            #return (b**a,one)
        #a = float(a)
        #if int(a) == a:
            #va = int(a)
        
        #if a == 0:
            #return (b**0,one)
        
        #un,c = _CompositeUnit(bunit._cpow(a))
        
        #return ((b**a)*c,un)
    
    def _ufunc(self,func,*args,**kwds):
        if self.linear:
            try:
                nl = next(a for a in args 
                          if isinstance(a,Quantity) and not a.unit.linear)
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
        if not isinstance(v,Unit):
            return Quantity._make(v,unit=self)
        
        if not v.linear:
            return v.__rmul__(self)
        
        if v is one:
            return self
        
        return _CompositeUnit(self.units + v.units)
    
    def __rmul__(self,v):
        if not isinstance(v,Unit):
            return Quantity._make(v,unit=self)
        
        if not v.linear:
            return v.__mul__(self)
        
        if v is one:
            return self
        
        return _CompositeUnit(self.units + v.units)
        
    def __truediv__(self,v):
        if v is one:
            return self
        
        if not isinstance(v,Unit):
            return Quantity._make(1/v,unit=self)
            
        if not v.linear:
            return v.__rtruediv__(self)
            
        vi = [(e[0],-e[1]) for e in v.units]
        return _CompositeUnit(self.units + vi)

    def __rtruediv__(self,v):
        return Quantity._make(1/v,unit=self)
        
        if v is one or v == 1:
            return _CompositeUnit([(e[0],-e[1]) for e in self.units])
        
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
        
    def _convtree(self):
        if not self._have_conversion:
            self._find_conv()
        Unit._convtree(self)
        
    def _convs(self,d,c,h):
        if not self._have_conversion:
            self._find_conv()
        Unit._convs(self,d,c,h)

    
class _One(Unit):
    """
    The only instance of this class should be the unit `one`.
    """
    def __mul__(self,v):
        if isinstance(v,Unit):            
            return v
        return Quantity._make(v)
        
    def __truediv__(self,v):
        if isinstance(v,Unit):            
            return v**-1
        return Quantity._make(1/v)

    def __rtruediv__(self,v):
        if isinstance(v,Unit):            
            return v
        return Quantity._make(v)
        
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
    
    @property
    def is_dimensionless(self):
        return True
    
with Unit._builtin():
    one = _One('1','',add_symbol=False)
    

class Quantity(PrettyPrinter):
    """
    Instances of this class represent a quantity with a value and a unit.
    The behavior of Quantity instances under mathematical operations with
    other Quanitity object or numerical values depends on the unit. E.g. 
    an interger of float may be added to `Quantity(1,unit='%')` but not to 
    `Quantity(1,unit='m/s')`.  For operations involving only linear units, the
    units will be automatically converted to facilitate the operation, e.g.
    `Quantity(1,unit='psi')` may be added to `Quantity(1,unit='psi')` but not
    `Quantity(1,unit='db(uPa)'.  Manual unit conversions can be realized by
    calling the `Quantity.convert` method, or in place by setting the 
    `Quantity.unit` property.
    
    Quantity instances may be created directly or by multiplying or dividing 
    a number by a `Unit` instance: `Quantity(2.2,unit='cm')` is equivalent to
    `2.2 * unit('cm')`.

    Parameters
    ----------

    value: number like including `ummy`
        the value of the Quantity

    unit:  `str` or `Unit`
        The `Unit` instance representing the unit of the Quantity or a string
        that references the `Unit` instance.
    """
    
    splonk_func_ret = False
    
    _arraytype = None
    
    @classmethod
    def _make(cls,value,unit=one):
        if (isinstance(value,np.ndarray) or isinstance(value,list) or
            isinstance(value,tuple)) and cls._arraytype is not None:
            return cls._arraytype(value,unit=unit)
        return cls(value,unit=unit)
            
    
    @staticmethod
    def _add_unit_sp(fmt,unit):
        if unit is None or unit is one:
            return ''
            
        if unit == '':
            return ''
            
        if unit.startswith('\t'):
            unit = unit[1:]
        else:
            if fmt == 'latex':
                unit = r'\:' + unit
            elif fmt == 'html':
                unit = '&nbsp;' + unit
            else:
                unit = ' ' + unit
            
        return unit
    
    def __init__(self,value,unit=one):
        self._value = value
        if isinstance(unit,Unit):
            self._unit = unit
        else:
            self._unit = Unit.unit(unit)
        self._old = None
        self.autoconvert = False
        
    @property
    def value(self):
        return self._value
           
    @property
    def unit(self):
        """
        Gets or sets the unit for the Quantity.
        
        If this property is set, a unit conversion will be performed.  The value 
        it is set to may be a string, `None`, a `Unit` object, or the integer 1.
        Both 1 and `None` will be interpreted as the Unit instance `one`. A
        `NoUnitConversionFoundError` will be raised if the unit conversion is
        not possible.
        
        Example
        -------
            
        >>> x = Quantity(0.001,unit='V')
        >>> x
        0.001 V
        >>> x.unit = 'uV'
        >>> x
        1000.0 uV
        """
        return self._unit
    @unit.setter
    def unit(self,unit):     
        #Always convert back to the original units so we don't accumulate
        #rounding errors.
        if self._old is None:
            self._old = (self.value,self._unit)
        else:
            self.value,self._unit = self._old
            
        unit = Unit.unit(unit)
        value = self._unit.convert(self.value,unit)
        
        self._value = value
        self._unit = unit
        return
        
    def convert(self,unit):
        """
        Returns a copy of the Quantity with converted units.  A 
        `NoUnitConversionFoundError` will be raised if the unit conversion is
        not possible.
        
        
        Parameters
        ----------
        unit:  `str` or `Unit`
            The unit for the `x` value and if `uunit` is `None`, the
            uncertainty.  It must be string, None, a `Unit` object, or the
            integer 1.  Both 1 and `None` will be interpreted as the Unit
            instance `one`.

        uunit `str`, `Unit` or None, optional
            The unit for the uncertainty `U`.  If this is `None` then `U`
            will have the same units as `x`.  The default is `None`.
        """
        unit = Unit.unit(unit)
        value = self.unit.convert(self.value,unit)
        return type(self)(value,unit=unit)
        
    def reduce_unit(self):
        """
        Cancels factors in a Quantity's unit when possible.  This modifies the
        calling gummy and returns `None`.
        
        Example
        -------
        
        >>> x = Quantity(5,unit='mm/m')
        >>> x.reduce_unit()
        >>> x
        0.005
        """
        un = self._unit.mulr(one)[0]
        self._unit = un
    
    @property
    def c(self):
        """
        This read-only property is used as a conversion flag during calculations.
        When an arithmetic operation is carried out between two Quantaties with
        different units, a unit conversion on one of the input quantities may be
        required to complete the calculation.  Attach this flag to the unit that
        you prefer be converted.

        Examples
        --------
        
        >>> a = Quantity(1,unit='cm')
        >>> b = Quantity(2,unit='mm')
        >>> a + b
        1.2 cm
        >>> a.c + b
        12 mm
        >>> a + b.c
        1.2 cm
        >>> a*b
        0.2 cm**2
        >>>a.c*b
        20 mm**2
        """
        c = Quantity(self.value,self.unit)
        c.autoconvert = True
        return c
 
    def tostring(self,fmt=None,**kwds):
        unit = self._unit.tostring(fmt=fmt,**kwds)
        unit = self._add_unit_sp(fmt,unit)
        if isinstance(self.value,PrettyPrinter):
            value = self.value.tostring(fmt=fmt,**kwds)
        else:
            value = str(self.value)
        return value + unit
    
    def copy(self,tofloat=False):
        if tofloat:
            try:
                return type(self)(self.value.tofloat(),unit=self.unit)
            except:
                return type(self)(float(self._value),unit=self._unit)
            
        try:
            return type(self)(self.value.copy(),unit=self.unit)
        except:
            return type(self)(self.value,self.unit)
        
    def tofloat(self):
        return self.copy(tofloat=True)
    
    def totuple(self):
        return (self.value,self.unit)
    
    def splonk(self):
        """
        returns self.value if self.unit is one else returns self
        """
        if self.unit is one:
            try:
                return self.value.splonk()
            except AttributeError:
                return self.value
        return self
        
    def _bop(self,v,f):
        if issubclass(type(v),type(self)):
            make =  v._make
        else:
            make =  self._make
            
        if isinstance(v,Quantity):
            vunit = v._unit
            v = v.value
            aconv = self.autoconvert
        else:
            vunit = one
            aconv = False
            
        r,runit = f(self.value,vunit,v,aconv)
            
        return make(r,unit=runit)
        
    def __add__(self, v):
        return self._bop(v,self.unit._add)
                
    def __radd__(self, v):
        return self._bop(v,self.unit._radd)
    
    def __sub__(self, v):
        return self._bop(v,self.unit._sub)
                
    def __rsub__(self, v):
        return self._bop(v,self.unit._rsub)
    
    def __mul__(self, v):
        return self._bop(v,self.unit._mul)
                
    def __rmul__(self, v):
        return self._bop(v,self.unit._rmul)
    
    def __truediv__(self, v):
        return self._bop(v,self.unit._truediv)
                
    def __rtruediv__(self, v):
        return self._bop(v,self.unit._rtruediv)
    
    def __floordiv__(self, v):
        return self._bop(v,self.unit._floordiv)
                
    def __rfloordiv__(self, v):
        return self._bop(v,self.unit._rfloordiv)
    
    def __pow__(self, v):
        return self._bop(v,self.unit._pow)
                
    def __rpow__(self, v):
        return self._bop(v,self.unit._rpow)
    
    def __mod__(self, v):
        return self._bop(v,self.unit._mod)
                
    def __rmod__(self, v):
        return self._bop(v,self.unit._rmod)                
        
    def __neg__(self):
        r,runit = self.unit._neg(self.value)
        return self._make(r,unit=runit)
        
    def __pos__(self):
        r,runit = self.unit._pos(self.value)
        return self._make(r,unit=runit)
        
    def __abs__(self):
        r,runit = self.unit._abs(self.value)
        return self._make(r,unit=runit)
    
    def _cmp(self,v,f):
        if isinstance(v,Quantity):
            if self._unit is not v.unit:
                try:
                    v = v.convert(self.unit)
                    return f(self.value,v.value)
                except:
                    raise IncompatibleUnitsError('Quantities with incompatible units cannot be compared ')
        else:
            try:
                s = self.convert(one)
                return f(s.value,v)
            except:
                raise IncompatibleUnitsError('Quantities with incompatible units cannot be compared')
        
    def __eq__(self, v):
        return self._cmp(v,lambda x,y: x == y)
    
    def __ne__(self, v):
        return self._cmp(v,lambda x,y: x != y)
        
    def __lt__(self, v):
        return self._cmp(v,lambda x,y: x < y)
    
    def __le__(self, v):
       return self._cmp(v,lambda x,y: x <= y)
        
    def __gt__(self, v):
        return self._cmp(v,lambda x,y: x > y)
        
    def __ge__(self, v):
        return self._cmp(v,lambda x,y: x >= y)
    
    def _ufunc(self,func,*args,**kwds):
        # handles numpy functions applied to Quantity arguments
        
        s = self
        for a in args:
            if issubclass(type(a),type(s)):
                s = a
        make = s._make
        
        try:
            x,xunit = self.unit._ufunc(func,*args,**kwds)
            return make(x,unit=xunit)
        except NotImplementedError:
            try:
                x = [a if not isinstance(a,Quantity) else a.convert(one).value 
                     for a in args]
                if self.splonk_func_ret:
                    return func(*x,**kwds)
                return make(func(*x,**kwds))
            except NoUnitConversionFoundError:
                raise NotImplementedError()
      
    def __array_ufunc__(self,ufunc,method,*args,**kwds):
        # handles numpy ufunc's applied to Quantity arguments
        if method != '__call__':
            return None
        
        return self._ufunc(ufunc,*args,**kwds)
    
    def __array_function__(self,func,method,*args,**kwds):        
        return self._ufunc(func,*args,**kwds)
    
    #def __float__(self):
        #return float(self.value)
    
    #def __int__(self):
        #return int(self.value)

    #def __complex__(self):
        #return complex(self.value)
    
    @property
    def real(self):
        return self.copy()
    
    @property
    def imag(self):
        return type(self)(0,unit=self.unit)


class QuantityArray(Quantity):
    """
    A subclass of Quantity.  Instance of this class represent a list, tuple, 
    or numpy array of values all with the same unit.  Elements of the array 
    are returned as `Quantity` instances.  Instances of this class can be
    created directly or by multiplying a list, tuple, or numpy array by a
    `Unit` instance.
    
    Parameters
    ----------

    value: `list`, `tuple` or `ndarray`
        the value of the Quantity

    unit:  `str` or `Unit`
        The `Unit` instance representing the unit of the Quantity or a string
        that references the `Unit` instance.
    """
    
    _elementtype = Quantity
    
    def __len__(self):
        return len(self.value)
    
    def __getitem__(self,key):
        return self._elementtype(self.value[key],unit=self.unit)

    def __contains__(self,quantity):
        if isinstance(quantity,self._elementtype):
            if self._unit is not quantity.unit:
                try:
                    quantity = quantity.convert(self.unit)
                    return quantity.value in self.value
                except:
                    return False
        else:
            try:
                s = self.convert(one)
                return quantity in s.value
            except:
                return False
            
    def __iter__(self):
        for e in self.value:
            yield self._elementtype(e,unit=self.unit)
            
    def __array__(self):
        return type(self)(np.asarray(self.value),unit=self.unit)
    
    
Quantity._arraytype = QuantityArray
        

