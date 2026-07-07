# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:32:16 2026

@author: Parksh
"""
import numpy as np
from abc import ABCMeta

from .printing import PrettyPrinter,MetaPrettyPrinter
from .abc import AbcQuantity,AbcQuantityArray,AbcUnit
from .exceptions import IncompatibleUnitsError
from .util import _isscalar

class MetaQuantity(MetaPrettyPrinter,ABCMeta):
    pass

class Quantity(PrettyPrinter,AbcQuantity,metaclass=MetaQuantity):
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

    value: numeric (including `ummy`)
        the value of the Quantity

    unit:  `str` or `Unit`
        The `Unit` instance representing the unit of the Quantity or a string
        that references the `Unit` instance.
    """
    
    splonk_func_ret = False
    
    _arraytype = None
    
    @classmethod
    def _make(cls,value,unit=1):
        if _isscalar(value):
            return cls(value,unit=unit)
        if cls._arraytype is None:
            if isinstance(value,np.ndarray):
                return np.array([cls._make(v,unit=unit) for v in value])
            return tuple(cls._make(v,unit=unit) for v in value)
        else:
            return cls._arraytype(value,unit=unit)
    
    @staticmethod
    def _add_unit_sp(fmt,unit):
        if unit is None or unit == 1:
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
    
    def __init__(self,value,unit=1):
        try:
            hash(value)
        except TypeError:
            raise TypeError('value is unhashable')
        self._value = value
        if unit == 1:
            self._unit = 1
        else:
            from ._unit import Unit
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
        if not isinstance(self._unit,AbcUnit):
            from ._unit import Unit
            self._unit = Unit.unit(self._unit)
        return self._unit
    @unit.setter
    def unit(self,unit):     
        #Always convert back to the original units so we don't accumulate
        #rounding errors.
        if self.unit == 1 and self._unit == 1:
            return
        
        from ._unit import Unit
        if self._old is None:
            self._old = (self.value,self.unit)
        else:
            self._value,self._unit = self._old
            
        unit = Unit.unit(unit)
        value = self._unit.convert(self.value,unit)
        
        self._value = value
        self._unit = unit
        return
    
    @property
    def unit_is_one(self):
        return self._unit == 1
        
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
        if unit == 1 and self._unit == 1:
            return type(self)(self.value)
        
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
        from ._unit import one
        un,f = self._unit._mul_cancel(one)
        self._unit = un
        self._value *= f
    
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
        """
        tostring(fmt='unicode')
        
        returns a string representation of the Quantity.  fmt may be "unicode",
        "html","latex" or "ascii".  The default is "unicode".
        """
        if self.unit_is_one:
            unit = ''
        else:
            unit = self._unit.tostring(fmt=fmt,**kwds)
            unit = self._add_unit_sp(fmt,unit)
        if isinstance(self.value,PrettyPrinter):
            value = self.value.tostring(fmt=fmt,**kwds)
        else:
            value = str(self.value)
        return value + unit
    
    def copy(self,totype=None):
        """
        copy(tofloat=False)
        
        returns a copy of self.  If tofloat is True, the self.value will be
        converted to float.  The default is False.
        """
        if totype is float:
            try:
                return type(self)(self.value.tofloat(),unit=self.unit)
            except:
                return type(self)(float(self._value),unit=self._unit)
        elif totype is not None:
            return type(self)(totype(self._value),unit=self._unit)
            
        try:
            return type(self)(self.value.copy(),unit=self.unit)
        except:
            return type(self)(self.value,self.unit)
        
    def tofloat(self):
        """
        returns a copy of self with value float(self.value) equivalent to 
        copy(tofloat=True)
        """
        return self.copy(totype=float)
    
    def totuple(self):
        """
        returns the tuple (self.value,self.unit)
        """
        return (self.value,self.unit)
    
    def tobaseunit(self):
        """
        Returns a Quantity equal to self converted to unit self.unit.base
        """
        return self.convert(self.unit.base)
    
    def splonk(self):
        """
        returns self.value if self.unit is one else returns self
        """
        if self.unit == 1:
            try:
                return self.value.splonk()
            except AttributeError:
                return self.value
        return self
        
    def _bop(self,v,f):
        from ._unit import one
        if issubclass(type(v),type(self)):
            make =  v._make
        else:
            make =  self._make
            
        if isinstance(v,AbcQuantity):
            vunit = v.unit
            v = v.value
            aconv = self.autoconvert
        else:
            vunit = one
            aconv = False
            
        r,runit = f(self.value,vunit,v,aconv)
            
        return make(r,unit=runit)
        
    def __add__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value + v)
            if v.unit_is_one:
                return type(self)(self.value + v.value)
            
        return self._bop(v,self.unit._add)
                
    def __radd__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v + self.value)
            if v.unit_is_one:
                return type(self)(v.value + self.value)
            
        return self._bop(v,self.unit._radd)
    
    def __sub__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value - v)
            if v.unit_is_one:
                return type(self)(self.value - v.value)
            
        return self._bop(v,self.unit._sub)
                
    def __rsub__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v - self.value)
            if v.unit_is_one:
                return type(self)(v.value - self.value)
            
        return self._bop(v,self.unit._rsub)
    
    def __mul__(self, v):
        if isinstance(v,AbcUnit):
            x,u = v._rmul(1,self.unit,self.value,self.autoconvert)
            return type(self)(x,unit=u)
        
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value*v)
            if v.unit_is_one:
                return type(self)(self.value*v.value)
            
        return self._bop(v,self.unit._mul)
                
    def __rmul__(self, v):
        if isinstance(v,AbcUnit):
            x,u = v._mul(1,self.unit,self.value,self.autoconvert)
            return type(self)(x,unit=u)
        
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v*self.value)
            if v.unit_is_one:
                return type(self)(v.value*self.value)
            
        return self._bop(v,self.unit._rmul)
    
    def __truediv__(self, v):
        if isinstance(v,AbcUnit):
            x,u = v._rtruediv(1,self.unit,self.value,self.autoconvert)
            return type(self)(x,unit=u)
        
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value/v)
            if v.unit_is_one:
                return type(self)(self.value/v.value)
            
        return self._bop(v,self.unit._truediv)
                
    def __rtruediv__(self, v):
        if isinstance(v,AbcUnit):
            x,u = v._truediv(1,self.unit,self.value,self.autoconvert)
            return type(self)(x,unit=u)
        
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v/self.value)
            if v.unit_is_one:
                return type(self)(v.value/self.value)
            
        return self._bop(v,self.unit._rtruediv)
    
    def __floordiv__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value//v)
            if v.unit_is_one:
                return type(self)(v.value//self.value)
            
        return self._bop(v,self.unit._floordiv)
                
    def __rfloordiv__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v//self.value)
            if v.unit_is_one:
                return type(self)(v.value//self.value)
            
        return self._bop(v,self.unit._rfloordiv)
    
    def __pow__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value**v)
            if v.unit_is_one:
                return type(self)(self.value**v.value)
            
        return self._bop(v,self.unit._pow)
                
    def __rpow__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v**self.value)
            if v.unit_is_one:
                return type(self)(v.value**self.value)
            
        return self._bop(v,self.unit._rpow)
    
    def __mod__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(self.value%v)
            if v.unit_is_one:
                return type(self)(self.value%v.value)
            
        return self._bop(v,self.unit._mod)
                
    def __rmod__(self, v):
        if self.unit_is_one:
            if not isinstance(v,AbcQuantity):
                return type(self)(v%self.value)
            if v.unit_is_one:
                return type(self)(v.value%self.value)
            
        return self._bop(v,self.unit._rmod)                
        
    def __neg__(self):
        if self.unit_is_one:
            return type(self)(-self.value)
        
        r,runit = self.unit._neg(self.value)
        return self._make(r,unit=runit)
        
    def __pos__(self):
        if self.unit_is_one:
            return type(self)(self.value)
        
        r,runit = self.unit._pos(self.value)
        return self._make(r,unit=runit)
        
    def __abs__(self):
        if self.unit_is_one:
            return type(self)(np.abs(self.value))
        
        r,runit = self.unit._abs(self.value)
        return self._make(r,unit=runit)
    
    def _cmp(self,v,f):
        if isinstance(v,Quantity):
            if self._unit is not v.unit:
                try:
                    v = v.convert(self.unit)
                except IncompatibleUnitsError:
                    raise IncompatibleUnitsError('Quantities with incompatible units cannot be ordered')
            return f(self.value,v.value)
        else:
            try:
                s = self.convert(1)
                return f(s.value,v)
            except IncompatibleUnitsError:
                raise IncompatibleUnitsError('Quantities with incompatible units cannot be ordered')
        
    def __eq__(self, v):
        try:
            return self._cmp(v,lambda x,y: x == y)
        except IncompatibleUnitsError:
            return False
    
    def __ne__(self, v):
        try:
            return self._cmp(v,lambda x,y: x != y)
        except IncompatibleUnitsError:
            return True
        
    def __lt__(self, v):
        return self._cmp(v,lambda x,y: x < y)
    
    def __le__(self, v):
       return self._cmp(v,lambda x,y: x <= y)
        
    def __gt__(self, v):
        return self._cmp(v,lambda x,y: x > y)
        
    def __ge__(self, v):
        return self._cmp(v,lambda x,y: x >= y)
    
    def _ufunc(self,func,*args,**kwds):
        b = np.broadcast(*args)
        if b.shape == ():
            return self._iufunc(func,*args,**kwds)
            
        ret = np.array([self._iufunc(func,*a,**kwds) for a in b])
        
        shape = b.shape
        if isinstance(ret[0],np.ndarray):
            shape += ret[0].shape
        
        return np.reshape(ret,shape)
    
    def _iufunc(self,func,*args,**kwds):
        # handles numpy functions applied to Quantity arguments
        
        s = self
        af = list(args)
        uone = True
        for i,a in enumerate(args):
            if issubclass(type(a),type(s)):
                s = a
            if isinstance(a,Quantity):
                af[i] = a.value
                if not a.unit_is_one:
                    uone = False
        make = s._make
        
        if uone:
            return make(func(*af,**kwds))
        
        try:
            x,xunit = self.unit._ufunc(func,*args,**kwds)
            return make(x,unit=xunit)
        except NotImplementedError:
            try:
                x = [a if not isinstance(a,Quantity) else a.convert(1).value 
                     for a in args]
                if self.splonk_func_ret:
                    return func(*x,**kwds)
                return make(func(*x,**kwds))
            except IncompatibleUnitsError:
                raise NotImplementedError()
      
    def __array_ufunc__(self,ufunc,method,*args,**kwds):
        # handles numpy ufunc's applied to Quantity arguments
        if method != '__call__':
            return None
        
        return self._ufunc(ufunc,*args,**kwds)
    
    def __array_function__(self,func,method,*args,**kwds):   
        #return self._ufunc(func,*args,**kwds)
        return self._ufunc(func,*args[0],**args[1])
    
    def __float__(self):
        s = self.convert(1)
        return float(s.value)

    def __complex__(self):
        s = self.convert(1)
        return complex(s.value)
    
    def __bool__(self):
        return self != 0
    
    @property
    def real(self):
        if self.unit.linear:
            s = self
        else:
            s = self.tobaseunit()
         
        try:
            v = s.value.real
        except:
            v = float(s.value)
        return type(self)(v,unit=s.unit)
     
    @property
    def imag(self):
        if self.unit.linear:
            s = self
        else:
            s = self.tobaseunit()
        
        try:
            v = s.value.imag
        except:
            v = 0
        return type(self)(v,unit=s.unit)
    
    def __hash__(self):
        r = self.tobaseunit().totuple()
        if r[1] == 1:
            return hash(r[0])
        return hash(r)
    
    def conjugate(self):
        try:
            v = self.value.conjugate()
        except:
            v = complex(float(self.value))
        return type(self)(v,unit=self.unit)
    
class QuantityArray(Quantity,AbcQuantityArray):
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
                s = self.convert(1)
                return quantity in s.value
            except:
                return False
            
    def __iter__(self):
        for e in self.value:
            yield self._elementtype(e,unit=self.unit)
            
    def __array__(self):
        return type(self)(np.asarray(self.value),unit=self.unit)
    
    
Quantity._arraytype = QuantityArray

        
