# -*- coding: utf-8 -*-

# module gummy

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

"""
The gummy object defined here represents the core functionality of this 
package, though the code in this module only integrates the unit conversion 
machinery from the .unit module with the other functions and provides pretty
printing functionality.  The first-order uncertainty propagation code resides 
in the .ummy module.  The Monte-Carlo propagation machinery is defined in the
.distributions module and is integrated with the ummy object in the nummy
module.  The gummy object, in turn, inherits from the nummy object.
"""

import numpy as np
from .ummy import ummy,immy,_isscalar,_format_exp,_to_decimal,_decimal_str
from .nummy import nummy,get_name
from .exceptions import IncompatibleUnitsError
from .unit import Unit,one,Quantity,MetaQuantity,MFraction
from .distributions import Distribution,MultivariateDistribution
from .pmethod import _Pmthd,loc_from_k
from math import isnan, isinf
from numbers import Integral,Rational,Real,Complex
from decimal import Decimal,localcontext,InvalidOperation


try:
    import mpmath as mp
except:
    mp = None

def _ku(k,u):
    try:
        return k*u
    except:
        try:
            if isinstance(u,Integral):
                raise
            return type(u)(k)*u # in case gummy.u is a decimal.Decimal 
        except:
            return float(k)*float(u)
        
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
    
def _p_to_str(fmt,p):
    if fmt == 'latex':
        pct = '\\%'
    else:
        pct = '%'
        
    x = 100*Decimal(p)
    if x <= 90:
        return '{:.0f}'.format(x) + pct
    
    y = 100 - x
    y = round(y,1-y.adjusted()-1)
    y = round(y,1-y.adjusted()-1)
    x = x.quantize(y)
    return '{:f}'.format(x) + pct
    
def _k_to_str(k):
    return '{:.1f}'.format(k)
    
def _dof_to_str(dof,fmt=None):
    if isinf(dof) or dof > 99:
        if fmt == 'html':
            return '&infin;'
        if fmt == 'latex':
            return r'\infty'
        if fmt == 'ascii':
            return 'inf'
        return '\u221E'
        
    if isinstance(dof,Integral):
        return str(dof)
    return '{:.1f}'.format(dof)
    
def _set_covariance_matrix(gummys, matrix):
    nummys = [g.value for g in gummys]
    nummy._set_covariance_matrix(nummys, matrix)
    for g in gummys:
        g._set_U(None,None)
        
def _set_correlation_matrix(gummys, matrix):
    nummys = [g.value for g in gummys]
    nummy._set_correlation_matrix(nummys, matrix)
        

class MetaGummy(MetaQuantity):
    # A metaclass to define some "classproperties" for gummy
    
    @property
    def cimethod(cls):
        """
        str in {'shortest', 'symmetric'}
        
        Get or set the method for calculating the confidence interval from 
        Monte-Carlo data.  If this property is set at the class level, it will
        change the default `cimethod` value for new gummys but will not affect
        gummys that have already been created.
        
        Can be set either to the string 'shortest' or the string 'symmetric'.
        This property gets or sets the method for calculating confidence
        intervals from Monte-Carlo data.  
        
        If it is set to 'shortest', the confidence interval will be taken to be 
        the shortest interval that includes the desired fraction of the probability 
        distribution.  
        
        If it is set to 'symmetric', then the confidence interval will be set so 
        that, for n Monte-Carlo samples and a coverage probability of `p`, then
        `n`*(1-`p`)/2 samples lie below the lower limit of the confidence interval
        and the same number of samples lie above the upper limit of the confidence 
        interval.
        """
        return nummy._cimethod
    @cimethod.setter
    def cimethod(cls,value):
        value = value.lower().strip()
        if value not in ['shortest','symmetric']:
            raise ValueError('cimethod ' + str(value) + ' is not recognized')
        nummy._cimethod = value
        
    @property
    def bayesian(cls):
        """
        `bool`
        
        Read/write at the class level, but read-only at the instance level.
        The default value is `False`; this should only be changed once at the
        beginning of the session.  This property affects how the level of 
        confidence `p` (sometimes called coverage probability) of an expanded
        uncertainty is related to the coverage factor `k` for a gummy based on
        data with finite degrees of freedom.

        Standard uncertainties are often based on the standard deviation of a set
        of measurements (and the assumption that these measurements are drawn
        from a normally distributed population).  Traditionally (e.g. the GUM
        2008 edition) the standard uncertainty is taken to be the standard
        deviation of the mean (s/sqrt(n), where s is the sample standard deviation
        and n is the number of measurements).  However there is some "extra
        uncertainty" because the sample standard devation not exactly equal to
        the population standard deviation.  This is taken into account by using
        a Student's t distribution to calculate the expanded uncertainty.  However
        it has been pointed out, by those who advocate a Bayesian point of view,
        that the probability distribution for the measurand here is best described
        by a shifted and scaled Student's t distribution.  So the standard
        uncertainty should be the standard deviation of this distribution which
        is s*sqrt{(n-1)/[n*(n-3)]}.  Thus

        u(bayesian) = sqrt[dof/(dof - 2)]*u(traditional)

        where dof = n - 1 and the "extra uncertainty" is incorporated directly
        into the standard uncertainty.
        
        Example
        -------
        >>> gummy.bayesian = True
        >>> g = gummy(1,0.03,dof=5)
        >>> g.bayesian
        True
        """
        return nummy._bayesian
    @bayesian.setter
    def bayesian(cls,v):
        nummy._bayesian = bool(v)
        
    @property
    def style(cls):
        """
        Get or set the default display style for new gummys.  This is a string
        with one of the following values:
        
        - "pm" or "+-" gives, e.g. in ascii format "(1.00 +/- 0.12)e-12 cm"
        - "pmi" or "+-i" gives e.g. "1.00e-12 cm +/- 1.2e-13 cm"
        - "concise" or "()" gives e.g. "1.00(12)e-12 cm"
        - "ueq" or "u equals" gives e.g. "1.00e-12 cm with u = 1.2e-13 cm"
        - "x" or "x only" gives e.g. "1.00e-12 cm"
        - "xf" gives e.g. "1.00e-12"
        - "u" or "u only" gives "1.2e-13 cm"
        - "uf" gives "1.2e-13"
        
        The following styles display a mean and confidence interval based on
        data from a Monte-Carlo simulation.  If no simulated data is available
        "no simulated data" is returned.
        
        - "pmsim" gives e.g. "(1.01 + 0.11 - 0.13)e-12 cm"
        - "pmsimi" gives e.g. "1.01e-12 cm + 1.1e-13 cm - 1.3e-11 cm"
        - "mcisim" gives e.g. "mean = 1.01e-12 cm, confidence interval = 
          [8.8e-13 cm, 1.13e-12 cm]"
        - "cisim" gives e.g "[8.8e-13 cm, 1.13e-12 cm]"
        - "usim" gives the standard deviation e.g. "1.2e-13 cm"
        - "ufsim" gives the standard deviation e.g. "1.2e-13"
        
        Note if `uunit` has been defined for the gummy instance, then concise style
        and pm style are not valid and the display will default to pmi style.
        
        The style can be set either at the class or the instance level.
        
        To control whether the coverage factor (k), the level of confidence (p),
        the degrees of freedom (dof), and name are displayed as part of the 
        uncertainty statement set the following variables to True or False (at 
        either the class or instance level):
        
        - `show_k`
        - `show_p`
        - `show_dof`
        - `show_name`
        
        `show_k`, `show_p`, and `show_dof` can also be set to `None`, which will
        allow the gummy object to decide whether to display the corresponding
        values.
        
        The number of significant digits in the uncertainty to be displayed is
        set with the `nsig` attribute.  This can be set either at the class level
        or the instance level.  When the gummy module is loaded `nsig` is set to 2.
        
        Scientific notation can be turned on or off by setting the `sci_notation`
        attribute to `True` or `False`.  If `sci_notation` is set to `None`, scientific notation
        will be used if the exponent of gummy `x` value is greater than than the
        value of the `sci_notation_high` attribute or less than the `sci_notation_low`
        attribute.  When the gummy module is loaded, `sci_notation` is `None`,
        `sci_notation_high` is 7 and `sci_notation_low` is -3.  These attributes
        can be set either at the class or instance level.
        
        Setting the `solidus` attribute to `True` uses a solidus (forward slash)
        to separate units with positive and negative exponents e.g. "m/s" and when
        `solidus` is set to `False` negative exponents are used e.g. "m s**-1".
        Setting the `mulsep` attribute to `True` inserts a dot or * between units while setting
        this to `False` a space separates units. When the gummy module is loaded,
        both `solidus` and `mulsep` are set to `False`.  They may be set at either
        the instance or class level.
        """
        return gummy._style
    @style.setter
    def style(cls,value):
        gummy._style = gummy._get_style(value)
        
    @property
    def cmp_k(cls):
        """
        Get or set the coverage factor for comparisons between gummys.
        Setting this property sets the `cmp_p` property to `None`.
        """
        return gummy._cmp_k
    @cmp_k.setter
    def cmp_k(cls,v):
        if v <= 0:
            raise ValueError('cmp_k <= 0')
        gummy._cmp_k = v
        gummy._cmp_p = None
        
    # The default level of confidence for new ummies
    @property
    def cmp_p(cls):
        """
        Get or set the probability level to use when comparing gummys.
        Setting this property sets the `cmp_k` property to `None`.
        """
        return gummy._cmp_p
    @cmp_p.setter
    def cmp_p(cls,v):
        if v <= 0 or v >= 1:
            raise ValueError('cmp_p is not in the interval (0,1)')
        gummy._cmp_p = v
        gummy._cmp_k = None
        
    @property
    def p_method(cls):
        """`str` in {'loc', 'cp', 'gauss', 'ccp', 'chebyshev'}
        
        This sets the default `p_method` attribute for newly created gummys which
        determines how the coverage factor is calculated from a given level of 
        confidence or coverage probability `p`.
        
        If `p_method` = 'loc', then the uncertainty is assumed to be represented
        by a normal probability distribution if `dof` = `float('inf')` and shifted
        and scaled Student's t distribution otherwise.  If `p_method` = 'gauss' or
        'cp' then the Gauss inequality is used, and if `p_method` = 'chebyshev' or
        'ccp' then the Chebyshev inequality is used.  For `p` = 0.95 and
        `dof` = `float('inf')`, `p_method` = 'loc' gives `k` = 2.0, while
        `p_method` = 'gauss' gives `k` = 3.0 and p_method = 'chebyshev' gives `k` = 4.5.
        """
        return gummy._p_method.method
    @p_method.setter
    def p_method(cls,v):
        if v is None:
            v = 'loc'
        gummy._p_method = _Pmthd(v)
        
    @property
    def max_dof(cls):
        """
        `int`
        
        Gets or sets the maximum finite value for dof.  Any dof larger than
        max_dof will be rounded to float('inf').
        """
        return ummy.max_dof
    @max_dof.setter
    def max_dof(cls,v):
        ummy.max_dof = v
        
    @property
    def nsig(cls):
        """
        `int`
        
        Gets or sets the number of significant digits in the uncertainty to 
        display.
        """
        return ummy.nsig
    @nsig.setter
    def nsig(cls,v):
        ummy.nsig = v
        
    @property
    def thousand_spaces(cls):
        """
        `bool`
        
        Gets or sets a bool value that determines whether to insert a space to
        group digits in x.
        """
        return ummy.thousand_spaces
    @thousand_spaces.setter
    def thousand_spaces(cls,v):
        ummy.thousand_spaces = bool(v)
        
    @property
    def sci_notation(cls):
        """
        `bool` or `None`
        
        Setting this to True or False forces or prevents scientific notation
        in the display of the value.  If this is set to `None` thr scientific
        notation will be used if the x value is below `10**sci_notation_low` or
        above `10**sci_notation_high`.
        """
        return ummy.sci_notation
    @sci_notation.setter
    def sci_notation(cls,v):
        if v is not None:
            v = bool(v)
        ummy.sci_notation = v
        
    @property
    def sci_notation_high(cls):
        """
        see the sci_notation property
        """
        return ummy.sci_notation_high
    @sci_notation_high.setter
    def sci_notation_high(cls,v):
        ummy.sci_notation_high = int(v)
        
    @property
    def sci_notation_low(cls):
        """
        see the sci_notation property
        """
        return ummy.sci_notation_low
    @sci_notation_low.setter
    def sci_notation_low(cls,v):
        ummy.sci_notation_low = int(v)
    
    @property
    def rounding_u(cls):
        """
        `bool`
        
        If this is set to `True` then the uncertainty of the ummy or gummy 
        includes a contribution to account for the finite resolution of the
        x-value.
        """
        return ummy.rounding_u
    @rounding_u.setter
    def rounding_u(cls,v):
        ummy.rounding_u = bool(v)
        
    @property
    def max_digits(cls):
        """
        `int`
        
        Gets or sets the maximum number of digits in x to display.
        """
        return ummy.max_digits
    @max_digits.setter
    def max_digits(cls,v):
        v = int(v)
        if v < 0:
            raise ValueError('max_digits must be >= 0')
        ummy.max_digits = v
        
    @property
    def solidus(cls):
        return Unit.solidus
    @solidus.setter
    def solidus(cls,v):
        Unit.solidus = bool(v)
        
    @property
    def mulsep(cls):
        return Unit.mulsep
    @mulsep.setter
    def mulsep(cls,v):
        Unit.mulsep = bool(v)
        
    
class gummy(Quantity,metaclass=MetaGummy):
    """
    A gummy object represents a numerical value with an uncertainty and (or) a
    unit.  They can be used in place of float values in Python expressions and 
    the uncertainty can be propagated with with both first-order and 
    Monte-Carlo methods.  Correlations between gummys are tracked.
    
    
    Parameters
    ----------
    x:  `float`, `Distribution` or `gummy`
        Either a float representing the value of the gummy or a sub-class of 
        `Distribution` that represents the probability distribution of the gummy.  
        If `x` is a Distribution neither `u` nor `dof` should be specified.    
        If `x` is a gummy then all other parameters are ignored and ``gummy(x)``
        is equivalent to ``x.copy(formatting=True)``.

    u:  `float` >= 0, optional,
        A number representing the uncertainty in `x`.  By default `u` 
        is taken to be the standard ("1-sigma") uncertainty.  But if `k` or `p`
        are specified then `u` is taken to be the corresponding expanded 
        uncertainty. The default value for `u` is 0.

    unit:  `str` or `Unit`, optional
        The units of the gummy.  This may be a Unit object or, more commonly, a 
        string that references a Unit object, e.g. ``gummy(1,unit='kg')`` or
        ``gummy(3.67,0.22,unit='m/s')``.  The default is `one`.

    dof:  `float` or `int` > 0, optional
        The number of degrees of freedom upon which the uncertainty is based.  
        The default is ``float('inf')``.

    k:  float > 0 or None, optional
        The coverage factor if `u` is an expanded uncertainty.  The value of the 
        `u` parameter is divided by the coverage factor to get the  standard 
        uncertainty for the new gummy.  If the parameter `p` is specified then
        the coverage factor is calculated using `p` and the value of the `k` 
        parameter is ignored.  The default for `k` is 1.

    p:  `float` between 0 and 1 or `None`, optional
       The level of confidence or probability  level if `u` is an expanded 
       uncertainty.  If this parameter is specified, then the coverage factor 
       is calculated using `dof` and the method specified
       by the `p_method` parameter.  The standard uncertainty for the gummy is 
       then set to `u` divided by the coverage factor.

    p_method:  {'loc', 'cp', 'gauss', 'ccp', 'chebyshev', `None`}, optional
        If the `p` parameter is specified, `p_method `sets the method that is 
        used to calculate the coverage factor.  If `p_method` is None then the 
        value of gummy.p_method is used (which in turn has a default of 'loc').
        
        If `p_method` = 'loc', then the uncertainty is assumed to be represented
        by a normal probability distribution if `dof` = float('inf') and shifted
        and scaled Student's t distribution otherwise.  If `p_method` = 'gauss' 
        or 'cp' then the Gauss inequality is used, and if `p_method` = 'chebyshev'
        or 'ccp' then the Chebyshev inequality is used.  For `p` = 0.95 and 
        `dof` = float('inf'), `p_method` = 'loc' gives `k` = 2.0, while 
        `p_method` = 'gauss' gives `k` = 3.0 and `p_method` = 'chebyshev' gives 
        `k` = 4.5.

    uunit: `str`, `Unit` or `None`, optional
        This represents the units of `u`.  It may be  a unit with the same 
        dimension as the unit parameter, e.g. a measurement result of 3 m with 
        an uncertainty of 1 mm can be represented by 
        ``gummy(3,1,unit='m',uunit='mm')`` The uunit parameter can also be a
        dimensionless unit if u represents a relative uncertainty, e.g. the gummy
        above can be equivalently represented by ``gummy(3,0.1,unit='m',uunit='%')``.
        If this is set to None, then the units of `u` are taken to be the same as
        those of `x` (as given by the unit parameter).  The default is `None`.
            
    utype:  `str` or `None`, optional
        An arbitrary string value labeling the uncertainty type.  When 
        a calculation is performed with gummys, the combined uncertainty of 
        effective degrees of freedom from one particular uncertainty type can be 
        found in the calculation result with the ufrom and doffrom methods.  E.g.
        you can create a set of gummys with uncertainties assigned either utype 
        "A" or utype "B", insert them into a measurement equation and find the 
        combined utype "A" uncertainty.  The default is `None`.
            
    name:  `str` or `None`, optional
        An arbitrary string naming the gummy.  The name is used when displaying 
        the gummy value and serves no other function.  The default is `None`.
    """
    
    # The following three attributes control whether k, p, and dof is displayed
    # along with the value and uncertainty.  They may be set to None (automatic 
    # based on the k, p, or dof settings), True, or False.
    show_k = None
    show_p = None
    show_dof = None
    show_name = True
    
    _solidus = None 
    _mulsep = None 
    
    slashaxis = True
                    
    _cmp_k = None
    _cmp_p = 0.95
    _ubreakdown = None
    _p_method = _Pmthd('loc')
    _Ubr = None
    
    _style = 'concise' # see the MeteGummy style property
    
    exception_on_fmt_error = False  # if False, an exception while trying to 
                                    # print the will cause _ret_something() to 
                                    # be called.  Can be set to True for debugging
                                    
    _arraytype = None
    
    def __init__(self,x,u=0,unit=one,dof=float('inf'),k=1,p=None,uunit=None,
                 utype=None,name=None):
        self._old = None
        self.autoconvert = False
        
        if isinstance(x,gummy):
            self._value = nummy(x.value)
            self._value._fp = self._get_p
            self._unit = x._unit
            self._U = self._value.u
            self._value._name = x._value._name
            self._k = 1
            self._pm = None
            self._set_k = True
            return
        
        if isinstance(x,Quantity):
            unit = x.unit
            x = x.value

        if unit is not one:
            unit = Unit.unit(unit)
        self._unit = unit
        
        if isinstance(x,ummy):
            self._value = nummy(x)
            self._value._fp = self._get_p
            self._U = self._value.u
            self._k = 1
            self._pm = None
            self._set_k = True
            return
        
        if p is not None:
            p = float(p)
            self._k = self._p_method.fptok(p,dof,gummy.bayesian)
            self._pm = p
            self._set_k = False
        else:
            if isinstance(k,Integral):
                k = int(k)
            else:
                k = float(k)
            self._k = k
            self._pm = None
            self._set_k = True
                
        if isinstance(x,Distribution):
            self._value = nummy(x,utype=utype,name=name)
            if uunit is None:
                self._U = _ku(self._k,self._value.u)
            else:
                uunit = Unit.unit(uunit)
                if uunit is unit:
                    uunit = None
                self._U = None
                self._set_U(self._k,uunit)
            return
        
        if isinstance(u,gummy):
            uunit = u.unit
            u = u.x
        elif isinstance(u,Quantity):
            uunit = u.unit
            u = u.value
        
        if uunit is not None:
            if u != 0:
                uunit = Unit.unit(uunit)
                if uunit is unit:
                    uunit = None
            else:
                uunit = None

        if uunit is not None and u != 0:
            U = Quantity(u,unit=uunit)
            # Now try to find the uncertainty in the x units
            if not unit.linear:
                u = unit.from_uunit(u,uunit)
            elif unit.is_dimensionless:
                if not uunit.is_dimensionless:
                    raise IncompatibleUnitsError('no conversion found for unit ' + str(uunit) + ' to one')
                if uunit is one:
                    u = U.convert(unit).value
                else:
                    u = abs(x)*U.convert(one).value
            else:
                try:
                    u = U.convert(unit).value
                except IncompatibleUnitsError:
                    # If no conversion was found for uunit to unit, see
                    # if unit can be converted to one.  In this case the u
                    # passed to the intializer was a relative uncertainty.
                    u = abs(x)*U.convert(one).value
                    
        if self._k != 1:
            try:
                u = u/self._k
            except TypeError:
                if isinstance(u,Integral):
                    u = u/MFraction(self._k)
                else:
                    u = u/type(u)(self._k)

        self._value = nummy(x,u=u,dof=dof,utype=utype,name=name)
        self._value._fp = self._get_p
        
        self._U = None
        self._set_U(self._k,uunit)
                        
    @property
    def x(self):
        """
        Gets the gummy's value.  Usually this is the mean of the probability
        distribution.  This property is read-only and returns a float.
        """
        return self._value._x
        
    @property
    def u(self):
        """
        Gets the gummy standard uncertainty in the units given by the `units`
        property (not the `uunits` property).  The standard uncertainty is
        sometimes called the "1-sigma" uncertainty.  The This property is read-only
        and returns a float.
        """
        return self._value.u
    
    @property
    def dof(self):
        """
        `float`, read-only

        Returns the number or degrees of freedom that the uncertainty of the 
        gummy is based on.  If the gummy was created as the result of an 
        operation between two or more other gummys, then the dof is the effective
        number of degrees of freedom calculated using the Welch-Satterthwaite 
        approximation.
        """

        return self.value.dof
    
    @property
    def isindependent(self):
        """
        `bool`, read-only

        Returns `True` if the gummy is an independent variable.  That is the 
        ummy has u > 0 and was not correlated with any other gummy's when it was 
        created or is perfectly correlated or anti-correlated (correlation 
        coefficeint 1 or -1) with such an ummy.'
        """
        
        return self.value.isindependent
        
    @property
    def U(self):
        """
        Gets the expanded uncertainty.  This property is read-only but changing 
        the coverage factor (setting the `k` property), the level of confidence
        (setting the `p` property), or changing the units of the uncertainty
        (setting the `uunit` property) will change the `U` return value.  This
        property returns a float.
        
        The expanded uncertainty `U` (see the `U` property) is related to the
        standard uncertainty u (see the u property) by `U` = `k`*`u`.  The coverage
        factor `k` can be set directly or the desired level of confidence `p`
        (see the `p` property) can be set an k will be calculated based on
        either a normal distribution if the number degrees of freedom (see the 
        dof property) for the uncertainty is infinite or a Student's t
        distribution otherwise.
        
        Setting this property will change the value of the `p` property.

        Examples
        --------
        >>> g = gummy(1,u=0.1)  # Setting k = 2 gives U = 2*u
        >>> g.k = 2
        >>> g.U
        0.2

        `U` may be expressed in different units from `x` by setting the uunit
        property.

        >>> g = gummy(2,0.001,unit='m')
        >>> g.uunit = 'mm'
        >>> g.U
        1
        
        U can also be expressed as a relative uncertainty:
        
        >>> g.uunit = '%'
        >>> g.U
        0.05
        
        The expanded uncertainty is used when printing the gummy:
        
        >>> g
        2.0000 m +/- 0.050%
        """
        if self._Ubr is not None:
            if isinstance(self._Ubr[0],Quantity):
                return [i._value for i in self._Ubr]
            else:
                return self._Ubr
            
        if isinstance(self._U,Quantity):
            return self._U.value
        return self._U
    
    def _set_U(self,k=None,unit=None):
        self._U = self._iset_U(k=k,unit=unit)
        
        if self._ubreakdown is not None:
            self._Ubr = [self._iset_U(k=k,unit=unit,u=float(self.ufrom(t))) for 
                       t in self._ubreakdown]
        else:
            self._Ubr = None
    
    def _iset_U(self,k=None,unit=None,u=None):
        # this is called when k, p, unit, or uunit is changed to get the new
        # value of _U
        
        if u is None:
            u = self.u
            
        if unit is None:
            if isinstance(self._U,Quantity):
                unit = self._U.unit
            
        if u == 0:
            if unit is None:
                return 0
            else:
                return Quantity(0,unit=unit)
        
        if isinf(u) or isnan(u):
            if unit is None:
                return u
            else:
                return Quantity(u,unit=unit)
                
        if k is None:
            k = self.k
        
        if unit is None or unit is self.unit:
            return _ku(k,u)
                
        else:
            if self._unit.linear:
                try:
                    if self.unit.is_dimensionless:
                        raise IncompatibleUnitsError()
                    return Quantity(_ku(k,u),unit=self.unit).convert(unit)
                    
                except IncompatibleUnitsError:
                    try:
                        try:
                            r = abs(_ku(k,u)/self.x)
                        except TypeError:
                            if isinstance(self.x,Integral):
                                r = abs(MFraction.fromnum(_ku(k,u))/self.x)
                            else:
                                r = abs(type(self.x)(_ku(k,u))/self.x)
                        return Quantity(r).convert(unit)
                    except ZeroDivisionError:
                        if not Unit.unit(unit).is_dimensionless:
                            raise IncompatibleUnitsError('no conversion found from unit ' + str(unit) + ' to one')
                        return Quantity(float('inf'),unit=unit)
            else:
                return Quantity(self.unit.to_uunit(_ku(k,u),unit),unit)

                    
    @property
    def Usim(self):
        """
        Gets the expanded uncertainty plus and minus components calculated
        from Monte-Carlo data at the level of confidence given by the `p` property.
        The return value is a tuple equivalent to ``(cisim[1] - xsim,xsim - cisim[0])``,
        see the `cisym` and `xsym` properties.  If no simulated data is available a
        `NoSimulatedDataError` will be raised when this property is called; see the
        sim and simulate methods.  `Usym` is ready-only, but changing the `p` or
        `k` properties will affect `Usym`.
        """
        if not isinstance(self._U,Quantity):
            return self._value.Usim
            
        if self.uunit_is_rel:
            x = self.xsim
            ci = self.cisim
            try:
                a = (ci[1] - x)/abs(x)
                b = (x - ci[0])/abs(x)
            except ZeroDivisionError:
                return (float('inf'),float('inf'))
            unit = self._U.unit
            return (one.convert(a,unit),one.convert(b,unit))
        else:
            unit = self._U.unit
            x = float(self.unit.convert(self.xsim,unit))
            ci1 = float(self.unit.convert(self.cisim[1],unit))
            ci0 = float(self.unit.convert(self.cisim[0],unit))
            return (ci1 - x, x - ci0)
        
        
    @property
    def xsim(self):
        """
        Gets the mean value of the Monte-Carlo data.  If no simulated data is 
        available a `NoSimulatedDataError` will be raised when this property is
        called; see the sim and simulate methods.  This property is read-only.
        """
        return self._value.xsim
        
    @property
    def usim(self):
        """
        Gets the standard deviation of the Monte-Carlo data.  If no simulated 
        data is available a `NoSimulatedDataError` will be raised when this property
        is called; see the sim and simulate methods.  This property is read-only.
        """
        return self._value.usim
        
    @property
    def cisim(self):
        """
        Gets the confidence interval using the Monte-Carlo data with the 
        level of confidence as set with the `p` property.  A tuple is returned
        with the first element giving the lower limit and the second element 
        giving the upper limit of the confidence interval.  If the `cimethod` property
        is set to "shortest", the confidence interval with the shortest length 
        will be found and if the `cimethod` is set to "symmetric" the probabilities
        outside the interval on either side will be equal.  If no simulated 
        data is available a `NoSimulatedDataError` will be raised when this property
        is called; see the `sim` and `simulate` methods.  This property is read-only,
        but changing the `p`, `k`, or `cimethod` properties will affect the return
        value.
        """
        return self._value.cisim
    
    @property
    def cimethod(self):
        """
        `str` in {'shortest', 'symmetric'}), default is 'shortest'
        
        Get or set the method for calculating the confidence interval from 
        Monte-Carlo data.  If this property is set at the class level, it will
        change the default `cimethod` value for new gummys but will no affect
        gummys that have already been created.
        
        Can be set either to the string 'shortest' or the string 'symmetric'.
        This property gets or sets the method for calculating confidence
        intervals from Monte-Carlo data.  
        
        If it is set to 'shortest', the confidence interval will be taken to be 
        the shortest interval that includes the desired fraction of the probability 
        distribution.  
        
        If it is set to 'symmetric', then the confidence interval will be set so 
        that, for n Monte-Carlo samples and a coverage probability of `p`, then
        `n`*(1-`p`)/2 samples lie below the lower limit of the confidence interval
        and the same number of samples lie above the upper limit of the confidence 
        interval.
        """
        
        return self.value._cimethod
    @cimethod.setter
    def cimethod(self,value):
        value = value.lower().strip()
        if value not in ['shortest','symmetric']:
            raise ValueError('cimethod ' + str(value) + ' is not recognized')
        self.value._cimethod = value
        
    @property
    def simdata(self):
        """
        `numpy.ndarray`, read-only

        Returns an array containing the Monte-Carlo simulation data.  A 
        `NoSimulatedDataError` is raised if no Monte-Carlo data is available.
        """
        return self.value.simdata
    
    @property
    def simsorted(self):
        """
        `numpy.ndarray`, read-only

        Returns a sorted array containing the Monte-Carlo simulation data.  A 
        `NoSimulatedDataError` is raised if no Monte-Carlo data is available.
        """
        return self.value.simsorted
    
    @property
    def distribution(self):
        """
        read-only

        Returns ths `Distribution` instance associated with the gummy.
        """
        return self.value.distribution
    
    @property
    def ksim(self):
        """
        read-only

        Returns ``0.5*(gummy.Usim[0] + gummy.Usim[1])/gummy.usim``
        """
        return self.value.ksim
    
    @property
    def independent(self):
        """
        `bool`, read-only

        Returns `False` if the owning gummy was created from a operation involving
        other gummys or has zero uncertainty and `True` otherwise.
        """
        return self.value.isindependent
            
    @property
    def name(self):
        """
        gets or sets an optional name for the gummy, may be `str`, `None` or a
        length four `tuple` of `str`.  If name is set to a length four tuple
        the elements are, in order, the unicode name, the html name, the latex
        name and the ASCII name.  Getting this property returns only the unicode
        name not the full tuple.  Use the `get_name` method to get the name
        in html, latex or ASCII format.
        """
        return self.value.name
    @name.setter
    def name(self,v):
        self.value.name = v
        
    def get_name(self,fmt='unicode',norm=None):
        """

        Parameters
        ----------
        fmt: The format, must be a `str` in {'unicode','html','latex','ascii'}.
              The default is 'unicode'.
        norm: An optional function which returns the name in nomral text.  This
              function is applied to the name before it is returned if
              fmt = 'latex', name has been set to single string, and the name
              is more than one character long.  The default is 
              '\\text{' + name + '}'.

        Returns
        -------
        `str`, the name in the requested format.

        """
        return self.value.get_name(fmt,norm)
    
    @property
    def unit(self):
        """
        Gets or sets the unit for `x `and, if the `uunit` attribute is
        `None`, the units for the uncertainty.
        
        If this property is set, a unit conversion will be performed.  The value 
        it is set to may be a string, `None`, a `Unit` object, or the integer 1.
        Both 1 and `None` will be interpreted as the Unit instance `one`. A
        `NoUnitConversionFoundError` will be raised if the unit conversion is
        not possible.
        
        Example
        -------
            
        >>> x = gummy(0.001,unit='V')
        >>> x
        0.001 V
        >>> x.unit = 'uV'
        >>> x
        1000.0 uV
        """

        return self._unit
    @unit.setter
    def unit(self,u):
        Quantity.unit.fset(self,u)
        self._U = None
        self._set_U()
            
    @property
    def uunit(self):
        """
        Gets or sets the units for the expanded uncertainty (see the `U`
        property).  This property may be set to a unit with the same dimension
        as the `unit` property or to a dimensionless unit such at "%" or "ppm"
        in which case `U` will be a relative uncertainty.  Setting `uunit` to `None`
        puts `U` in the units given by the `units` property.
        
        Examples
        --------

        `U` may be expressed in different units from x by setting the `uunit`
        property.
        
        >>> g = gummy(2,0.001,unit='m')
        >>> g.uunit = 'mm'
        >>> g.U
        1
        
        U can also be expressed as a relative uncertainty:
        
        >>> g.uunit = '%'
        >>> g.U
        0.05
        
        Setting uunit to None changes the units of U back to the original units:
        >>> g.uunit = None
        >>> g.U
        0.001
        
        More examples:
        >>> g = gummy(0.001,0.0000012,unit='V')
        >>> g
        (0.001 000 0 +/- 0.000 001 2) V
        >>> g.unit = 'uV'
        >>> g
        (1000.0 +/- 1.2) uV
        >>> g.uunit = '%'
        >>> g.unit = 'mV'
        >>> g
        1.0000 mV +/- 0.12%
        >>> g.unit = 'uV'
        >>> g
        1000.0 uV +/- 0.12%
        >>> g.uunit = None
        >>> g
        (1000.0 +/- 1.2) uV
        """
        if isinstance(self._U,Quantity):
            return self._U.unit
        else:
            return None
    @uunit.setter
    def uunit(self,unit):
        if self.u == 0:
            return
        if unit is None:
            unit = self._unit
        else:
            unit = Unit.unit(unit)
        self._U = None
        self._set_U(None,unit)
        
    @property
    def uunit_is_rel(self):
        """
        Returns True if gummy.U is a relative uncertainty and False otherwise.
        This property is read-only.
        """
        if not isinstance(self._U,Quantity):
            return False
        try:
            if self._U.unit is one:
                return False
            self._U.unit.convert(self._U,one)
            return True
        except:
            return False
        
    @property
    def k(self):
        """
        Get or set the coverage factor; must be > 0.
        
        The expanded uncertainty `U` (see the `U` property) is related to
        the standard uncertainty `u` (see the `u` property) by `U` = `k`*`u`.
        The coverage factor `k` can be set directly or the desired level
        of confidence `p` (see the `p` property) can be set an `k` will be
        calculated based on either a normal distribution if the number
        degrees of freedom (see the `dof` property) for the uncertainty is
        infinite or a Student's t distribution otherwise.
        
        Setting this property will change the value of the `p` property.
        
        Examples
        --------
        
        Setting `k` = 2 gives `U` = 2*`u`:
        
        >>> g = gummy(1,u=0.1)
        >>> g.k = 2
        >>> g.U
        0.2
        
        The `p` property has been set based on the value of `k`:
        
        >>> g.p
        0.95449973610364158
        
        Changing the `p` property will change the value of the `k` property:
        
        >>> g.p = 0.9973
        >>> g.k
        2.9999769927034015
        """
        if self.u == 0:
            return None
        return self._k
    @k.setter
    def k(self,v):
        if v <= 0:
            raise ValueError('k <= 0')
        self._k = v
        self._pm = None
        self._set_k = True
        self._set_U(self._k,None)
        
    def _get_p(self):
        #if self.u == 0:
            #return 1
        if self._pm is None:
            self._pm = self._p_method.fktop(self._k,self.dof,self.bayesian)
        if self._pm < 0:
                return 0
        return self._pm
        
    @property
    def p(self):
        """
        Get or set the level of confidence; must be in the interval (0,1).
        
        The expanded uncertainty `U` (see the `U` property) is related to
        the standard uncertainty `u` (see the `u` property) by `U` = `k`*`u`.
        The coverage factor `k` can be set directly or the desired level
        of confidence `p` (see the `p` property) can be set an `k` will be
        calculated based on either a normal distribution if the number
        degrees of freedom (see the `dof` property) for the uncertainty is
        infinite or a Student's t distribution otherwise.
        
        Setting this property will change the value of the "k" property.
        
        Examples
        --------

        Setting `k` = 2 gives `U` = 2*`u`:
        
        >>> g = gummy(1,u=0.1)
        >>> g.k = 2
        >>> g.U
        0.2
        
        The `p` property has been set based on the value of `k`:
        
        >>> g.p
        0.95449973610364158
        
        Changing the `p` property will change the value of the `k` property:
        
        >>> g.p = 0.9973
        >>> g.k
        2.9999769927034015
        """
        return self._get_p()
    @p.setter
    def p(self,v):
        if isinstance(v,str) and v.strip().lower() == 'ssd':
            self.p = loc_from_k(1)
            return
        if v <= 0 or v >= 1:
            raise ValueError('p is not in the interval (0,1)')
        self._pm = v
        self._k = self._p_method.fptok(v,self.dof,self.bayesian)
        self._set_k = False
        self._set_U(self._k,None)
        
    def correlation(self,gummy):
        """
        Returns the correlation coefficient between `self` and `g`.
        """
        if isinstance(gummy,Quantity):
            gummy = gummy.value
        return self.value.correlation(gummy)
    
    def covariance(self,gummy):
        """
        Returns the covariance between `self` and `g`.
        """
        if isinstance(gummy,Quantity):
            gummy = gummy.value
        return self.value.covariance(gummy)
    
    @staticmethod
    def correlation_matrix(gummys):
        """
        Returns the correlation matrix of a list or array of gummys.
        """
        m = [g.value if isinstance(g,Quantity) else g for g in gummys]
        return [[b.correlation(a) if isinstance(b,ummy) else 0 
                for b in m] for a in m]
        
    @staticmethod
    def covariance_matrix(gummys):
        """
        Returns the variance-covariance matrix of a list or array of gummys.
        """
        m = [g.value if isinstance(g,Quantity) else g for g in gummys]
        return [[b.covariance(a) if isinstance(b,ummy) else 0 
                for b in m] for a in m]
    
    def correlation_sim(self,gummy):
        """
        Returns the correlation coefficient, calculated from Monte-Carlo data, 
        between the owning gummy and the gummy `g`.
        
        See the method `gummy.correlation(g)` for the corresponding result based
        on first order error propagation.
        """
        if isinstance(gummy,Quantity):
            gummy = gummy.value
        return self.value.correlation_sim(gummy)
    
    def covariance_sim(self,gummy):
        """
        Returns the covariance, calculated from Monte-Carlo data, between the 
        owning gummy and the gummy `g`.
        
        See the method `gummy.covariance(g)` for the corresponding result based
        on first order error propagation.
        """
        if isinstance(gummy,Quantity):
            gummy = gummy.value
        return self.value.covariance_sim(gummy)
    
    @staticmethod
    def correlation_matrix_sim(gummys):
        """
        The staticmethod takes a list of gummys an returns the correlation
        matrix calculated from Monte-Carlo data.  The return value is numpy 
        ndarray.
        
        See the method `gummy.correlation_matrix(gummys)` for the corresponding
        result based on first order error propagation.
        """
        m = [g.value if isinstance(g,Quantity) else g for g in gummys]
        return [[b.correlation_sim(a) if isinstance(b,nummy) else 0 
                for b in m] for a in m]
    
    @staticmethod
    def covariance_matrix_sim(gummys):
        """
        The staticmethod takes a list of gummys an returns the variance-covariance
        matrix calculated from Monte-Carlo data.  The return value is numpy
        ndarray.
        
        See the method gummy.covariance_matrix(gummys) for the corresponding
        result based on first order error propagation.
        """
        m = [g.value if isinstance(g,Quantity) else g for g in gummys]
        return [[b.covariance_sim(a) if isinstance(b,nummy) else 0 
                for b in m] for a in m]
    
    @property
    def finfo(self):
        return self.value.finfo
    
    @property
    def utype(self):
        """
        `str`, `None` or a list containing strings and possibly `None`
    
        An arbitrary string value labeling the uncertainty type or or a 
        list of types if the gummy was constructed from independent
        variables with different utypes.
        """
        return self.value.utype
        
    def ufrom(self,x):
        """
        Gets the standard uncertainty contributed from particular gummys
        or utypes if all other free variables are held fixed.
        
        Parameters
        ----------
        x:  `gummy`, `str`, or array_like
            A gummy, a string referencing a utype or a list containing
            gummys and strings.
            
        Returns
        -------
        `float`
        
        Example
        -------
        >>>  a = gummy(1.2,0.2,utype='A')
        >>>  b = gummy(3.2,0.5,utype='A')
        >>>  c = gummy(0.9,0.2,utype='B')
        >>>  d = a + b + c
        >>>  d.ufrom('A')
        0.53851648071345048
        """
        return self.value.ufrom(x)
    
    def doffrom(self,x):
        """
        Gets the degrees of freedom contributed from particular gummys or
        utypes if all other free variables are held fixed.  Caution:  any
        correlations in the calculations can cause errors in dof calculations.
        
        Parameters
        ----------
        x:  `gummy`, `str`, or array_like
            A gummy, a string referencing a utype or a list containing
            gummys and strings.

        Returns
        -------
        `float`

        Example
        -------
        >>>  a = gummy(1.2,0.2,dof=5,utype='A')
        >>>  b = gummy(3.2,0.5,dof=7,utype='A')
        >>>  c = gummy(0.9,0.2,utype='B')
        >>>  d = a + b + c
        >>>  d.doffrom('A')
        9.0932962619709627
        """
        return self.value.doffrom(x)
    
    def ufromsim(self,x):
        """
        Gets the standard deviation of the Monte-Carlo data only allowing the
        independent variables in `x` to vary.  Independent istributions not in 
        `x` are held fixed.  `sim` or `simulate` must be called to generate
        Monte-Carlo data before calling this method.
        
        Parameters
        ----------
        x:  `gummy`, `str`, or array_like
            A gummy, a string referencing a utype or a list containing
            gummys and strings.
            
        Returns
        -------
        `float`
        """
        if _isscalar(x):
            x = [x]
        x = [i.value if isinstance(i,Quantity) else i for i in x]
        
        return self.value.ufromsim(x)
    
    def datafrom(self,x,save=True):
        """
        Recomputes the Monte-Carlo `simdata `with only the varaibles in `x` 
        allowed to vary.  `sim` or `simulate` must be called to generate
        Monte-Carlo data before calling this method.  This method cannot be 
        called with save == `True` from from a gummy representing an independent 
        variable (that is from a gummy not created by by mathematical operations
        between two or more other gummy's).
        
        Parameters
        ----------
        ufrom:  list containing `gummy` or `str`
            all independent gummys not in the list or having a utype
            not in the list are held fixed at their `.x` value
            
        save:  If `save` is `True` the recomputed data is stored in the `simdata`
            attribute and `None` is returned.  If `save` is `False` then the
            recomputed data is returned and the `simdata` attribute is not
            overwritten.
        
        Returns
        -------
        'numpy.array' if `save` is `False`, otherwise returns `None`

        Raises
        ------
        `NoSimulatedDataError`:
            if no simulated data is available from a call to
            `Distribution.simulate`.
        `RuntimeError`:  
            if this method is called from an independent `gummy` 
        """
        if _isscalar(x):
            x = [x]
        x = [i.value if isinstance(i,Quantity) else x for i in x]
        
        return self.value.datafrom(x,save=save)
        
    @property
    def style(self):
        """
        Get or set the default display style for new gummys.  This is a string
        with one of the following values:
        
        - "pm" or "+-" gives, e.g. in ascii format "(1.00 +/- 0.12)e-12 cm"
        - "pmi" or "+-i" gives e.g. "1.00e-12 cm +/- 1.2e-13 cm"
        - "concise" or "()" gives e.g. "1.00(12)e-12 cm"
        - "ueq" or "u equals" gives e.g. "1.00e-12 cm with u = 1.2e-13 cm"
        - "x" or "x only" gives e.g. "1.00e-12 cm"
        - "xf" gives e.g. "1.00e-12"
        - "u" or "u only" gives "1.2e-13 cm"
        - "uf" gives "1.2e-13"
        
        The following styles display a mean and confidence interval based on
        data from a Monte-Carlo simulation.  If no simulated data is available
        "no simulated data" is returned.
        
        - "pmsim" gives e.g. "(1.01 + 0.11 - 0.13)e-12 cm"
        - "pmsimi" gives e.g. "1.01e-12 cm + 1.1e-13 cm - 1.3e-11 cm"
        - "mcisim" gives e.g. "mean = 1.01e-12 cm, confidence interval = 
          [8.8e-13 cm, 1.13e-12 cm]"
        - "cisim" gives e.g "[8.8e-13 cm, 1.13e-12 cm]"
        - "usim" gives the standard deviation e.g. "1.2e-13 cm"
        - "ufsim" gives the standard deviation e.g. "1.2e-13"
        
        Note if `uunit` has been defined for the gummy instance, then concise style
        and pm style are not valid and the display will default to pmi style.
        
        The style can be set either at the class or the instance level.
        
        To control whether the coverage factor (k), the level of confidence (p),
        the degrees of freedom (dof), and name are displayed as part of the 
        uncertainty statement set the following variables to True or False (at 
        either the class or instance level):
        
        - `show_k`
        - `show_p`
        - `show_dof`
        - `show_name`
        
        `show_k`, `show_p`, and `show_dof` can also be set to `None`, which will
        allow the gummy object to decide whether to display the corresponding
        values.
        
        The number of significant digits in the uncertainty to be displayed is
        set with the `nsig` attribute.  This can be set either at the class level
        or the instance level.  When the gummy module is loaded `nsig` is set to 2.
        
        Scientific notation can be turned on or off by setting the `sci_notation`
        attribute to `True` or `False`.  If `sci_notation` is set to `None`, scientific notation
        will be used if the exponent of gummy `x` value is greater than than the
        value of the `sci_notation_high` attribute or less than the `sci_notation_low`
        attribute.  When the gummy module is loaded, `sci_notation` is `None`,
        `sci_notation_high` is 7 and `sci_notation_low` is -3.  These attributes
        can be set either at the class or instance level.
        
        Setting the `solidus` attribute to `True` uses a solidus (forward slash)
        to separate units with positive and negative exponents e.g. "m/s" and when
        `solidus` is set to `False` negative exponents are used e.g. "m s**-1".
        Setting the `mulsep` attribute to `True` inserts a dot or * between units while setting
        this to `False` a space separates units. When the gummy module is loaded,
        both `solidus` and `mulsep` are set to `False`.  They may be set at either
        the instance or class level.
        """
        return self._style
    @style.setter
    def style(self,value):
        self._style = gummy._get_style(value)
        
    @staticmethod
    def _get_style(text):
        text = text.strip().lower()
        text = text.replace(' ','')
        if text == 'pm' or text == 'plusminus' or text == '+-' or text == '+/-':
            return 'pm'
        if text == 'pmi' or text == 'plusminusi' or text == '+-i' or text == '+/-i':
            return 'pmi'
        if text == 'ueq' or text == 'uequal' or text == 'uequals':
            return 'ueq'
        if text == 'concise' or text == '()':
            return 'concise'
        if text == 'concisef' or text == '()f':
            return 'concisef'
        if text == 'x' or text == 'xonly':
            return 'x'
        if text == 'u' or text == 'uonly':
            return 'u'
        if text == 'pmsim' or text == 'plusminussim' or text == '+-sim' or text == '+/-sim' or text == 'simpm' or text == 'simplusminus' or text == 'sim+-' or text == 'sim+/-':
            return 'pmsim'
        if text == 'pmsimi' or text == 'plusminussimi' or text == '+-simi' or text == '+/-simi' or text == 'simpmi' or text == 'simplusminusi' or text == 'sim+-i' or text == 'sim+/-i':
            return 'pmsimi'
        if text == 'cisim' or text == 'confidenceintervalsim' or text == 'simci' or text == 'simconfidenceinterval':
            return 'cisim'
        if text == 'mcisim' or text == 'meanconfidenceintervalsim' or text == 'simmci' or text == 'simmeanconfidenceinterval':
            return 'mcisim'
        if text in ['usim','ufsim','xsim','xfsim','xf','uf','xunit','uunit']:
            return text
        raise ValueError('style ' + str(text) + ' is not recognized')
        
    @property
    def solidus(self):
        if self._solidus is None:
            return Unit.solidus
        return self._solidus
    @solidus.setter
    def solidus(self,v):
        self._solidus = bool(v)
        
    @property
    def mulsep(self):
        if self._mulsep is None:
            return Unit.mulsep
    @mulsep.setter
    def mulsep(self,v):
        self._mulsep = bool(v)
        return self._mulsep
        
    @property
    def bayesian(self):
        """
        `bool`

        Read/write at the class level, but read-only at the instance level.
        The default value is `False`; this should only be changed once at the
        beginning of the session.  This property affects how the level of
        confidence `p` (sometimes called coverage probability) of an expanded
        uncertainty is related to the coverage factor `k` for a gummy based on
        data with finite degrees of freedom.

        Standard uncertainties are often based on the standard deviation of a set
        of measurements (and the assumption that these measurements are drawn
        from a normally distributed population).  Traditionally (e.g. the GUM
        2008 edition) the standard uncertainty is taken to be the standard
        deviation of the mean (s/sqrt(n), where s is the sample standard deviation
        and n is the number of measurements).  However there is some "extra
        uncertainty" because the sample standard devation not exactly equal to
        the population standard deviation.  This is taken into account by using
        a Student's t distribution to calculate the expanded uncertainty.  However
        it has been pointed out, by those who advocate a Bayesian point of view,
        that the probability distribution for the measurand here is best described
        by a shifted and scaled Student's t distribution.  So the standard
        uncertainty should be the standard deviation of this distribution which
        is s*sqrt{(n-1)/[n*(n-3)]}.  Thus

        u(bayesian) = sqrt[dof/(dof - 2)]*u(traditional)

        where dof = n - 1 and the "extra uncertainty" is incorporated directly
        into the standard uncertainty.

        Example
        -------
        >>> gummy.bayesian = True
        >>> g = gummy(1,0.03,dof=5)
        >>> g.bayesian
        True
        """
        return self.value._bayesian
        
    @property
    def nsig(self):
        """
        `int`
        
        Gets or sets the number of significant digits in the uncertainty to 
        display.
        """
        return self.value.nsig
    @nsig.setter
    def nsig(self,v):
        self.value.nsig = v
        
    @property
    def thousand_spaces(self):
        """
        `bool`
        
        Gets or sets a bool value that determines whether to insert a space to
        group digits in x.
        """
        return self.value.thousand_spaces
    @thousand_spaces.setter
    def thousand_spaces(self,v):
        self.value.thousand_spaces = bool(v)
        
    @property
    def sci_notation(self):
        """
        `bool` or `None`
        
        Setting this to True or False forces or prevents scientific notation
        in the display of the value.  If this is set to `None` thr scientific
        notation will be used if the x value is below `10**sci_notation_low` or
        above `10**sci_notation_high`.
        """
        return self.value.sci_notation
    @sci_notation.setter
    def sci_notation(self,v):
        self.value.sci_notation = v
        
    @property
    def sci_notation_high(self):
        """
        see the sci_notation property
        """
        return self.value.sci_notation_high
    @sci_notation_high.setter
    def sci_notation_high(self,v):
        self.value.sci_notation_high = v
        
    @property
    def sci_notation_low(self):
        """
        see the sci_notation property
        """
        return self.value.sci_notation_low
    @sci_notation_low.setter
    def sci_notation_low(self,v):
        self.value.sci_notation_low = v
        
    @property
    def max_digits(self):
        """
        `int`
        
        Gets or sets the maximum number of digits in x to display.
        """
        return self.value.max_digits
    @max_digits.setter
    def max_digits(self,v):
        v = int(v)
        if v < 0:
            raise ValueError('max_digits must be >= 0')
        self.value.max_digits = v
        
            
    def copy(self,formatting=True,totype=None):
        """
        Returns a copy of the gummy.  If the `formatting` parameter is
        `True` the display formatting information will be copied and if
        `False` the display formatting will be set to the default for a
        new gummy.  The default for `formatting` is `True`.  If 'totype`
        is defined the x and u properties will be converted to type `totype`
        before copying.
        """
        
        r = type(self)(self._value.copy(formatting=formatting,totype=totype),
                       unit = self._unit)
        r._old = self._old
                
        if formatting:
            if self._style != type(r)._style:
                r._style = self._style
            if self.show_p != type(r).show_p:
                r.show_p = self.show_p
            if self.show_k != type(r).show_k:
                r.show_k = self.show_k
            if self.show_dof != type(r).show_dof:
                r.show_dof = self.show_dof
            if self.mulsep != type(r).mulsep:
                r.mulsep = self.mulsep
            if self.solidus != type(r).solidus:
                r.solidus = self.solidus
            r._k = self._k
            r._pm = self._pm
            r._set_k = self._set_k
            if totype is not None:
                r._set_U(unit=self.uunit)
            else:
                r._U = self._U
        else:
            r._U = self.u
            r._k = 1
            r._pm = None
            r._set_k = True
        
        return r
        
    def graft(self,unit):
        """
        Returns a copy of the gummy with different units but the same `x` and
        `u` values.  This is different from ``gummy.convert(unit)`` in that
        ``gummy.convert(unit)`` changes the `x` and `u `values to express the
        same quantity  in different units while `gummy.graft(unit)` simply
        tacks on a different unit to the same numerical values.
        
        Parameters
        ----------
        unit:  `str` or `Unit`
            The unit for the `x` value and if `uunit` is `None`, the
            uncertainty.  It must be string, None, a `Unit` object, or the
            integer 1.  Both 1 and `None` will be interpreted as the Unit
            instance `one`.
        """     
        return self*Unit.unit(unit)/self.unit
        
    @staticmethod
    def simulate(gummys,n=100000,ufrom=None):
        """
        Generates Monte-Carlo data for one or more gummys.  Calling this method
        erases previously generated Monte-Carlo data for all gummys.  See also
        the `gummy.sim()` method to generate data for one gummy only.
        
        Parameters
        ----------
        n:  `int` > 0, optional
            The number of samples to generate.  The default value is 100000.
            
        gummys: A list or array of `gummy` for which to generate the Monte-Carlo
            data.

        ufrom: `None`, `gummy`, `str` or array_like
            If this is not `None`, then only the gummys referenced here will be
            allowed to vary, and all other gummys will be held fixed at their
            mean values.  This can be a gummy, a string referencing a utype or
            a list containing  gummys and strings.  The default value is `None`.
        """
        if ufrom is not None:
            if isinstance(ufrom,(Quantity,ummy,str,Distribution)):
                ufrom = [ufrom]
            ufrom = [g.value if isinstance(g,Quantity) else g for g in ufrom]
        gummys = [g.value if isinstance(g,Quantity) else g for g in gummys]
        return nummy.simulate(gummys,n=n,ufrom=ufrom)
    
    def sim(self,n=100000,ufrom=None):
        """
        Generates Monte-Carlo data for this gummy (as well as for any gummys 
        that this gummy depends on).  Calling this method erases previously 
        generated Monte-Carlo data for all gummys, so use the 
        `gummy.simulate()` staticmethod if you need Monte-Carlo data for
        several gummys simultaneously.
        
        Parameters
        ----------
        n:  `int` > 0, optional
            The number of samples to generate.  The default value is 100000.
            
        ufrom: `None`, `gummy`, `str` or array_like
            If this is not `None`, then only the gummys referenced here will be
            allowed to vary, and all other gummys will be held fixed at their
            mean values.  This can be a gummy, a string referencing a utype or
            a list containing  gummys and strings.  The default value is `None`.
        """
        return gummy.simulate([self],n,ufrom)
    
    @staticmethod
    def clear_all():
        """
        Clears Monte-Carlo data from all existing gummys.
        """
        Distribution.clear_all()
        
    def clear(self):
        """
        Clears the gummys Monte-Carlo data.
        """
        self.value.clear()
        
    @classmethod
    def _plotlabel(cls,g,symbol=None,exponent=0,math=None,norm=None,
                   slashaxis=None):
        # This returns a string that is used for an axis label for gummy.hist() and
        # gummy.covplot().
        
        if math is None:
            math = cls.latex_math_plot
        if norm is None:
            norm = cls.latex_norm_plot
        if slashaxis is None:
            if isinstance(g,nummy):
                slashaxis = g.slashaxis
            else:
                slashaxis = gummy.slashaxis
            
        if g is None:
            name = ''
        else:
            try:
                name = g.get_name(fmt='latex',norm=norm)
                if name is None:
                    name = norm('value')
            except:
                name = str(g)
                if len(name) > 1:
                    name = norm(name.strip())
            
        if symbol is None and isinstance(g,gummy) and g.unit is not one:
            unit = g.unit
            if slashaxis:
                s = False
            else:
                s = g.solidus
            symbol = unit.tostring(fmt='latex',solidus=s)
        if symbol == '':
            symbol = None
            
        xl = name
        if symbol is not None:
            if xl != '':
                if slashaxis:
                    if exponent != 0:
                        xl += '\\,/\\,10^{' + str(exponent) + '} ' + symbol
                    else:
                        xl += '\\,/\\,' + symbol
                else:
                    if exponent != 0:
                        xl += '(10^{' + str(exponent) + '} ' + symbol + ')'
                    else:
                        xl += '(' + symbol + ')'
            else:
                xl += symbol
        if xl != '':
            return math(xl)
        else:
            return None
        
    def hist(self,title=None,xlabel=None,p=None,show_p=True,title_style=None,
             mean_marker=True,mean_marker_options={},ci_marker=True,
             ci_marker_options={},hold=False,math=None,norm=None,**plot_options):
        """
        Plots a histogram of the Monte-Carlo data for the gummy.  Before calling
        this method gummy.sim or gummy.simulate must be called to generate the
        Monte-Carlo data.
        
        Parameters
        ----------
        all parameters are optional
        
        title:  `str` or `None`, optional
            A title for the plot.  If this is omitted or set to `None` then
            a title will be generated using the gummy name (if it has one)
            and the mean value and confidence interval.  The title will also
            give the standard deviation of the date.  The formatting of the
            auto-generated title depends on the value of the `title_style`
            parameter.
            
        xlabel:  `str` or `None`, optional
            A label for the horizontal axis of the plot.  If this is omitted
            or set to None then a label will be generated using the name and
            unit of the gummy.  If xlabel is None and the gummy has no name
            and a unit of one, then the horizontal axis will not be labeled.
            
        p:  `float` between 0 and 1 or `None`
            The probability for the confidence interval (as printed in the
            title and indicated by the `ci_markers`). If this is `None` then
            the value of the `gummy.p` class property is used. The default
            is `None`.
            
        show_p:  `bool`, optional
            Whether or not to show the level of confidence in the title if the
            title is auto-generated.  The default value is `True`
            
        title_style:  {'pmsim','pmsimi','cisim','mcisym','xsim','xfsim','usim','ufsim'}, optional
            The style for displaying the value in the title. See the
            `gummy.style` property for details.  It this is None or omitted
            then the value of the `gummy.style` class property is used.
            
        mean_marker:  `bool`, optional
            Whether or not to display a vertical line at the mean value The
            default is `True`.
           
        mean_marker_options:  `dict`, optional
            A dictionary containing keywords to be passed to the `pyplot.axvline`
            method which draws the mean marker.  For example setting this to
            ``{'color'='r','linewidth'=4}`` makes the mean marker red and with
            a thickness of four points.
           
        ci_marker:  `bool`, optional
            Whether or not to display vertical lines at the upper and lower
            limits of the confidence interval.  The default is `True`.
           
        ci_marker_options:  `dict`, optional
            A dictionary containing keywords to be passed to the `pyplot.axvline`
            method which draws the confidence interval markers.
            
        hold:  `bool`, optional
            If this is `False` ``pyplot.show()`` is called before this method exits.
            If it is `True` ``pyplot.show()`` is not called.  The default is
            `False`.
           
        plot_options:
             These are optional keyword arguments that are passed to the `pyplot.hist`
             method.  For example bins=50 overrides the default number of bins (100).
             For other options see the `pyplot.hist` documentation.
        """
        import matplotlib.pyplot as plt
                    
        g = self.copy(formatting=True)
        if p is not None and p != self.p:
            g.p = p
            
        if math is None:
            math = type(self).latex_math_plot
        if norm is None:
            norm = type(self).latex_norm_plot
        if xlabel is None:
            xlabel = gummy._plotlabel(self,math=math)
        if title is None:
            if title_style is None:
                title_style = self.style
            if title_style not in ['pmsim','pmsimi','cisim','mcisym','xsim','xfsim','usim','ufsim']:
                title_style = 'pmsim'
            else:
                title_style = self.style
            g.style = title_style
            g.show_p = show_p
            g.show_k = False
            g.show_dof = False
            if g.xsim >= g.cisim[0] and g.xsim <= g.cisim[1] and not isinf(g.xsim):
                title0 = g.tostring(fmt='latex',norm=norm)
                title0 = math(title0)
            else:
                title0a = norm('mean') + '=' 
                title0a += g.tostring(fmt='latex',style='xsim',norm=norm)
                title0b = norm('confidence interval') + '='
                title0b += g.tostring(fmt='latex',style='cisim',norm=norm)
                title0 = math(title0a)+ '\n' + math(title0b)
            title1 = norm('standard deviation') + ' = '
            title1  += self.tostring(fmt='latex',style='usim',norm=norm)
            title1 = math(title1)
            title = title0 + '\n' + title1
            
        if xlabel is not None:
            plot_options['xlabel'] = xlabel
        if title is not None:
            plot_options['title'] = title
            
        ci = None
        if 'range' not in plot_options and ci_marker:
            ci = g.cisim
            a = (ci[1] - ci[0])/3
            plot_options['range'] = (ci[0] - a,ci[1] + a)
        
        self.value.hist(hold=True,**plot_options)
        
        if mean_marker:
            if 'linewidth' not in mean_marker_options and 'lw' not in mean_marker_options:
                mean_marker_options['linewidth'] = 1
            if 'color' not in mean_marker_options:
                mean_marker_options['color'] = 'k'
            if 'zorder' not in mean_marker_options:
                mean_marker_options['zorder'] = 3
                
            xsim = g.xsim
            plt.axvline(xsim,**mean_marker_options)
        
        if ci_marker:
            if 'linewidth' not in mean_marker_options and 'lw' not in ci_marker_options:
                ci_marker_options['linewidth'] = 1
            if 'color' not in ci_marker_options:
                ci_marker_options['color'] = '0.5'
            if 'zorder' not in ci_marker_options:
                ci_marker_options['zorder'] = 3
                
            if ci is None:
                ci = g.cisim
            plt.axvline(ci[0],**ci_marker_options)
            plt.axvline(ci[1],**ci_marker_options)
        
        if not hold:
            plt.show()
        
    @staticmethod
    def covplot(x,y,title=None,xlabel=None,ylabel=None,mean_marker=False,
                mean_marker_options={},hold=False,math=None,**plot_options):
        """
        Creates scatter plot showing the covariance between two gummys.
        
        Parameters
        ----------
        x:  `gummy`
            The gummy to plot on the horizontal axis.
        
        y:  `gummy`
            The gummy to plot on the vertical axis.
        
        title:  `str` or `None`, optional
            A title for the plot.  If this is omitted or set
            to None then the correlation will be displayed as the title.
           
        xlabel:  `str` or `None`, optional
            A label for the horizontal axis.  If this os
            omitted or None then that axis will be labeled either "x" or with
            the `x` gummy's unit.
           
        ylabel:  `str` or `None`, optional
            A label for the vertical axis.  If this os
            omitted or None then that axis will be labeled either "y" or with
            the `y` gummy's unit.
           
        mean_marker:  `bool`, optional
            Whether or not to display line markers at the mean
            values of `x` and `y`.  The default is `False`.
           
        mean_marker_options:  `dict`, optional
            A dictionary of options to be passed to the
            `pyplot.axvline` and `pyplot.axhline` methods that draw the `mean_marker`.
           
        hold:  `bool`, optional
            If this is `False` then ``pyplot.show()`` is called before this method
            exits.  If it is `True` ``pyplot.show()`` is not called.  The default is
            `False`.
           
        plot_options:  These are optional keyword arguments that are passed to
             the `pyplot.plot` method.  For example ``ms=0.1`` decreases the size of the
             dots in the plot.
        """
        import matplotlib.pyplot as plt
        
        if math is None:
            math = nummy.latex_math_plot
        if xlabel is None:
            xlabel =  gummy._plotlabel(x,math=math)
        if ylabel is None:
            ylabel =  gummy._plotlabel(y,math=math)
            
        nummy.covplot(x,y,title=title,xlabel=xlabel,ylabel=ylabel,
                                   hold=True,**plot_options)
        
        if mean_marker:
            if 'linewidth' not in mean_marker_options and 'lw' not in mean_marker_options:
                mean_marker_options['linewidth'] = 1
            if 'color' not in mean_marker_options:
                mean_marker_options['color'] = '0.5'
            if 'zorder' not in mean_marker_options:
                mean_marker_options['zorder'] = 3
            plt.axvline(x.xsim,**mean_marker_options)
            plt.axhline(y.xsim,**mean_marker_options)
            
        if not hold:
            plt.show()
            
    def toummy(self):
        """
        returns a `Quantity` with an `ummy` value.
        """
        return Quantity(self.value.toummy(),unit=self.unit)
    
    def splonk(self):
        """
        splonks the gummy
        """
        return Quantity(self.value.toummy().splonk(),unit=self.unit).splonk()
        
    @property
    def ubreakdown(self):
        """
        `list` of `str` or `None`, default is `None`
        
        If this is set to a list containing strings referencing utypes then 
        when the gummy is printed, the uncertainty from each utype will be 
        displayed separately.
        
        Examples
        --------
        >>>  a = gummy(1.2,0.2,utype='A')
        >>>  b = gummy(3.2,0.5,utype='A')
        >>>  c = gummy(0.9,0.2,utype='B')
        >>>  d = a + b + c
        >>>  d
        5.30(57)
        >>>  d.ubreakdown = ['A','B']
        >>>  d
        5.30(54)(20)
        >>>  d.style = 'ueq'
        >>>  d
        5.30 with u(A) = 0.54 and u(B) = 0.20
        """
        return self._ubreakdown
    @ubreakdown.setter
    def ubreakdown(self,v):
        if v is None:
            self._ubreakdown = None
        else:
            try:
                self._ubreakdown = list(v)
            except TypeError:
                self._ubreakdown = [v]
                
        try:
            self._set_U()
        except:
            self._ubreakdown = None
            self._Ubr = None
            raise
    
    def budget(self,xlist,**kwds):
        """
        Returns a Budget object that can be used to display an uncertainty
        budget table listing the the contributions of the gummys in `xlist`
        to the total uncertainty in the calling gummy.

        To display the table use the `Budget.html()` or `Budget.latex()`
        methods  in a console or notebook that supports this type of output
        or the python  built-in function to get a unicode table.

        The `Budget.tohtml()` and `Budget.tolatex()` methods can be used
        to get strings with the html or latex code.

        The `Budget.df` property can be used to retrieve a pandas `DataFrame`
        with the table.  Also `Budget.df_str`, `Budget.df_html` and
        `Budget.df_latex` return DataFrames with formatted strings as entries
        rather than numerical values.

        Parameters
        ----------
        xlist:   array_like of `gummy`
            The independent variables.  Warnings will be generated if the
            gummys in this list over determine `self` (that is if not all
            variables in this list can be treated as independent variables)
            or under determine `self` (that is if some variables
            contributing to the uncertainty in `self` are missing).
            
        uunit:  `str` or Unit, optional
            Unit to use to express the uncertainties.  This useful if you
            wish to express all uncertainties  as relative uncertainty unit
            (e.g. %).

        k, p: `float`, optional
            k or p values for the expanded uncertainty; do not specify both `k`
            and `p`; if neither are specified the the `k` and `p` values of
            `self` are used

        sort:  `bool`, optional
            Whether or not to sort the gummys in `xlist` by significance.

        columns: `list` of `str` or `None`, optional
            Allows the user to select the columns (and ordering of the columns)
            for display.  The available columns are:

            "component" or "name":  the names of the gummy, displayed by default

            "description":  description given in the description parameter list,
                displayed by default if the description parameter is not None

            "unit":  the unit of the gummy, displayed by default

            "value":  the x value of the gummy, displayed by default

            "u" or "uncertainty":  The uncertainty of the gummy.  This is the
                standard uncertainty except possible in the last row where an
                expanded uncertainty is displayed.  This column is displayed by
                default.

            "dof":  the degrees of freedom for the uncertainty, displayed by default
                if any uncertainty has finite degrees of freedom

            "type":  the uncertainty type, displayed by default if any gummy has
                a type defined

            "s" or "significance":  the sensitivity coefficient (below) multiplied
                by the standard uncertainty, displayed by default

            "d", "derivative" or "partial":  the partial derivative of the y gummy
                with resect to the gummy in that row

            "c" or "sensitivity coefficient":  the absolute value of "d"

            "custom":  value given in the custom parameter list, displayed by default
                if the custom parameter is not None

            The columns displayed can also be set with the `columns` property.

        column_names:  `dict` or `None`, optional
            Names to display as column headers, if this is `None` then the default
            names are used.  The dictionary should use as keys any of the column
            names listed above in the columns parameter description and as values
            the desired heading for this column.  The column names can also be set
            with the `column_names` property.

        show_subtotals: `bool`, optional
            If any uncertainty types are defined, the combined standard uncertainty
            for each type is displayed in the table.  This can also be changed by
            setting the `show_subtotals` attribute.  The default is `True`.

        show_expanded_u:  `bool` or `None`, optional
            Whether or not to display the expanded uncertainty in the last row.  If
            this is `None`, then the expanded uncertainty is displayed if `
            `self.k != 1``. This can also be changed by setting the `show_expanded_u`
            attribute.

        show_s:  `bool`, optional
            Whether or not to show the significance column.  This is ignored if the
            `columns` parameter is not `None`.  The default can be changed by setting
            the attribute class attribute `Budget.show_s`.  The default value is `True`.

        show_d:  `bool`, optional
            Whether or not to show the partial derivatives column.  This is ignored
            if columns is not None.  The default can be  changed by setting the
            attribute class attribute `Budget.show_d`.  The default is `True`.

        show_c:  `bool`, optional
            Whether or not to show the sensitivity coefficient column.  This is ignored
            if columns is not None.  The default can be changed by setting the
            attribute class attribute `Budget.show_c`.  Teh default value is `False`

        units_on_values:  `bool` or `None`
            If this is `True`, units are shown in the value and u columns and if `False`
             the units are in a separate column.  If `None` then the units are in a
             separate column unless  any gummy in `xlist` or `y` has a `uunit` defined.

        sim:  `bool`, optional
            If True, the combined uncertainty and partial derivatives will be calculated
            using Monte-Carlo data.  The default is `False`

        css:  `str` or `None`, optional
            A css header to be used when displaying the table in HTML format.  If this is
            `None` then `Budget.default_css` will be used.

        description:  `list` of `str` or `None`, optional
            An optional column of descriptions to be printed in the table.  This should
            be a description for `y` then for each `x`, and followed, optionally, by
            subtotal and expanded uncertainty descriptions.

        description_math_mode:  `bool`, optional
            If this is `False`, then when using a LaTeX format, the description is put
            in normal text mode rather than math mode.  The default is `False`

        custom:  `list` of `str` or `None`, optional
            An optional column of  additional information to be printed in the table.
            This should be a value for `y` then for each `x`, and followed, optionally,
            by subtotal and expanded uncertainty values.

        custom_heading:  `str` or `None`, optional
            A heading for the custom column.

        custom_math_mode:  `bool`
            If this is False, then when using a LaTeX format, the custom value is put
            in normal text mode rather than math mode.  The default value is `False`.

        solidus, mulsep, slashaxis:
            see gummy.solidus, gummy.mulsep and gummy.slashaxis
        """
        from .budget import Budget
        return Budget(self,xlist,**kwds)
        
    def tohtml(self,**kwds):
        """
        Returns a string with the value as an html fragment.  All parameters
        are optional.  Any parameters that are not `None override the corresponding
        attributes of the calling gummy.  This is equivalent to the gummy.tostring
        method with the `fmt` parameter set to 'html'.
        """
        return super().tohtml(**kwds)

    def tolatex(self,**kwds):
        """
        Returns a string with the value as an LaTeX fragment.  It is assumed
        that LaTeX is in math mode.  All parameters are optional.  Any parameters
        that are not None override the corresponding attributes of the calling
        gummy.  This is equivalent to the `gummy.tostring` method with the `fmt`
        parameter set to 'latex'.
        """
        return super().tolatex(**kwds)

    def toascii(self,**kwds):
        """
        Returns a string with the value formatted so that only ASCII characters
        are used.  All parameters are optional.  Any parameters that are not
        `None` override the corresponding attributes of the calling gummy.  This
        is equivalent to the `gummy.tostring` method with the `fmt` parameter set
        'ascii'.
        """
        return super().toascii(**kwds)

    def latex(self,**kwds):
        """
        Prints the gummy using LaTeX if this method is executed from an Ipython
        console or from a Jupyter or Ipython notebook.  All parameters are
        optional.  Any parameters that are not `None` override the corresponding
        attributes of the calling gummy.
        """
        super().latex(**kwds)

    def html(self,**kwds):
        """`
        Prints the gummy using HTML if this method is executed from an Ipython
        console or from a Jupyter or Ipython notebook.  All parameters are
        optional.  Any parameters that are not `None` override the corresponding
        attributes of the calling gummy.
        """
        super().html(**kwds)

    def unicode(self,**kwds):
        """
        Prints the gummy.  All parameters are optional.  Any parameters that are
        not `None` override the corresponding attributes of the calling gummy.
        """
        super().unicode(**kwds)

    def ascii(self,**kwds):
        """
        Prints the gummy formatted using only ASCII characters.  All parameters
        are optional.  Any parameters that are not `None` override the corresponding
        attributes of the calling gummy.
        """
        super().ascii(**kwds)
                    
    def _retsomthng(self,fmt='unicode'):
        # if tostring encounters and exception, return as much as possible
        if self.exception_on_fmt_error:
            raise
        if self.name is None:
            try:
                return(str(self.x) + '{' + str(self.u) + '}' + '??')
            except:
                try:
                    return(str(self.x) + '(??)')
                except:
                    try:
                        return('??{' + str(self.u) + '}')
                    except:
                        return('??')
        else:
            try:
                return(str(self.get_name()).strip() + ' = ' + str(self.x) + '{' + str(self.u) + '}' + '??')
            except:
                try:
                    return(str(self.get_name()).strip() + ' = ' + str(self.x) + '{??}')
                except:
                    try:
                        return(str(self.get_name()).strip() + ' = ??{' + str(self.u) + '}')
                    except:
                        try:
                            return(str(self.get_name()).strip() + ' = ' + '??')
                        except:
                            return('??')
        
    def tostring(self,fmt=None,style=None,k=None,p=None,show_k=None,
                 show_p=None,show_dof=None,show_name=None,name=None,
                 norm=None,raw=False,nsig=None,solidus=None,
                 mulsep=None,**kwds):
        """
        Returns a string displaying the value of the gummy in the desired format.
        The `fmt` parameter is a string with the value in {"unicode","latex",
        "html","ascii"} or `None`.  `fmt` will default to 'ascii' if `self.printer`
        is 'ascii' or 'unicode' otherwise.  Any other parameters that are not
        `None` override the corresponding attributes of `self`.
        """
        try:
            if fmt is None:
                if self.printer == 'ascii':
                    fmt = 'ascii'
                else:
                    fmt = 'unicode'
                    
            if norm is None:
                norm = type(self).latex_norm
            if k is not None and k != self._k:
                if p is not None:
                    raise ValueError('k and p may not both be specified')
                g = self.copy(formatting=True)
                g.k = k
                return g.tostring(fmt=fmt,style=style,show_k=show_k,
                                  show_p=show_p,show_dof=show_dof,
                                  norm=norm,raw=raw,nsig=nsig)
            if p is not None and p != self.p:
                g = self.copy(formatting=True)
                g.p = p
                return g.tostring(fmt=fmt,style=style,show_k=show_k,
                                  show_p=show_p,show_dof=show_dof,
                                  norm=norm,raw=raw,nsig=nsig)
            if fmt is None or fmt == '':
                fmt = 'unicode'
            else:
                fmt = fmt.strip().lower()
                
            if style is not None:
                style = gummy._get_style(style)
            
            if self.x is np.ma.masked:
                return '--'        
            
            if fmt not in {'html','latex','utf-8','unicode','ascii'}:
                raise ValueError('format ' + str(fmt) + ' is unrecognized.')
                
            if hasattr(self.unit,'_format_xu'):
                v = self._unit._format_xu(self,fmt,style,norm,nsig,solidus=solidus,mulsep=mulsep)
            else:
                v = self._format_xu(fmt,style,norm,nsig,solidus=solidus,mulsep=mulsep)

            if raw:
                return v
                
            style = v[0]
            
            if style == 'override':
                # Allow a non-linear unit to override the regular formatting.
                return v[1]
                
            if show_name is None:
                show_name = self.show_name
            if show_name:
                if name is None:
                    name = self.get_name(fmt=fmt,norm=norm)
                if name is None:
                    name = ''
                else:
                    name += ' = '
            else:
                name = ''
                
            if style == 'x' or style == 'xsim':
                txt = v[1][0] + v[1][1] + v[1][2]
                txt = txt.strip()
                return name + txt
            elif style == 'xf' or style == 'xfsim':
                txt = v[1][0] + v[1][1]
                txt = txt.strip()
                return txt
            elif style == 'xunit':
                return v[1][2]
            elif style == 'u' or style == 'usim':
                txt = v[2][0] + v[2][1] + v[2][2]
                txt = txt.strip()
                return txt
            elif style == 'uf' or style == 'ufsim':
                txt = v[2][0] + v[2][1]
                txt = txt.strip()
                return txt
            elif style == 'uunit':
                if len(v) < 3 or len(v[2]) < 3 or v[2][2] is None:
                    return ''
                return v[2][2]
            elif style == 'ueq':
                if self._ubreakdown is None or len(self._ubreakdown) != len(v) - 2:
                    if self.k == 1:
                        if fmt == 'html':
                            pm = ' with <i>u</i>&nbsp;=&nbsp;'
                        elif fmt == 'latex':
                            pm = norm(' with ') + 'u = '
                        else:
                            pm = ' with u = '
                    else:
                        if fmt == 'html':
                            pm = ' with <i>U</i>&nbsp;=&nbsp;'
                        elif fmt == 'latex':
                            pm = norm(' with ') + 'U = '
                        else:
                            pm = ' with U = '
                    txt = v[1][0] + v[1][1] + v[1][2] + pm + v[2][0] + v[2][1] + v[2][2]
                
                else:
                    txt = v[1][0] + v[1][1] + v[1][2]
                    for i,t in enumerate(v[2:]):
                        b = self._ubreakdown[i]
                        if not isinstance(b,str):
                            if isinstance(b,gummy):
                                nm = b.get_name(fmt=fmt,norm=norm)
                                if nm is not None:
                                    b = nm
                                else:
                                    b = str(i+1)
                            else:
                                b = str(i+1)
                        if fmt == 'html':
                            if i == 0:
                                pm = ' with <i>u</i><sub>'+b+'</sub>&nbsp;=&nbsp;'
                            elif i == len(v) - 3:
                                pm = ' and <i>u</i><sub>'+b+'</sub>&nbsp;=&nbsp;'
                            else:
                                pm = ', <i>u</i><sub>'+b+'</sub>&nbsp;=&nbsp;'
                        elif fmt == 'latex':
                            if i == 0:
                                pm = norm(' with ') + 'u_{' + norm(b) + '} = '
                            elif i == len(v) - 3:
                                pm = norm(' and ') + 'u_{' + norm(b) + '} = '
                            else:
                                pm = norm(', ') + 'u_{' + norm(b) + '} = '
                        else:
                            if i == 0:
                                pm = ' with u(' + b + ') = '
                            elif i == len(v) - 3:
                                pm = ' and u(' + b + ') = '
                            else:
                                pm = ', u(' + b + ') = '
                        txt += pm + t[0] + t[1] + t[2]
            elif style == 'pm':
                if fmt == 'html':
                    pm = ' &plusmn; '
                elif fmt == 'latex':
                    pm = r' \,\pm\, '
                elif fmt == 'ascii':
                    pm = ' +/- '
                else:
                    pm = ' \u00B1 '
                    
                txt = v[1][0]
                for t in v[2:]:
                    txt += pm + t[0]
                if v[1][1] != '' or v[2][2] != '':
                    txt = '(' + txt + ')' + v[1][1] + v[1][2]
            elif style == 'pmsim':
                txt = v[1][0] + ' + ' + v[2][0] + ' - ' + v[3][0]
                if v[1][1] != '' or v[2][2] != '':
                    txt = '(' + txt + ')' + v[1][1] + v[1][2]
            elif style == 'pmi':
                if fmt == 'html':
                    pm = ' &plusmn; '
                elif fmt == 'latex':
                    pm = r' \,\pm\, '
                elif fmt == 'ascii':
                    pm = ' +/- '
                else:
                    pm = ' \u00B1 '
                txt = v[1][0] + v[1][1] + v[1][2] 
                for t in v[2:]:
                    txt += pm + t[0] + t[1] + t[2]
            elif style == 'pmsimi':
                txt = v[1][0] + v[1][1] + v[1][2] 
                txt += ' + ' + v[2][0] + v[2][1] + v[2][2]
                txt += ' - ' + v[3][0] + v[3][1] + v[3][2]
            elif style == 'cisim':
                txt = '[' + v[2][0] + v[2][1] + v[2][2]
                txt += ', ' + v[3][0] + v[3][1] + v[3][2] + ']'
            elif style == 'mcisim':
                if name == '':
                    if fmt == 'latex':
                        txt = norm('mean') + ' = '
                    else:
                        txt = 'mean = '
                else:
                    txt = ''
                txt += v[1][0] + v[1][1] + v[1][2]
                if fmt == 'latex':
                    txt += norm(', confidence interval') + ' = '
                else:
                    txt += ', confidence interval = '
                txt += '[' + v[2][0] + v[2][1] + v[2][2]
                txt += ', ' + v[3][0] + v[3][1] + v[3][2] + ']'
            elif style == 'concisef':
                txt = v[1][0]
                for t in v[2:]:
                    txt += '(' + t[0] + ')'
                return txt + v[1][1]
            else:
                # concise style:
                txt = v[1][0]
                for t in v[2:]:
                    txt += '(' + t[0] + ')'
                txt += v[1][1] + v[1][2]
                
            if style not in ['x','xf','u','uf','xsim','xfsim','usim','ufsim']:
                if show_k is None:
                    show_k = self.show_k
                    
                if style in ['pmsim','pmsimi','cisim','mcisim']:
                    setk = False
                else:
                    setk = self._set_k
                    
                if show_k is None:
                    show_k = (self.k != 1 and setk)
                else:
                    show_k = bool(show_k)
    
                if show_p is None:
                    show_p = self.show_p
                if show_p is None:
                    show_p = (self.k != 1 and not self._set_k)
                else:
                    show_p = bool(show_p)
                if show_p:
                    if self.p is None:
                        show_p = False
                    
                if show_dof is None:
                    show_dof = self.show_dof
                if show_dof is None:
                    show_dof = (self.dof < 10)
                else:
                    show_dof = bool(show_dof)
                    
                if style == 'ueq':
                    ns = 1
                else:
                    ns = 0
                if show_k:
                    ns += 1
                        
                    if style == 'ueq':
                        if show_p or show_dof:
                            itxt = ', '
                        else:
                            itxt = ' and '
                    else:
                        itxt = ' with '
                        
                    if fmt == 'latex':
                        itxt = norm(itxt)
                        
                    if style in ['pmsim','pmsimi','cisim','mcisim']:
                        k = self.ksim
                    else:
                        k = self.k
                    if fmt == 'html':
                        itxt += '<i>k</i>&nbsp;=&nbsp;' + _k_to_str(k)
                    else:
                        itxt += 'k = ' + _k_to_str(k)
                    
                    txt += itxt
                        
                if show_p:
                    ns += 1
                    p = self.p
                        
                    if show_k and show_dof:
                        itxt0 = ','
                    elif show_k and style == 'ueq':
                        itxt0 = ' and'
                    elif show_k or style == 'ueq':
                        itxt0 = ' and'
                    else:
                        itxt0 = ' with'
                        
                    ltxt = _p_to_str(fmt,p)
                    if ltxt[0] == '1' or ltxt[0] == '8':
                        itxt0 += ' an ' 
                        itxt1 = ' ' + self._p_method.text
                    else:
                        itxt0 += ' a ' 
                        itxt1 = ' ' + self._p_method.text
                    if fmt == 'latex':
                        itxt = norm(itxt0) + ltxt + norm(itxt1)
                    else:
                        itxt = itxt0 + ltxt + itxt1
                        
                    txt += itxt
                    
                if show_dof:                    
                    if (style == 'ueq') + show_k + show_p > 1:
                        itxt = ' and '
                    elif ns == 1:
                        itxt = ' and '
                    else:
                        itxt = ' with '
                        
                    if fmt == 'html':
                        itxt += '<i>&nu;</i>&nbsp;=&nbsp;' + _dof_to_str(self.dof,fmt)
                    elif fmt == 'latex':
                        itxt = norm(itxt)
                        itxt += r'\nu = ' + _dof_to_str(self.dof,fmt)
                    else:
                        itxt += _dof_to_str(self.dof) + ' degrees of freedom'
                    
                    txt += itxt
                    
            txt = txt.strip()
    
            if fmt == 'html':
                txt = '<span>' + txt + '</span>'
                
            return name + txt
        except:
            return self._retsomthng(fmt=fmt)
    
    def _format_xu(self,fmt,style,norm,nsig,xsig=None,solidus=None,mulsep=None):
        if style is None:
            style = self.style
            
        if nsig is None:
            nsig = self.nsig
        if nsig > self.max_digits - 2:
            nsig = self.max_digits - 2
        if nsig < 1:
            nsig = 1
            
        if solidus is None:
            solidus = self.solidus
            
        if mulsep is None:
            mulsep = self.mulsep 
            
        if style in  ['u','uf']:
            if isinstance(self._U,Quantity):
                un = self._U.unit
                u = self._U.value
            else:
                un = self.unit
                u = self._U
            if style == 'uf':
                un = 1
            uret =  gummy(u,unit=un)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1]
            return (style,('','',''),uret)
            
        if style in ['pmsim','cisim','pmsimi','mcisim','xsim','xfsim','usim','ufsim']:
            sim = True
            if self.simdata is None:
                if fmt == 'latex':
                    return ('override',norm('no simulated data'))
                return ('override','no simulated data')
            u = self.usim
        else:
            sim = False
            u = self.u
            
        if style in  ['usim','ufsim']:
            if isinstance(self._U,Quantity):
                un = self._U.unit
            else:
                un = self.unit
            u = self.usim
            if style == 'ufsim':
                un = 1
            uret =  gummy(u,unit=un)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1]
            return (style,('','',''),uret)
        
        xsym = _add_unit_sp(fmt,self.unit.tostring(fmt=fmt,solidus=solidus,mulsep=mulsep,strip=False))

        if sim:
            x = self.xsim
        else:
            x = self.x
            
        xd = _to_decimal(x,max_digits=self.max_digits)
        with localcontext(prec=nsig + 2):
            ud = _to_decimal(u)*_to_decimal(self._k)
        
        elip = ''
        uexp = 0
        xexp = 0
        dp = 0
        with localcontext(prec=self.max_digits):
            if ud == 0 or not ud.is_finite():
                if style != 'xf':
                    style = 'x'
                    
                if xd.is_finite():
                    nm = None
                    if xsig is None:
                        xsig = self.max_digits
                        if self.finfo is not None and self.finfo.precision > 0 and self.finfo.precision <= xsig:
                            xsig = self.finfo.precision
                            nm = True
                    else:
                        nm = False
                        
                    xd = round(xd,xsig-xd.adjusted()-1)
                    xd = round(xd,xsig-xd.adjusted()-1) # round again in case 9 round up to 10 in the last line
                    
                    if xd == 0:
                        xexp = 0
                    else:
                        xexp = xd.adjusted()
                    
                    if isinstance(x,Rational) and not isinstance(x,Integral) and x != xd:
                        fstr = str(x.numerator) + '/' + str(x.denominator)
                        if len(fstr) < self.max_digits + 1:
                            return (style,(fstr,'',xsym))
                        
                    if (nm is None and x == xd) or nm:
                        xd = xd.normalize()
                    elif nm is None:
                        elip = '...'
                
            else:
                ud = round(ud,nsig-ud.adjusted()-1)
                ud = round(ud,nsig-ud.adjusted()-1) 
                uexp = ud.adjusted()
    
                if xd.is_finite():
                    try:
                        if xsig is None:
                            xd = xd.quantize(ud)
                        else:
                            xd = round(xd,xsig-xd.adjusted()-1)
                    except InvalidOperation:
                        xd = round(xd,self.max_digits-xd.adjusted()-1)
                        xexp = xd.adjusted()
                        elip = '...'
                        if style not in ['x','xf','xsim','xfsim','ueq','u equals']:
                            style = 'pmi'
                            
                    if abs(xd) < abs(ud) and style not in ['x','xf','xsim','xfsim']:
                        dp = xd.adjusted() - uexp
                        xexp = uexp
                    else:
                        xexp = xd.adjusted()
                else:
                    if style in ['concise','pm']:
                        style = 'pmi'
                    elif style == 'pmsim':
                        style == 'pmsimi'
                    
        scin = False
        if self.sci_notation is None:
            if xexp > self.sci_notation_high or xexp < self.sci_notation_low:
                scin = True
        else:
            scin = self.sci_notation
        if xexp > self.max_digits - 1:
            scin = True
        elif xexp < 0 and len(xd.as_tuple().digits) + self.sci_notation_low > self.max_digits - 1:
            scin = True
        elif style == 'concise' and xd.is_finite() and (xd.as_tuple()[2] > 0 or len(xd.as_tuple()[1]) < nsig):
            scin = True
            
        if not scin:
            xexp = 0
            
        if xexp == 0:
            fstr =_decimal_str(xd,dalign=0,fmt=fmt)
        else:
            fstr = _decimal_str(xd,dplace=dp,fmt=fmt)
            
        xret = (fstr + elip,_format_exp(fmt,xexp),xsym)
                    
        if style in ['x','xf','xsim','xfsim']:
            return (style,xret)
                
        ugummy = isinstance(self._U,Quantity) or (self._Ubr is not None and isinstance(self._Ubr[0],Quantity))
        if ugummy:
            if sim:
                usm = self.Usim
                uret0 = gummy(usm[0],unit=self._U.unit)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1]
                uret1 = gummy(usm[1],unit=self._U.unit)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1]
                style = 'pmsimi'
                return (style,xret,uret0,uret1)
            else:
                if self._Ubr is None:
                    uret = [gummy(self._U)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1]]
                else:
                    uret = [gummy(i)._format_xu(fmt,'x',norm,nsig,xsig=nsig,solidus=solidus,mulsep=mulsep)[1] for i in self._Ubr]
                if style == 'pm' or style == 'concise':
                    style = 'pmi'
                return tuple([style,xret] + uret)

        with localcontext(prec=nsig + 2):
            if sim:
                u = abs(_to_decimal(self.cisim[1] - self.cisim[0]))
            else:
                u = abs(_to_decimal(self._U))
            
            if u.is_finite():
                if xd.is_finite():
                    uu = u.quantize(xd)
                    if uu == 0:
                        u = round(u,nsig-u.adjusted()-1)
                        u = round(u,nsig-u.adjusted()-1)
                    else:
                        u = uu
                else:
                    u = round(u,nsig-u.adjusted()-1)
                    u = round(u,nsig-u.adjusted()-1)
            
            if self._ubreakdown is None:
                ub = [u]
            else:
                ub = []
                for b in self._Ubr:
                    bb = abs(_to_decimal(b))
                    if bb.is_finite():
                        if xd.is_finite():
                            bb = bb.quantize(xd)
                        else:
                            bb = round(bb,nsig-bb.adjusted()-1)
                            bb = round(bb,nsig-bb.adjusted()-1)
                    ub.append(bb)

            
        uret = []
        if style == 'concise':
            for ue in ub:
                utxt = _decimal_str(ue,dplace=None,fmt=fmt)
                uret.append((utxt,'',xsym))
        elif style == 'pm':
            for ue in ub:
                utxt = _decimal_str(ue,fmt=fmt,dalign=xexp)
                uret.append((utxt,_format_exp(fmt,ue),xsym))
        else:
            if style == 'pmsim' or u == 0 or not u.is_finite():
                uexp = xexp
            else:
                uexp = u.adjusted()
                
            if self.sci_notation is None:
                if uexp > self.sci_notation_high or uexp < self.sci_notation_low:
                    scin = True
            else:
                scin = self.sci_notation
            if uexp > self.max_digits - 1:
                scin = True
            elif uexp < 0 and nsig + self.sci_notation_low > self.max_digits - 1:
                scin = True
                
            if not scin:
                uexp = 0
                
            if style in ['ueq','pmi']:
                if uexp == 0:
                    for ue in ub:
                        utxt = _decimal_str(ue,fmt=fmt,dalign=0)
                        uret.append((utxt,_format_exp(fmt,uexp),xsym))
                else:
                    for ue in ub:
                        utxt = _decimal_str(ue,fmt=fmt,dalign=uexp)
                        uret.append((utxt,_format_exp(fmt,uexp),xsym))
            elif style in ['pmsim','pmsimi']:
                u0 = _to_decimal(self.Usim[0]).quantize(u)
                u1 = _to_decimal(self.Usim[1]).quantize(u)
                if uexp == 0:
                    utxt0 = _decimal_str(u0,fmt=fmt,dalign=0)
                    utxt1 = _decimal_str(u1,fmt=fmt,dalign=0)
                    if utxt0 == utxt1:
                        if style == 'pmsim':
                            style = 'pm'
                        else:
                            style = 'pmi'
                        uret.append((utxt0,'',xsym))
                    else:
                        uret.append((utxt0,'',xsym))
                        uret.append((utxt1,'',xsym))
                else:
                    utxt0 = _decimal_str(u0,fmt=fmt,dalign=uexp)
                    utxt1 = _decimal_str(u1,fmt=fmt,dalign=uexp)
                    uetxt = _format_exp(fmt,uexp)
                    if utxt0 == utxt1:
                        if style == 'pmsim':
                            style = 'pm'
                        else:
                            style = 'pmi'
                        uret.append((utxt0,uetxt,xsym))
                    else:
                        uret.append((utxt0,uetxt,xsym))
                        uret.append((utxt1,uetxt,xsym))
                
            elif style in ['cisim','mcisim']:
                for ci in self.cisim:
                    cin = _to_decimal(ci).quantize(u)
                    dp = 0
                    if cin.is_finite():
                        xnexp = cin.adjusted()
                        if cin == 0:
                            xnexp += nsig - 1
                            dp = -nsig + 1
                    else:
                        xnexp = 0
                    if self.sci_notation is None:
                        scin = (xnexp > self.sci_notation_high or xnexp < self.sci_notation_low)
                    else:
                        scin = self.sci_notation
                    if cin == 0 and xnexp > nsig-1:
                        scin = True
                    if scin:
                        cin = _decimal_str(cin,fmt=fmt,dplace=dp)
                        xentxt = _format_exp(fmt,xnexp)
                        uret.append((cin,xentxt,xsym))                 
                    else:
                        uret.append(( _decimal_str(cin,fmt=fmt,dalign=0),'',xsym))
                    
        return tuple([style,xret]+uret)
        
    @classmethod
    def create(cls,x,u=0,unit=one,dof=float('inf'),k=1,p=None,uunit=None,
               utype=None,name=None,correlation_matrix=None,
               covariance_matrix=None):
        """
        A class method that creates a list of correlated gummys.
        
        Parameters
        ----------
        x:
            Either a list of floats corresponding to the x-value of each
            gummy or an and instance of a MultivariateDistribution sub-class.

        u, unit, dof, k, p, uunit, utype, and name:
            Lists that correspond to the parameters in the gummy
            initializer (with the i-th value  in each list passed to
            the initializer for the i-th gummy). With the exception of
            the "name" parameter, these may also be a single value with
            this same value is to passed to each initializer.

        correlation_matrix:
            A list or array to be used as the correlation matrix of the
            gummys. This is optional and must be set to the default value
            of `None` if the covariance_matrix is specified.

        covariance_matrix:
            A list or array to be used as the variance-covariance matrix
            of the gummys.  If the covariance matrix is`  specified the `u`
            parameter will be ignored  This parameter is optional and must
            be set to the default value of `None` if the `correlation_matrix`
            is specified.  If both the `correlation_matrix`  and the
            `covariance_matrix` are `None` (or omitted) then the gummys
            will be uncorrelated.
                
        Returns
        -------
        a list of gummys
        
        Notes
        -----
        This package does not implement a multivariate Student's t
        distribution that has differing degrees of freedom for each component.
        So if if the elements of dof are finite and not all the same and
        either a correlation_matrix or a covariance_matrix is defined, the
        joint distribution for Monte-Carlo calculations (but not first-
        order calculations) will default to a multivariate normal distribution.
        """
        
        if isinstance(x,MultivariateDistribution):
            ret = nummy.create(x,u=u,dof=dof,name=name,
                       correlation_matrix=correlation_matrix,
                       covariance_matrix=covariance_matrix)
            ret = [gummy(r) for r in ret]
                
            if unit is not one:
                if isinstance(unit,str) or isinstance(unit,Unit):
                    unit = [unit]*len(x)
                for i,r in enumerate(ret):
                    r._unit = Unit.unit(unit[i])
                
            if p is not None:
                if _isscalar(p):
                    p = [p]*len(x)
                for i,r in enumerate(ret):
                    r.p = p[i]
            else:
                if _isscalar(k):
                    k = [k]*len(x)
                for i,r in enumerate(ret):
                    r.k = k[i]
                    
            if uunit is not None:
                if isinstance(uunit,str) or isinstance(uunit,Unit):
                    uunit = [uunit]*len(x)
                for i,r in enumerate(ret):
                    r.uunit = uunit[i]
                
            return ret
                
        if correlation_matrix is not None or covariance_matrix is not None:
            if any([isinstance(v,Distribution) for v in x]):
                raise TypeError('Distribtuion instances may not be used in x if a correlation_matrix nor a covariance_matrix is defined')
            if dof is not None and not _isscalar(dof):
                raise TypeError('dof cannot be set individually of a correlation of covariance matrix is specified')
            if utype is not None and not isinstance(utype,str):
                raise TypeError('utype cannot be set individually of a correlation of covariance matrix is specified')


        if covariance_matrix is not None:
            if correlation_matrix is not None:
                raise TypeError('correlation_matrix and covariance_matrix cannot both be specified')
            if u != 0 and u is not None:
                raise TypeError('covariance_matrix and u cannot both be specified')
                
        n = len(x)
        
        if _isscalar(u):
            u = [u]*n
            
        if isinstance(unit,str) or isinstance(unit,Unit):
            unit = [unit]*n
            
        if dof is None:
            dof = [float('inf')]*n
        elif _isscalar(dof):
            dof = [dof]*n
            
        if p is None:
            p = [None]*n
        elif _isscalar(p):
            p = [p]*n
            
        if k is None:
            k = [1]*n
        elif _isscalar(k):
            k = [k]*n
            
        if uunit is None:
            uunit = [None]*n
        elif isinstance(uunit,str) or isinstance(uunit,Unit):
            uunit = [uunit]*n
            
        if utype is None:
            utype = [None]*n
        elif isinstance(utype,str):
            utype = [utype]*n
            
        if name is None:
            name = [None]*n
            
        ret = [cls(x[i],u=u[i],unit=unit[i],dof=dof[i],k=k[i],p=p[i],
                   uunit=uunit[i],utype=utype[i],name=name[i]) for i in range(n)]
            
        if correlation_matrix is not None or covariance_matrix is not None:
            x = [r.x for r in ret]
            u = [r.u for r in ret]
            nret = nummy.create(x,u=u,dof=dof[0],utype=utype[0],name=name,
                                correlation_matrix=correlation_matrix,
                                covariance_matrix=covariance_matrix)
            
            for n,r in zip(nret,ret):
                r._value = n
                if r.uunit is None:
                    r._U = _ku(r._k,r._value.u)
                else:
                    r._set_U(r._k,r.uunit)
            
        return ret
    
    @classmethod
    def apply(cls,function,derivative,*args):
        """
        A classmethod that applies a function to one or more gummy or jummy 
        objects propagating the uncertainty.
        
        Parameters
        ----------
        function: `function`
              The the function to be applied. For `gummy.apply`, 'function'
              should take one or more float arguments and return a float value 
              or float array.  For `jummy.apply`, 'function' may also take and
              return complex values.

        derivative:  `function`
              The name of a second function that gives the derivatives
              with respect to the arguments of `function`.  `derivative` should
              take an equal number of arguments as `function`.  If `function`
              takes one argument `derivative` should return a float and if
              `function` takes more than one argument then `derivative` should
              return a tuple, list or array of floats that contains the derivatives
              with respect to each argument.  In the case of `jummy.apply`, the
              derivatives with respect to each argument may be real or complex
              values, in which case `function` is assumed to be holomorphic.  Or
              the derivative may be a 2 x 2 matrix of the form:

                              [[ du/dx, du/dy ],
                               [ dv/dx, dv/dy ]]

             where function(x + j*y) = u + j*v.

        *args:  `gummy`, `jummy`, or `float`
              One or more arguments to which `function` will be applied.  These
              arguments need not all be `Dfunc` objects; arguments  such as
              floats will be taken to be constants with no uncertainty.
              They may also be numpy ndarrays in which case the usual numpy
              broadcasting rules apply.
              
        Returns
        -------
        `gummy`, `jummy`:
            If none of the arguments are `gummy` or `jummy`
            then the return value is the same type as the return value of `function`.
            Otherwise `gummy.apply` returns a `gummy` and `jummy.apply` returns either a
            `gummy` or a `jummy` depending on whether `function` has a float or
            a complex return value.
            
        
        Examples
        --------
            
        >>> import numpy as np
        >>> x = gummy(0.678,u=0.077)
        >>> gummy.apply(np.sin,np.cos,x)
        0.627 +/- 0.060
        
        >>> x = gummy(1.22,u=0.44)
        >>> y = gummy(3.44,u=0.67)
        >>> def dhypot(x,y):
        ...     return (x1/sqrt(x1**2 + x2**2),x2/np.sqrt(x1**2 + x2**2))
        >>> gummy.apply(np.hypot,dhypot,x,y)
        3.65 +/- 0.65
        """
        
        return cls._apply(function,derivative,*args)
    
    @classmethod
    def _apply(cls,function,derivative,*args,fxdx=None):
        
        if fxdx is None:
            args = [a.convert(one).value if isinstance(a,Quantity) else a for a in args]
            x = [a.x if isinstance(a,ummy) else a for a in args]
            fx = function(*x)
            d = derivative(*x)
        else:
            fx,d,x = fxdx
        
        if not _isscalar(fx):
            return [cls.apply(lambda *y: function(*y)[i],
                               lambda *y: derivative(*y)[i],
                               *args,fxdx=(fx[i],d[i],x)) 
                    for i in range(len(fx))]
        r = nummy._apply(function,derivative,*args,fxdx=(fx,d,x))
        return cls(r)
        
    @classmethod
    def napply(cls,function,*args,fxx=None):
        """
        gummy.napply(function, arg1, arg2, ...) and
        jummy.napply(function, arg1, arg2, ...)
        
        A classmethod that applies a function to one or more gummy or jummy 
        objects propagating the uncertainty.  This method is similar to apply 
        except that the derivatives are computed numerically so a derivative 
        function does not need to be supplied.
        
        Parameters
        ----------
        function: `function`
            The the function to be applied. For `gummy.apply`, 'function'
            should take one or more float arguments and return a float value
            r float array.  For `jummy.apply`, 'function' may also take and
            return complex values.

        *args:  `gummy`, `jummy`, or `float`
              One or more arguments to which `function` will be applied.  These
              arguments need not all be `Dfunc` objects; arguments  such as
              floats will be taken to be constants with no uncertainty.
              They may also be numpy ndarrays in which case the usual numpy
              broadcasting rules apply.

        Returns
        -------
        `gummy`, `jummy`:
            If none of the arguments are `gummy` or `jummy`
            then the return value is the same type as the return value of `function`.
            Otherwise `gummy.apply` returns a `gummy` and `jummy.apply` returns either a
            `gummy` or a `jummy` depending on whether `function` has a float or
            a complex return value.
            
        
        Examples
        --------
            
        >>> import numpy as np
        >>> x = gummy(0.678,u=0.077)
        >>> gummy.napply(np.sin,x)
        0.627 +/- 0.060
        
        >>> x = gummy(1.22,u=0.44)
        >>> y = gummy(3.44,u=0.67)
        >>> gummy.napply(np.hypot,x,y)
        3.65 +/- 0.65
        """
        
        return cls._napply(function,*args)
        
    @classmethod
    def _napply(cls,function,*args,fxx=None):
        if fxx is None:
            args = [a.convert(one).value if isinstance(a,Quantity) else a for a in args]
            x = [a.x if isinstance(a,ummy) else a for a in args]
            fx = function(*x)
        else:
            fx,x = fxx
            
        if not _isscalar(fx):
            return [cls.napply(lambda *y: function(*y)[i],*args,fxx=(fx[i],x)) 
                    for i in range(len(fx))]
        return cls(nummy._napply(function,*args,fxx=fxx))
    
    def __array_ufunc__(self,ufunc,method,*args,**kwds):
        if method != '__call__':
            return None
        
        if any(isinstance(a,np.ndarray) for a in args):
            args = [np.array(a) if isinstance(a,gummy) else a for a in args]
            return ufunc(*args)
        
        return super().__array_ufunc__(ufunc,method,*args,**kwds)
                
    def _nprnd(self,f):
        ret = self._value._nprnd(f)
        ret._unit = self._unit
        return ret
    
    def __add__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) + v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self) + v
        
        return super().__add__(v)
    
    def __radd__(self,v):
        if isinstance(v,np.ndarray):
            return v + np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v + jummy(self)
        
        return super().__radd__(v)
    
    def __sub__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) - v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self) - v
        
        return super().__sub__(v)
    
    def __rsub__(self,v):
        if isinstance(v,np.ndarray):
            return v - np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v - jummy(self)
        
        return super().__rsub__(v)
    
    def __mul__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)*v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self)*v
        
        return super().__mul__(v)
    
    def __rmul__(self,v):
        if isinstance(v,np.ndarray):
            return v*np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v*jummy(self)
        
        return super().__rmul__(v)
    
    def __truediv__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)/v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self)/v
        
        return super().__truediv__(v)
    
    def __rtruediv__(self,v):
        if isinstance(v,np.ndarray):
            return v/np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v/jummy(self)
        
        return super().__rtruediv__(v)
        
    def __pow__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self)**v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self)**v
        
        return super().__pow__(v)
    
    def __rpow__(self,v):
        if isinstance(v,np.ndarray):
            return v**np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v**jummy(self)
        
        return super().__rpow__(v)
        
    def __floordiv__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) // v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self) // v
        
        return super().__floordiv__(v)
        
    def __rfloordiv__(self,v):
        if isinstance(v,np.ndarray):
            return v // np.array(self)
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return v // jummy(self)
        
        return super().__rfloordiv__(v)
        
    def __mod__(self,v):
        if isinstance(v,np.ndarray):
            return np.array(self) % v
        
        if isinstance(v,Complex) and not isinstance(v,Real):
            return jummy(self) % v
        
        return super().__mod__(v)
    
    def __rmod__(self,v):
        if isinstance(v,np.ndarray):
            return v % np.array(self)
        
        return super().__rmod__(v)
    
    @property
    def real(self):
        """
        returns a copy of the gummy
        """
        return self.copy(formatting=False)
    
    @property
    def imag(self):
        if self.unit.linear:
            s = self
        else:
            s = self.tobaseunit()
        
        return type(self)(0,unit=s.unit)
    
    def conjugate(self):
        """
        returns a copy of the gummy
        """
        return self.copy(formatting=False)
    
    def angle(self):
        if self.x >= 0:
            return type(self)(0)
        else:
            return type(self)(np.pi)
        

class jummy(immy):
    
    show_name = True
    
    _element_type = gummy
    
    def __init__(self,real=None,imag=None,r=None,phi=None,cov=None,name=None):
        """
        A jummy object represents a complex valued quantity with `gummy`
        real and imaginary components.
        
        Parameters
        ----------
        real, imag, r, phi:  `float` or `gummy`
            The value may be specified in  either cartesian coordinates
            using `real` and `imag` or polar coordinates with `r` and `phi`.
            The pair `real`, `imag` or `r`, `phi` may both be `gummy` or
            both be `float`.  If they are `float` then `cov` may
            also be specified.
                
        cov:  2 x 2 array_like of `float`, optional
            The variance-covariance matrix for either the pair `real`,
            `imag` or the pair `r`, `phi`.
            
        name:  `str`
            An optional name for the jummy.
        """
        self.name = name
        super().__init__(real=real,imag=imag,r=r,phi=phi,cov=cov)
        if self._ridef:
            try:
                if self._real.unit is not self._imag.unit:
                    # check that the units are compatible
                    self._real.convert(self._imag.unit)
            except:
                raise IncompatibleUnitsError('the real an imaginary parts must have compatible units')
        else:
            if not self._phi.unit.is_dimensionless:
                raise IncompatibleUnitsError('phi must have dimensonless units')
            
    def tostring(self,fmt='unicode',norm=None,show_name=None,name=None,
                 style=None,**kwds):
        if show_name is None:
            show_name = self.show_name
        if show_name:
            if name is None:
                name = self.get_name(fmt)
            if name is None:
                name = ''
            else:
                name += ' = '
        else:
            name = ''
            
        if style is None:
            style = self.style
        if style == 'polar':
            r = self.r.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
            
            if self.phi.unit.linear:
                if self.phi.x < 0:
                    i = abs(self.phi).tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
                    sign = '-'
                else:
                    i = self.phi.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
                    sign = ''
            else:
                i = self.phi.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
                if i.startswith('-'):
                    i = i[1:]
                    sign = '-'
                else:
                    sign = ''
                    
            if fmt == 'html':
                if sign == '':
                    sign ='&thinsp;'
                ret = name + r + '&thinsp;&middot;&thinsp;<i>e</i><sup>'
                ret += sign + '<i>' + self._imag_symbol + '</i>&thinsp;'
                ret += i + '</sup>'
            elif fmt == 'latex':
                ret = name + r + '\\,\\cdot\\,e^{' + sign + self._imag_symbol + i + '}'
            else:
                ret = name + r + ' exp(' + sign + self._imag_symbol + i + ')'
            return ret
                
        r = self.real.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
        if self.imag.unit.linear and self.imag.value == 0:
            return r
        if self.real.unit.linear and self.real.value == 0:
            r = ''
        
        if self.imag.unit.linear:
            if self.imag.x < 0:
                i = abs(self.imag).tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
                if r == '':
                    sign = '-'
                else:
                    sign = ' - '
            else:
                i = self.imag.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
                if r == '':
                    sign = ''
                else:
                    sign = ' + '
        else:
            i = self.imag.tostring(fmt=fmt,style='concise',k=1,norm=norm,**kwds)
            if i.startswith('-'):
                i = i[1:]
                if r == '':
                    sign = ' - '
                else:
                    sign = '-'
            else:
                if r == '':
                    sign = ''
                else:
                    sign = ' + '
            
        if fmt == 'html':
            i = '<i>' + self._imag_symbol + '</i>&thinsp;' + i
        elif fmt == 'latex':
            i = self._imag_symbol + '\\,' + i
        else:
            i = self._imag_symbol + i
        
        return name + r + sign + i

    def toimmy(self):
        """
        returns an immy representation of the jummy
        """
        return immy(real=self.real.toummy(),imag=self.imag.toummy())
    
    def splonk(self):
        """
        splonks the jummy
        """
        return self.toummy().splonk()
    
    @property
    def name(self):
        if self._name is None:
            return None
        if isinstance(self._name,str):
            return self._name
        return self._name[0]
    @name.setter
    def name(self,v):
        if v is None:
            self._name = None
            return
        elif isinstance(v,str):
            self._name = v.strip()
            return
        
        try:
            if len(v) != 4:
                raise ValueError('the name must be a string or a length 4 tuple or str')
        except TypeError:
            raise ValueError('the name must be a string or a length 4 tuple of str')
            
        try:
            n = v[0].strip()
            self._name = tuple([n if e is None else e.strip() for e in v])
        except AttributeError:
            raise ValueError('the name must be a string or a length 4 tuple of str')
            
    def get_name(self,fmt='unicode',norm=None):
        if norm is None:
            norm = type(self).latex_norm
        return get_name(self._name,fmt,norm)
    