.. _hand_made_doc:

MetroloPy
=========

Most of the functionality in the MetroloPy package lies in the
gummy_ object.

.. _gummy:

class gummy
===========

class metrolopy.\ **gummy(x, u=0, unit=one, dof=float('inf'), k=1, p=None, 
p_method=None, uunit=None, utype=None, name=None)**

A gummy object represents a numerical value with an uncertainty and (or)
a unit. They can be used in place of float values in Python expressions
and the uncertainty can be propagated with with both first-order and
Monte-Carlo methods. Correlations between gummys are tracked. gummys
must be real valued; the jummy_ class is an extension of the
gummy class to complex numbers.

gummy parameters
----------------
(all parameters except *x* are optional)

-  **x**: (float or Distribution_ ) Either a
   number representing the value of the gummy or a
   Distribution_ instance that represents the
   probability distribution of the gummy. If *x* is a
   Distribution_ any value given for *u* or *dof*
   will be ignored. The same Distribution_ instance
   should never be used as the *x* parameter for two different gummys,
   e.g. to define a gummy with a :ref:`uniform distribution<UniformDist>`, 
   use::

       >>> g = gummy(UniformDist(center = 2.25, half_width = 0.12))

   Note that if *x* is a float value, the uncertainty distribution is
   assumed to be either a :ref:`Normal distribution<NormalDist>` or a
   :ref:`shifted and scaled Students's t distribution<TDist>` depending on
   the value of the *dof* parameter.

-  **u**: (float >= 0) A number representing the uncertainty in *x*.
   By default *u* is taken to be the standard ("1-sigma") uncertainty, 
   however if *k* or *p* are specified then *u* is taken to be the
   corresponding expanded uncertainty. Also by default, when the *uunit*
   parameter is ommitted, *u* is assumed to have the same
   units as *x*. Note that the :ref:`U<upper-u>` property of the newly
   created gummy will be equal to the value of the *u* parameter. The
   :ref:`u<lower-u>` property of the newly created gummy will only be equal to 
   the *u* parameter if the *k*, *p* and *uunit* parameters are omitted.

-  **unit**: (str or :ref:`Unit<class-Unit>` ) The units of the gummy. This
   may be a :ref:`Unit<class-Unit>` object or, more commonly, a string that
   references a Unit object, e.g. ``gummy(1,unit='kg')`` or
   ``gummy(3.67,u=0.22,unit='m/s')``. The default is the Unit instance
   *one*. Spaces or asterisks ('*') may be used to represent unit
   multiplication, a forward slash ('/') division and double asterisks
   ('**') powers, e.g ``gummy(3.67,u=0.22,unit='m s**-1')``.  If a unit
   name contains a space, '\*' or '/' character then the name must be
   enclosed in square brackets, e.g: ``'[light year]'``. 
   
-  **dof**: (float or int > 0) The number of degrees of freedom upon
   which the uncertainty is based. The default is ``float('inf')``.

-  **k**: (float > 0) The coverage factor for *u*. The value of the
   *u* parameter is divided by the coverage factor to get the standard
   uncertainty for the new gummy. If the paramter *p* is specified then
   the coverage factor is calculated using *p* and the value of the *k*
   parameter is ignored. The default for *k* is 1. For example::

       >>> g = gummy(3.24,u=0.22,k=2)
       >>> g
       3.24(22) with k = 2.0
       >>> g.U
       0.22
       >>> g.u
       0.11    

-  **p**: (float between 0 and 1 or None) The level of
   confidence or probability level for *u*. If this parameter is
   specified, then the coverage factor *k* is calculated using *dof* and
   the method specified by the *p_method* parameter. The standard
   uncertainty for the gummy is then set to *u* divided by the coverage
   factor *k*. For example::

       >>> g = gummy(3.24,u=0.22,p=0.95)
       >>> g
       3.24(22) with a 95% level of confidence
       >>>g.k
       1.959963984540054
       >>> g.U
       0.22
       >>> g.u
       0.11224696052342387

.. _p_method:

-  **p_method**: ({'loc', 'cp', 'gauss',
   'ccp', 'chebyshev', None}) If the *p* parameter is
   specified, *p_method* sets the method that is used to calculate the
   coverage factor. If p_method is omitted or set to 'loc', then the
   uncertainty is assumed to be represented by a normal probability
   distribution if *dof* = float('inf') and shifted and scaled Student's
   t distribution otherwise. If *p_method* = 'gauss' or 'cp' then the
   Guass inequality is used, and if *p_method* = 'chebyshev' or 'ccp'
   then the Chebyshev inequality is used. For *p* = 0.95 and *dof* =
   float('inf'), *p_method* = 'loc' gives *k* = 2.0, while *p_method*
   = 'gauss' gives *k* = 3.0 and *p_method* = 'chebyshev' gives *k* =
   4.5.
   
.. _parameter-uunit:

-  **uunit**: (str or :ref:`Unit<class-Unit>` ) This represents the units of
   *u*. It may be a unit with the same dimension as the *unit*
   parameter, e.g. a measurement result of 3 m with an uncertainty of 1
   mm can be represented by ``gummy(3,0.001,unit='m')`` or equivalently
   ``gummy(3,1,unit='m',uunit='mm')`` The *uunit* parameter can also be
   a dimensionless unit if *u* represents a relative uncertainty, e.g.
   the gummy above can be also represented by
   ``gummy(3,0.1,unit='m',uunit='%')``. If *uunit* is set to None, then
   the units of *u* are taken to be the same as those of x (as given by
   the *unit* parameter). The default is None. A
   ``NoUnitConversionFoundError`` exception will be generated if *uunit*
   is not dimensionless and no conversion exists between *uunit* and
   *unit*.

.. _utype:

-  **utype**: (str) An arbitrary string value labeling the
   uncertainty type. When a calculation is performed with gummys, the
   combined uncertainty of effective degrees of freedom from one
   particular uncertainty type can be found in the calculation result
   with the ufrom and doffrom methods. E.g. you can create a set of
   gummys with uncertainties assigned either utype "A" or utype "B",
   insert them into a measurement equation and find the combined utype
   "A" uncertainty.

.. _parameter-name:

-  **name**: (str) An arbitrary string naming the gummy. The name is
   used when displaying the gummy value or labeling plot axes and serves
   no other function.

basic gummy properties
----------------------

.. _x:

-  **x**: (read-only) Gets the value of the gummy. This property is
   read-only, but changing the :ref:`unit<property-unit>` property will
   change *x*.

.. _lower-u:

-  **u**: (read-only) Gets the standard uncertainty of the gummy in
   the units set by the :ref:`unit<property-unit>` property. Note that
   setting the uunit_ property only affects the value
   of the *U* property and not the *u* property.

.. _const:

-  **const**: (read-only) Return True if *u* == 0 and False otherwise.

.. _dof:

-  **dof**: (read-only) Gets the degrees of freedom associated with *u*.

.. _property-utype:

-  **utype**: (read-only) Gets the the uncertainty type. See the
   utype_ parameter.

.. _upper-u:

-  **U**: (read-only) gets the "expanded" uncertainty, this will
   depend of the values of the :ref:`unit<property-unit>`, *k*, *p*, and
   *p_method* properties, also *U* may be expressed in different units
   from *x* and *u* and may be a relative uncertainty, see the
   uunit_ property.

.. _k:

-  **k**: gets or sets the coverage factor for the expanded
   uncertainty *U*, *U* = *k*\ \*\ *u*, setting the *p* property will
   change the value of this parameter

.. _p:

-  **p**: gets or sets the level of confidence for the expanded
   uncertainty, changing this property will change the *k* property,
   the relation between the value of this property and the property *k*
   is defined by the *p_method* property
   
.. _property-p-method:

-  **p_method**: see the p_method_ parameter

.. _name:

-  **name**: An arbitrary string naming the gummy. The name is used when
   displaying the gummy value or labeling plot axes and serves no other
   function.

-  **ubreakdown**: (list of str or None, default is
   None) If this is set to a list containing strings referencing
   :ref:`utypes<utype>` then when the gummy is printed, the
   uncertainty from each utype will be displayed separately.
   Example::

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

-  **independent**: (read-only) Returns ``False`` if the gummy was
   created from a mathematical operation involving other gummys and
   ``True`` otherwise.

-  **real**: returns a ``self``

-  **imag**: returns a zero value gummy

-  **bayesian**: (bool) Read/write at the class level, but read-only
   at the instance level. The default value is ``False``; this property
   should only be changed once at the beginning of the session. This
   property affects how a standard uncertainty based on data with finite
   degrees of freedom is defined and thus how the level of confidence
   p_ (sometimes called coverage probability) of an
   expanded uncertainty is related to the coverage factor
   k_. Standard uncertainties are often based on the
   standard deviation of a set of measurements (and the assumption that
   these measurements are drawn from a normally distributed population).
   Traditionally (e.g. the GUM 2008 edition) the standard uncertainty is
   taken to be the standard deviation of the mean (*s*/sqrt(*n*), where *s* is
   the sample standard deviation and *n* is the number of measurements).
   However there is some "extra uncertainty" because the sample standard
   deviation does not exactly equal to the population standard deviation.
   This is taken into account by using a Student's *t* distribution to
   calculate the expanded uncertainty. However it has been pointed out,
   by those who advocate a Bayesian point of view, that the probability
   distribution for the measurand here is best described by a shifted
   and scaled Student's *t* distribution. So the standard uncertainty
   should be the standard deviation of the Student's *t* distribution which is
   *s*\*sqrt{(*n*-1)/[*n*\*(*n*-3)]}. Thus the relationship between the Bayesian
   and traditional standard uncertainty definitions is::

       u(bayesian) = [dof/(dof - 2)]*u(traditional)

   where *dof* = *n* - 1 and the "extra uncertainty" of the traditional method
   is incorporated directly into the standard uncertainty.
   
-  **rounding_u**: (bool)  Set at the class level.  If this is set to
   ``True`` then uncertainty is added to account for floating point rounding
   errors.   An uncertainty proportional to the machine epsilon is added to
   the uncertainty whenever a gummy is created with a floating point data type. 
   Then this uncertainty is propagated like any other uncertainty. This 
   can give some idea of the magnitude of the floating point errors, but is 
   not a substitute for a full numerical error analysis.  The default value
   is ``False``.

basic gummy methods
-------------------

-  **correlation(g)**: Returns the correlation coefficient between
   ``self`` and *g*.

-  **covariance(g)**: Returns the covariance between
   ``self`` and *g*.

-  static **correlation_matrix(gummys)**: Returns the correlation matrix of a
   list or array of gummys. The return value is a numpy.ndarray.

-  static **covariance_matrix(gummys)**: Returns the variance-covariance
   matrix of a list or array of gummys. The return value is a
   numpy.ndarray.

-  **copy(formatting=True)**: Returns a copy of the gummy. The copy will
   be have a correlation coefficient of 1 with the original gummy. If
   the *formatting* parameter is ``True`` the display formatting
   information will be copied and if ``False`` the display formatting
   will be set to the default for a new gummy.

-  **ufrom(x)**: Gets the standard uncertainty contributed from
   particular gummys or utype_ if all other
   free variables are held fixed. The parameter *x* may be a string
   referencing a utype or a list containing gummys and (or) strings.
   This method returns a float. Example::

       >>>  a = gummy(1.2,0.2,utype='A')
       >>>  b = gummy(3.2,0.5,utype='A')
       >>>  c = gummy(0.9,0.2,utype='B')
       >>>  d = a + b + c
       >>>  d.ufrom('A')
       0.53851648071345048   

-  **doffrom(x)**: Gets the effective degrees of freedom contributed
   from particular gummys or utype_ if all
   other free variables are held fixed. The parameter *x* may be a
   string referencing a utype or a list containing gummys and (or)
   strings. This method returns a float. Example::

       >>>  a = gummy(1.2,0.2,dof=5,utype='A')
       >>>  b = gummy(3.2,0.5,dof=7,utype='A')
       >>>  c = gummy(0.9,0.2,utype='B')
       >>>  d = a + b + c
       >>>  d.doffrom('A')
       9.0932962619709627

.. _create:

-  classmethod **create(x, u=None, unit=None, dof=None, k=None, p=None, uunit=None,
   utype=None, name=None,correlation\_matrix=None,
   covariance\_matrix=None)**: Creates a list of
   correlated gummys.

   **create parameters** (only *x* is required, all others are optional):

   -  **x**: Either a list of floats corresponding to the x-value of
      each gummy or an instance of a MultivariateDistribution sub-class.

   -  **u**, **unit**, **dof**, **k**, **p**, **uunit**, **utype**, and
      **name**: Lists that correspond to the parameters in the gummy
      initializer (with the i-th value in each list passed to the
      initializer for the i-th gummy). With the exception of the "name"
      parameter, these may also be a single value with this same value
      is to passed to each initializer.

   -  **correlation_matrix**: An list or array to be used as the
      correlation matrix of the gummys. This is optional and must be set
      to the default value of None if the covariance\_matrix is
      specified.

   -  **covariance_matrix**: An list or array to be used as the
      variance-covariance matrix of the ummys. If the covariance matrix
      is specified the u parameter will be ignored This parameter is
      optional and must be set to the default value of None if the
      correlation\_matrix is specified. If both the correlation\_matrix
      and the covariance\_matrix are None (or omitted) then the gummys
      will be uncorrelated.

   **create returns**: a list of gummys

   Note: This package does not implement a multivariate Students's *t*
   distribution that has differing degrees of freedom for each
   component. So if if the elements of dof are finite and not all the
   same and either a correlation\_matrix or a covariance\_matrix is
   defined, the joint distribution for Monte-Carlo calculations (but not
   first-order calculations) will default to a multivariate normal
   distribution.

.. _budget:

-  **gummy.budget(xlist, uunit=None, k=None, p=None, sort=True,
   columns=None, column\_names=None, xnames=None, yname=None,
   show\_subtotals=True, show\_expanded\_u=None, description=None,
   description\_math\_mode=False, custom=None, custom\_heading=None,
   custom\_math\_mode=False, css=None, solidus=None, mulsep=None,
   show\_s=True, show\_d=False, show\_c=False, units\_on\_values=None,
   sim=False)**: Returns a Budget object that can be used to display an
   uncertainty budget table listing the the contributions of the gummys
   in *xlist* to the total uncertainty in the calling gummy ``self``.

   To display the table use the Budget.html() or Budget.latex() methods
   in a console or notebook that supports this type of output or the
   python ``print`` function to get a unicode table.

   The Budget.tohtml() and Budget.tolatex() methods can be used to get
   strings with the html or latex code.

   The Budget.df property can be used to retrieve a pandas DataFrame
   with the table. Also Budget.df\_str, Budget.df\_html and
   Budget.df\_latex return DataFrames with formatted strings as entries
   rather than numerical values.

   **metrolopy.gummy.budget parameters** (*xlist* is required, all others are 
   optional):

   -  **xlist**: (list of gummy) The independent variables.
      Warnings will be generated if the gummys in this list over
      determine ``self`` (that is if not all variables in this list can
      be treated as independent variables) or under determine ``self``
      (that is if some variables contributing to the uncertainty in
      ``self`` are missing).

   -  **uunit**: (str or :ref:`Unit<class-Unit>`, default is None) Unit to use
      to express the uncertainties. This useful if you wish to express
      all uncertainties as relative uncertainty unit (e.g. %).

   -  **k** and **p**: (float or None, default is None) *k*
      or *p* values for the expanded uncertainty of the total combined
      uncertainty; specify either *k* or *p* and not both; if neither
      are specified the the *k* and *p* values of ``self`` are used.

   -  **sort**: (bool, default is ``True``) Whether or not to sort
      the gummys in *xlist* by significance.

   -  **columns**: (list of str or None) Allows the user to select the
      columns (and ordering of the columns) for display. The available
      columns are:

      -  "component" or "name": the names of the gummy, displayed by
         default

      -  "description": description given in the description parameter
         list, displayed by default if the description parameter is not
         None

      -  "unit": the unit of the gummy, displayed by default

      -  "value": the x value of the gummy, displayed by default

      -  "u" or "uncertainty": The uncertainty of the gummy. This is the
         standard uncertainty except possible in the last row where an
         expanded uncertainty is displayed. This column is displayed by
         default.

      -  "dof": the degrees of freedom for the uncertainty, displayed by
         default if any uncertainty has finite degrees of freedom

      -  "type": the uncertainty type, displayed by default if any gummy
         has a utype_ defined

      -  "s" or "significance": the sensitivity coefficient (below)
         multiplied by the standard uncertainty, displayed by default

      -  "d", "derivative" or "partial": the partial derivative of
         ``self`` with respect to the gummy in that row

      -  "c" or "sensitivity coefficient": the absolute value of "d"

      -  "custom": value given in the custom parameter list, displayed
         by default if the custom parameter is not None

   The columns displayed can also be set with the Budget.columns
   property.

   -  **column\_names**: (dict or None) Names to display as
      column headers, if this is None then the default names are used.
      The dictionary should use as keys any of the column names listed
      above in the columns parameter description and as values the
      desired heading for this column. The column names can also be set
      with the Budget.column\_names property.

   -  **show\_subtotals**: (bool, default is ``False``) If any
      uncertainty types are defined, the combined standard uncertainty
      for each type is displayed in the table. This can also be changed
      by setting the Budget.show\_subtotals attribute.

   -  **show\_expanded\_u**: (bool or None, default is None)
      Whether or not to display the expanded uncertainty in the last
      row. If this is None, then the expaned uncertianty is displayed if
      ``self.k != 1``. This can also be changed by setting the
      Budget.show\_expanded\_u attribute.

   -  **show\_s**: (bool, default is ``True``) Whether or not to
      show the significance column. This is ignored if columns is not
      None. The default can be changed by setting the attribute class
      attribute Budget.show\_s.

   -  **show\_d**: (bool, default is ``True``) Whether or not to
      show the partial derivatives column. This is ignored if columns is
      not None. The default can be changed by setting the attribute
      class attribute Budget.show\_d.

   -  **show\_c**: (bool, default is ``False``) Whether or not to
      show the sensitivity coefficient column. This is ignored if
      columns is not None. The default can be changed by setting the
      attribute class attribute Budget.show\_c.

   -  **units\_on\_values**: (bool or None, default is
      None): If this is ``True``, units are shown in the value and u
      columns and if ``False`` the units are in a separate column. If
      None then the units are in a separate column unless ``self`` or
      any gummy in *xlist* has a uunit defined.

   -  **sim**: (bool, default is ``False``): If ``True``, the
      combined uncertainty and partial derivatives will be calculated
      using Monte-Carlo data.

   -  **css**: (str or None, defualt is None) css header to
      be used when displaying the table in HTML format. If this is None
      then Budget.default\_css will be used.

   -  **description**: (list of str or None, default is
      None) An optional column of descriptions to be printed in the
      table. This should be a description for ``self`` then for each
      *x*, and followed, optionally, by subtotal and expanded
      uncertainty descriptions.

   -  **description\_math\_mode**: (bool, default is ``False``) If
      this is ``False``, then when using a LaTeX format, the description
      is put in normal text mode rather than math mode.

   -  **custom**: (list of str or None, default is None)
      An optional column of additional information to be printed in the
      table. This should be a value for ``self`` then for each *x*, and
      followed, optionally, by subtotal and expanded uncertainty values.

   -  **custom\_heading**: (str or None, default is None) A
      heading for the custom column.

   -  **custom\_math\_mode**: (bool, default is ``False``) If this
      is ``False``, then when using a LaTeX format, the custom value is
      put in normal text mode rather than math mode.

   -  **solidus** and **mulsep**: Affects unit formatting, see the gummy
      attributes solidus_ and mulsep_

.. _conjugate:

-  **conjugate()**: returns a copy of ``self``

.. _angle:

-  **angle()**: returns ``gummy(numpy.pi)`` if ``self.x`` < 0 and
   ``gummy(0)`` otherwise

arithmetic operations and functions involving gummys
----------------------------------------------------

The standard Python arithmetic operations are allowed between gummys and
between gummys and floats or integers: addition, subtraction,
multiplication, division, floor division, exponentiation, modulus,
absolute value, and negation. These operations are allowed with
complex types, with the result a jummy_ rather than a
gummy instance. For addition and subtraction, the units must be
compatible (the of units of the two operands do not need to be the same,
but a conversion must exist between the units, see also the
c_ property). Exponents must be dimensionless (that is a
conversion from the exponent unit to the unit *one* must exist) and if
the exponent has an uncertainty, the base must be dimensionless.
Nonlinear units such as the decibel and the degree Celsius affect the
behavior of gummys under certain operations.

The gummy module installs a number of common mathematical
functions_ that can be applied directly to dimensionless
gummys, e.g::

    >>> import gummy as uc
    >>> g = uc.gummy(0.123,0.022) 
    >>> uc.sin(g)
    0.851(12)

For numpy version 1.13 or later, many numpy functions can be applied
directly to dimensionless gummys, e.g::

    >>> import numpy as np
    >>> import gummy as uc
    >>> g = uc.gummy(0.123,0.022) 
    >>> np.cos(g)
    -0.525(19)
        

The two class methods immediately below may also be used to apply an
arbitrary numerical function to one or more gummys.

gummy methods for applying numerical functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _apply:

-  classmethod **apply(function, derivative, arg1, arg2, ...)**:
   Applies a function to one or more dimensionless gummy objects
   propagating the uncertainty.

   **apply Parameters**:

   -  **function**: The the function to be applied. This must be a
      Python function that takes one more float arguments and return a
      float value or float array.

   -  **derivative**: The name of a second function that gives the
      derivatives with respect to the arguments of *function*.
      *derivative* should take an equal number of arguments as
      *function*. If *function* takes one argument *derivative* should
      return a float and if *function* takes more than one argument then
      derivative should return a tuple, list or array of floats that
      contains the derivatives with respect to each argument.

   -  **arg1, arg2, ...**: One or more arguments to which *function*
      will be applied. These arguments need not all be gummys objects;
      arguments such as floats will be taken to be constants with no
      uncertainty. They may also be numpy ndarrays in which case the
      usual numpy broadcasting rules apply. All gummy arguments must be
      dimensionless, that there must exist a conversion to the unit
      *one*.

   -  **return value**: If none of the arguments arg1, arg2, ... are
      gummy then the return value is *function* directly applied to the
      arguments. Otherwise the return value is a gummy.

   Examples::

       >>> import numpy as np
       >>> x = gummy(0.678,u=0.077)
       >>> gummy.apply(np.sin,np.cos,x)
       0.627(60)

       >>> x = gummy(1.22,u=0.44)
       >>> y = gummy(3.44,u=0.67)
       >>> def dhypot(x,y):
       ...     return (x1/sqrt(x1**2 + x2**2),x2/np.sqrt(x1**2 + x2**2))
       >>> gummy.apply(np.hypot,dhypot,x,y)
       3.65(65)

.. _napply:

-  classmethod **napply(function, derivative, arg1, arg2, ...)**:
   Applies a function to one or more dimensionless gummy objects
   propagating the uncertainty. This method is similar to apply except
   that the derivatives are computed numerically so a derivative
   function does not need to be supplied.

  **napply parameters**:

   -  **function**: The the function to be applied. This must be a
      Python function that takes one more float arguments and return a
      float value or float array.

   -  **derivative**: The name of a second function that gives the
      derivatives with respect to the arguments of *function*.
      *derivative* should take an equal number of arguments as
      *function*. If *function* takes one argument *derivative* should
      return a float and if *function* takes more than one argument then
      derivative should return a tuple, list or array of floats that
      contains the derivatives with respect to each argument.

   -  **arg1, arg2, ...**: One or more arguments to which *function*
      will be applied. These arguments need not all be gummys objects;
      arguments such as floats will be taken to be constants with no
      uncertainty. They may also be numpy ndarrays in which case the
      usual numpy broadcasting rules apply. All gummy arguments must be
      dimensionless, that there must exist a conversion to the unit
      *one*.

   -  **return value**: If none of the arguments arg1, arg2, ... are
      gummy then the return value is *function* directly applied to the
      arguments. Otherwise the return value is a gummy.

   Examples::

       >>> import numpy as np
       >>> x = gummy(0.678,u=0.077)
       >>> gummy.napply(np.sin,x)
       0.627(60)

       >>> x = gummy(1.22,u=0.44)
       >>> y = gummy(3.44,u=0.67)
       >>> gummy.napply(np.hypot,x,y)
       3.65(65)      

gummy properties and methods related to units and unit conversion
-----------------------------------------------------------------

Units are represented by instances of the :ref:`Unit<class-Unit>` class or
sub-classes, however the user rarely needs to interact directly with
these objects as strings can be used in place of :ref:`Unit<class-Unit>`
objects in all properties and methods dealing with units. It is,
however, straight forward for to create custom units.

gummy properties related to units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _property-unit:

-  **unit**: Gets or sets the units for the values of the
   x_ and :ref:`u<lower-u>` properties, also sets the
   units for :ref:`U<upper-u>` if *uunit* is None. When setting this
   property either a :ref:`Unit<class-Unit>` object or a string referencing a
   :ref:`Unit<class-Unit>` object may be used. A ``NoUnitConversionFoundError``
   exception will be generated if no conversion exists between the
   original unit and the new unit. Example::

       >>> g = gummy(1,unit='J')
       >>> g.unit = 'erg'
       
   Spaces or asterisks ('*') may be used to represent unit
   multiplication, a forward slash ('/') division and double asterisks
   ('**') powers, e.g 'm/s' or 'm s**-1'.  If a unit
   name contains a space, '\*' or '/' character then the name must be
   enclosed in square brackets, e.g: '[light year]'.

.. _uunit:

-  **uunit**: gets or sets the units for :ref:`U<upper-u>`. Setting
   *uunit* to None puts :ref:`U<upper-u>` in the same units as
   x_. If *uunit* is a dimensionless unit (e.g. *one*,
   '%', 'ppm' or 'um/m') then :ref:`U<upper-u>` is a
   relative uncertainty. When setting this property either a
   :ref:`Unit<class-Unit>` object or a string referencing a Unit object may be
   used.

.. _unit_is_rel:

-  **uunit_is_rel**: Returns ``True`` if the :ref:`U<upper-u>`
   property will return a relative uncertainty and ``False`` otherwise.

.. _c:

-  **c**: This read-only property is used as a conversion flag during
   calculations. When an arithmetic operation is carried out between two
   gummys with different units, a unit conversion on one of the input
   quantities may be required to complete the calculation. Attach this
   flag to the unit that you prefer be converted if you do not which
   gummy to make the choice. For example::

       >>> a = gummy(1,u=0.01,unit='cm')
       >>> b = gummy(2,u=0.2,unit='mm')
       >>> a + b
       1.200(22) cm
       >>> a.c + b
       12.00(22) mm
       >>> a + b.c
       1.200(22) cm

gummy methods related to units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **convert(unit,uunit=None)**: Returns a copy of the original gummy
   with converted units. *unit* and *uunit* are either a strings or Unit
   instances for the new units for the *x* and *U* properties of the new
   gummy respectively.

-  **graft(unit,uunit=None)**: Returns a copy of the original gummy,
   keeping the same numerical value of the *x* and *U* properties of the
   original gummy, but with new units *unit* and *uunit* respectively.

-  **reduce\_unit()**: Cancels factors in a gummy's unit when possible.
   This modifies the calling gummy and returns None. For example
   example::

       >>> g = gummy(5,unit='mm/m')
       >>> g
       5 mm/m
       >>> g.reduce_unit()
       >>> g
       0.005

.. _search_units:

unit search function
~~~~~~~~~~~~~~~~~~~~

The following function is not part of the gummy class but is useful when
dealing with units.

-  metrolopy.\ **search\_units(search=None, fmt=None, show\_all=False, units=
   None, prnt=True)**: Prints a list of all units that match the search
   terms. If this function is called with no arguments, then a list of
   all loaded units is printed.

   **search\_units parameters**:

   -  **search**: (str) A space separated list of search terms to case
      insensitively match. If this is omitted or set equal to None a
      list of all loaded units will be printed. The default is None.

   -  **fmt**:
      ({'html', 'latex', 'unicode', 'ascii', None}, optional)
      The output format. If None, then the gummy.printer value is used.
      If latex output is selected, Markdown is actually used with the
      unit symbols and conversion displayed using inline LaTeX.

   -  **show\_all**: (bool,optional) If ``True`` units are shown
      with each prefix listed on a separate line (e.g. the millisecond
      and the microsecond are listed in addition to the second) and
      interval units are shown. If ``False`` only the base unit is
      shown. The default is ``False``.

   -  **units**: (list of str or :ref:`Unit<class-Unit>`,optional) A list of
      units to print. If this parameter is specified the values of the
      search and show\_all parameters are ignored.

   -  **prnt**: (bool,optional) If this is ``True``, the results are
      printed. If it is ``False`` the results are returned as a string.
      The default is True.

gummy properties and methods related to Monte-Carlo simulation
--------------------------------------------------------------

Gummy allows uncertainty propagation using Monte-Carlo simulation in
addition to first order error propagation. Before using many of the
properties and methods in this section, Monte-Carlo data must be
generated by calling the sim_ method (to generate
Monte-Carlo data for one gummy) or the simulate_
static method (to generate Monte-Carlo data for one or more gummys).
Note that these methods will erase all previous Monte-Carlo data from
all gummys before generating new data. So if you want data available for
multiple gummys use the simulate_ static method
rather than the sim_ method. A ``NoSimulatedDataError``
exception will be raised if no simulated data is available when a
property is accessed or a method called that needs Monte-Carlo data.

gummy properties related to Monte-Carlo simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **xsim**: (read-only) Gets the mean of the simulated data.

-  **usim**: (read-only) Gets the standard deviation of the simulated
   data.

-  **cisim**: (read-only) Gets a tuple giving the lower followed by the
   upper bound of the confidence interval calculated from the simulated
   data. The confidence level for the interval is equal to the
   p_ property of the gummy. See the *cimethod* property
   for details on how the interval is calculated.

-  **cimethod**: ({'shortest', 'symmetric'}) 
   Gets or sets the method for calculating the confidence
   interval from the Monte-Carlo data. It can be set either to the
   string 'shortest' or the string 'symmetric' (the default is
   'shortest'). The 'shortest' confidence interval will be taken
   to be the shortest interval that includes the desired fraction of the
   probability distribution. If the confidence interval is
   'symmetric', then it will be set so that, for *n* Monte-Carlo
   samples and a coverage probability of *p*, *n*\ \*(1-\ *p*)/2 samples
   lie below the lower limit of the confidence interval and the same
   number of samples lie above the upper limit of the confidence
   interval. This property can be set at the class or instance level.

-  **Usim**: (read-only) Gets the tuple: (*xsim* - *cisim*\ [0], *xsim*
   + *cisim*\ [1])

-  **ksim**: (read-only) Gets 0.5\*(\ *Usim*\ [0] + *Usim*\ [1])/*usim*.

-  **simdata**: (read-only) Gets a ``numpy.ndarray`` containing the
   simulated data.

-  **simsorted**: (read-only) Gets a ``numpy.ndarray`` containing the
   simulated data sorted from smallest to largest.

-  **distribution**: (read-only) Gets the
   Distribution_ sub-class instance representing the
   gummy. If the gummy was created as a result of mathematical
   operations involving other gummys, then the distribution will be an
   instance of the ``Convolution`` sub-class of
   Distribution_. For other gummys the distribution
   can be specified using the *x* parameter in the
   initializer. If *x* is not specified as a
   Distribution_, the distribution will be taken as
   either the NormalDist_ or TDist_
   sub-classes of Distribution_.

gummy methods related to Monte-Carlo simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sim:

-  **sim(n=100000,ufrom=None)**: Generates Monte-Carlo data for the
   calling gummy. Calling this method erases previously generated
   Monte-Carlo data for all gummys, so use the gummy.simulate
   staticmethod if you need Monte-Carlo data for several gummys
   simultaneously. *n* is the number of samples to generate. The default
   value for *n* is 100000. The *ufrom* keyword can be used to
   separately analyze various contributions to the uncertainty. If
   *ufrom* is specified, only the gummys referenced by *ufrom* will be
   allowed to vary and all other gummys will be held fixed. *ufrom* can
   be a gummy, a string referencing a *utype*, or a list containing
   gummys and (or) string references to utypes.

.. _simulate:

-  static **simulate(gummys,n=100000,ufrom=None)**: Generates
   Monte-Carlo data for one or more gummys references by the parameter
   *gummys*. *gummys* can be can be a gummy, a string referencing a
   *utype*, or a list containing gummys and (or) string references to
   utypes. Calling this method erases previously generated Monte-Carlo
   data for all gummys. The *n* and *ufrom* parameters are the same as
   for the .sim method.

-  static **set_seed(seed)**: Sets the seed for the
   ``numpy.random.RandomState`` object shared by all ``Distribution``
   instances.

-  **hist(title=None, xlabel=None, p=None, show_p=True, title_style=None,
   mean_marker=True, mean_marker_options={}, ci_marker=True,
   ci_marker_options={}, hold=False, \*\*plot_options)**:
   Plots a histogram of the Monte-Carlo data for the gummy. Before
   calling this method either the .sim or .simulate method must be
   called to generate theMonte-Carlo data. 
   
   **hist parameters** (all are optional):

   -  **title**: (str or None) A title for the plot. If this is omitted or
      set to None then a title will be generated using the gummy name
      (if it has one) and the mean value and confidence interval. The
      title will also give the standard deviation of the date. The
      formatting of the auto-generated title depends on the value of the
      title_style parameter.

   -  **xlabel**: (str or None) A label for the horizontal axis of the plot.
      If this is omitted or set to None then a label will be generated
      using the name and unit of the gummy. If xlabel is None and the
      gummy has no name and a unit of one, then the horizontal axis will
      not be labeled.

   -  **p**: (float between 0 and 1 or None) The probability for the
      confidence interval (as printed in the title and indicated by the
      ci_markers). If this is none then the value of the gummy.p
      property is used. The default is None.

   -  **show_p**: (bool) Whether or not to show the level of confidence in
      the title if the title is auto-generated.

   -  **title_style**: (str in
      {'pmsim','pmsimi','cisim','mcisym','xsim','xfsim',
      'usim','ufsim'}) The style for displaying the value in the title.
      See the gummy.style property for details. It this is None or
      omitted then the value of the gummy.style property is used.

   -  **mean_marker**: (bool) Whether or not to display a vertical line at
      the mean value (as given by gummy.xsim). The default is True.

   -  **mean_marker_options**: (dict) A dictionary containing keywords to
      be passed to the pyplot.axvline method which draws the mean marker.
      For example setting this to {'color'='r','linewidth'=4} makes the
      mean marker red and with thickness of four points.

   -  **ci_marker**: (bool) Whether or not to display vertical lines at the
      upper and lower limits of the confidence interval. The default is
      True.

   -  **ci_marker_options**: (dict) A dictionary containing keywords to be
      passed to the pyplot.axvline method which draws the confidence
      interval markers.

   -  **hold**: (bool) If this is False pyplot.show() is called before this
      method exits. If it is True pyplot.show() is not called. The
      default is False.

   -  **plot_options**: These are optional keyword arguments that are
      passed to the pyplot.hist method. For example bins=50 overrides the
      default number of bins (100). For other options see the pyplot.hist
      documentation.

.. _covplot:

-  static **covplot(x, y, title=None, xlabel=None, ylabel=None, mean_marker=False,
   mean_marker_options={}, hold=False, \*\*plot_options)**: 
   Creates scatter plot showing the covariance between
   two gummys.

   **covplot paramters** (all but *x* and *y* are optional):

   -  **x**: (gummy) The gummy to plot on the horizontal axis.

   -  **y**: (gummy) The gummy to plot on the vertical axis.

   -  **title**: (str or None) A title for the plot. If this is
      omitted or set to None then the correlation will be displayed as
      the title.

   -  **xlabel**: (str or None) A label for the horizontal axis.
      If this is omitted or None then that axis will be labeled either
      "x" or with the *x* gummy's unit.

   -  **ylabel**: (str or None) A label for the vertical axis. If
      this is omitted or None then that axis will be labeled either "y"
      or with the *y* gummy's unit.

   -  **mean_marker**: (bool) Whether or not to display line markers at
      the mean values of *x* and *y*. The default is False.

   -  **mean_marker_options**: (dict) A dictionary of options to be
      passed to the pyplot.axvline and pyplot.axhline methods that draw
      the mean\_marker.

   -  **hold**: (bool) If this is False pyplot.show() is called before
      this method exits. If it is True pyplot.show() is not called. The
      default is False.

   -  **plot_options**: These are optional keyword arguments that are
      passed to the pyplot.plot method. For example ms=0.1 decreases the
      size of the dots in the plot.

gummy properties and methods related to display and formatting
--------------------------------------------------------------

gummy formatting properties and attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of these properties and attributes can be set at either the class or
instance level.

-  **style**: Get or set the default display style for new gummys. This
   is a string with one of the following values:

   -  ``'pm'`` or ``'+-'`` may give, e.g. in ascii format
      ``(1.00 +/- 0.12)e-12 cm``
   -  ``'pmi'`` or ``'+-i'`` gives e.g. ``1.00e-12 cm +/- 1.2e-13 cm``
   -  ``'concise'`` or ``'()'`` gives e.g. ``1.00(12)e-12 cm``
   -  ``'ueq'`` or ``'u equals'`` gives e.g.
      ``1.00e-12 cm with u = 1.2e-13 cm``
   -  ``'x'`` or ``'x only'`` gives e.g. ``1.00e-12 cm``
   -  ``'xf'`` gives e.g. ``1.00e-12``
   -  ``'u'`` or ``'u only'`` gives ``1.2e-13 cm``
   -  ``'uf'`` gives ``1.2e-13``

   The following styles display a mean and confidence interval based on
   data from a Monte-Carlo simulation (if no simulated data is available
   the string ``'no simulated data'`` is returned):

   -  ``'pmsim'`` gives e.g. ``(1.01 + 0.11 - 0.13)e-12 cm``
   -  ``'pmsimi'`` gives e.g. ``1.01e-12 cm + 1.1e-13 cm - 1.3e-11 cm``
   -  ``'mcisim'`` gives e.g.
      ``mean = 1.01e-12 cm, confidence interval = [8.8e-13 cm, 1.13e-12 cm]``
   -  ``'cisim'`` gives e.g ``[8.8e-13 cm, 1.13e-12 cm]``
   -  ``'usim'`` gives the standard deviation e.g. ``1.2e-13 cm``
   -  ``'ufsim'`` gives the standard deviation e.g. ``1.2e-13``

   Note if uunit has been defined for the gummy instance, then concise
   style and pm style are not valid and the display will default to pmi
   style.

   The style can be set either at the class or the instance level.

-  **nsig**: (int) This number of significant figures to display for
   the uncertainty. The default is 2. This can be set at the class or
   instance level.  If nsig is set to -1, then the uncertainty will not
   be explicitly printed but the x value will be rounded so that the 
   expanded uncertainty U is between 1 and 10 counts in the last
   digit.

-  **show\_k**: (bool or None) Determines whether the coverage
   factor (*k*) is printed along with the value and uncertianty. If this
   is set to None (the default value), then *k* will be displayed if
   the k_ property or *k* parameter
   has been set to a value other than 1.

-  **show\_p**: (bool or None) Determines whether the level of
   confidence for the uncertainty (*p*) is printed along with the value
   and uncertainty. If this is set to None (the default value), then
   *p* will be displayed if the p_ property or
   *p* parameter has been set.

-  **show\_dof**: (bool or None) Determines whether the
   effective number of degrees of freedom is printed along with the
   value and uncertainty. If this is set to None (the default
   value), then *dof* will be displayed if the value of the
   dof_ property is less than 10.

-  **show\_name**: (bool) If this is ``True`` and the
   name_ of the gummy is not None then the gummy will
   be displayed as ``[*name*] = ...``

-  **sci\_notation**: (bool or None) Determines whether the
   scientific notation will be used. If this is set to None (the
   default value) then scientific notation will be used if the value of
   the x_ property has an exponent greater than the
   *sci\_notation\_high* attribute or an exponent less than
   *sci\_notation\_low* attribute.

-  **sci\_notation\_high**: (int) See the *sci\_notation* attribute,
   default is 7.

-  **sci\_notation\_low**: (int) See the *sci\_notation* attribute,
   default is -3.

.. _solidus:

-  **solidus**: (bool) Determines whether unit are displayed with
   negative exponents (if *solidus* is ``False``) or if a solidus
   (back-slash) is used to separate units in the numerator and
   denominator of a composite unit (when *solidus* is ``True``, this is
   the default).

.. _mulsep:

-  **mulsep**: (bool) If ``True`` an asterix, '*', or dot,
   '·', is used to separate units in a composite unit, and a space is
   used if *mulsep* is ``False``. The default is ``False``.

.. _slashaxis:

-  **slashaxis**: (bool) Affects plots axes labeling where units are
   used. If ``True`` e.g. an axis may be labeled 'Time / s' while if
   ``False`` it would be labeled 'Time (s)'.

.. _printer:

-  **printer**: Get or set the prefered display printer. This is a
   string with one of the following values: 'any', 'latex',
   'html', 'unicode', 'ascii', or 'any_but_latex'.
   'any' will usually pick html or latex output when running in an
   IPython console or Jupyter notebook and unicode otherwise.
   "any\_but\_latex" will usually pick html when running in an IPython
   console or Jupyter notebook and unicode otherwise. 'latex' and 'html'
   are only available when running under IPython. If these printers are
   not available the display will default to 'unicode'.

gummy display methods
~~~~~~~~~~~~~~~~~~~~~

-  **latex(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Prints the gummy using LaTeX if this method is
   executed from a latex capable Ipython console or from a Jupyter or
   Ipython notebook. All parameters are optional. Any parameters that
   are not None override the corresponding attributes of the calling
   gummy.

-  **html(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Prints the gummy using HTML if this method is
   executed from an Ipython console or from a Jupyter or Ipython
   notebook. All parameters are optional. Any parameters that are not
   None override the corresponding attributes of the calling gummy.

-  **unicode(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Prints the gummy using unicode characters. Any
   parameters that are not None override the corresponding attributes of
   the calling gummy. If ``g`` is a gummy, then ``g.unicode()`` is equivalent to
   ``print(g.tounicode(...))``.

-  **ascii(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Prints the gummy using unicode characters. Any
   parameters that are not None override the corresponding attributes of
   the calling gummy. If ``g`` is a gummy, then ``g.ascii()`` equivalent to
   ``print(g.toascii(...))``.

-  **tostring(fmt=None, style=None, k=None, p=None, show\_k=None,
   show\_p=None, show\_dof=None, show\_name=None, nsig=None,
   solidus=None, mulsep=None)**: Returns a string displaying the value
   of the gummy in the desired format. The fmt parameter is a string
   with the value in {"unicode","latex","html","ascii"} or None. fmt
   will default to 'ascii' if self.printer is'ascii' or 'unicode'
   otherwise. Any other parameters that are not None override the
   corresponding attributes of self.

-  **tohtml(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Returns a string with the value as an html fragment.
   All parameters are optional. Any parameters that are not None
   override the corresponding attributes of the calling gummy. This is
   equivalent to the gummy.tostring method with the *fmt* parameter set
   to 'html'.

-  **tolatex(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Returns a string with the value as an LaTeX fragment.
   It is assumed that LaTeX is in math mode. All parameters are
   optional. Any parameters that are not None override the corresponding
   attributes of the calling gummy. This is equivalent to the
   gummy.tostring method with the *fmt* parameter set to 'latex'.

-  **toascii(style=None, k=None, p=None, show\_k=None, show\_p=None,
   show\_dof=None, show\_name=None, nsig=None, solidus=None,
   mulsep=None)**: Returns a string with the value formatted so that
   only ASCII characters are used. All parameters are optional. Any
   parameters that are not None override the corresponding attributes of
   the calling gummy. This is equivalent to the gummy.tostring method
   with the *fmt* parameter set 'ascii'.

.. _Distribution:

class Distribution and sub-classes
==================================

The ``Distribution`` class is the abstract base class for objects which
represent the the probability distributions that the gummy Monte-Carlo
samples are drawn from. Instances of these ``Distributions`` can used as
the *x* parameter when creating gummys. The
distributions below are built into the gummy package and custom
distributions can also be defined by the user.

.. _arcsindist:

- class metrolopy.\ **ArcSinDist(center=None, half_width=None,
  lower_limit=None, upper_limit=None)**: Arcsin distribution, specify 
  either *center* and *half_width* or *lower_limit* and *upper_limit*. 

.. _binomialdist:

- class metrolopy. **BinomialDist(n, p)**: Binomial distribution with number of
  trials *n* and success probability *p*.

.. _convolution:

- class metrolopy.\ **Convolution(func, d1, d2, ...)**: Normally this Distribution
  is not created directly, but is the result of mathematical operations
  involving gummys. This sub-class represents distributions resulting from
  applying *func* to *d1*, *d2*, ... The
  function *func* takes an the same number of scalar arguments as there
  are *d1*, *d2*, ... parameters and returns a scalar. *d1*, *d2*, ... can
  be either instances of ``Distribution`` subclasses or scalar values.
    
.. _curvlineartrapdist:

- class metrolopy.\ **CurvlinearTrapDist(center=None, half_width=None, 
  limit_half_range=None, lower_limit=None, upper_limit=None)**: 
  Curvlinear trapezoidal
  distribution, *limit_half_range* is required. Also either *center* and
  *half_width* or *lower_limit* and *upper_limit* are required. This is
  intended to represent a variable that follows a uniform distribution but
  where the upper and lower limits are not exactly known and may vary by
  up to the *limit_half_range* from the given lower and upper limit
  values. 

.. _exponentialdist:

- class metrolopy.\ **ExponentialDist(scale=None, rate=None)**:
  Exponential distribution with probability density function::
  
      f(x;rate) = rate*exp(-rate*x). 
      
  Specify either *scale* or *rate* (*scale* = 1/*rate*), but not both. 

.. _gammadist:

- class metrolopy.\ **GammaDist(shape, scale)**: Gamma distribution with the
  *shape* and *scale* parameters. 
    
.. laplacedist:

- class metrolopy.\ **LaplaceDist(x, scale)**: Laplace distribution with location
  parameter *x* and *scale* parameter. 

.. _lognormaldiat:

- class metrolopy.\ **LogNormalDist(mu=None,sigma=None)**: Log-normal distribution
  where the logrithm of the random variable has mean *mu* and standard 
  deviation *sigma*. 

.. _multinormaldist:

- class metrolopy.\ **MultiNormalDist(mean, cov)**:
  Multivariate normal distribution. The parameter *mean* is a list of mean
  values for each dimension and *cov* is the variance-covariance matrix.

.. _mutitdist:

- class metrolopy.\ **MultiTDist(mean, cov, dof)**: Multivariate shifted and
  scaled Students's *t* distribution. The parameter *mean* is a list of mean
  values for each dimension and *cov* is the variance-covariance matrix.
  The parameter *dof* is the number of degrees of freedom and must be
  scalar; all dimensions must have the same number of degrees of freedom.
  
.. _multivariatedist:

- class metrolopy.\ **MultvariateDistribution(nd)**: Abstract base class for
  mulit-variate distributions. *nd* is the number of dimensions. 

.. _NormalDist:

- class metrolopy.\ **NormalDist(x, s)**: Normal distribution with mean *x* and
  standard deviation *s*. 

.. _poissondist:

- class metrolopy.\ **PoissonDist(lam)**: Poisson
  distribution with rate parameter *lam*. 

.. _TDist:

- class metrolopy.\ **TDist(x, s, dof)**: Shifted and scaled Students's *t*
  distribution with degrees of freedom *dof*, mean *x*, and scale factor *s*. 

.. _trapezoidaldist:

- class metrolopy.\ **TrapezoidalDist(lower_limit, upper_limit, top_to_base_ratio)**:
  Trapezoidal distribution 
    
.. _triangulardist:

- class metrolopy.\ **TriangularDist(mode, half_width=None, left_width=None,
  right_width=None, lower_limit=None, upper_limit=None)**:
  Triangular distribution. For a symmetric distribution specify
  *half_width*, otherwise specify two, and only two, of the parameters
  *left_width*, *right_width*, *lower_limit*, *upper_limit*.

.. _UniformDist:

- class metrolopy.\ **UniformDist(center=None, half_width=None,
  lower_limit=None, upper_limit=None)**:
  A uniform distribution. Specify two, and only two, of the parameters
  *center*, *half_width*, *lower_limit* and *upper_limit* . 
  
.. _weibulldist:

- class metrolopy.\ **WeibullDist(shape, scale)**: Weibull distribution with
  *shape* and *scale* parameters.

custom distributions
--------------------

Custom distributions can be implemented by creating a class that
inherits from the ``Distribution`` class and implements the following
methods:

-  **random(self, n=None)**: Return a numpy array of *n* values
   drawn from the distribution. If *n* is None then a single scalar
   value should be returned. Preferably use, as a random number
   generator, the numpy ``RandomState`` object accessed with the
   ``Distribution.random_state`` static method.

-  **x(self)**: A scalar "center" of the distribution. This is used to
   get the x_ value of a gummy defined with the
   distribution.

-  **u(self)**: A scalar "standard uncertainty" of the distribution
   (usually the standard deviation). This is used to get the
   :ref:`u<lower-u>` value of a gummy defined with the distribution.

for example::

    class ChiSquaredDist(Distribution):
        def __init__(self,dof):
            self.dof = dof

        def random(self,n=None):
            return Distribution.random_state().chisquare(self.dof,n)

        def x(self):
            return self.dof

        def u(self):
            return 2*self.dof           

custom multi-variate distributions
----------------------------------

To create a multi-variate distribution, inherit from the
``MultivariateDistribution`` class and define the following methods:

-  **simulate(self, n)**: Return a numpy array of *n* samples drawn
   from the distribution. Preferably use, as a random number generator,
   the numpy ``RandomState`` object accessed with the
   ``Distribution.random_state()`` static method.  For a distribution with
   number of dimensions *nd*, the shape of the returned array must be
   (*nd*, *n*).

-  **x(self)**: a list or array with the "center" of the distribution
   for each dimension (usually the mean of the distribution). This is
   used to get the x_ values of gummys defined with the
   distribution.

-  **u(self)**: A list or array with "standard uncertainty" of the
   distribution (usually the standard deviation) for each dimension.
   This is used to get the :ref:`u<lower-u>` values of gummys defined
   with the distribution.

and the following read-only property

-  **cov**: (read-only property) Returns the variance-covariance matrix

The ``__init__`` function must also call the
``MultivariateDistribution`` ``__init__`` with the number of
dimensions nd of the distribution e.g. ``super().__init__(nd)``.

for example::

    class DirechletDist(MultivariateDistribution):
        def __init__(self,alpha):
            self.alpha = alpha

            a0 = sum(alpha)
            b = (a0**2*(a0 + 1))
            self._x = [a/a0 for a in alpha]
            self._u = [a*(a0 - a)/b for a in alpha]
            self._cov = [[ai*(a0 - ai)/b if i == j else -ai*aj/b 
                          for i,ai in enumerate(alpha)] 
                          for j,aj in enumerate(alpha)]

            super().__init__(len(alpha))

        def _simulate(self,n):
            self.simdata = Distribution.random_state().dirichlet(self.alpha,n).T

        def x(self):
            return self._x

        def u(self):
            return self._u

        @property
        def cov(self):
            return self._cov

.. _functions:

built-in mathematical functions
===============================

The following functions are installed as part of the gummy package can
take gummy, jummy, or float arguments. Arguments must be
dimensionless for transcendental functions.

metrolopy.\ **absolute(x)**: equivalent to abs(x)

metrolopy.\ **add(x1,x2)**: equivalent to z1 + x2

metrolopy.\ **angle(x)**: returns the complex argument of x

metrolopy.\ **arccos(x)**: inverse cosine of x

metrolopy.\ **arccosh(x)**: inverse hyperbolic cosine of x

metrolopy.\ **arcsin(x)**: inverse sine of x

metrolopy.\ **arcsinh(x)**: inverse hyperbolic sine of x

metrolopy.\ **arctan(x)**: inverse tangent of x

metrolopy.\ **arctanh(x)**: inverse hyperbolic tangent of x

metrolopy.\ **arctan2(x,y)**: arctan of x/y and giving the correct quadrant

metrolopy.\ **araound(x,n=0)**: x rounded to n digits

metrolopy.\ **cbrt(x)**: cube root of x

metrolopy.\ **ceil(x)**: ceiling of x

metrolopy.\ **conj(x)**: complex conjugate of x

metrolopy.\ **cos(x)**: cosine of x

metrolopy.\ **cosh(x)**: hyperbolic cosine of x

metrolopy.\ **cross(\*args,\*\*kwds)**: alias for numpy.cross

metrolopy.\ **cumpord(\*args,\*\*kwds)**: alias for numpy.cumprod

metrolopy.\ **cumsum(\*args,\*\*kwds)**: alias for numpy.cumsum

metrolopy.\ **diff(\*args,\*\*kwds)**: alias for numpy.diff

metrolopy.\ **divide(x1,x2)**: equivalent to x1 / x2

metrolopy.\ **divmod(x1,x2)**: returns (x1 // x2, x1 % x2)

metrolopy.\ **ediff1d(\*args,\*\*kwds)**: alias for numpy.ediff1d

metrolopy.\ **exp(x)**: natural exponential function

metrolopy.\ **expm1(x)**: exp(x) - 1

metrolopy.\ **exp2(x)**: exponential function with base 2

metrolopy.\ **fabs(x)**: equivalent to the built-in python function abs

metrolopy.\ **fix(x)**: returns x rounded towards zero

metrolopy.\ **floor(x)**: floor of x

metrolopy.\ **floor_divide(x1,x2)**: equivalent to x1 // x2

metrolopy.\ **gradient(\*args,\*\*kwds)**: alias for numpy.gradient

metrolopy.\ **heaviside(x,h0)**: Heavyside function: 0 for x < 0, h0 at x
== 0, and 1 for x > 0

metrolopy.\ **imag(x)**: returns the imaginary part of x

metrolopy.\ **log(x)**: natural logarithm of x

metrolopy.\ **logaddexp(x,y)**: loge(ex + ey)

metrolopy.\ **logaddexp2(x,y)**: log2(2x + 2y)

metrolopy.\ **log1p(x)**: natural logarithm of 1 + x

metrolopy.\ **log2(x)**: logarithm base 2 of x

metrolopy.\ **log10(x)**: logarithm base 10 of x

metrolopy.\ **mod(x)**: x1 % x2

metrolopy.\ **modf(x1,x2)**: returns (x1 % 1, x1 // 1)

metrolopy.\ **multiply(x1,x2)**: equivalent to x1 \* x2

metrolopy.\ **negative(x)**: equivalent to -x

metrolopy.\ **power(x1,x2)**: equivalent to x1\*\*x2

metrolopy.\ **prod(\*args,\*\*kwds)**: alias for numpy.prod

metrolopy.\ **real(x)**: returns the real part of x

metrolopy.\ **reciprocal(x)**: equivalent to 1/x

metrolopy.\ **remainder(x)**: x1 % x2

metrolopy.\ **rint(x)**: x rounded to the nearest integer value

metrolopy.\ **sin(x)**: sine of x

metrolopy.\ **sign(x)**: sign function, -1 for x < 0, 0 for x == 0, and 1
for x > 0

metrolopy.\ **sinh(x)**: hyperbolic sine of x

metrolopy.\ **square(x)**: square of x

metrolopy.\ **sqrt(x)**: the square root of x

metrolopy.\ **subtract(x1,x2)**: equivalent to x1 - x2

metrolopy.\ **sum(\*args,\*\*kwds)**: alias for numpy.sum

metrolopy.\ **tan(x)**: tangent of x

metrolopy.\ **tanh(x)**: hyperbolic tangent of x

metrolopy.\ **true\_divide(x1,x2)**: equivalent to x1 / x2

metrolopy.\ **trunc(x)**: x rounded towards zero

.. _class-Unit:

class Unit
==========

The gummy class uses ``Unit`` instances to represent physical units. A
number of units are loaded with the gummy package. See the
search_units_ function to get a list of all
available units. Custom units can also be defined by creating instances
of the ``Unit`` class or a sub-class. Though you
can assign the instance to a variable, this is not necessary since units
can be accessed using string names. E.g. we can define:

::

        >>> Unit('weird meter','wm',conversion=Conversion('m',0.9144),add_symbol=True)
        >>> gummy(3.3,unit='wm')
        3.3 wm
        

class metrolopy.\ **Unit(name, symbol, conversion=None,
short_name=None, add_symbol=False,
html_symbol=None, latex_symbol=None,
ascii_symbol=None, description=None, order = -1)**:
Creating an instance of this class creates a representation of a
physical unit and adds it to the unit library. Units already in the unit
library or derived units made up of other units in the unit library can
be accessed by passing a text string with the unit name or symbol to the
static method :ref:`unit<unit-unit>` . (In most cases you
do not need to call the :ref:`unit<unit-unit>` method directly; a gummy object
will automatically call this method when a gummy property or method is
passed a string that references a unit).

Unit parameters
---------------

-  **name**: (str) The name of the unit. The name can be used access
   the unit with the :ref:`unit<unit-unit>` method, but
   note that if you define a Unit with an identical name to a previously
   defined unit then the older name will be shadowed.

-  **symbol**: (str) A unicode symbol used when displaying the unit.
   If the *add_symbol* parameter is set to ``True``, then this symbol
   can also be used to access the unit with the
   :ref:`unit<unit-unit>` method. A gummy is normally
   printed with a space between the value and the unit, however this
   space is removed if the symbol string starts with a tab character.

-  **conversion** (``Conversion`` instance, optional, default is
   None) A ``Conversion`` instance representing the conversion to
   another unit. The conversion takes as its first argument the other
   unit and as the second argument the conversion factor (float or gummy)
   to the other unit, e.g.::

       >>> Unit('inch','in')
       >>> Unit('foot','ft',conversion=Conversion('in',12))
       >>> Unit('yard','yd',conversion=Conversion('ft',3))

   Circular conversion chains must be avoided. This will generate a
   ``CircularUnitConversionError`` exception::

       >>> Unit('inch','in',conversion=Conversion('yd',1/36))
       >>> Unit('foot','ft',conversion=Conversion('in',12))
       >>> Unit('yard','yd',conversion=Conversion('ft',3))

   The exception will not be raised until a unit conversion is attempted
   using one of these units and not immediately after the units are
   defined. An equivalent and allowable way of defining the first set of
   units above is::

       >>> Unit('inch','in')
       >>> Unit('foot','ft',conversion=Conversion('in',12))
       >>> Unit('yard','yd',conversion=Conversion('in',36))

-  **short_name**: (str, optional, default is None) An
   alternate short name for the unit that can be used with the
   :ref:`unit<unit-unit>` method.

-  **add_symbol**: (bool, optional, default is ``False``) If this
   is ``True``, the symbol, in addition to the *name* and *short_name*
   can be used to look up the with the
   :ref:`unit<unit-unit>` method.

-  **html_symbol, latex_symbol, ascii_symbol**: (str, optional,
   default is None): html, latex, and ascii versions of the symbol
   if they are different from the unicode representation of the symbol.
   A gummy is normally printed with a space between the value and the
   unit, however this space is removed if the symbol string starts with
   a tab character.

-  **description**: (str, optional, default is None) A description
   of the unit. Words used in the description can be searched using the
   search_units_ function.

-  **order**: (int, default is -1) When composite units are printed, the 
   symbols with the lowest *order* value will be printed to the left (unless
   this behavior is overridden with the :ref:`reorder<unit-reorder>` method).

Unit static methods
-------------------

.. _unit-unit:

-  static **unit(text, exception=True)**: This method is called whenever a 
   string referencing a Unit is 
   passed to a gummy property or method.  This method finds and returns a Unit
   instance from the Unit library. The parameter *text* may be a string representing
   the unit. The string can contain the name, short name or (if the unit
   was created with *add_symbol* set to ``True``) the symbol of the unit.
   This parameter may also be a combination of names and/or symbols of 
   several different units.
   Spaces or the character '\*' represent multiplication, the character
   '/' represents division and the string '\*\*' represents the power
   operator.  For example txt can be: ``'kg m**2/s'`` or equivalently
   ``'kilogram*metre*metre*second**-1'`` or ``'(kg/s)*m**2'``. If a unit
   name contains a space, '\*' or '/' character then the name must be
   enclosed in square brackets, e.g: ``'[light year]'``. If *text* is a
   Unit instance, then that instance is returned.  If *text* is the integer
   1 or the string ``'1'``, then the instance *one* is returned.
   If the parameter *exception* is True a ``UnitNotFoundError`` or
   ``UnitLibError`` is raised if a unit is not found that matches
   *text*. If the parameter *exception* is False and a unit is not found, then this
   method returns ``None`` without raising an exception. The default is
   ``True``.  

.. _unit-reorder:

-  static **reorder(txt)**: This changes the order in which the symbols of
   composite derived units are printed. For example::

       >>> print(Unit.unit('ft lb'))
       ft lb
       >>> print(Unit.unit('lb ft'))  #This is the same unit as above and prints identically
       ft lb
       >>> Unit.reorder('lb ft')  #Now the order will be changed when the unit is displayed
       >>> print(Unit.unit('ft lb'))  
       lb ft

-  static **alias(alias, unit)**: Creates an alias (an alternate name) that can
   be used to reference a Unit. The parameter *alias* is a string
   containing the new alias. The parameter *unit* is a string
   referencing the Unit (or the Unit instance itself) that will be
   assigned the alias.
   
   
Unit properties
---------------

-  **aliases**: (read-only) Returns a set of the unshadowed aliases of the unit. 
   
-  **shadowed_aliases**: (read-only) Returns a set of any aliases that have
   been shadowed by other unit definitions.
   
-  **is_dimensionless**: (read-only) Returns `True` if a conversion exists 
   between `self` and the Unit instance `one`, and `False` if not.

-  **units**: (read-only) Returns a list of the constituent units and their 
   exponents, e.g. for kg m**2/s *units* would return [(kg, 1), (m, 2), (s, -1)].
   

Unit sub-classes
----------------

For examples of unit definitions using the following Unit sub-classes
see the siunit.py and usunits.py modules in the gummy package

-  class metrolopy.\ **PrefixedUnit**: Creates a set of units with SI prefixes
   (..., kilo, mega, giga, ...).  For example::
   
       PrefixedUnit('metre','m',additional_names=('meter',),add_symbol=True,
                    order=1,description='SI unit of length',
                    base_description='SI base unit of length')
                    
       PrefixedUnit('inch','in',Conversion('m',0.0254),prefixes=['micro'],
                    add_symbol=True,description='unit of length')

   The definition above for the metre generates a set of units using all of the
   SI prefixes.  But the prefixes key word is used with the inch definition so 
   that only the inch and microinch are generated.

-  class metrolopy.\ **BinaryPrefixedUnit**: Creates a set of unit with binary
   prefixes (..., kibi, mebi, gibi, ...)

-  class metrolopy.\ **NonlinearUnit**: Abstract base class for ``LogUnit`` and
   ``OffsetUnit``

-  class metrolopy.\ **LogUnit**: Generates logrithmic units (e.g. decibel or 
   neper).  For example::
   
       LogUnit('decibel sound pressure level','dB',
               LogConversion(gummy(20,unit='uPa'),20,10,np.log10),
               short_name='dB(SPL)',add_symbol=False,
               description='sound pressure level in air')
               
   The conversion is defined here with::
   
       LogConversion(reference, multiplier, log_base, log_func, offset=0) 
   
   so that the conversion to the LogUnit from the *reference* unit is given by::
   
       multiplier*log_func(x/reference) + offset
       
   and the conversion back is given by::
   
       reference*log_base**(x - offset)/multiplier
       
   gummys with a LogUnit may only be added or subtrated from gummys with that
   same unit.  gummys with a LogUnit may only be multiplied or divided by 
   gummys with a linear unit.

-  class metrolopy.\ **OffsetUnit**: Generated units with an offset origin 
   (degree Celsius or degree Fahrenheit).  For example::
   
       OffsetUnit('degree Fahrenheit', '\u00B0F', OffsetConversion('degR',459.67),
                latex_symbol='^{\circ}F' ,ascii_symbol='degF', add_symbol=True, 
                description='unit of temperature')
                
   The conversion is defined with::
   
       OffsetConversion(unit, offset)
       
   where *unit* must be linear unit (with not offset origin) that differs from
   the OffsetUnit only by the offset.  In addition to the OffsetUnit, an
   _IntervalUnit is generated which has name equal to the OffsetUnit name with
   ' interval' appended and '-i' appended to the short name or symbol alias. 
   The _IntervalUnit appears when OffsetUnits are subtracted or when an 
   OffsetUnit is used in a _CompositeUnit.  A gummy with an OffsetUnit may be 
   added to another quanitity only if that quaitiy is a gummy with the 
   corresponding _IntervalUnit.
   
-  class metrolopy.\ **_CompositeUnit**:  Represents a derived unit made up 
   of several Unit (or Unit sub-class) instances combined.
   
-  class metrolopy.\ **_One**:  The instance of this class *one* represents
   the number 1 and is the default unit for a gummy.  Dimensionless units are
   defined as any Unit where a conversion to *one* exists.

.. _jummy:

class jummy
=====================

class
metrolopy.\ **jummy(real=None,imag=None,r=None,phi=None,cov=None,unit=one**

A jummy object represents a complex valued quantity with gummy real and
imaginary components.

jummy parameters
----------------

-  **real,imag,r,phi**: (float or gummy) The value may be specified
   in either Cartesian coordinates using the *real* and *imag*
   parameters or polar coordinates with the *r* and *phi* parameters.
   The pair *real*, *imag* or *r*, *phi* may both be gummy or both be
   float. If they are float then *cov* and *unit* may also be
   specified.

-  **cov**: (2 x 2 list, ``tuple`` or ``numpy.ndarray`` of
   float) The variance-covariance matrix for either the pair *real*,
   *imag* or the pair *r*, *phi*.

-  **unit**: (str or ``Unit`` or list or ``tuple``, or
   ``numpy.ndarray`` of length 2 of str or ``Unit``) Units for
   *real*, *imag* or *r*, *phi*. In the case that *real* and *imag* are
   specified with different units, there must exist a conversion between
   the two units. Units for *phi* must be dimensionless.

jummy properties
----------------

-  **x**: (read-only) returns ``complex(jummy.real.x,jummy.imag.x)``

-  **cov**: (read-only) returns the variance-covariance matrix between
   ``jummy.real`` and ``jummy.imag``

-  **real**: (read-only) a gummy representing the real part of the value

-  **imag**: (read-only) a gummy representing the imaginary part of the
   value

-  **unit**: Gets or sets the units of jummy.real and jummy.imag. If the
   units of jummy.real are different from jummy.imag then a tuple of
   Unit and length 2 is returned. Otherwise a :ref:`Unit<class-Unit>` instance
   is returned.

jummy methods
-------------

In addition to the methods below, the `gummy class display
methods <#gummy-display-methods>`__ can also be used with the jummy
class

-  **conjugate()**: returns the (jummy valued) complex conjugate

-  **angle()** returns a gummy representing Arg(jummy)

-  **copy(self,formatting=True)**: Returns a copy of the jummy. If the
   *formatting* parameter is ``True`` the display formatting information
   will be copied and if ``False`` the display formatting will be set to
   the default for a new jummy.

-  classmethod **apply(function, derivative, arg1, arg2, ...)**: A classmethod that
   applies a function to one or more jummy objects propagating the
   uncertainty.

   **apply parameters**:

    -  **function**: The the function to be applied. The function should
       take one or more float or ``complex`` arguments and return a
       float or ``complex`` value.
    
    -  **derivative**: The name of a second function that gives the
       derivatives with respect to the arguments of *function*. *derivative*
       should take an equal number of arguments as *function*. If *function*
       takes one argument *derivative* should return a float and if
       *function* takes more than one argument then derivative should return
       a ``tuple``, list or ``numpy.ndarray`` of float that contains
       the derivatives with respect to each argument. The derivatives with
       respect to each argument may be real or complex values, in which case
       *function* is assumed to be holomorphic. Or the derivative may be a 2
       x 2 matrix of the form::
    
                                 [[ du/dx, du/dy ],
                                  [ dv/dx, dv/dy ]]
    
       where function(x + j\ *y) = u + j*\ v.
    
    -  arg1, arg2, ...\*\*: One or more arguments to which *function* will
       be applied. These arguments need not all be jummy objects; arguments
       such as floats will be taken to be constants with no uncertainty.
       They may also be numpy.ndarrays in which case the usual numpy
       broadcasting rules apply.

    **apply returns**:  If none of the arguments *arg1*, *arg2*, ... are gummy 
    or jummy then the
    return value is the same type as the return value of *function*.
    Otherwise apply returns either a gummy or a jummy depending on whether
    *function* has a float or a complex return value.

-  classmethod **napply(function, arg1, arg2, ...)**: A classmethod that applies a
   function to one or more jummy objects propagating the uncertainty.
   This method is similar to ``jummy.apply`` except that the derivatives
   are computed numerically so a derivative function does not need to be
   supplied.

curve fitting
=============

.. _Fit:

class Fit
---------

class metrolopy.\ **Fit(x, y=None, f=None, p0=None, ux=None, uy=None,
sigma_is_known=True, xunit=None, yunit=None, solver=None,
maxiter=None, nprop=False, \*\*keywords)**

A class for performing non-linear fitting. The fit function may be
passed in the parameter *f* or may be specified by sub-classing ``Fit``
and overriding the f_ method.

Fit parameters
~~~~~~~~~~~~~~

All parameters except x are optional

-  **x**: The x-coordinates of the data. This is a list or
   numpy.ndarray of floats or gummys (all point must be of the same
   type, floats and gummys may not be mixed). The x-coordinates may be
   one dimensional or may be multi-dimensional. For d-dimensional
   coordinates with (with N total data points) this parameter should be
   of the form::

              [[x1[1], x1[2], ... , x1[N]],
               [x2[1], x2[2], ... , x2[N]],
               .
               .
               .
               [xd[1], xd[2], ... , xd[N]]]

   If gummys are given, then the must be dimensionless unless the
   get_puints method is implemented in a subclass.

-  **y**: The y-coordinates of the data (shape and type requirements are
   the same as for the x-coordinates). This may be omitted only if the
   odr solver is used.

-  **f**: The fit function. For d dimensional x-coordinates and k fit
   parameters it should be of the form f(x1,x2,...,xd,p1,p2,...,pk) and
   return a float or (if *y* is multi-dimensional) a list or array of
   floats. This parameter is required unless the f method is overridden
   in a subclass.

-  **p0**: (list or numpy.ndarray of float) The inital
   values for the fit parameters. This parameter is required unless the
   get_p0 method is overridden in a subclass.

-  **ux**: (float or list or numpy.ndarray of float)
   Uncertainty in the x values. This should not be specified if the *x*
   argument contains gummys. If *ux* is specified then only the odr
   solver may be used. The default is None.

-  **uy**: (float or list or numpy.ndarray of float)
   Uncertainty in the y values. This should not be specified if the y
   argument contains gummys. The default is None.

-  **sigma_is_known**: (bool) If this is ``True`` then any
   uncertainties in the data (either as gummys in the *x* or *y* values
   or in the *ux* or *uy* parameters) are used to calculate the
   uncertainties in the fit. Otherwise, the uncertainties are based on
   the standard deviation of the residuals and the uncertainties in the
   data are used only for weighting the data points. This parameter is
   ignored if *nprop* is True.

-  **xunits**, **yunits** (str, default None) units for the x
   and y coordinates. These should not be specified if the *x* and *y*
   parameters contain gummys. These may only be specified if the
   get\_punits method is overridden in a subclass.

-  **solver**: ({'nls' , 'odr', None}) If
   *solver* = 'nls' then scipy.optimize.leastsq is used to perform the
   fit. If *solver* = 'odr' then scipy.odr is used. 'nls' may not be
   used if the y-coordinate is None or multi-dimensional or if there is
   uncertainty in the x-coordinates. If this is None, then 'nls' will be
   used when possible. Any keyword parameters not recognized by ``Fit``
   will be passed to the solver.

-  **maxiter**: (int) The maximum number of iterations that the
   solver may use. f this is None or omitted then the default value for
   the solver will be used.

-  **nprop**: (bool, default ``False``) If this is ``True`` then
   uncertainties in the fit will be numerically calculated by varying
   each data point. This will not work if there are more than a few data
   points or if the fit is not very stable. If this is ``False`` than
   the covariance matrix generated by the solver will be used to
   calculate the uncertainties.

-  **other keywords**: Any additional keyword parameters will be passed
   to the solver.

Fit properties
~~~~~~~~~~~~~~

-  **p**: (read-only, list of gummy) The fitted values for the fit
   function parameters as correlated gummys

-  **pf**: (read-only, list of float) The fitted values for the
   fit function parameters as floats

-  **res**: (read-only, numpy.ndarray of float ) the fit residuals

-  **s**: (read-only,float) the standard deviation (or, when there
   are uncertainties for the input data, the square root of the reduced
   chi-squared) of the residuals

-  **cov**: (read-only, numpy.ndarray) the covariance matrix
   generated by the solver

-  **fit_output**: (read-only) the raw output of the solver

-  **x**: (read-only, numpy.ndarray of float or of gummy)
   numpy array of the x-coordinates of the data.

-  **xf**: (read-only, numpy.ndarray of float) numpy array of the
   x-coordinates of the data as floats

-  **xdim**: (read-only,int) the dimension of the x-coordinates

-  **ux**: (read-only, float, numpy.ndarray of float or
   None): uncertainties in the x-coordinates

-  **y**: (read-only, numpy.ndarray of float or of gummy)
   numpy array of the y-coordinates of the data.

-  **yf**: (read-only, numpy.ndarray of float) numpy array of the
   yx-coordinates of the data as floats

-  **ydim**: (read-only,int) the dimension of the y-coordinates

-  **uy**: (read-only, float, numpy.ndarray of float or
   None): uncertainties in they-coordinates

-  **count**: (read-only, int) the number of data points

-  **p0**: (read-only, list of float) The initial values for the
   fit function parameters

-  **solver**: (read-only, str) the solver used

-  **punits**: (read-only, list of :ref:`Unit<class-Unit>`) the units of
   the fit parameters

-  **nparam**: (read-only, int) the number of fit parameters

Fit methods
~~~~~~~~~~~

-  **ypred(x1,x2,...)**: Takes xdim floats and returns a gummy
   representing the predicted value at that x-coordinate.

-  **ypredf(x1,x2,...)**: Takes xdim floats and returns a float giving
   the predicted value at that x-coordinate.

-  **plot(data\_format='ko', data\_options={}, show\_data=True,
   error\_bars=True, error\_bar\_k=1, fit\_format='k-', fit\_options={},
   show\_fit=True, cik=None, cip=None, ciformat='g-', cioptions={},
   clk=None, clp=None, clformat='r-', cloptions = {}, xmin=None,
   xmax=None, xlabel=None, ylabel=None, hold=False,
   plot\_points=None)**: plots the data (only available if *x* and *y* are
   one-dimensional)

   **plot parameters** (all parameters are optional):
   
   -  **data\_format**: (str) The format string passed to pyplot.plot
      or pyplot.errorbar when plotting the data points. The default is
      'ko'.

   -  **data\_options**: (dict) A dictionary containg key words that
      are passed to pyplot.plot or pyplot.errorbar when plotting the data
      points.

   -  **show\_data**: (bool) Whether or not to plot the data points.
      The default is ``True``.

   -  **error\_bars**: (bool) Whether or not to plot error bars on
      the data points (if uncertainty values were defined for the data).
      The default is ``True``.

   -  **error\_bar\_k**: (float or int) Coverage factor for the
      error bars. The length of the error bars are determined by
      multiplying the standard uncertainty for each data point by this
      quantity. The default value is 1.

   -  **fit\_format**: (str) The format string passed to pyplot.plot
      or pyplot.errorbar when plotting the fitted curve. The default is
      'k-'.

   -  **fit\_options**: (dict) A dictionary containg key words that
      are passed to pyplot.plot or pyplot.errorbar when plotting the
      fitted curve.

   -  **show\_fit**: (bool) Whether or not to plot the fitted curve.
      The default is ``True``.

   -  **xmin** and **xmax**: (float) The lower and upper limits of
      the fitted, confidence interval and control limit curves. If this
      is None, the limits are equal to x1 +/- (x2 -
      x1)\*Fit.over\_plot where x1 is the x value of the first data
      point, x2 is the x value of the last data point and Fit.over\_plot
      is an attribute of the Fit object with default value 0.05.

   -  **xlabel** and **ylabel**: (str) Labels for the x and y axes.
      If units are defined for the x or y axes, the unit symbol will be
      added to the end of the labels defined here. If these are set to
      None, then the values of the ``Fit.xlabel`` and ``Fit.ylabel``
      attributes will be used. The default is None.

   -  **plot\_points**: (int) The number of points to use in each
      curve when plotting the fit, confidence interval, and control
      limit curves. If this is set to None, then the value of the
      Fit.plot\_points attribute will be used, which has a default value
      of 100.

   -  **hold**: (bool) If hold is ``False`` then ``pyplot.show()`` is
      executed just before this function returns.

   -  **cik**: (float or None) Coverage factor for the
      uncertainty bands in the plot. If *cik* and *cip* are None (
      the default values) then uncertainty bands will not be shown. Do
      not specify both *cik* and *cip*.

   -  **cip**: (float or None) Confidence level for the
      uncertainty bands in the plot. If *cik* and *cip* are None (
      the default values) then uncertainty bands will not be shown. Do
      not specify both *cik* and *cip*.

   -  **ciformat**: (str, default is 'g-') Format string passes
      to the pyplot.plot command that plots the uncertainty bands.

   -  **cioptions**: (dict) Keywork options passed to the pyplot.plot
      command that plots the uncertainty bands.

   -  **clk**,\ **clp**,\ **clformat**, and **cloptions**: Control limit
      options, same as above for the uncertainty bands. The control
      limit band is the control limit coverage factor multiplied by the
      RSS of the fit uncertainty and the standard deviation of the
      residuals.

Fit abstract methods
~~~~~~~~~~~~~~~~~~~~

These methods may be overridden when sub-classing ``Fit``

.. _f:

-  **f(self,x1,x2,...,xd,p1,p2,...,pk)**: The fit function. The
   function to fit. It must either have signature f(self,x,p1,p2,...,pn)
   where there are p1 to pn are the n fit parameters and the independent
   variable *x* has one dimension, or f(self,x1,x2,...,xm,p1,p2,...,pn)
   where the independent variable x has m dimensions at each
   observation. *f* should return either a float or a 1-d array of
   floats depending on the dimension of the response variable *y*.

-  **jac(self,x1,x2,...,xk,p1,p2,...,pk)**: The Jacobian. This method
   may optionally be overridden in a derived class. If not overridden,
   this method throws a NotImplementedError() the derivatives will be
   calculated numerically.

   It must have the same signature as the *f* method and return a list
   of derivatives of the form::

           [df/dx1,df/dx2,...,df/dp1,df/dp2,...] 

   if f returns a scalar or::

           [[df1/dx1,df1/dx2,...,df1/dp1,df1/dp2,...],
            [df2/dx1,df2/dx2,...,df2/dp1,df2/dp2,...],...]

   if f returns a 1-d array [f1,f2,...].

-  **get\_p0(self)**: Returns an initial guess for the fit parameters
   based on the input *x* and *y* data. This is not required, but if it
   is not implemented then the p0 parameter is a required parameter for
   the **init** method.

-  **get\_punits(self)**: Returns a list units for the fit parameters.
   This is not required, but if it is not implemented then only float
   values or dimensionless gummys may be as the x and y parameters and
   the xunit and yunit parameters to the **init** method may not be
   used.

-  **funicode(self), flatex(self), fhtml(self)**: Returns a str
   containing unicode, latex, and html representations of the fit
   function.

sub-classes of Fit for some common functions
--------------------------------------------

-  metrolopy.\ **PolyFit(x, y, deg=1, p0=None, ux=None, uy=None,
   sigma\_is\_known=True, xunit=None, yunit=None, solver=None,
   maxiter=None, nprop=False, \*\*keywords)**:  Fits the *x*, *y* data to a 
   polynomial. In addition to the parameters
   for Fit_, ``PolyFit`` takes the parameter *deg* which is
   the degree of the polynomial. The solver parameter can take the
   string ``'ols'`` in addition to the 'odr' and 'nls'
   solvers defined by class ``Fit``. A linear
   fit algorithm, ols, will be used by default if *x* and *y* are one
   dimensional and there is no uncertainty in the x-values. The odr
   solver must be used if there is uncertainty in the x-values or if the
   y-coordinates are multi-dimensional. By the nonlinear least squares
   solver, nls, will be used by if *x* is multi-dimensional. Initial
   values *p0* may be specified if the nls or odr solvers are used, but
   are not required. Both the *x* and *y* data may have units.

-  metrolopy.\ **SinFit(x, y, p0=None, ux=None, uy=None,
   sigma\_is\_known=True, xunit=None, yunit=None, solver=None,
   maxiter=None, nprop=False, \*\*keywords)**: Fits the x,y data to a
   function of the form::
   
       p[0]*sin(p[1]*x + p[2]) + p[3]
       
   See the Fit_ class for the parameters, properties and methods of
   this class.  This class is pretty good at guessing initial parameters.
   
-  metrolopy.\ **ExpFit(x, y, p0=None, ux=None, uy=None,
   sigma\_is\_known=True, xunit=None, yunit=None, solver=None,
   maxiter=None, nprop=False, \*\*keywords)**: Fits the x,y data to a
   function of the form::
   
        p[0]*exp(x/p[1]) + p[2]
        
   See the Fit_ class for the parameters, properties and methods of
   this class.

-  metrolopy.\ **DoubleExpFit(x, y, p0=None, ux=None, uy=None,
   sigma\_is\_known=True, xunit=None, yunit=None, solver=None,
   maxiter=None, nprop=False, \*\*keywords)**: Fits the x,y data to a
   function of the form::
   
       p[0]*exp(x/p[1]) + p[2]*exp(x/p[3]) + p[4]

   See the Fit_ class for the parameters, properties and methods of 
   this class.

-  metrolopy.\ **OneOverTFit(x, y, p0=None, ux=None, uy=None,
   sigma\_is\_known=True, xunit=None, yunit=None, solver=None,
   maxiter=None, nprop=False, \*\*keywords)**: Fits the x,y data to a
   function of the form::
   
       p[0]/x + p[1]
       
   See the Fit_ class for the parameters, properties and methods of 
   this class.



