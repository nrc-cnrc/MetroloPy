# MetroloPy

tools for dealing with physical quantities:  uncertainty propagation and unit conversion

---

MetroloPy is a pure python package and requires Python 3.5 or later and the SciPy stack (NumPy, SciPy and Pandas).  It looks best in a Jupyter Notebook.

Install MetroloPy with `pip install metrolopy`  or 
`conda install -c conda-forge metrolopy`.

Physical quantities can then be represented in Python as `gummy` objects with an uncertainty and (or) a unit:

<pre><code>&gt;&gt;&gt; import metrolopy as uc
&gt;&gt;&gt; a = uc.gummy(1.2345,u=0.0234,unit='cm')
&gt;&gt;&gt; a
1.234(23) cm

&gt;&gt;&gt; b = uc.gummy(3.034,u=0.174,unit='mm')
&gt;&gt;&gt; f = uc.gummy(uc.UniformDist(center=0.9345,half_width=0.096),unit='N')
&gt;&gt;&gt; p = f/(a*b)
&gt;&gt;&gt; p
2.50(21) N/cm<sup>2</sup>

&gt;&gt;&gt; p.unit = 'kPa'
&gt;&gt;&gt; p.uunit = '%'
&gt;&gt;&gt; p
25.0 kPa &plusmn; 8.5%
</code></pre>

MetroloPy can do much more including Monte-Carlo uncertainty propagation, generating uncertainty budget tables, and curve fitting.  It can also handle expanded uncertainties, degrees of freedom, correlated quantities, and complex valued quantities. See:

* [a tutorial](https://nrc-cnrc.github.io/MetroloPy/_build/html/_static/tutorial.html) (or  <a href="https://nrc-cnrc.github.io/MetroloPy/_build/html/_downloads/tutorial.ipynb" download> download the tutorial as Jupyter notebook</a>)
* [the documentation](https://nrc-cnrc.github.io/MetroloPy/)
* [the issues page on GitHub](https://github.com/nrc-cnrc/Metrolopy/issues)
* [a list of the units built into MetroloPy](https://nrc-cnrc.github.io/MetroloPy/_static/units.html)
* [a list of the physical constants built into MetroloPy](https://nrc-cnrc.github.io/MetroloPy/_static/constants.html)

## new in version 0.6.0

* A constant library has been added with physical constants that can be accessed
  by name or alias with the `constant` function.  The `search_constants` function 
  with no argument gives a listing of all built-in constants.  Each constant 
  definition includes any correlations with other constants.

* The `Quantity` class has been added to represent a general numerical value
  multiplied by a unit and the `unit` function has been added to retrieve
  `Unit` instances from the unit library by name or alias.  `Unit` instances 
  can now be multiplied and divided by other `Unit` instances to produce
  composite units, can be multiplied and divided by numbers to produce 
  `Quantity` instances or multiply or divide `Quantity` instances.  The 
  `gummy` class is now a subclass of `Quantity` with a `nummy` value rather 
  than a subclass of `nummy`.  A `QuantityArray` class has been introduced
  to represent an array of values all with the same unit.  Multiplying a `Unit`
  instance by a list, tuple, or numpy array produces a `QuantityArray` instance.

* The `immy` class has been introduced as an `ummy` valued counterpart of the 
  `jummy` class for representing complex values with uncertainties.  `immy` 
  and `jummy` values can now be displayed in a polar representation in addition 
  to a cartesian representation.  `immy` and `jummy` .r and .phi properties 
  have been added to access the magnitude and argument of the values as a 
  complement to the .real and .imag properties.



