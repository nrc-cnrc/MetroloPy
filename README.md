
# MetroloPy

tools for dealing with physical quantities:  uncertainty propagation and unit conversion

---

MetroloPy is a pure python package and requires Python 3 and the SciPy stack (NumPy, SciPy, Pandas, and IPython).  It looks best in a Jupyter Notebook.

Install MetroloPy with pip:

```
$ pip install metrolopy
```

Physical quantities can then be represented in Python as `gummy` objects:

```
>>> import metrolopy as uc
>>> a = uc.gummy(1.2345,u=0.0234,unit='cm')
>>> a
1.234(23) cm

>>> b = uc.gummy(3.034,u=0.174,unit='mm')
>>> f = uc.gummy(uc.UniformDist(center=0.9345,half_width=0.096),unit='N')
>>> p = f/(a*b)
>>> p
2.50(21) N/cm²

>>> p.unit = 'kPa'
>>> p.uunit = '%'
>>> p
25.0 kPa ± 8.5%
```

MetroloPy can do much more including Monte-Carlo uncertainty propagation, generating uncertainty budget tables, and curve fitting.  It can also handle expanded uncertainties, degrees of freedom, correlated quantities, and complex valued quantities. See:

* [a tutorial](https://nrc-cnrc.github.io/MetroloPy/_static/tutorial.html) (or  <a href="https://nrc-cnrc.github.io/MetroloPy/_downloads/tutorial.ipynb" download> download the tutorial as Jupyter notebook</a>)
* [the documentation](https://nrc-cnrc.github.io/MetroloPy/)
