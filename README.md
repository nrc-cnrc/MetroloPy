# MetroloPy

tools for dealing with physical quantities:  uncertainty propagation and unit conversion

---

MetroloPy is a pure python package and requires Python 3 and the SciPy stack (NumPy, SciPy, Pandas, and IPython).  It looks best in a Jupyter Notebook.

Install MetroloPy with pip:

```
$ pip install metrolopy
```

Physical quantities can then be represented in Python as `gummy` objects:

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
