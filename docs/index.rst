.. metrolopy documentation master file, created by
   sphinx-quickstart on Fri Feb 22 16:15:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
===================   
MetroloPy, the docs
===================

tools for dealing with physical quantities:  uncertainty propagation and unit 
conversion


getting started
===============

MetroloPy is a pure python package and requires Python 3 and the 
`SciPy <https://www.scipy.org/>`_ stack (NumPy, SciPy, Pandas, and IPython).  
It looks best in a `Jupyter notebook <https://jupyter.org/>`_.

Install MetroloPy with ``pip install metrolopy``  or 
``conda install -c conda-forge metrolopy``.

Physical quantities can then be represented in Python as `gummy`_ objects
with an uncertainty and (or) a unit:

.. _gummy: hand_made_doc.html#class-gummy

.. raw:: html

    <div class="highlight-default"><div class="highlight"><pre><span></span><span class="gp">&gt;&gt;&gt; </span><span class="kn">import</span> <span class="nn">metrolopy</span> <span class="k">as</span> <span class="nn">uc</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">a</span> <span class="o">=</span> <span class="n">uc</span><span class="o">.</span><span class="n">gummy</span><span class="p">(</span><span class="mf">1.2345</span><span class="p">,</span><span class="n">u</span><span class="o">=</span><span class="mf">0.0234</span><span class="p">,</span><span class="n">unit</span><span class="o">=</span><span class="s1">&#39;cm&#39;</span><span class="p">)</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">a</span>
    <span class="go">1.234(23) cm</span>
    
    <span class="gp">&gt;&gt;&gt; </span><span class="n">b</span> <span class="o">=</span> <span class="n">uc</span><span class="o">.</span><span class="n">gummy</span><span class="p">(</span><span class="mf">3.034</span><span class="p">,</span><span class="n">u</span><span class="o">=</span><span class="mf">0.174</span><span class="p">,</span><span class="n">unit</span><span class="o">=</span><span class="s1">&#39;mm&#39;</span><span class="p">)</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">f</span> <span class="o">=</span> <span class="n">uc</span><span class="o">.</span><span class="n">gummy</span><span class="p">(</span><span class="n">uc</span><span class="o">.</span><span class="n">UniformDist</span><span class="p">(</span><span class="n">center</span><span class="o">=</span><span class="mf">0.9345</span><span class="p">,</span><span class="n">half_width</span><span class="o">=</span><span class="mf">0.096</span><span class="p">),</span><span class="n">unit</span><span class="o">=</span><span class="s1">&#39;N&#39;</span><span class="p">)</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">p</span> <span class="o">=</span> <span class="n">f</span><span class="o">/</span><span class="p">(</span><span class="n">a</span><span class="o">*</span><span class="n">b</span><span class="p">)</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">p</span>
    <span class="go">2.50(21) N/cm<sup>2</sup></span>
    
    <span class="gp">&gt;&gt;&gt; </span><span class="n">p</span><span class="o">.</span><span class="n">unit</span> <span class="o">=</span> <span class="s1">&#39;kPa&#39;</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">p</span><span class="o">.</span><span class="n">uunit</span> <span class="o">=</span> <span class="s1">&#39;%&#39;</span>
    <span class="gp">&gt;&gt;&gt; </span><span class="n">p</span>
    <span class="go">25.0 kPa &plusmn; 8.5%</span>
    </pre></div>
    </div>
    
MetroloPy can do much more including `Monte-Carlo uncertainty propagation`_, 
generating `uncertainty budget tables`_, and `curve fitting`_.  It can also handle 
`expanded uncertainties`_, `degrees of freedom`_, `correlated quantities`_, and 
`complex valued quantities`_.

..  _Monte-Carlo uncertainty propagation: _static/tutorial.html#Monte-Carlo-uncertainty-propagation

.. _uncertainty budget tables: _static/tutorial.html#an-uncertainty-budget

.. _curve fitting: _static/tutorial.html#curve-fitting

.. _expanded uncertainties: _static/tutorial.html#expanded-uncertainty

.. _degrees of freedom: _static/tutorial.html#degrees-of-freedom-and-uncertainty-types

.. _correlated quantities: _static/tutorial.html#creating-correlated-gummys

.. _complex valued quantities: hand_made_doc.html#jummy

further reading
===============

* `a tutorial <_static/tutorial.html>`_ (or  :download:`download the tutorial as a Jupyter notebook <tutorial.ipynb>`)
* :ref:`the API documentation <hand_made_doc>` (see the auto-generated docs below for an alternative set of documentation)
* `a list of the measurement units built into MetroloPy <_static/units.html>`_
* `a list of the physical constants built into MetroloPy <https://nrc-cnrc.github.io/MetroloPy/_static/constants.html>`_
* :ref:`package development and to do <todo>`
* `the issues page on GitHub <https://github.com/nrc-cnrc/Metrolopy/issues>`_
* `the source code on GitHub <https://github.com/nrc-cnrc/Metrolopy/>`_


the auto-generated docs (indices and tables)
============================================

These pages were automatically generated from the doc strings in the source
code.  They are slightly more comprehensive but perhaps slightly more confusing
than the handmade API docs referenced in the further reading section above.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Where did the name for the gummy class come from?
=================================================

The name comes from the 
`JGCM/ISO Guide to the Expression of Uncertainty in Measurement <https://www.bipm.org/en/publications/guides/gum.html>`_ 
which is also known as the GUM.  The GUM is the international standard used by
metrology institutes and calibration labs for expressing and propagating 
measurement uncertainty.  The gummy object implements many of the 
recommendations outlined in the GUM.

Other references of note are the 
`draft 9th edition of the SI Brochure <https://www.bipm.org/en/measurement-units/rev-si/#communication>`_ 
which contains the definitions of the SI units, and 
`NIST Special Publication 1038 <https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1038.pdf>`_ 
where many of the US customary units are defined.  Some physical constants
used for unit definitions are from the 
`2014 CODATA recommended values <http://www.codata.org/committees-and-groups/fundamental-physical-constants>`_
and the 
`IAU 2009 system of astronomical constants <http://maia.usno.navy.mil/NSFA/IAU2009_consts.html>`_.


version history
===============

* Version 0.5.0, built 26 March 2019, is the first public release.
* Version 0.5.1, built 2 April 2019, fixed a major bug that generated negative 
  uncertainties in some cases and fixed some other minor bugs.  Improved support 
  for fraction.Fraction and mpmath.mpf values.
* Version 0.5.2, built 5 April 2019, fixed a major bug that propagated
  uncertainty incorrectly if a gummy was created with an uncertainty set with
  an integer data type.  Fixed several other minor bugs.
* Version 0.5.3, built 10 April 2019, minor bug fixes.
* Version 0.5.4, built 15 April 2019, minor bug fixes.
* Version 0.5.5, built 7 May 2020, minor bug fixes.
* Version 0.5.6, built 24 September 2020, minor bug fixes.
* Version 0.5.7, built 26 September 2020, minor change to setup.py.
* Verison 0.6.0, built 19 October 2020, added the `Quantity` and `immy` classes 
  as well as a library of physical constants.
* Version 0.6.1 built 19 October 2020, bug fixes
* Version 0.6.1 built 10 July 2021, bug fixes


author
======

Harold Parks, parksh@nrc.ca, National Research Council Canada


license
=======

MetroloPy is distributed as free software under the terms of the 
`GNU General Public License version 3 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_. 
In practice, this license imposes 
no restriction on using MetroloPy. However, if you want to further convey 
verbatim or modified versions of the code you must do so under the same license 
terms. Please contact NRC if you wish to license MetroloPy under different 
terms.
