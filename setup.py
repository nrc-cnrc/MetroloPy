# -*- coding: utf-8 -*-
"""
Created on Wed May 10 09:36:41 2017

@author: h parks
"""

from setuptools import setup

__version__ = None
with open('metrolopy/version.py') as f:
    exec(f.read())
    
with open("README.md", "r") as fh:
    long_description = fh.read()
            
setup(name = 'metrolopy',
      version = __version__,
      description = 'tools for dealing with measured quantities: uncertainty propagation and unit conversion',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author = 'Harold Parks, National Research Council Canada',
      author_email = 'parksh@nrc.ca',
      url = 'http://nrc-cnrc.github.io/MetroloPy/',
      packages = ['metrolopy','metrolopy.tests'],
      package_data = {'metrolopy':['license.txt']},
      setup_requires = ['numpy'],
      install_requires=['scipy','matplotlib','pandas'],
      extras_require = {'pretty':['IPython']},
      zip_safe = True,
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Framework :: Jupyter",
        "Framework :: IPython",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Physics"
        ]   
      )


    
    
