# -*- coding: utf-8 -*-

# module printing

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
This module contains code to support pretty printing of gummys using LaTex and 
HTML.
"""

try:
    from IPython.display import display,Math,HTML,Markdown
    ipython_installed = True
except ImportError:
    ipython_installed = False


def print_html(text):
    if ipython_installed:
        display(HTML(text))
    else:
        raise NotImplementedError('ipython must be installed to print in HTML format')
        
def print_latex(text):
    if ipython_installed:
        display(Math(text))
    else:
        raise NotImplementedError('ipython must be installed to print in HTML format')
        
def print_markdown(text):
    if ipython_installed:
        display(Markdown(text))
    else:
        raise NotImplementedError('ipython must be installed to print in Markdown format')
        

        
    
# Called from _repr_latex_ or Ipyhton.display.Math and gummy.latex_math 
# Not sure if this is necessary; may depend on the version of Ipython.
def _latex_math(text):
    return '$ ' + text + ' $'
 
def _latex_norm(text):
    return '\\text{' + text + '}'

# Same as above for text used in plot titles and labels
def _latex_norm_plot(txt):
    txt = txt.replace(' ',' \\ ')
    return '\\mathrm{' + txt + '}'


def _find_printer(value):
    if value is None:
        return 'any_but_latex'
    value = value.lower().strip()
    if value not in ['latex','html','ascii','unicode','any','any_but_latex']:
        raise ValueError('printer ' + str(value) + ' is not recognized')
    return value
    
class MetaPrettyPrinter(type):
    @property
    def printer(self):
        """
        Get or set the prefered display printer.  This is a string with one of 
        the following values:
        
        "any", "latex", "html", "unicode", "ascii", or "any_but_latex"
        
        "any" will usually pick html or latex output when running in an IPython 
        console or Jupyter notebook and unicode otherwise.
        
        "any_but_latex" will usually pick html when running in an IPython 
        console or Jupyter notebook and unicode otherwise.
        
        "latex" and "html" are only available when running under IPython.  If 
        these printers are not available the display will default to "unicode".
        """
        return PrettyPrinter._printer
    @printer.setter
    def printer(self,value):
        PrettyPrinter._printer = _find_printer(value)
        
class PrettyPrinter(metaclass=MetaPrettyPrinter):
    # Base class that can be inherited to provide HTML and LaTeX ouput to an
    # ipython console or Jupyter notebook.
    
    # The inheriting class must define a tostring(self,fmt='unicode',**kwds)
    # method that accepts fmt values in ['html',latex','unicode','ascii'].
    
    break_on_printing_error = False
    
    # Set latex_math to None in the inherited class if the output of the latex
    # and _repr_latex should not be in LaTeX math mode.
    latex_math = _latex_math
    
    latex_norm = _latex_norm
    latex_math_plot = _latex_math
    latex_norm_plot = _latex_norm_plot
    
    _printer = 'any_but_latex'
            
    @property
    def printer(self):
        """
        Get or set the preferred display printer.  This is a string with one of
        the following values:
        
        "any", "latex", "html", "unicode", "ascii", or "any_but_latex"
        
        "any" will usually pick html or latex output when running in an IPython 
        console or Jupyter notebook and unicode otherwise.
        
        "any_but_latex" will usually pick html when running in an IPython 
        console or Jupyter notebook and unicode otherwise.
        
        "latex" and "html" are only available when running under IPython.  If 
        these printers are not available the display will default to "unicode".
        """
        return self._printer
    @printer.setter
    def printer(self,value):
        self._printer = _find_printer(value)
        
    def _repr_html_(self):
        if self._printer not in ['any','html','any_but_latex']:
            return None

        txt = self.tostring(fmt='html')

        return txt
            
    def _repr_latex_(self):
        if self._printer not in ['any','latex']:
            return None
        
        txt = self.tostring(fmt='latex')

        if self.latex_math is None:   
            return txt
        return type(self).latex_math(txt)
    
    def __str__(self):
        return self.tostring(fmt='unicode')
            
    def __repr__(self):
        if self. _printer == 'ascii':
            return self.tostring(fmt='ascii')
        return self.tostring(fmt='unicode')
            
    def __format__(self,fmt):
        return self.tostring(fmt=fmt)
    
    def tohtml(self,**kwds):
        """
        Returns a string representing this object formatted for html; equivalent
        to PrettyPrinter.tostring(fmt='html',**kwds).  See the tostring method.
        """
        return self.tostring(fmt='html',**kwds)
    
    def tolatex(self,**kwds):
        """
        Returns a string representing this object formatted for LaTeX; equivalent
        to PrettyPrinter.tostring(fmt='latex',**kwds).  See the tostring method.
        """
        return self.tostring(fmt='latex',**kwds)
    
    def tounicode(self,**kwds):
        """
        Returns a string representing this object; equivalent to 
        PrettyPrinter.tostring(fmt='unicode',**kwds) and PrettyPrinter.__str__().  
        See the tostring method.
        """
        return self.tostring(fmt='unicode',**kwds)
    
    def toascii(self,**kwds):
        """
        Returns a string representing this object formatted using only ASCII
        characters; equivalent to PrettyPrinter.tostring(fmt='ascii',**kwds).  
        See the tostring method.
        """
        return self.tostring(fmt='ascii',**kwds)
    
    def latex(self,math=None,**kwds):
        """
        Prints a representation of the object using LaTeX formatting if this method
        is called from an IPython console or Juptyer notebook. See the tostring
        method.
        """
        if math is None:
            math = type(self).latex_math
        
        if math is None:
            print_latex(self.tolatex(**kwds))
        else:
            print_latex(math(self.tolatex(**kwds)))
        
    def html(self,**kwds):
        """
        Prints a representation of the object using HTML formatting if this method
        is called from an IPython console or Juptyer notebook. See the tostring 
        method.
        """
        print_html(self.tohtml(**kwds))
        
    def unicode(self,**kwds):
        """
        Prints a representation of the object. Equivalent to 
        print(cls.tostring(fmt='unicode')).  See the tostring method.
        """
        print(self.tounicode(**kwds))
        
    def ascii(self,**kwds):
        """
        Prints a representation of the object using only ASCII characters.
        Equivalent to print(cls.tostring(fmt='ascii')).  See the tostring
        method.
        """
        print(self.toascii(**kwds))
    
    
def set_printer(value): 
    """
    Sets the preferred default display printer.  This is a string with one of
    the following values:
    
    "any", "latex", "html", "unicode", or "ascii"
    
    "any" will usually pick html or latex output when running in an IPython 
    console or Jupyter notebook and unicode otherwise.
    
    "latex" and "html" are only available when running under IPython.  If 
    these printers are not available the display will default to "unicode".
    """
    if value is None:
        PrettyPrinter._printer = 'any_but_latex'
        return
    value = value.lower().strip()
    if value not in ['latex','html','ascii','unicode','any','any_but_latex']:
        raise ValueError('printer ' + str(value) + ' is not recognized')
    PrettyPrinter._printer = value
