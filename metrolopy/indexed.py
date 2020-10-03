# -*- coding: utf-8 -*-

# module indexed

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
indexed
"""

from importlib import import_module


class _Builtin:
    """
    See the Indexed._builtin() class methos for the use of this class
    """
    def __init__(self,_cls):
        self._cls = _cls
        
    def  __enter__(self):
        self._cls._open_lib = self._cls._builtin_lib
    def __exit__(self,type, value, traceback):
        self._cls._open_lib = self._cls._lib

class Indexed:

    # to speed loading of the module, only wait to create instances until they
    # are needed.  If a instance cannot be found the modules below are loaded in 
    # reverse order.
    #_builtins_to_import = []
    
    # define these in a subclass at the class level:
    # _builtin_lib = {}
    # _lib = {}
    # _open_lib = _builtin_lib
    
    @classmethod
    def alias(cls,alias,inst):
        """
        Creates an alias that can be used to reference an instance.
        
        Parameters
        ----------
        alias:  `str`
            a string containing the new alias

        inst:  `str` or `Indexed` instance
            A string referencing the `Indexed` instance that will be assigned
            the alias or the `Indexed` instance its self.
        """
        inst = cls.get(inst)
        cls._open_lib[alias] = inst
        inst._add_alias(alias)
        
    @staticmethod
    def load(library_name):
        import_module(library_name)
    
    @classmethod
    def get(cls,name,exception=True):
        if isinstance(name,cls):
            return name
        c = cls._lib.get(name)
        if c is None:
            c = cls._builtin_lib.get(name)
            if c is None:
                while len(cls._builtins_to_import) > 0:
                    import_module(cls._builtins_to_import.pop(),cls.__module__)
                    c = cls._builtin_lib.get(name)
                    if c is not None:
                        break    
    
        if c is None and exception:
            cls._raise_not_found(name)
        return c
    
    @classmethod
    def _raise_not_found(cls,name):
        raise ValueError(cls.object_name + ' "' + str(name) + '" was not found')
    
    @classmethod
    def _builtin(cls):
        """
        Allows Indexed instances to be added to Indexed._builtin_lib rather
        than Indexed._lib using a "with" statement.
        
        Example
        -------
        
        >>> with Unit._builtin():
        ...     one = Unit('1','',add_symbol=False)
        
        """
        return _Builtin(cls)
    
    @staticmethod
    def format_latex(text):
        return text
    
    def __init__(self,name,symbol=None,short_name=None,add_symbol=False,
                 html_symbol=None,latex_symbol=None,ascii_symbol=None,
                 description=None):
        self.name = name
        self.description=description
        
        if symbol is None:
            self.symbol = name
        else:
            self.symbol = symbol
            
        if html_symbol is None:
            self.html_symbol = symbol
        else:
            self.html_symbol=html_symbol
        if latex_symbol is None:
            self.latex_symbol = self.format_latex(symbol)
        else:
            self.latex_symbol = self.format_latex(latex_symbol)
        if ascii_symbol is None:
            self.ascii_symbol=symbol
        else:
            self.ascii_symbol=ascii_symbol
            
        self._aliases = set()
        
        self._open_lib[name] = self
        if add_symbol:
            if symbol.strip() != name:
                self._open_lib[symbol.strip()] = self
                self._aliases.add(symbol.strip())
            if (ascii_symbol is not None and ascii_symbol.strip() != name and 
                   ascii_symbol.strip() != symbol.strip()):
                self._open_lib[ascii_symbol.strip()] = self
                self._aliases.add(ascii_symbol.strip())
        if short_name is not None:
            self._open_lib[short_name] = self
            self.short_name = short_name
            self._aliases.add(short_name)
        else:
            if add_symbol:
                if ascii_symbol is not None:
                    self.short_name = ascii_symbol
                else:
                    self.short_name = symbol
            else:
                self.short_name = name
        
    @property
    def aliases(self):
        """read-only
        
        Returns a set of the unshadowed aliases of the instannce.  To add aliases
        to the unit use the `Indexed.alias` static method.
        """
        return {a for a in self._aliases if (self.get(a,exception=False) is self and a != self.name)}
    
    @property
    def shadowed_aliases(self):
        """read-only
        
        Returns a set of the shadowed aliases of the instance.
        """
        sha = {a for a in self._aliases if self.get(a,exception=False) is not self}
        if self.get(self.name) is not self:
            sha.add(self.name)
        return sha
    
    def _add_alias(self,alias):
        self._aliases.add(alias)
        
    def tostring(self,fmt=None,**kwds):
        """
        Returns a string containing the symbol for the instance in the format
        given by the keyword `fmt` which may be set to a string the values 
        'html', 'latex', 'ascii' or 'unicode'.  
        """
        if fmt is None or fmt == 'unicode':
            return self.symbol
        if fmt is None:
            return self.symbol
        if fmt == 'html':
            return self.html_symbol
        if fmt == 'latex':
            return self.latex_symbol
        if fmt == 'ascii':
            return self.ascii_symbol
        raise ValueError('format ' + str(fmt) + ' is not recognized')

    