# -*- coding: utf-8 -*-

# module prefixedunit

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
Classes to automatically generate prefixed units from a base unit.
"""

from .unit import Unit,Conversion

class PrefixedUnit(Unit):    
    """PrefixedUnit(name,symbol,conversion=None,short_name=None,
                    additional_names=None,additional_short_names=None,
                    add_symbol=False,html_symbol=None,latex_symbol=None,
                    ascii_symbol=None,description=None,order=-1,prefixes=None,
                    base_prefix=None,base_description=None)
    
    Creating an instance of this class not only creates the base unit but also
    adds units with with a set of prefixes to the unit library.
    
    If the `prefixes` keyword is `None` then units with the standard SI prefixes
    will be generated. 
    """
    prefix_definitions = {
                'yocto':[1e-24,'y',None,None,None],
                'zepto':[1e-21,'z',None,None,None],
                'atto':[1e-18,'a',None,None,None],
                'femto':[1e-15,'f',None,None,None],
                'pico':[1e-12,'p',None,None,None],
                'nano':[1e-9,'n',None,None,None],
                'micro':[1e-6,'\u03BC','&mu;',None,'u'],
                'milli':[0.001,'m',None,None,None],
                'centi':[0.01,'c',None,None,None],
                'deci':[0.1,'d',None,None,None],
                'deca':[10,'da',None,None,None],
                'hecto':[100,'h',None,None,None],
                'kilo':[1000,'k',None,None,None],
                'mega':[1000000,'M',None,None,None],
                'giga':[1e9,'G',None,None,None],
                'tera':[1e12,'T',None,None,None],
                'peta':[1e15,'P',None,None,None],
                'exa':[1e18,'E',None,None,None],
                'zetta':[1e21,'Z',None,None,None],
                'yotta':[1e24,'Y',None,None,None],
                }
    
    @staticmethod
    def _add_prefix(pdef,symbol,index):
        if pdef is None:
            return symbol
        
        pfx = pdef[index]
        if pfx is None:
            pfx = pdef[1]
        
        if symbol.startswith('\t\t'):
            return '\t\t' + pfx + symbol[2:]
        if symbol.startswith('\t'):
            return '\t' + pfx + symbol[2:]
        
        return pfx + symbol
    
    def __init__(self,name,symbol,conversion=None,short_name=None,
                 additional_names=None,additional_short_names=None,
                 add_symbol=False,html_symbol=None,latex_symbol=None,
                 ascii_symbol=None,linear=True,description=None,order=-1,
                 prefixes=None,base_prefix=None,base_description=None,**kwds):
        self.linear=linear
        if base_description is not None:
            self.description=base_description
        else:
            self.description=description
            
        if symbol.startswith('\t\t'):
            raise ValueError('the symbol may not start with more than one tab character')
        if html_symbol is None:
            html_symbol = symbol
        if latex_symbol is None:
            latex_symbol = symbol
        if ascii_symbol is None:
            ascii_symbol = symbol
        
        #reset _used_units incase we are shadowing any unit definitions
        Unit._used_units = {} 
        
        self._aliases = set()
        
        if short_name is None:
            if add_symbol:
                if ascii_symbol is not None:
                    short_name = ascii_symbol
                else:
                    short_name = symbol
            else:
                short_name = name

        self.parent = kwds.get('parent')
                    
        if self.parent is None:
            if prefixes is None:
                self.prefixes = list(type(self).prefix_definitions.keys())
            else:
                self.prefixes = list(prefixes)
            pfx = base_prefix
            if pfx is not None:
                pdef = type(self).prefix_definitions[pfx]
            self.conversion = conversion
            self.order = order
            
            pfxs = list(self.prefixes)
            if base_prefix is not None:
                pfxs.append(None)
            for p in pfxs:
                type(self)(name,symbol,
                            conversion=conversion,
                            short_name=short_name,
                            additional_names=additional_names,
                            additional_short_names=additional_short_names,
                            add_symbol=add_symbol,
                            html_symbol=html_symbol,
                            latex_symbol=latex_symbol,
                            ascii_symbol=ascii_symbol,
                            linear=linear,
                            description=description,
                            base_prefix=base_prefix,
                            order=order,
                            parent=self,
                            pfx=p,
                            **kwds)
        else:
            pfx = kwds['pfx']
            if pfx is None:
                factor = 1/type(self).prefix_definitions[base_prefix][0]
            else:
                pdef = type(self).prefix_definitions[pfx]
                if base_prefix is None:
                    factor = pdef[0]
                else:
                    factor = pdef[0]/type(self).prefix_definitions[base_prefix][0]
            self.conversion = Conversion(self.parent,factor)
            self.order = order - 0.01
            self.prefixes = self.parent.prefixes                
        
        if pfx is None:
            self.name = name
            self.short_name = short_name
            self.symbol = symbol
            pdef = None
        else:
            self.name = pfx + name
            if pdef[4] is None:
                self.short_name = pdef[1] + short_name
            else:
                self.short_name = pdef[4] + short_name
            self.symbol = pdef[1] + symbol
            
        self.html_symbol = PrefixedUnit._add_prefix(pdef,html_symbol,2)
        self.latex_symbol = Unit.format_latex(PrefixedUnit._add_prefix(pdef,latex_symbol,3))
        self.ascii_symbol = PrefixedUnit._add_prefix(pdef,ascii_symbol,4)      
                
        Unit._open_lib[self.name] = self
        Unit._open_lib[self.short_name] = self
        if self.short_name != self.name:
            self._aliases.add(self.short_name)
        
        if add_symbol:
            if self.symbol.strip() != self.name:
                Unit._open_lib[self.symbol.strip()] = self
                self._aliases.add(self.symbol.strip())
            if (self.ascii_symbol is not None and 
                self.ascii_symbol.strip() != self.name and 
                self.ascii_symbol != self.symbol):
                Unit._open_lib[self.ascii_symbol.strip()] = self
                self._aliases.add(self.ascii_symbol.strip())
        if additional_names is not None:
            for n in additional_names:
                if pfx is None:
                    nm = n
                else:
                    nm = pfx + n
                Unit._open_lib[nm]= self
                self._aliases.add(nm)
        if additional_short_names is not None:
            for n in additional_short_names:
                if pfx is None:
                    nm = n
                else:
                    nm = pdef[1] + n
                Unit._open_lib[nm]= self
                self._aliases.add(nm)
                
class BinaryPrefixedUnit(PrefixedUnit):
    """PrefixedUnit(name,symbol,conversion=None,short_name=None,
                    additional_names=None,additional_short_names=None,
                    add_symbol=False,html_symbol=None,latex_symbol=None,
                    ascii_symbol=None,description=None,order=-1,prefixes=None,
                    base_prefix=None,base_description=None)
    
    Creating an instance of this class not only creates the base unit but also
    adds units with with a set of prefixes to the unit library.
    
    If the `prefixes` keyword is `None` then units with the following prefixes will
    be generated:  kibi, mebi, gibi, tebi, pebi, exbi, zibi, yobi, kilo, mega, 
    giga, tera, peta, exa, zetta and yotta.
    """
    prefix_definitions = {
                'kilo':[1000,'k',None,None,None],
                'mega':[1e6,'M',None,None,None],
                'giga':[1e9,'G',None,None,None],
                'tera':[1e12,'T',None,None,None],
                'peta':[1e15,'P',None,None,None],
                'exa':[1e18,'E',None,None,None],
                'zetta':[1e21,'Z',None,None,None],
                'yotta':[1e24,'Y',None,None,None],
                'kibi':[1024,'Ki',None,None,None],
                'mebi':[1024**2,'Mi',None,None,None],
                'gibi':[1024**3,'Gi',None,None,None],
                'tebi':[1024**4,'Ti',None,None,None],
                'pebi':[1024**5,'Pi',None,None,None],
                'exbi':[1024**6,'Ei',None,None,None],
                'zibi':[1024**7,'Zi',None,None,None],
                'yobi':[1024**8,'Yi',None,None,None]
                }
        