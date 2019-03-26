# -*- coding: utf-8 -*-

# module budget

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
This module defines the Budget object which displays a table listing the 
uncertainty contributions to a gummy.
"""

import numpy as np
from .gummy import gummy
from .unit import one
from .exceptions import BudgetWarning
from .printing import PrettyPrinter


import warnings
from collections import OrderedDict

class Budget(PrettyPrinter):
    """
    A class that facilitates the creation of uncertainty budget tables.
    
    To display the table use the `Budget.html` or `Budget.latex`
    methods  in a console or notebook that supports this type of output
    or the python  built-in function to get a unicode table.

    The `Budget.tohtml` and `Budget.tolatex` methods can be used
    to get strings with the html or latex code.

    The `Budget.df` property can be used to retrieve a pandas `DataFrame`
    with the table.  Also `Budget.df_str`, `Budget.df_html` and
    `Budget.df_latex` return DataFrames with formatted strings as entries
    rather than numerical values.

    Parameters
    ----------
    y:  `gummy`
        the dependant variable

    xlist:   array_like of `gummy`
        The independent variables.  Warnings will be generated if the
        gummys in this list over determine `y` (that is if not all
        variables in this list can be treated as independent variables)
        or under determine `y` (that is if some variables
        contributing to the uncertainty in `y` are missing).

    uunit:  `str` or Unit, optional
        Unit to use to express the uncertainties.  This useful if you
        wish to express all uncertainties  as relative uncertainty unit
        (e.g. %).

    k, p: `float`, optional
        k or p values for the expanded uncertainty; do not specify both `k`
        and `p`; if neither are specified the the `k` and `p` values of `y`
        are used

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
        this is `None`, then the expanded uncertainty is displayed if ``y.k != 1``.
        This can also be changed by setting the `show_expanded_u` attribute.

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
    
    default_css = """
.budget_table {
    display: table;
}
    
.budget_row {
    display: table-row;
}

.budget_header_row {
    display: table-row;
    font-weight: bold;
}

.budget_first_col_header_cell
{
    display: table-cell;
    text-align: center;
    border-bottom: solid;
    border-width: thin;
    padding-left: 5px;
    padding-right: 5px;
}

.budget_header_cell {
    display: table-cell;
    text-align: center;
    border-bottom: solid;
    border-width: thin;
    padding-left: 5px;
    padding-right: 5px;
    
}

.budget_first_col_x_cell {
    display: table-cell;
    text-align: center;
    padding-left: 5px;
    padding-right: 5px;
    min-width: 8ex;
}

.budget_x_cell {
    display: table-cell;
    text-align: center;
    padding-left: 5px;
    padding-right: 5px;
    min-width: 8ex;
}
"""
    
    @staticmethod
    def _formatsig(x,fmt=None):
        if x < 10**-9:
            return '0'
        if x >= 0.01:
            return '{:.2f}'.format(x)
        if x >= 0.001:
            return '{:.3f}'.format(x)
        return Budget._format_float(x,nsig=1,fmt=fmt)
        
    @staticmethod
    def _format_float(x,fmt=None,nsig=None):
        g = gummy(x,abs(x))
        g.style = 'xf'
        if nsig is not None:
            g.nsig = nsig
        return g.tostring(fmt)
    
    units_on_values = None
    show_s = True
    show_d = False
    show_c = True
    default_columns = None
    default_column_names = None
        
    def __init__(self,y,xlist,uunit=None,units_on_values=None,
                 sort=True,solidus=None,mulsep=None,slashaxis=None,columns=None,
                 column_names=None,xnames=None,yname=None,show_subtotals=True,
                 show_expanded_u=None,description=None,description_math_mode=False,
                 custom=None,custom_heading=None,custom_math_mode=False,
                 show_s=None,show_d=None,show_c=None,css=None,k=None,p=None,
                 sim=False):
        from pandas import DataFrame
        self._DataFrame = DataFrame
        
        if solidus is None:
            self.solidus = gummy.solidus
        else:
            self.solidus = solidus
            
        if mulsep is None:
            self.mulsep = gummy.mulsep
        else:
            self.mulsep = mulsep
            
        if slashaxis is None:
            self.slashaxis = gummy.slashaxis
        else:
            self.slashaxis = slashaxis
            
        if show_s is not None:
            self.show_s = show_s
        if show_d is not None:
            self.show_d = show_d
        if show_c is not None:
            self.show_c = show_c
        
        self.sim = sim
        
        got_dof = np.isfinite(y.dof) or any(np.isfinite(g.dof) for g in xlist)
        got_unit = (y.unit is not one) or any(g.unit is not one for g in xlist)
        got_uunit = (((y.uunit is not None) or any(g.uunit is not None for g in xlist))
                     and uunit is None)
        got_tag = any(g.utype is not None for g in xlist)
            
        self._y = y.copy(formatting=True)
        if k is not None or p is not None:
            if k is not None and p is not None:
                raise ValueError('k and p may not both be specified')
            if k is not None:
                self._y.k = k
            else:
                self._y.p = p
            
        if uunit is not None:
            self._y.uunit = uunit
        self._yu = y.copy(formatting=True)
        self._yu.k = 1
        
        self.nx = len(xlist)
        
        ind = list(range(len(xlist)))
        
        if xnames is not None and len(xnames) != len(xlist):
            raise ValueError('len(xnames) != len(xlist)')
        if xnames is not None:
            xnames = [str(nm).strip() for nm in xnames]
            hxnames = ['<span><i>' + str(nm).strip() + '</i></span>' 
                       if len(nm) <= 1 else str(nm).strip() for nm in xnames]
            lxnames = [str(nm).strip() if len(nm) <= 1 else self.latex_norm(nm) for nm in xnames]
            
        else:
            xnames = []
            hxnames = []
            lxnames = []
            for i,x in enumerate(xlist):
                if x.name is None:
                    xnames.append('x[' + str(i+1) + ']')
                    hxnames.append('<span><i>x</i><sub>' + str(i+1) + '</sub></span>')
                    lxnames.append('x_{' + str(i+1) + '}')
                else:
                    xnames.append(str(x.name).strip())
                    if not isinstance(x.name,str) or len(x.name) > 1:
                        hxnames.append(str(x.name).strip())
                        lxnames.append(type(self).latex_norm(str(x.name).strip()))
                    else:
                        hxnames.append('<span><i>' + str(x.name).strip() + '</i></span>')
                        lxnames.append(str(x.name).strip())
                        
            
        if yname is None:
            if y.name is None:
                yname = 'y'
            else:
                yname = y.name
        self._yname = str(yname).strip()
        if len(yname) <= 1:
            self._hyname = '<span><i>' + str(yname).strip() + '</i></span>'
            self._lyname = str(yname).strip()
        else:
            self._hyname = str(yname).strip()
            self._lyname = type(self).latex_norm(str(yname).strip())
            
        
        if units_on_values is not None:
            self.units_on_values = units_on_values
        
        self.show_subtotals = show_subtotals
        self.show_expanded_u = show_expanded_u
        self.css = css
        
        self.uunit = uunit
        
        if description is None:
             self.ydescription = None
             self.yedescription = None
             self.tdescription = None
             description = (self.nx+1)*['']
        else:
            if len(description) < self.nx + 1:
                raise ValueError('description must have a length of at least len(xlist) + 1')
            self.ydescription = description[0]
            if len(description) > self.nx + 1:
                self.yedescription = description[-1]
                if len(description) > self.nx + 1:
                    self.tdescription = description[(self.nx + 1):-1]
                else:
                    self.tdescription = None
            else:
                self.yedescription = None
                self.tdescription = None
        self.description_math_mode = description_math_mode
        if custom is None:
            self.ycustom = None
            self.yecustom = None
            self.tcustom = None
            custom = (self.nx+1)*['']
        else:
            if len(custom) < self.nx + 1:
                raise ValueError('custom must have a length of at least len(xlist) + 1')
            self.ycustom = custom[0]
            if len(custom) > self.nx + 1:
                self.yecustom = custom[-1]
                if len(custom) > self.nx + 1:
                    self.tcustom = custom[(self.nx + 1):-1]
                else:
                    self.tcustom = None
            else:
                self.yecustom = None
                self.tcustom = None
        self.custom_math_mode = custom_math_mode
        if custom_heading is None:
            self.custom_heading = ''
        else:
            self.custom_heading = custom_heading
            
        x = xlist
            
        x = [g.copy(formatting=True) for g in x]
        for g in x:
            g.k = 1
            if uunit is not None:
                g.uunit = uunit
            
        if self.sim:
            cm = gummy.correlation_matrix_sim(x)
            yu = y.usim
        else:
            cm = gummy.correlation_matrix(x)
            yu = y._u
        svd = np.linalg.svd(cm,compute_uv=False)
        for v in svd:
            if v < 1e-6:
                warnings.warn('the x variables are over determined; they cannot all be taken as independant variables',BudgetWarning,stacklevel=2)
                break
            
        b = [y.correlation(z) for z in x]
        s = np.linalg.lstsq(cm,b)[0]
        tu = 0
        
        if self.sim:
            d = [p*yu/q.usim for p,q in zip(s,x)]  
            for i in range(len(x)):
                for j in range(len(x)):
                    tu += d[i]*d[j]*x[i].correlation(x[j])*x[i].usim*x[j].usim
        else:
            d = [p*yu/q._u for p,q in zip(s,x)]  
            for i in range(len(x)):
                for j in range(len(x)):
                    tu += d[i]*d[j]*x[i].correlation(x[j])*x[i]._u*x[j]._u
        tu = np.sqrt(tu)
        
        if self.sim:
            mlim = 0.05
        else:
            mlim = 0.001
        if np.abs((yu - tu)/yu) > mlim:
            warnings.warn('a source of uncertainty seems to be missing from the x variables',BudgetWarning,stacklevel=2)
        
        if sim:
            fsigs = [np.abs(i*j.usim/yu) for i,j in zip(d,x)]
        else:
            fsigs = [np.abs(i*j._u/yu) for i,j in zip(d,x)]
        fsc = list(np.abs(d))
        
        sigs = []
        html_sigs = []
        latex_sigs = []
        html_d = []
        latex_d = []
        str_d = []
        sc = []
        html_sc = []
        latex_sc = []
        fdof = []
        dof = []
        html_dof = []
        latex_dof = []
     
        for i,xi in enumerate(x):
            sigs.append(Budget._formatsig(fsigs[i]))
            html_sigs.append(Budget._formatsig(fsigs[i],fmt='html'))
            latex_sigs.append(Budget._formatsig(fsigs[i],fmt='latex'))
            str_d.append(Budget._format_float(d[i]))
            html_d.append(Budget._format_float(d[i],fmt='html'))
            latex_d.append(Budget._format_float(d[i],fmt='latex'))
            sc.append(Budget._format_float(fsc[i]))
            html_sc.append(Budget._format_float(fsc[i],fmt='html'))
            latex_sc.append(Budget._format_float(fsc[i],fmt='latex'))
            fdof.append(x[i]._dof)
            dof.append(gummy._dof_to_str(x[i].dof))
            html_dof.append(gummy._dof_to_str(x[i].dof,fmt='html'))
            latex_dof.append(gummy._dof_to_str(x[i].dof,fmt='latex'))
        
        a = np.array([x,ind,xnames,hxnames,lxnames,dof,html_dof,latex_dof,fdof,
                      sigs,html_sigs,latex_sigs,fsigs,str_d,html_d,latex_d,d,sc,
                      html_sc,latex_sc,fsc,description[1:(self.nx + 1)],custom[1:(self.nx + 1)]]).T

        self._dfx = DataFrame(a,columns=['x','i','xnames','hxnames','lxnames',
                                         'dof','html_dof','latex_dof','fdof',
                                          's','html_s','latex_s','fs','d',
                                          'html_d','latex_d','fd','sc','html_sc',
                                          'latex_sc','fsc','description','custom'])
        if sort:
            self._dfx.sort_values(by=['fs'],inplace=True,ascending=False)
            
        stags = set([])
        for r in xlist:
            if r.utype is not None:
                stags.add(r.utype)
        self.type = OrderedDict()
        if stags:
            lst = list(stags)
            lst.sort()
            for t in lst:
                u = self._y.ufrom(t,sim=self.sim)
                g = gummy(self._y.x,u=u,unit=self._y.unit)
                if self._y.uunit is not None:
                    g.uunit = self._y.uunit
                self.type[t] = (u,self._y.doffrom(t),g)               

        if self.units_on_values is None:
           self.units_on_values = got_uunit
            
        if columns is None:
            if self.default_columns is not None:
                self.columns = self.default_columns
            else:
                cols = ['name']
                if self.ydescription is not None:
                    cols.append('description')
                if got_unit and not self.units_on_values:
                    cols.append('unit')
                cols.append('value')
                cols.append('u')
                if got_dof:
                    cols.append('dof')
                if got_tag:
                    cols.append('type')
                if self.show_d:
                    cols.append('d')
                if self.show_c:
                    cols.append('c')
                if self.show_s:
                    cols.append('s')
                if self.ycustom is not None:
                    cols.append('custom')
            
            self.columns = cols
        else:
            self.columns = columns
            
        if column_names is None and self.default_column_names is not None:
            self.column_names = self.default_column_names
        else:
            self.column_names = column_names
        
    @property
    def columns(self):
        """
        `list` of `str` or `None`

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
        """
        return self._columns
    @columns.setter
    def columns(self,cols):
        self._columns = []
        for c in cols:
            c = c.strip().lower()
            if c == 'name' or c == 'names' or c == 'component' or c == 'components':
                self._columns.append('component')
            elif c == 'unit' or c == 'units':
                self._columns.append('unit')
            elif c == 'value' or c == 'values':
                self._columns.append('value')
            elif c == 'u' or c == 'uncertainty' or c == 'uncertainties':
                self._columns.append('u')
            elif c == 'dof':
                self._columns.append('dof')
            elif c == 'type' or c == 'tag':
                self._columns.append('type')
            elif c == 's' or c == 'significance' or c == 'signifincances':
                self._columns.append('s')
            elif c == 'd' or c == 'derivative' or c == 'derivatives' or c == 'partial' or c == 'partials':
                self._columns.append('d')
            elif c == 'c' or c == 'sensitivity coefficient' or c == 'sensitivity coefficients':
                self._columns.append('c')
            elif c == 'description':
                self._columns.append('description')
            elif c == 'custom':
                self._columns.append('custom')
            else:
                raise ValueError('column name ' + str(c) + ' is unrecognized')
                
    @property
    def column_names(self):
        """
       `dict`

        Names to display as column headers.  The dictionary should use as keys
        any of the column names listed above in the columns parameter
        description and as values the desired heading for this column.
        """
        return self._column_names
    @column_names.setter
    def column_names(self,names):
        self._column_names = {}
        if names is None:
            return
        for c,v in names.items():
            c = c.strip().lower()
            v = v.strip()
            if c == 'name' or c == 'names' or c == 'component' or c == 'components':
                self._column_names['component'] = v
            elif c == 'unit' or c == 'units':
                self._column_names['unit'] = v
            elif c == 'value' or c == 'values':
                self._column_names['value'] = v
            elif c == 'u' or c == 'uncertainty' or c == 'uncertainties':
                self._column_names['u'] = v
            elif c == 'dof':
                self._column_names['dof'] = v
            elif c == 'type' or c == 'tag':
                self._column_names['type'] = v
            elif c == 's' or c == 'significance' or c == 'significances':
                self._column_names['s'] = v
            elif c == 'd' or c == 'derivative' or c == 'derivatives' or c == 'partial' or c == 'partials':
                self._column_names['d'] = v
            elif c == 'c' or c == 'sensitivity coefficient' or c == 'sensitivity coefficients':
                self._column_names['c'] = v
            elif c == 'description':
                self._columns.names['description'] = v
            elif c == 'custom':
                self._columns.names['custom'] = v
            else:
                raise ValueError('original column name ' + str(c) + ' is unrecognized')
                
    @property
    def k(self):
        """
        Gets or sets the `k` value for the expanded uncertainty of `y`.
        """
        return self._y.k
    @k.setter
    def k(self,value):
        self._y.k = value
        
    @property
    def p(self):
        """
        Gets or sets the `p` value for the expanded uncertainty of `y`.
        """
        return self._y.p
    @p.setter
    def p(self,value):
        self._y.p = value
            
    @property
    def df(self):
        """
        read-only

        Returns a Panda's DataFrame with the the budget table.
        """
        self._update()
        return self._df
        
    @property
    def df_str(self):
        """
        read-only

        Returns a Panda's DataFrame with the the budget table, with
        all values displayed as strings.
        """
        self._update()
        return self._df_str
        
    @property
    def df_html(self):
        """
        read-only

        Returns a Panda's DataFrame with the the budget table, with
        entries displayed using HTML.
        """
        self._update()
        return self._df_html
        
    @property
    def df_latex(self):
        """
        read-only

        Returns a Panda's DataFrame with the the budget table, with all
        entries displayed using LaTeX.
        """
        self._update()
        return self._df_latex
        
    def _update(self):        
        show_exu = self.show_expanded_u or (self.show_expanded_u is None and self._y.k != 1)
       
        self._df = self._DataFrame()
        self._df_str = self._DataFrame()
        self._df_html = self._DataFrame()
        self._df_latex = self._DataFrame()
                        
        for c in self._columns:
            if c == 'component':
                nms = list(self._dfx['xnames'])
                hnms = list(self._dfx['hxnames'])
                lnms = list(self._dfx['lxnames'] )
                if self.show_subtotals:
                    for t in self.type.keys():
                        nms.append('Combined type ' + str(t))
                        hnms.append('<span><i>u<sub>c</sub></i> type ' + str(t))
                        lnms.append('u_{c} \\text{ type ' + str(t) + '}')
                nms.append(self._yname)
                hnms.append(self._hyname)
                lnms.append(self._lyname)
                if show_exu:
                    if self.sim and not self._y._set_k:
                        k = self._y.ksim
                    else:
                        k = self._y.k
                    nms.append('Uc at k = ' + gummy._k_to_str(k))
                    hnms.append('<span><i>U<sub>c<sub></i> at <i>k</i> = ' + gummy._k_to_str(k) + '</span>')
                    lnms.append('U_c ' + type(self).latex_norm(' at ') + ' k = ' + gummy._k_to_str(k))
                    
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = 'Component'
                    cnml = type(self).latex_norm('Component')
                else:
                    cnml = cnm
                    
                self._df[cnm] = nms
                self._df_str[cnm] = nms
                self._df_html[cnm] = hnms
                self._df_latex[cnml] = lnms
                
            if c == 'unit':
                sym = []
                hsym = []
                lsym = []
                for r in self._dfx['x']:
                    if isinstance(r,gummy):
                        sym.append(r.unit.tostring(solidus=self.solidus,mulsep=self.mulsep))
                        hsym.append(r.unit.tostring(fmt='html',solidus=self.solidus,mulsep=self.mulsep))
                        lsym.append(r.unit.tostring(fmt='latex',solidus=self.solidus,mulsep=self.mulsep))
                    else:
                        sym.append('')
                        hsym.append('')
                        lsym.append('')
                        
                if isinstance(self._y,gummy):
                    yun = self._y.unit.tostring(solidus=self.solidus,mulsep=self.mulsep)
                    hyun = self._y.unit.tostring(fmt='html',solidus=self.solidus,mulsep=self.mulsep)
                    lyun = self._y.unit.tostring(fmt='latex',solidus=self.solidus,mulsep=self.mulsep)
                else:
                    yun = ''
                    hyun = ''
                    lyun = ''
                
                n = 1
                if self.show_subtotals:
                    n += len(self.type.values())
                            
                if show_exu:
                    n += 1

                sym += n*[yun]
                hsym += n*[hyun]
                lsym += n*[lyun]
                    
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = 'Unit'
                    cnml = type(self).latex_norm('Unit')
                self._df[cnm] = sym
                self._df_str[cnm] = sym
                self._df_html[cnm] = hsym
                self._df_latex[cnml] = lsym
                
            if c == 'value':
                if self.units_on_values:
                    if self.sim:
                        xst = 'xsim'
                    else:
                        xst = 'x'
                else:
                    if self.sim:
                        xst = 'xfsim'
                    else:
                        xst = 'xf'
                    
                fvalues = [g.x for g in self._dfx['x']]
                values = [g.tostring(style=xst,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                hvalues = [g.tostring(fmt='html',style=xst,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                lvalues = [g.tostring(fmt='latex',style=xst,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                
                if self.show_subtotals:
                    for t in self.type.values():
                        fvalues.append(None)
                        values.append('')
                        hvalues.append('')
                        lvalues.append('')
                        
                if self.sim:
                    fvalues.append(self._y.xsim)
                else:
                    fvalues.append(self._y.x)
                values.append(self._y.tostring(style=xst))
                hvalues.append(self._y.tostring(fmt='html',style=xst))
                lvalues.append(self._y.tostring(fmt='latex',style=xst))
                
                if show_exu:
                    fvalues.append(None)
                    values.append('')
                    hvalues.append('')
                    lvalues.append('')
                    
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = 'Value'
                    cnml = type(self).latex_norm('Value')
                self._df[cnm] = fvalues
                self._df_str[cnm] = values
                self._df_html[cnm] = hvalues
                self._df_latex[cnml] = lvalues
                
            if c == 'u':
                if self.units_on_values:
                    ust_g = 'u'
                    if self.sim:
                        ust = 'usim'
                    else:
                        ust = 'u'
                else:
                    ust_g = 'uf'
                    if self.sim:
                        ust = 'ufsim'
                    else:
                        ust = 'uf'
                    
                fu = [g.u for g in self._dfx['x']]
                u = [g.tostring(style=ust,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                hu = [g.tostring(fmt='html',style=ust,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                lu = [g.tostring(fmt='latex',style=ust,solidus=self.solidus,mulsep=self.mulsep) for g in self._dfx['x']]
                
                if self.show_subtotals:
                    for t in self.type.values():
                        fu.append(t[2].U)
                        u.append(t[2].tostring(style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                        hu.append(t[2].tostring(fmt='html',style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                        lu.append(t[2].tostring(fmt='latex',style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                            
                fu.append(self._yu.U)
                u.append(self._yu.tostring(style=ust,solidus=self.solidus,mulsep=self.mulsep))
                hu.append(self._yu.tostring(fmt='html',style=ust,solidus=self.solidus,mulsep=self.mulsep))
                lu.append(self._yu.tostring(fmt='latex',style=ust,solidus=self.solidus,mulsep=self.mulsep))
                if show_exu:
                    if self.sim:
                        g = gummy(self._y.x,self._y.ksim*self._y.usim,unit=self._y.unit)
                    else:
                        g = self._y
                    fu.append(g.U)
                    u.append(g.tostring(style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                    hu.append(g.tostring(fmt='html',style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                    lu.append(g.tostring(fmt='latex',style=ust_g,solidus=self.solidus,mulsep=self.mulsep))
                    
                cnm = cnmh = self.column_names.get(c)
                if cnm is None:
                    if self.uunit is None or self.units_on_values:
                        cnm = cnml = 'u'
                        cnmh = '<span><i>u</i></span>'
                    else:
                        if self.slashaxis:
                            cnm = 'u\u2009/\u2009' + self._y.uunit.tostring(solidus=self.solidus,mulsep=self.mulsep).strip()
                            cnml = 'u\\,/\\,' + self._y.uunit.tostring(fmt='latex',solidus=self.solidus,mulsep=self.mulsep).strip()
                            cnmh = '<span><i>u</i>&thinsp;/&thinsp;' + self._y.uunit.tostring(fmt='html',solidus=self.solidus,mulsep=self.mulsep).strip() + '</span>'
                        else:
                            cnm = 'u\u2009(' + self._y.uunit.tostring(solidus=self.solidus,mulsep=self.mulsep).strip() + ')'
                            cnml = 'u\\,(' + self._y.uunit.tostring(fmt='latex',solidus=self.solidus,mulsep=self.mulsep).strip() + ')'
                            cnmh = '<span><i>u</i>&thinsp;(' + self._y.uunit.tostring(fmt='html',solidus=self.solidus,mulsep=self.mulsep).strip() + ')</span>'
                self._df[cnm] = fu
                self._df_str[cnm] = u
                self._df_html[cnmh] = hu
                self._df_latex[cnml] = lu
                
            if c == 'dof':                
                fdof = list(self._dfx['fdof'])
                dof = list(self._dfx['dof'])
                hdof = list(self._dfx['html_dof'])
                ldof = list(self._dfx['latex_dof'] )
                if self.show_subtotals:
                    for t in self.type.values():
                        fdof.append(t[1])
                        dof.append(gummy._dof_to_str(t[1]))
                        hdof.append(gummy._dof_to_str(t[1],'html'))
                        ldof.append(gummy._dof_to_str(t[1],'latex'))
                fdof.append(self._y.dof)
                dof.append(gummy._dof_to_str(self._y.dof))
                hdof.append(gummy._dof_to_str(self._y.dof,fmt='html'))
                ldof.append(gummy._dof_to_str(self._y.dof,fmt='latex'))
                if show_exu:
                    fdof.append(None)
                    dof.append('')
                    hdof.append('')
                    ldof.append('')
                    
                cnm = self.column_names.get(c)
                if cnm is None:
                    self._df['DoF'] = fdof
                    self._df_str['DoF'] = dof
                    self._df_html['<span><i>&nu;<sub>eff</sub></i></span>'] = hdof
                    self._df_latex[r'\nu_{eff}'] = ldof
                else:
                    self._df[cnm] = fdof
                    self._df_str[cnm] = dof
                    self._df_html[cnm] = hdof
                    self._df_latex[cnm] = ldof
                
            if c == 'type':
                tags = []
                for r in self._dfx['x']:
                    if r.utype is not None:
                        tags.append(str(r.utype))
                    else:
                        tags.append('')
                        
                for t in self.type.keys():
                    if self.show_subtotals:
                        tags.append(str(t))                
                
                tags.append('')
                if show_exu:
                    tags.append('')
            
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = 'Type'
                    cnml = type(self).latex_norm('Type')
                self._df[cnm] = tags
                self._df_str[cnm] = tags
                self._df_html[cnm] = tags
                self._df_latex[cnml] = [type(self).latex_norm(t) for t in tags]
                
            if c == 's':                
                fs = list(self._dfx['fs'])
                s = list(self._dfx['s'])
                hs = list(self._dfx['html_s'])
                ls = list(self._dfx['latex_s'])
                if self.show_subtotals:
                    for v in self.type.values():
                        if self.sim:
                            yu = self._y.usim
                        else:
                            yu = self._y._u
                        fs.append(v[0]/yu)
                        s.append(Budget._formatsig(v[0]/yu))
                        hs.append(Budget._formatsig(v[0]/yu,'html'))
                        ls.append(Budget._formatsig(v[0]/yu,'latex'))
                fs.append(None)
                s.append('')
                hs.append('')
                ls.append('')
                if show_exu:
                    fs.append(None)
                    s.append('')
                    hs.append('')
                    ls.append('')
                    
                cnm = self.column_names.get(c)
                if cnm is None:
                    self._df['s'] = fs
                    self._df_str['s'] = s
                    self._df_html['<span><i>s</i></span>'] = hs
                    self._df_latex['s'] = ls
                else:
                    self._df[cnm] = fs
                    self._df_str[cnm] = s
                    self._df_html[cnm] = hs
                    self._df_latex[cnm] = ls
                
            if c == 'd':                
                fd = list(self._dfx['fd'])
                d = list(self._dfx['d'])
                hd = list(self._dfx['html_d'])
                ld = list(self._dfx['latex_d'] )
                if self.show_subtotals:
                    for v in self.type.values():
                        fd.append(None)
                        d.append('')
                        hd.append('')
                        ld.append('')
                fd.append(None)
                d.append('')
                hd.append('')
                ld.append('')
                if show_exu:
                    fd.append(None)
                    d.append('')
                    hd.append('')
                    ld.append('')
                    
                cnm = self.column_names.get(c)
                if cnm is None:
                    self._df['dy/dx'] = fd
                    self._df_str['dy/dx'] = d
                    self._df_html['<span>&part;<i>y</i>/&part;<i>x</i></span>'] = hd
                    self._df_latex[r'\frac{\partial y}{\partial x}'] = ld
                else:
                    self._df[cnm] = fd
                    self._df_str[cnm] = d
                    self._df_html[cnm] = hd
                    self._df_latex[cnm] = ld
                
            if c == 'c':                
                fc = list(self._dfx['fsc'])
                cc = list(self._dfx['sc'])
                hc = list(self._dfx['html_sc'])
                lc = list(self._dfx['latex_sc'])
                if self.show_subtotals:
                    for v in self.type.values():
                        fc.append(None)
                        cc.append('')
                        hc.append('')
                        lc.append('')
                fc.append(None)
                cc.append('')
                hc.append('')
                lc.append('')
                if show_exu:
                    fc.append(None)
                    cc.append('')
                    hc.append('')
                    lc.append('')
                    
                cnm = self.column_names.get(c)
                if cnm is None:
                    self._df['|dy/dx|'] = fc
                    self._df_str['|dy/dx|'] = cc
                    self._df_html['|<span>&part;<i>y</i>/&part;<i>x</i></span>|'] = hc
                    self._df_latex[r'\left\lvert\frac{\partial y}{\partial x} \right\rvert'] = lc
                else:
                    self._df[cnm] = fc
                    self._df_str[cnm] = cc
                    self._df_html[cnm] = hc
                    self._df_latex[cnm] = lc
                    
            if c == 'description':
                d = list(self._dfx['description'])
                if self.show_subtotals:
                    if (self.tdescription is not None and 
                        len(self.tdescription) == len(self.type.values())):
                        d += self.tdescription
                    else:
                        d += len(self.type.values())*['']
                d += [self.ydescription]
                if show_exu:
                    if self.yedescription is not None:
                        d += [self.yedescription]
                    else:
                        d += ['']
                
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = 'Description'
                    cnml = type(self).latex_norm('Description')
                self._df[cnm] = d
                self._df_str[cnm] = d
                self._df_html[cnm] = d
                if not self.description_math_mode:
                    self._df_latex[cnml] = [type(self).latex_norm(i) if i != '' 
                                            else '' for i in d]
                else:
                    self._df_latex[cnml] = d
                    
            if c == 'custom':
                d = list(self._dfx['custom'])
                if self.show_subtotals:
                    if (self.tcustom is not None and 
                        len(self.tcustom) == len(self.type.values())):
                        d += self.tcustom
                    else:
                        d += len(self.type.values())*['']
                d += [self.ycustom]
                if show_exu:
                    if self.yecustom is not None:
                        d += [self.yecustom]
                    else:
                        d += ['']
                
                cnm = cnml = self.column_names.get(c)
                if cnm is None:
                    cnm = self.custom_heading
                    if self.custom_heading == '':
                        cnml = ''
                    else:
                        cnml = type(self).latex_norm(self.custom_heading)
                self._df[cnm] = d
                self._df_str[cnm] = d
                self._df_html[cnm] = d
                if not self.custom_math_mode:
                    self._df_latex[cnml] = [type(self).latex_norm(i) if i != '' 
                                            else '' for i in d]
                else:
                    self._df_latex[cnml] = d
                        
        
    def tostring(self,fmt='unicode'):
        """
        Returns a string representation of the budget table

        Parameters
        ----------
        fmt: {'unicode','html','latex','ascii'}, optional
            encoding for the output.  The default is 'unicode'.
        """
        fmt = fmt.strip().lower()
        
        if fmt == 'unicode':
            return self.df_str.to_string(index=False)
        if fmt == 'html':
            return self._tohtml()
        if fmt == 'latex':
            return self._tolatex()
        if fmt == 'ascii':
            raise ValueError('fmt ascii is not available for Budget')
        raise ValueError('fmt ' + str(fmt) + ' is not valid')
        
    def _tolatex(self):
        self._update()
        nc = len(self._df_latex.columns)
        if self.show_subtotals:
            nt = len(self.type)
        else:
            nt = 0
        if self.nx == 0 or nc == 0:
            return ''
        ndf = len(self._df_latex.index)
        
        txt = r'\begin{array}{ '
        for i in range(nc):
            txt += 'c '
        txt += '}\n'
        
        aa = False
        for c in self._df_latex.columns:
            if aa:
                txt += ' & '
            else:
                aa = True
            txt += c
        txt += r' \\'
        txt += '\n'
        txt += r'\hline'
        txt += '\n'
       
        for i in range(ndf):
            aa = False
            for j in range(0,nc):
                v = self._df_latex.iloc[i,j]
                v = v.strip()
                if aa:
                    txt += ' & '
                else:
                    aa = True
                if v != '':
                    txt += v
            txt += r' \\'
            txt += '\n'
            if i == self.nx - 1 or i == self.nx+nt-1 or i == ndf-2:
                txt += r'\hline'
                txt += '\n'
        txt += '\end{array}'
        return txt

    def _tohtml(self):
        self._update()
        nx = self.nx
        nc = len(self._df_html.columns)
        if self.show_subtotals:
            nt = len(self.type)
        else:
            nt = 0
        if nx == 0 or nc == 0:
            return ''
        ndf = len(self._df_html.index)
            
        txt = '<div>\n<style>'
        if self.css is not None:
            txt += self.css
        else:
            txt += self.default_css
        txt += '</style>\n<div class="budget_table">\n    <div class="budget_header_row">\n'
        for c in self._df_html.columns:
            txt += '        <div class = "budget_header_cell">' + c + '</div>\n'
        txt += '    </div>\n'
      
        for i in range(ndf):
            txt += '    <div class="budget_row">\n'
            for j in range(0,nc):
                txt += '        <div class="'
                if i < nx - 1:
                    if j == 0:
                        txt += 'budget_first_col_x_cell'
                    else:
                        txt += 'budget_x_cell'
                elif i == nx - 1:
                    if j == 0:
                        txt += 'budget_first_col_header_cell'
                    else:
                        txt += 'budget_header_cell'
                elif i < (nx+nt-1):
                    if j == 0:
                        txt += 'budget_first_col_x_cell'
                    else:
                        txt +=  'budget_x_cell'
                elif i == (nx+nt-1):
                    if j == 0:
                        txt += 'budget_first_col_header_cell'
                    else:
                        txt += 'budget_header_cell'
                elif i == ndf-1:
                    if j == 0:
                        txt += 'budget_first_col_x_cell'
                    else:
                        txt += 'budget_x_cell'
                else:
                    if j == 0:
                        txt += 'budget_first_col_header_cell'
                    else:
                        txt += 'budget_header_cell'
                v = self._df_html.iloc[i,j]
                if v == '':
                    txt += '">  </div>\n'
                else:
                    txt += '">' + self._df_html.iloc[i,j] + '</div>\n'
            txt += '    </div>\n'
            
        txt += '</div>\n</div>'
        return txt