# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:08:11 2026

@author: Parksh
"""

import numpy as np
from .fit import Fit
from .distributions import Distribution

class DistFit(Fit):
    def __init__(self,x,cdf,p0,**kwds):
        
        if isinstance(cdf,Distribution):
            c = type(cdf)
        elif isinstance(cdf,type) and issubclass(cdf,Distribution):
            c = cdf
        else:
            c = None
            
        if c is not None:
            self.p_names = c.p_names
            self.p_names_html = c.p_names_html
            self.p_names_latex = c.p_names_latex
            self.p_names_ascii = c.p_names_ascii
            self.cdf = lambda *args:c.cdf(*args)
            self.Dist = c
        else:
            self.cdf = cdf
            self.p_names = None
            self.p_names_html = None
            self.p_names_latex = None
            self.p_names_ascii = None
            self.Dist = None
            
        x = np.sort(x)
        super().__init__(x,p0=p0,**kwds)
        
    def f(self,x,*p,**kwds):
        z = len(x)*np.diff(self.cdf(x,*p,**kwds),prepend=0)
        return np.log(z*np.exp(z)) - 1
    
    def plot(self,ylabel=None,xlabel=None,title=None,hold=False,
             plot_points=None,show_fit=True,show_data=True,xmin=None,xmax=None,fit_options={},
             fig_options={},sub_plot_options={},hist_options={}):
        import matplotlib.pyplot as plt
        
        if plot_points is None:
            plot_points = self.plot_points
            
        fig = plt.figure(**fig_options)
        ax = fig.add_subplot(**sub_plot_options)
            
        if show_fit:
            if xmin is None:
                p1 = self.xf[0]
            else:
                p1 = xmin
            if xmax is None:
                p2 = self.xf[-1]
            else:
                p2 = xmax
            r = p2 - p1
            if xmin is None:
                st = p1 - self.over_plot*r
            else:
                st = p1
            if xmax is None:
                en = self.xf[-1] + self.over_plot*r
            else:
                en = p2
            fx = np.linspace(st,en,plot_points)
            dx = (en - st)/plot_points
            fy = np.diff(self.cdf(fx,*self.pf))/dx
            fx = fx[1:]
            ax.plot(fx,fy,**fit_options)
        
        if show_data:
            if not 'bins' in hist_options:
                hist_options['bins'] = 100
            if not 'density' in hist_options:
                hist_options['density'] = True
            if not 'histtype' in hist_options:
                hist_options['histtype'] ='stepfilled'
            ax.hist(self.xf,**hist_options)
            
        if ylabel is not None:
            ax.ylabel(ylabel)
        if xlabel is not None:
            ax.xlabel(xlabel)
        if title is not None:
            ax.title(title)
        if not hold:
            plt.show()
            
    def plot_cdf(self,ylabel=None,xlabel=None,title=None,hold=False,
                 plot_points=None,show_fit=True,show_data=True,xmin=None,
                 xmax=None,fit_options={},fit_format='k-',
                 fig_options={},sub_plot_options={},data_options={},data_format='ko'):
        import matplotlib.pyplot as plt
        
        if plot_points is None:
            plot_points = self.plot_points
            
        fig = plt.figure(**fig_options)
        ax = fig.add_subplot(**sub_plot_options)
            
        if show_fit:
            if xmin is None:
                p1 = self.xf[0]
            else:
                p1 = xmin
            if xmax is None:
                p2 = self.xf[-1]
            else:
                p2 = xmax
            r = p2 - p1
            if xmin is None:
                st = p1 - self.over_plot*r
            else:
                st = p1
            if xmax is None:
                en = self.xf[-1] + self.over_plot*r
            else:
                en = p2
            fx = np.linspace(st,en,plot_points)
            fy = self.cdf(fx,*self.pf)
            if fit_format is None:
                ax.plot(fx,fy,**fit_options)
            else:
                ax.plot(fx,fy,fit_format,**fit_options)
        
        if show_data:
            if 'ms' not in data_options and 'markersize' not in data_options:
                if self.count > 100:
                    data_options['ms'] = 1
                elif self.count > 30:
                    data_options['ms'] = 2
                elif self.count > 20:
                    data_options['ms'] = 3
            if 'ls' not in data_options and 'linestyle' not in data_options:
                data_options['ls'] = 'none'
                    
            if data_format is None:
                ax.plot(self.xf,(np.arange(self.count) + 1)/self.count,**data_options)
            else:
                ax.plot(self.xf,(np.arange(self.count) + 1)/self.count,data_format,**data_options)
        
        if ylabel is not None:
            ax.ylabel(ylabel)
        if xlabel is not None:
            ax.xlabel(xlabel)
        if title is not None:
            ax.title(title)
        if not hold:
            plt.show()
        