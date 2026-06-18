# -*- coding: utf-8 -*-

# module __init__

# Copyright (C) 2026 National Research Council Canada
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

from ._gummy import gummy,jummy
from ._unit import Unit,Conversion,one,unit
from .quantity import Quantity,ArrayQuantity
from ._ummy import ummy,immy
from .dof import DoF
from .budget import Budget
from .prefixedunit import PrefixedUnit,BinaryPrefixedUnit
from .nonlinearunit import NonlinearUnit,NonlinearConversion
from .offsetunit import OffsetUnit,OffsetConversion
from .logunit import LogUnit,LogConversion
from .functions import (sin,cos,tan,arcsin,arccos,arctan,arctan2,sinh,cosh,
                        tanh,arcsinh,arccosh,arctanh,exp,exp2,expm1,log2,
                        log10,log1p,logaddexp,logaddexp2,sqrt,crt,square,add,
                        subtract,negative,multiply,divide,true_divide,
                        floor_divide,reciprocal,power,absolute,mod,remainder,
                        divmod,modf,angle,real,imag,conj,around,rint,fix,floor,
                        ceil,trunc,heaviside,sign,sum,prod,cumsum,cumprod,
                        diff,ediff1d,gradient,cross)
from .misc import (correlation_matrix,covariance_matrix,
                   correlation_matrix_sim,covariance_matrix_sim,
                   simulate,clear_all_sim,covplot)
from ._mean import (autocorrelation,n_eff,wmean,mean,sigma_trim,delta_diff,
                   delta_diff_mean,delta_sum,delta_sum_mean,mean_datetime)
from .distributions import (Distribution,Convolution,MultivariateDistribution,
                            MultivariateElement,NormalDist,TDist)
from .miscdistributions import (MultiNormalDist,MultiTElement,MultiTDist,
                                UniformDist,GammaDist,LaplaceDist,TriangularDist,
                                ExponentialDist,PoissonDist,BinomialDist,
                                CurvlinearTrapDist,TrapezoidalDist,ArcSinDist,
                                LogNormalDist,WeibullDist,AveragedFrom)
from .exceptions import  (UnitError,IncompatibleUnitsError,UnitLibError,
                          UnitNotFoundError,CircularUnitConversionError,
                          UnitLibNotFoundError,ConstantNotFoundError,
                          NoSimulatedDataError,UnitWarning,GummyWarning,
                          FitWarning,BudgetWarning,UncertiantyPrecisionWarning)
from .printing import set_printer
from .unitutils import search_units,shadowed_units,convert,search_units_result
from .constant import (GummyConstant,JummyConstant,constant,search_constants,
                       shadowed_constants,search_constants_result)

from .fit import Fit
from .polyfit import PolyFit
from .miscfit import ExpFit,DoubleExpFit,OneOverTFit,SinFit
from .distfit import DistFit

from .mfraction import MFraction
from .abc import (UncertainValue,UncertainComplexValue,AbcQuantity,
                  AbcQuantityArray,AbcUnit)