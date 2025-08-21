# -*- coding: utf-8 -*-

# module exceptions

# Copyright (C) 2025 National Research Council Canada
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


class UnitError(ValueError):
    """
    Base class for Unit exceptions.
    """
    pass

class IncompatibleUnitsError(UnitError):
    """
    This exception is raised when an operation or conversion is attempted with
    Quantity instances that have units that are incompatible for that operation.
    """
    pass

class UnitLibError(UnitError):
    """
    This exception is raised when the UnitLibrary cannot parse a unit string.
    """
    pass

class UnitNotFoundError(UnitLibError):
    pass

class CircularUnitConversionError(UnitError):
    pass

class UnitLibNotFoundError(UnitLibError):
    pass

class ConstantNotFoundError(ValueError):
    pass

class NoSimulatedDataError(Exception):
    pass

class UnitWarning(Warning):
    pass

class GummyWarning(Warning):
    pass

class FitWarning(Warning):
    pass

class BudgetWarning(Warning):
    pass

class UncertiantyPrecisionWarning(Warning):
    pass