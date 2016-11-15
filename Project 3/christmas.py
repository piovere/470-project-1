"""Tools to calculate flux in multi-region medium for G energy groups

Lots of awesome things going on here
"""

import numpy as np


FUEL_NODES = 10
REFLECTOR_NODES = 10
ENERGY_GROUPS = 4

class Material(object):
    """Material for a reactor
    """

    def __init__(self):
        self._scattering = 0
        self._absorption = 0
        self._fission = 0
        self._width = 0

    @property
    def scattering(self):
        """The scattering cross sections

        This should be a (g * g) matrix where the (i,j)th element
        is the cross section for scattering from energy group i to
        energy group j.
        """
        return self._scattering

    @scattering.setter
    def scattering(self, matrix):
        self._scattering = matrix

    @property
    def fission(self):
        """The fission cross sections

        This should be a length-g array where g is the number
        of energy groups.
        """
        return self._fission

    @fission.setter
    def fission(self, fission):
        self._fission = fission

    @property
    def absorption(self):
        """Absorption cross sections

        This should be a length-g array where g is the number
        of energy groups.
        """
        return self._absorption

    @absorption.setter
    def absorption(self, absorption):
        self._absorption = absorption

    @property
    def width(self):
        """Width

        The one-dimensional width of the material in the core.
        """
        return self._width

    @width.setter
    def width(self, width):
        self._width = width

    def matrix(self, nodes):
        """The matrix that operates on the flux in the material

        Will be of size (n * n) where n = nodes - 1
        """
        sections = nodes - 1
        return np.eye(sections)
