""" Calculate the matrix for removal terms for G group diffusion
"""


from __future__ import division
import numpy as np


def matrix(scattering, absorption, diffusion, width, nodes):
    """ matrix for left side of diffusion equation
    """
    mat = np.eye(nodes - 1)

    dr = width / (nodes - 1)
    ft = -diffusion / dr ** 2

    # Non-boundary terms
    for j in range(1, nodes - 2):
        mat[j][j-1] = ft * (1 - 1 / (2 * j))
        mat[j][j] = -ft + absorption + np.sum(scattering)
        mat[j][j+1] = ft * (1 + 1 / (2 * j))

    # Reflective boundary condition at node 0
    # mat[0][0] = -2 * ft + absorption + np.sum(scattering)
    mat[0][0] = 2. * ft + absorption
    # mat[0][1] = 2 * ft
    mat[0][1] = 2. * ft

    # flux at rightmost = 0
    mat[nodes - 2][nodes - 2] = -2 * ft + absorption + np.sum(scattering)
    mat[nodes - 2][nodes - 3] = ft * (1 - 1 / (2 * (nodes - 1)))

    # print(mat)
    # print(type(mat))
    # print(mat.shape)
    return mat
