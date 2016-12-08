"""Calculate the eigenvalue (k) and eigenvector (flux) for a given matrix

"""


from __future__ import division
import numpy as np
import numpy.linalg as la
from .config import MIN_ERROR, MAX_ITERATIONS
from .matrix import matrix


def eigenfluxes(width, nodes, source, min_error=MIN_ERROR):
    """Calculate flux and k-value for given material
    """
    # Start with material matrix
    scattering = np.array([0., 0., 0., 0.])
    diffusion = 9.
    absorption = 0.1532
    mat = matrix(scattering, absorption, diffusion, width, nodes)

    # Seed flux and k values
    flux = np.ones(mat.shape[0])
    k = 1.

    # seed initial error values
    k_error = 100.
    flux_error = 100.

    num_iterations = 0

    # iteratively solve for flux
    try:
        while (k_error > MIN_ERROR or flux_error > MIN_ERROR) and num_iterations < MAX_ITERATIONS:
            num_iterations += 1

            if num_iterations % 10000 == 0:
                print(num_iterations)

            flux, flux_error, k, k_error = iterate(mat, flux, k, source)
    except ValueError as e:
        print('I\'m in the exception block')
        print(k_error)
        print(k)
        print(flux_error)
        raise e

    return {
        'k': k,
        'flux': flux,
        'k error': k_error,
        'flux error': flux_error,
        'iterations': num_iterations
    }


def iterate(mat, flux, k, source=None):
    """ One flux solver iteration
    """
    k_old = k
    flux_old = flux

    # npf = 2.43
    # fission = 0.1570 / npf
    npf = 0
    fission = 0

    b_vector = rightside(flux, fission, 0., npf, k, source)

    flux = la.solve(mat, b_vector)

    m_new = la.norm(flux)
    m_old = la.norm(flux_old)

    k = k_old * m_new / m_old

    flux_error = np.amax(abs(flux - flux_old) / flux)
    k_error = abs(k - k_old) / k

    return flux, flux_error, k, k_error


def rightside(flux, fission, scattering, neutrons_per_fission, k, sources=None):
    """

    :param flux: numpy.ndarray
    :param fission: float
    :param scattering: float
    :param neutrons_per_fission: float
    :param k: float
    param sources: numpy.ndarray
    """
    res = np.zeros_like(flux)
    if sources is None:
        sources = np.zeros_like(flux)
    fis = flux * fission * neutrons_per_fission / k
    scat = np.zeros_like(flux) * scattering

    res = sources + fis + scat
    # res[0] = res[0] / 2

    # print(res)
    return res
