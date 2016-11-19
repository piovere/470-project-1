"""Calculate the eigenvalue (k) and eigenvector (flux) for a given matrix

"""

from babyreactor.material import Material
import numpy as np
import numpy.linalg as la


MIN_ERROR = 0.001
MAX_ITERATIONS = 1e7

def eigenflux(material, width, nodes, source, min_error=MIN_ERROR, reflector=None):
    """Calculate flux and k-value for given material
    """
    # Start with material matrix
    a = material.matrix(nodes=nodes, width=width)

    # Seed flux and k values
    flux = np.ones(len(a[0]))
    k = 1.

    # seed initial error values
    k_error = 100.
    flux_error = 100.

    num_iterations = 0

    # iteratively solve for flux
    while (k_error > MIN_ERROR or flux_error > MIN_ERROR) and num_iterations < MAX_ITERATIONS:
        # keep the values from the last iteration
        flux_old = flux
        k_old = k

        flux = np.multiply(source, flux) / k

        flux = la.solve(a, flux)

        flux_norm = la.norm(flux)
        flux_old_norm = la.norm(flux_old)

        k = flux_norm / flux_old_norm

        k_error = abs(k - k_old) / k
        flux_error = np.sum(abs(flux - flux_old))

        # normalize flux for the next iteration
        flux = flux / la.norm(flux)

    return {
        'k': k,
        'flux': flux,
        'k error': k_error,
        'flux error': flux_error,
        'iterations': num_iterations
    }
