""" Quick dirty test script
"""


import numpy as np
from babyreactor.flux import eigenfluxes
from babyreactor.matrix import matrix


ABSORPTION = 0.15
SCATTERING = np.array([0.15, 0.15, 0.15, 0.15])
DIFFUSION = 9.
mat = matrix(SCATTERING, ABSORPTION, DIFFUSION, 10, 5)

print(eigenfluxes(150., 20, np.zeros((19))))
