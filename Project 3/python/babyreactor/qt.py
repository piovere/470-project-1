""" Quick dirty test script
"""


import numpy as np
from babyreactor.matrix import matrix


ABSORPTION = 0.15
SCATTERING = np.array([0.15, 0.15, 0.15, 0.15])
DIFFUSION = 9.
print(matrix(SCATTERING, ABSORPTION, DIFFUSION, 10, 5))
