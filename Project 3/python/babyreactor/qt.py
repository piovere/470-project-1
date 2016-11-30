from babyreactor.matrix import matrix
import numpy as np


ABSORPTION = 0.15
SCATTERING = np.array([0.15, 0.15, 0.15, 0.15])
DIFFUSION = 9.
matrix(SCATTERING, ABSORPTION, DIFFUSION, 10, 5)
