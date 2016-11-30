from babyreactor.matrix import matrix
import numpy as np


absorption = 0.15
scattering = np.array([0.15, 0.15, 0.15, 0.15])
diffusion = 9.
matrix(scattering, absorption, diffusion, 10, 5)