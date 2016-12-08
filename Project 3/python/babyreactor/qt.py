""" Quick dirty test script
"""


import numpy as np
import matplotlib.pyplot as plt
from babyreactor.flux import eigenfluxes
from babyreactor.matrix import matrix


ABSORPTION = 0.15
SCATTERING = np.array([0.15, 0.15, 0.15, 0.15])
DIFFUSION = 9.
# mat = matrix(SCATTERING, ABSORPTION, DIFFUSION, 10, 5)

nodes = 10
source = np.zeros((nodes - 1))
source[0] = 10.
data = eigenfluxes(150., nodes, source)
plt.plot(data['flux'])
plt.yscale('log')
plt.text(nodes - 1, data['flux'][-1], data['iterations'])
