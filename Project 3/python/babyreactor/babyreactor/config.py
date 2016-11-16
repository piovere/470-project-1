"""Configuration values

In case you did not receive prompts, or need some sensible defaults
"""

import numpy as np


FUEL_NODES = 10
REFLECTOR_NODES = 10
ENERGY_GROUPS = 4

reflector_scattering = np.array(
    [
        [0.37045, 0.04152, 0.00001, 0.00000],
        [0.00000, 0.98285, 0.07459, 0.01371],
        [0.00000, 0.00000, 0.76110, 0.31856],
        [0.00000, 0.00000, 0.00085, 1.96607]
    ]
)

reflector_fission = np.array(
    [0., 0., 0., 0.]
)

reflector_absorption = np.array(
    [0.00051, 0.00354, 0.01581, 0.04637]
)

reflector_transport = np.array(
    [0.20608, 0.60215, 0.56830, 1.21110]
)

fuel_scattering = np.array(
    [
        []
    ]
)
