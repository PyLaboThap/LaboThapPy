
import numpy as np
from scipy.interpolate import interp1d

def cordier_line():

    spec_speed = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                            2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30])

    spec_diam = np.array([9.66, 6.98, 5.65, 4.7, 4.05, 3.57, 3.25,
                          3.01, 1.89, 1.57, 1.40, 1.26, 1.19, 1.11,
                          1.04, 0.97, 0.93, 0.64, 0.48])

    line = interp1d(spec_speed,spec_diam)

    return line
