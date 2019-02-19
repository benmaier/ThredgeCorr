import os 
import numpy as np
from math import log10

#=========== SIMULATION DETAILS ========
projectname = "ThredgeCorr"
basename = "degree_dist_3"

seed = -1                          
N_measurements = 1

rhos = np.array([0.0, 0.025, 0.05, 0.1, 0.25, 0.5])
ks = np.array([4, 64, 256, 1024])
ns = 200 * 8**np.arange(4,dtype=int)

external_parameters = [
                        ( 'k', ks ),
                        ( 'n', ns ),
                        ( 'rho', rhos ),
                      ]
internal_parameters = [
                      ]
standard_parameters = [
                      ]

only_save_times = False

