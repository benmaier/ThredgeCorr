import os 
import numpy as np
from math import log10

#=========== SIMULATION DETAILS ========
projectname = "ThredgeCorr"
basename = "small_degree_dist"

seed = -1
N_measurements = 1

rhos = np.array([0.0, 0.05, 0.15, 0.5])
ks = np.array([100])
ns = np.array([1e5])

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

#============== QUEUE ==================
queue = "SGE"
memory = "1G"
priority = 0

#============ CLUSTER SETTINGS ============
username = "bfmaier"
server = "groot0.biologie.hu-berlin.de"
useratserver = username + u'@' + server

shell = "/bin/bash"
pythonpath = "/opt/python36/bin/python3.6"
name = basename + "_meas_" + str(N_measurements)
serverpath = "/home/"+username +"/"+ projectname + "/" + name 
resultpath = serverpath + "/results"

#=======================================
localpath = os.path.join(os.getcwd(),"results_"+name)
n_local_cpus = 3

#========================
git_repos = [
                #( "/home/"+username+"/repo/", pythonpath + " setup.py install --user" )
            ]

if __name__ == "__main__":
    print(ns, ks, rhos)
