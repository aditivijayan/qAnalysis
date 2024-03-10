import numpy as np
import math
import scipy.integrate as integrate
import h5py

def getdomain(file):
    infile = open(file)
    lines = infile.readlines()
    dom_range = np.zeros((2,3))
    ncell = np.zeros(3)
    dom_min = [0.0,0.0,0.0]
    dom_min[0] = float(lines[3].split()[2])
    dom_min[1] = float(lines[3].split()[3])
    dom_min[2] = float(lines[3].split()[4])
    
    dom_max = [0.0,0.0,0.0]
    dom_max[0] = float(lines[4].split()[2])
    dom_max[1] = float(lines[4].split()[3])
    dom_max[2] = float(lines[4].split()[4])

    ncell[0]=int(lines[15].split()[2])
    ncell[1]=int(lines[15].split()[3])
    ncell[2]=int(lines[15].split()[4])
    
    return dom_min, dom_max, ncell