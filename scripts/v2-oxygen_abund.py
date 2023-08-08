#!/usr/bin/env python
# coding: utf-8

# In[35]:

##This version of the code uses file names as parsed arguments.s
import yt
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import math
import matplotlib.colors as mcolors
import matplotlib.cm
import h5py
import scipy
import scipy.interpolate as interpolate
import os, glob
import time as ostime
from multiprocessing import Process, Queue, Pool
import multiprocessing
from itertools import repeat
import datetime
import argparse

# In[2]:


plt.rcParams['font.size']=22
plt.rcParams['axes.linewidth']=2
plt.rcParams['xtick.major.size']=10
plt.rcParams['xtick.minor.size']=5
plt.rcParams['xtick.major.width']=2
plt.rcParams['xtick.minor.width']=1
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.major.size']=10
plt.rcParams['ytick.minor.size']=5
plt.rcParams['ytick.major.width']=2
plt.rcParams['ytick.minor.width']=1
plt.rcParams['ytick.direction']='in'

pc    = 3.018e18
kpc   = 1.e3 * pc 
Msun  = 2.e33
yr_to_sec =3.15e7
Myr   = 1.e6 * yr_to_sec
kmps  = 1.e5
boltzmann_constant_cgs    = 1.38e-16
hydrogen_mass_cgs = 1.67e-24
gamma = 5./3.

# In[3]:

start_time = ostime.time()
##############################################################
file = '/scratch/pawsey0807/aditi/quokka_devbranch/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackle = h5py.File(file)
array = grackle['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-6, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

bins = 100
egas_arr = np.logspace(-16., -5., bins)
nH_arr   = np.logspace(-6.0, 4.0, int(bins))
logrho_arr = np.log10(nH_arr[:-1])
logEgas_arr = np.log10(egas_arr[:-1])

T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))
i=0
for egas in egas_arr:
    j=0
    for nH in nH_arr:
        C = (gamma - 1.) * egas / (boltzmann_constant_cgs*nH)
        minT = C*np.amin(table)
        maxT = C*np.amax(table)
        def func(T):
            mu = interpolate.interp2d(table_temp, table_nH, table,\
                              kind='linear', copy=True, bounds_error=False, fill_value=None)
            return C*mu(T,nH)[0] - T

        T[i,j] = scipy.optimize.toms748(func, minT, maxT)
        j+=1
    i+=1

##############################################################

# In[4]:
build = 'Re2pc'
run = 'Start28Myr/'
lev = 0

# build = 'MetalGrad/'
# run = '4pcBase-8.0Eta/'
# lev = 1

home  = "/scratch/pawsey0807/aditi/quokka_devbranch/quokka/sims/"
parser = argparse.ArgumentParser(description='Optional app description')
# parser.add_argument('--build', type=str, help='Build of the run')
# parser.add_argument('--run', type=str, help='Run name')
parser.add_argument('--file', type=str, help='Filename to be analysed')
# parser.add_argument('--lev', type=str, help='Level to be analysed', default=0)
args = parser.parse_args()
filename = args.file
# lev   = args.lev
# build = args.build
# run   = args.run


# In[5]:


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


# In[6]:


class Data:
    fac = 1
    lev = 0
    file = ''
    dom_min = [0.0, 0.0, 0.0]
    def getData(file):
        ds   = yt.load(file)
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)
        density = np.array(data['gasDensity'])
        time = ds.current_time.to('Myr')
        Egas = np.array(data["gasInternalEnergy"])
        return density, Egas, time
    
    
if(lev==0):
    fac = 1
else:
    fac = 2 * lev
Data.fac = fac


infile = os.path.join(home, build, run,  'metal_uniform_512.in')
dom_min, dom_max, ncells = getdomain(infile)
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))

outhome = '/scratch/pawsey0807/aditi/quokka_devbranch/Analysis/data/'

output_dir = os.path.join(outhome, build,run,'Oxygen')
# if not os.path.isfile(output_dir):
#     os.makedirs(output_dir)

outputfile_name = os.path.join(outhome,build,run,'Oxygen/outflow_' + filename.split('plt')[1] + '.h5')

if((os.path.isfile(outputfile_name))):
    print('File already exists!', filename)
    quit()
        
dom_min, dom_max, ncells = getdomain(infile)

xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))

dx = (dom_max[0]- dom_min[0])/(Data.fac*int(ncells[0]))
dy = (dom_max[1]- dom_min[1])/(Data.fac*int(ncells[1]))
dz = (dom_max[2]- dom_min[2])/(Data.fac*int(ncells[2]))

Lx = (dom_max[0]- dom_min[0])
Ly = (dom_max[1]- dom_min[1])
Lz = (dom_max[2]- dom_min[2])

inputfile = os.path.join(home, build, run, filename)

ds   = yt.load(filename)

data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)
timestep = ds.current_time.to('Myr')

rho = np.array(data['gasDensity'])
rhoZ = np.array(data['scalar_2'])
pz = np.array(data['z-GasMomentum'])
vz = pz/rho

dA   = np.full(rho.shape, dx * dy)  
dVol   = np.full(rho.shape, dx * dy * dz)  

net_mass_outflow_rate   = rho * vz * dA
net_scalar_outflow_rate = rhoZ *  vz * dA

total_mass_outflow_rate   = np.sum(net_mass_outflow_rate  , axis=(0,1))  
total_scalar_outflow_rate = np.sum(net_scalar_outflow_rate, axis=(0,1))

total_scalar_value        = np.sum(rhoZ * dVol, axis=(0,1) )
total_mass                = np.sum(rho  * dVol, axis=(0,1))

hfo = h5py.File(outputfile_name, 'w')
        
hfo.create_dataset('TotalOutflowRate'            , data=total_mass_outflow_rate)
hfo.create_dataset('TotalScalarOutflowRate'            , data=total_scalar_outflow_rate)

hfo.create_dataset('TotalScalarValue'            , data=total_scalar_value)
hfo.create_dataset('TotalMass'            , data=total_mass)

hfo.create_dataset('Zrange'  , data=zrange)
hfo.create_dataset('Timestep', data=timestep)
hfo.close()

print("--------Written file------->>",filename)