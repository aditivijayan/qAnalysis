import yt
import matplotlib.pyplot as plt
import numpy as np
import scipy
import matplotlib.colors as mcolors
import matplotlib.cm
import glob
import h5py
import os
import time as ostime
import scipy.interpolate as interpolate
from multiprocessing import Process, Queue
import multiprocessing
import matplotlib.cm     as cm
from scipy import optimize
from scipy import stats
import argparse

pc    = 3.018e18
kpc   = 1.e3 * pc 
Msun  = 2.e33
yr_to_sec =3.15e7
Myr   = 1.e6 * yr_to_sec
kmps  = 1.e5
boltzmann_constant_cgs    = 1.38e-16
hydrogen_mass_cgs = 1.67e-24
gamma = 5./3.

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

# build = 'MetalGrad/'
# run = '4pcBase-4.0Eta/'
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
outputfile_name = os.path.join(outhome,build,run,'MetAverage/met_phase_' + filename.split('plt')[1] + '.h5')

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
egas0 = np.array(data['gasInternalEnergy'])
dA   = np.full(rho.shape, dx * dy)  
dVol   = np.full(rho.shape, dx * dy * dz)  

rho0 = rho/hydrogen_mass_cgs


logrho_arr = np.log10(nH_arr[:-1])
logrho = np.log10(rho0)
delta_rho = logrho_arr[1] - logrho_arr[0]
idxrho = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

logEgas_arr = np.log10(egas_arr[:-1])
logEgas = np.log10(egas0)
delta_egas = logEgas_arr[1] - logEgas_arr[0]
idxegas = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

Temp = (1.-wgt_rho) * (1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
            wgt_rho  *    wgt_egas  * T[tuple(idxegas+1), tuple(idxrho+1)] +\
        (1. -wgt_rho) *    wgt_egas  * T[tuple(idxegas+1), tuple(idxrho)]   +\
            wgt_rho  * (1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  


neutral = (Temp<=7105.0)
wim     = (Temp>7105.0) & (Temp<=2.e4)
whim    = (Temp>2.e4)   & (Temp<5.e5)
hot     = (Temp>5.e5)

factor = 1.

met_neutral = np.sum(np.ma.array(rhoZ*dVol*factor, mask=~neutral), axis=(0,1)).data
met_wim     = np.sum(np.ma.array(rhoZ*dVol*factor, mask=~wim), axis=(0,1)).data
met_whim    = np.sum(np.ma.array(rhoZ*dVol*factor, mask=~whim), axis=(0,1)).data
met_hot     = np.sum(np.ma.array(rhoZ*dVol*factor, mask=~hot), axis=(0,1)).data

mass_neutral = np.sum(np.ma.array(rho*dVol*factor, mask=~neutral), axis=(0,1)).data
mass_wim     = np.sum(np.ma.array(rho*dVol*factor, mask=~wim), axis=(0,1)).data
mass_whim    = np.sum(np.ma.array(rho*dVol*factor, mask=~whim), axis=(0,1)).data
mass_hot     = np.sum(np.ma.array(rho*dVol*factor, mask=~hot), axis=(0,1)).data

outhome = '/scratch/pawsey0807/aditi/quokka_devbranch/Analysis/data/'
outputfile_name = os.path.join(outhome,build,run, 'MetAverage/met_phase_'+ filename.split('plt')[1] + '.h5')

hfo = h5py.File(outputfile_name, 'w')
hfo.create_dataset('MetNeutral', data=met_neutral)
hfo.create_dataset('MetWIM'    , data=met_wim)
hfo.create_dataset('MetWHIM'   , data=met_whim)
hfo.create_dataset('MetHot'    , data=met_hot)

hfo.create_dataset('MassNeutral', data=mass_neutral)
hfo.create_dataset('MassWIM'    , data=mass_wim)
hfo.create_dataset('MassWHIM'   , data=mass_whim)
hfo.create_dataset('MassHot'    , data=mass_hot)

hfo.create_dataset('Zrange'  , data=zrange)
hfo.create_dataset('Timestep', data=timestep)
hfo.close()

print("--------Written file------->>",outputfile_name)