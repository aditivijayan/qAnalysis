#!/usr/bin/env python
# coding: utf-8

# In[35]:


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


# In[3]:

start_time = ostime.time()

file = '/scratch/pawsey0807/aditi/quokka_devbranch/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackle = h5py.File(file)
array = grackle['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-6, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])


# In[4]:

pc    = 3.018e18
kpc   = 1.e3 * pc 
Msun  = 2.e33
yr_to_sec =3.15e7
Myr   = 1.e6 * yr_to_sec
kmps  = 1.e5
boltzmann_constant_cgs    = 1.38e-16
hydrogen_mass_cgs = 1.67e-24
gamma = 5./3.


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
    
    
home  = "/scratch/pawsey0807/aditi/quokka_devbranch/quokka/sims/"
build = "Re2pc/"
# run   = "FakeAMR/"
# run   = "29Myr/"
run = 'Start28Myr/'
lev =0

# build = "MetalGrad/"
# run = "4pcBase-8.0Eta/"
# lev = 1

# build = "HighSigGas/"
# run = "Run/"
# lev = 0

path = home + build + run
os.chdir(path)
list_file = glob.glob("plt*")


print("Analysing -- ", len(list_file), ' --files.')

if(lev==0):
    fac = 1
else:
    fac = 2 * lev
Data.fac = fac
infile = home +  build + run + 'metal_uniform_512.in'
dom_min, dom_max, ncells = getdomain(infile)
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))


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

# temperature_table = interpolate.RectBivariateSpline(egas_arr, nH_arr, T)



def getPhaseOutflowRate(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        f = item[0]
        infile = item[1]
        
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
        
        inhome = '/scratch/pawsey0807/aditi/quokka_devbranch/quokka/sims/'
        runname = infile.split('/')[-3] + '/' + infile.split('/')[-2] +  '/' 
        inputfile = inhome + runname + f
        
        ds   = yt.load(inputfile)
        
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)
        timestep = ds.current_time.to('Myr')
        
        
        rho = np.array(data['gasDensity'])
        rhoZ = np.array(data['scalar_2'])
        pz = np.array(data['z-GasMomentum'])
        vz = pz/rho
        Egas = np.array(data['gasInternalEnergy'])
        
        pz = np.array(data['z-GasMomentum'])
        vz = pz/rho
        dA   = np.full(rho.shape, dx * dy)  
        dVol   = np.full(rho.shape, dx * dy * dz)  
        egas0 = Egas
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


        net_mass_outflow_rate   = rho * vz * dA
    
        ZOinit = 8.6e-3
        MO = 1. * Msun
        rhoOxy_inj = MO * rhoZ/1.e3 #1.e3 is the normalization of passive_sclar
        
        net_mass_outflow_rate   = rho * vz * dA
        net_oxygen_outflow_rate_inj = rhoOxy_inj *  vz * dA

        hot = (Temp>1.e6)
        warm = (Temp<2.e4)
        X, Y, Z = np.meshgrid(xrange, yrange, zrange)
        # zout  = (vz*zrange>0.0)
        zout  = (np.abs(vz)>0.0)
        hout  = hot #* zout
        wout  = warm #* zout
        
        #Mass Outflow rate#
        total_mass_outflow_rate = np.sum(np.ma.array(net_mass_outflow_rate, mask=~zout), axis=(0,1)).data  
        hot_mass_outflow_rate   = np.sum(np.ma.array(net_mass_outflow_rate, mask=~hout), axis=(0,1)).data 
        warm_mass_outflow_rate  = np.sum(np.ma.array(net_mass_outflow_rate, mask=~wout), axis=(0,1)).data 
        
        
        #Inj Oxygen Mass Outflow rate#
        total_metal_outflow_rate  = np.sum(np.ma.array(net_oxygen_outflow_rate_inj, mask=~zout), axis=(0,1)).data 
        hot_metal_outflow_rate    = np.sum(np.ma.array(net_oxygen_outflow_rate_inj, mask=~hout), axis=(0,1)).data 
        warm_metal_outflow_rate   = np.sum(np.ma.array(net_oxygen_outflow_rate_inj, mask=~wout), axis=(0,1)).data 
        
        
        outhome = '/scratch/pawsey0807/aditi/quokka_devbranch/Analysis/'
        # runname = infile.split('/')[8] + '/' + infile.split('/')[9] + '/'
        outputfile_name = outhome +  runname + 'PhaseOutflowRates/phase_outflow_' + f.split('plt')[1] + '.h5'
        
        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('WarmOutflowRate'  , data=warm_mass_outflow_rate)
        hfo.create_dataset('HotOutflowRate'            , data=hot_mass_outflow_rate)
        hfo.create_dataset('TotalOutflowRate'            , data=total_mass_outflow_rate)

        hfo.create_dataset('WarmMetOutflowRate'  , data=warm_metal_outflow_rate)
        hfo.create_dataset('HotMetOutflowRate'            , data=hot_metal_outflow_rate)
        hfo.create_dataset('TotalMetOutflowRate'            , data=total_metal_outflow_rate)

        hfo.create_dataset('Zrange'  , data=zrange)
        hfo.create_dataset('Timestep', data=timestep)
        hfo.close()
      
        print("--------Written file------->>",f)


queue      = Queue()
listfile = list_file[4:14]
num = len(listfile)
infile_list = [infile]*num
num_workers = os.cpu_count()
print('No. of workers=', num_workers)
# cProfile.run('getTgasfromEgas(queue)')
the_pool = [Process(target=getPhaseOutflowRate, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()
count_exist=0
count_new=0
for i in range(num):
    outhome = '/scratch/pawsey0807/aditi/quokka_devbranch/Analysis/'
    runname = infile_list[i].split('/')[8] + '/' + infile_list[i].split('/')[9] + '/'
    outputfile_name = outhome + build + run + 'PhaseOutflowRates/phase_outflow_' + listfile[i].split('plt')[1] + '.h5'
    file_time = datetime.datetime.fromtimestamp(os.path.getmtime(path + list_file[i]))
    # time_last = datetime.datetime.fromtimestamp(os.path.getmtime(path + 'phase_outflow_1035000.h5'))
    # if(file_time<time_last):
    # print(listfile[i], infile_list[i])
    if((os.path.isfile(outputfile_name))):
        print('File--', outputfile_name, ' exists already!')
        count_exist+=1     
    else:
        queue.put((listfile[i],infile_list[i]))
        count_new+=1

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

count = 0 
# for i in range(num):
#     outhome = '/scratch/pawsey0807/aditi/quokka_devbranch/Analysis'
#     runname = infile_list[i].split('/')[8] + '/' + infile_list[i].split('/')[9] + '/'
#     outputfile_name = outhome +  runname + 'PhaseOutflowRates/phase_outflow_' + listfile[i].split('plt')[1] + '.h5'
#     if((os.path.isfile(outputfile_name))): count+=1
        

print("Total time= %s seconds ---" % (ostime.time() - start_time))
print('Total files analysed = ', count)
print("Files new=", count_new)
print("Files old=", count_exist)
