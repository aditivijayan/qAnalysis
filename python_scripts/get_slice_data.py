#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import argparse

grackledata = h5py.File(grackle)
array = grackledata['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-10, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
args = parser.parse_args()

input_folder = args.input_folder

i=0
bins = 100
egas_arr = np.logspace(-24., -5., bins)
nH_arr   = np.logspace(-6.0, 4.0, int(bins))
T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))


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

    
data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R4/')

def makeSlice(queue):
    f = 'plt7720000/'
    inputfile = os.path.join(data_path, f)
    infile   = os.path.join(data_path, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)

    zrange = np.linspace(dom_min[2], dom_max[2], (int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (int(ncells[1])))

    dx = (dom_max[0]- dom_min[0])/(int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(int(ncells[2]))

    Lx = (dom_max[0]- dom_min[0])
    Ly = (dom_max[1]- dom_min[1])
    Lz = (dom_max[2]- dom_min[2])

    plane = int(ncells[1]/2)
    lev = 0

    ds   = yt.load(inputfile)
    timestep = ds.current_time.to('Myr')
    data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions)
    print("0.This is okay!")
    rho_gas = np.array(data['gasDensity'])
    print("1.This is okay!")
    eint    = np.array(data['gasInternalEnergy'])
    print("2.This is okay!")
    vz = np.array(data['z-GasMomentum'])/rho_gas
    print("3.This is okay!")
    vx = np.array(data['x-GasMomentum'])/rho_gas
    vy = np.array(data['y-GasMomentum'])/rho_gas


    eint[eint<0.0] = 1.e-28
    egas0=eint
    density = rho_gas
    cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
    rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs


    logrho_arr = np.log10(nH_arr[:-1])
    logrho     = np.log10(rho0)
    delta_rho  = logrho_arr[1] - logrho_arr[0]
    idxrho     = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas     = np.log10(egas0)
    delta_egas  = logEgas_arr[1] - logEgas_arr[0]
    idxegas     = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
            wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
        (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
            wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  


    if not os.path.exists(output_folder):
        print(output_folder)
        os.makedirs(output_folder)

    outputfile_name = os.path.join(data_path, 'slice-' + f.split('plt')[1] + '.h5')
    hfo = h5py.File(outputfile_name, 'w')
    hfo.create_dataset('Zrange'       , data=zrange)
    hfo.create_dataset('Xrange'       , data=xrange)
    hfo.create_dataset('Timestep'       , data=timestep)
    hfo.create_dataset('Energy'       , data=egas[:,plane,:])
    hfo.create_dataset('Rho'       , data=rho_gas[:,plane,:])
    hfo.create_dataset('Temperature'       , data=temp[:,plane,:])
    hfo.create_dataset('Vz'       , data=vz[:,plane,:])

    hfo.close()
    print("--------Written file------->>", inputfile)

queue      = Queue()
start_time = ostime.time()
listfile = [data_path]
num_workers = os.cpu_count()


the_pool = [Process(target=makeSlice, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()

for i in range(1):
    queue.put((listfile[i]))

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

print("Total time= %s seconds ---" % (ostime.time() - start_time))



