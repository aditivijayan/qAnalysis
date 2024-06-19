#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import argparse

# file = '/g/data/jh2/av5889/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackledata = h5py.File(grackle)
array = grackledata['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-10, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
parser.add_argument('input_folder', type=str, help='Path to input folder containing plt files')
parser.add_argument('overwrite',  action='store_true', default=0, help='Overwrite existing files, default=0')
args = parser.parse_args()

input_folder = args.input_folder
overwrite = args.overwrite

i=0
bins = 100
egas_arr = np.logspace(-22., -5., bins)
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

    
# temperature_table = interpolate.RectBivariateSpline(egas_arr, nH_arr, T)

data_path = os.path.join(scratch, input_folder)
output_folder = os.path.join(h5_path, input_folder , 'TRhoHisto/')
os.chdir(data_path)
list_file = glob.glob("plt*")

def makeDTHisto(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        print('Reached here!')
        f = item[0]
        infile = item[1]
        
        infile   = os.path.join(data_path, 'metal_uniform.in')
        dom_min, dom_max, ncells = getdomain(infile)
        fac = 1
        
        zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
        
        dx = (dom_max[0]- dom_min[0])/(1*int(ncells[0]))
        dy = (dom_max[1]- dom_min[1])/(1*int(ncells[1]))
        dz = (dom_max[2]- dom_min[2])/(1*int(ncells[2]))
        dVol = dx*dy*dz
        
       
        
        inputfile = os.path.join(data_path, f)
        ds   = yt.load(inputfile)
        lev=0
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)
        timestep = ds.current_time.to('Myr')

        rho = np.array(data['gasDensity'])
        Egas = np.array(data['gasInternalEnergy'])
        
        nH    = rho/hydrogen_mass_cgs
        mass = rho * dVol

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

        temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
                wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
            (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
                wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  
        
        bins = 201

        (dlow, dhigh)   = (1.e-6, 1.e2)
        (tlow, thigh)   = (10., 2.e10)

        zout = (np.abs(zrange)>kpc)
        temp_outflow = temp[:,:, zout]
        nH_outflow = nH[:,:, zout]
        mass_outflow = mass[:,:, zout]
        
        dedges = np.logspace(np.log10(dlow), np.log10(dhigh), bins)
        tedges = np.logspace(np.log10(tlow), np.log10(thigh), bins)

        htot, binx, biny = np.histogram2d(nH.reshape(-1,1)[:,0], temp.reshape(-1,1)[:,0], \
                                    bins=[dedges, tedges], weights=mass.reshape(-1,1)[:,0])

        hout, binx, biny = np.histogram2d(nH_outflow.reshape(-1,1)[:,0], temp_outflow.reshape(-1,1)[:,0], \
                                    bins=[dedges, tedges], weights=mass_outflow.reshape(-1,1)[:,0])
        
        

        
        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
        
        outputfile_name =os.path.join(output_folder, 'histo_' + f.split('plt')[1] + '.h5') 

        if not((os.path.exists(outputfile_name) and overwrite==0)):
            hfo = h5py.File(outputfile_name, 'w')
            hfo.create_dataset('DensBins'       , data=binx)
            hfo.create_dataset('TempBins'       , data=biny)
            hfo.create_dataset('TotalMass'       , data=htot)
            hfo.create_dataset('OutflowMass'       , data=hout)
            hfo.create_dataset('Timestep', data=timestep)
            hfo.close()
            print("--------Written file------->>",f)

queue      = Queue()
start_time = ostime.time()
listfile = list_file
num = len(listfile)
infile   = os.path.join(data_path, 'metal_uniform.in')
infile_list = [infile]*num
num_workers = os.cpu_count()


the_pool = [Process(target=makeDTHisto, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()

for i in range(num):
    queue.put((listfile[i],infile_list[i]))

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

print("Total time= %s seconds ---" % (ostime.time() - start_time))