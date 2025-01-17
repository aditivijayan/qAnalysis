#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import argparse

#file = '/g/data/jh2/av5889/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackledata = h5py.File(grackle)
array = grackledata['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-10, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
parser.add_argument('--overwrite',  action='store_true', default=1, help='Overwrite existing files, default=0')
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

data_path = os.path.join(scratch, 'sims', input_folder)
output_folder = os.path.join(h5_path, input_folder , 'PhaseOutflowRates/')
os.chdir(data_path)
list_file = glob.glob("plt*")

def getPhaseOutflowRate(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        
        f = item[0]
        infile = item[1]
        
        infile   = os.path.join(data_path, 'metal_uniform.in')
        dom_min, dom_max, ncells = getdomain(infile)
        fac = 1
        zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
        xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
        yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))
        
        dx = (dom_max[0]- dom_min[0])/(1*int(ncells[0]))
        dy = (dom_max[1]- dom_min[1])/(1*int(ncells[1]))
        dz = (dom_max[2]- dom_min[2])/(1*int(ncells[2]))
        
        Lx = (dom_max[0]- dom_min[0])
        Ly = (dom_max[1]- dom_min[1])
        Lz = (dom_max[2]- dom_min[2])
        
        plane = int(ncells[1]*fac/2)
        lev = 0
    
        inputfile = os.path.join(data_path, f)
        ds   = yt.load(inputfile)
        timestep = ds.current_time.to('Myr')
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)

        rho_gas = np.array(data['gasDensity'])
        Eint    = np.array(data['gasInternalEnergy'])
        vz = np.array(data['z-GasMomentum'])/rho_gas
        vx = np.array(data['x-GasMomentum'])/rho_gas
        vy = np.array(data['y-GasMomentum'])/rho_gas
        rhoZ = np.array(data['scalar_0'])
        
        egas0=Eint
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

        dA   = np.full(rho_gas.shape, dx * dy)  
        dVol   = np.full(rho_gas.shape, dx * dy * dz) 
        net_mass_outflow_rate   = rho_gas * vz * dA
        
        ZOinit = 8.6e-3
        MO = 1. * Msun
        rhoOxy_inj = MO * rhoZ/1.e3 #1.e3 is the normalization of passive_sclar
        
        net_mass_outflow_rate   = rho_gas * vz * dA
        net_oxygen_outflow_rate_inj = rhoOxy_inj *  vz * dA

        hot = (temp>1.e6)
        warm = (temp<2.e4)

        zout  = (vz*zrange>0.0)
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
        
        total_scalar_value        = np.sum(rhoZ * dVol, axis=(0,1) )
        total_mass                = np.sum(rho_gas  * dVol, axis=(0,1))

        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
        
        outputfile_name =os.path.join(output_folder, 'phase_totoutflow_' + f.split('plt')[1] + '.h5')

        if not((os.path.exists(outputfile_name) and overwrite==0)):
            hfo = h5py.File(outputfile_name, 'w')
            hfo.create_dataset('WarmOutflowRate'  , data=warm_mass_outflow_rate)
            hfo.create_dataset('HotOutflowRate'            , data=hot_mass_outflow_rate)
            hfo.create_dataset('TotalOutflowRate'            , data=total_mass_outflow_rate)

            hfo.create_dataset('WarmMetOutflowRate'  , data=warm_metal_outflow_rate)
            hfo.create_dataset('HotMetOutflowRate'            , data=hot_metal_outflow_rate)
            hfo.create_dataset('TotalMetOutflowRate'            , data=total_metal_outflow_rate)

            hfo.create_dataset('TotalScalarValue'            , data=total_scalar_value)
            hfo.create_dataset('TotalMass'            , data=total_mass)

            hfo.create_dataset('Zrange'  , data=zrange)
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


the_pool = [Process(target=getPhaseOutflowRate, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()

for i in range(num):
    queue.put((listfile[i],infile_list[i]))

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

print("Total time= %s seconds ---" % (ostime.time() - start_time))
