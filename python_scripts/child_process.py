
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import sys

<<<<<<< HEAD
def process_file(foldername, filename):
    
    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    

    ds   = yt.load(inputfile)
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
    timestep = ds.current_time.to('Myr')
    egas    = np.array(data['gasEnergy'])
    eint    = np.array(data['gasInternalEnergy'])
    pvz = np.array(data['z-GasMomentum'])
    
    eout_rate = egas *  pvz * dx * dy * sign
    eth_rate  = eint *  pvz * dx * dy * sign

    energy_outflow_rate_z = np.sum(eout_rate, axis=(0,1))
    energy_outflow_rate = (energy_outflow_rate_z + energy_outflow_rate_z[::-1])/2.

    eint_outflow_z = np.sum(eth_rate, axis=(0,1))
    eint_outflow_rate = (eint_outflow_z + eint_outflow_z[::-1])/2

    output_folder = os.path.join(foldername , 'Eloading/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    outputfile_name =os.path.join(foldername, 'eloading_' + filename.split('plt')[1] + '.h5') 
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:
=======
def energy_loading(foldername, filename):
    h5path = '/gpfs/alpine2/proj-shared/ast196/aditi/qAnalysis/h5Data/'
    output_folder = os.path.join(h5path, foldername , 'Eloading/')
    print(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    outputfile_name = os.path.join(output_folder, 'eloading_' + filename.split('plt')[1] + '.h5')

    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    
    else:
    
        inputfile = os.path.join(foldername, filename)
        infile   = os.path.join(foldername, 'metal_uniform.in')

        dom_min, dom_max, ncells = getdomain(infile)
        fac = 1
        zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    
        sign = zrange/np.abs(zrange)    
        dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
        dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    

        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
        timestep = ds.current_time.to('Myr')
        rho_gas = np.array(data['gasDensity'])
        egas    = np.array(data['gasEnergy'])
        eint    = np.array(data['gasInternalEnergy'])
        vz = np.array(data['z-GasMomentum'])/rho_gas
    
        eout_rate = egas *  vz * dx * dy 
        eth_rate  = eint *  vz * dx * dy 

        energy_outflow_rate = np.sum(eout_rate, axis=(0,1))

        eint_outflow_rate = np.sum(eth_rate, axis=(0,1))


>>>>>>> 01b9b8a3f65e9dda89745545c6ced5d8351be271
        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Timestep'       , data=timestep)
        hfo.create_dataset('EnergyOutflowRate'       , data=energy_outflow_rate)
        hfo.create_dataset('ThEnergyOutflowRate'       , data=eint_outflow_rate)
        
        hfo.close()
        print("--------Written file------->>", inputfile)

<<<<<<< HEAD
def phase_mass(foldername, filename):

    file = '/g/data/jh2/av5889/freshquokka/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
    grackle = h5py.File(file)
    array = grackle['CoolingRates/Primordial/MMW'][()]
    table = array[:,0,:]
    table_nH   = np.logspace(-10., 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])

    i=0
    bins = 200
    egas_arr = np.logspace(-25., -5., bins)
    nH_arr   = np.logspace(-7.0, 4.0, int(bins))
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


    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    dV = dx * dy * dz

    ds   = yt.load(inputfile)
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
    timestep = ds.current_time.to('Myr')
    rho_gas = np.array(data['gasDensity'])
    eint    = np.array(data['gasInternalEnergy'])
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

    warm = temp<2.e4
    hot = temp>1.e6 

    mass = rho_gas * dV
    mass_warm = rho_gas[warm] * dV
    mass_hot  = rho_gas[hot] * dV
    vol = np.full(temp.shape, dx*dy)
    vol_warm = vol[warm]
    vol_hot = vol[hot]

    warm_mass_z  = np.sum(np.ma.array(mass, mask=~warm), axis=(0,1)).data 
    hot_mass_z  = np.sum(np.ma.array(mass, mask=~hot), axis=(0,1)).data 
    warm_vol_z = np.sum(np.ma.array(vol, mask=~warm), axis=(0,1)).data 
    hot_vol_z = np.sum(np.ma.array(vol, mask=~hot), axis=(0,1)).data 

    output_folder = os.path.join(foldername , 'Phases/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    outputfile_name =os.path.join(foldername, 'phase_' + filename.split('plt')[1] + '.h5') 
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:
        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Timestep'       , data=timestep)
        hfo.create_dataset('WarmMass'       , data=warm_mass_z)
        hfo.create_dataset('HotMass'       , data=hot_mass_z)
        hfo.create_dataset('WarmVol'       , data=warm_vol_z)
        hfo.create_dataset('HotVol'       , data=hot_vol_z)
        
        hfo.close()
        print("--------Written file------->>", inputfile)



=======
>>>>>>> 01b9b8a3f65e9dda89745545c6ced5d8351be271
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No file path provided.")
    else:
        foldername = sys.argv[1]
        filename = sys.argv[2]
        print(sys.argv)
<<<<<<< HEAD
        phase_mass(foldername, filename)
=======
        energy_loading(foldername, filename)
>>>>>>> 01b9b8a3f65e9dda89745545c6ced5d8351be271
