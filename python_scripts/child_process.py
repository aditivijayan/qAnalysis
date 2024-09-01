
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import sys

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


        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Timestep'       , data=timestep)
        hfo.create_dataset('EnergyOutflowRate'       , data=energy_outflow_rate)
        hfo.create_dataset('ThEnergyOutflowRate'       , data=eint_outflow_rate)
        
        hfo.close()
        print("--------Written file------->>", inputfile)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No file path provided.")
    else:
        foldername = sys.argv[1]
        filename = sys.argv[2]
        print(sys.argv)
        energy_loading(foldername, filename)
