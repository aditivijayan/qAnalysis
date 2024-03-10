#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import argparse

file = '/g/data/jh2/av5889/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackle = h5py.File(file)
array = grackle['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-6, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
parser.add_argument('input_folder', type=str, help='Path to input folder containing plt files')
args = parser.parse_args()

input_folder = args.input_folder

i=0
bins = 100
egas_arr = np.logspace(-16., -5., bins)
nH_arr   = np.logspace(-6.0, 4.0, int(bins))
T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))


i=0
bins = 100
egas_arr = np.logspace(-16., -5., bins)
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

data_path = os.path.join(scratch, 'sims/' + input_folder)
output_folder = os.path.join(fig_path, input_folder + '/TotalMass/')
os.chdir(data_path)
list_file = glob.glob("plt*")

def makeSlicePlot(queue):
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
        
        dx = xrange[1] - xrange[0]
        dy = yrange[1] - yrange[0]
        dz = zrange[1] - zrange[0]
        dV = dx * dy * dz

        lev = 0
    
        inputfile = os.path.join(data_path, f)
        ds   = yt.load(inputfile)
        timestep = ds.current_time.to('Myr')
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)

        rho_gas = np.array(data['gasDensity'])
        
        total_mass = np.sum(rho_gas*dV)
        time = ds.current_time.to('Myr')

        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
        
        outputfile_name =os.path.join(output_folder, 'total_mass_' + f.split('plt')[1] + '.txt')
        # plt.savefig(outputfile_name, bbox_inches='tight', dpi=160)
        np.savetxt(outputfile_name, np.column_stack([time, total_mass]))
        
        
queue      = Queue()
start_time = ostime.time()
listfile = list_file
num = len(listfile)
infile   = os.path.join(data_path, 'metal_uniform.in')
infile_list = [infile]*num
num_workers = os.cpu_count()


the_pool = [Process(target=makeSlicePlot, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()

for i in range(num):
    queue.put((listfile[i],infile_list[i]))

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

print("Total time= %s seconds ---" % (ostime.time() - start_time))