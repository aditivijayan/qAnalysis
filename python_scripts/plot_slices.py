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
parser.add_argument('--overwrite',  action='store_true', default=1, help='Overwrite existing files, default=0')
args = parser.parse_args()

input_folder = args.input_folder
overwrite = args.overwrite


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

    
# temperature_table = interpolate.RectBivariateSpline(egas_arr, nH_arr, T)

data_path = os.path.join(scratch, 'sims/', input_folder)
output_folder = os.path.join(fig_path, input_folder + 'Slice/')
print(output_folder)
os.chdir(data_path)
print(data_path)
list_file = glob.glob("plt*")
if not os.path.exists(output_folder):
    print(output_folder)
    os.makedirs(output_folder)
    print("Created output folder!\n")
 
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
        eint    = np.array(data['gasInternalEnergy'])
        vz = np.array(data['z-GasMomentum'])/rho_gas
        #vx = np.array(data['x-GasMomentum'])/rho_gas
        #vy = np.array(data['y-GasMomentum'])/rho_gas
        
        plane = (int)(ncells[1]/2)

        egas0=eint
        rho_slice = rho_gas[:,plane,:]
        eint_slice = eint[:,plane,:]

        density = rho_slice
        egas0=eint_slice
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

        
        print("Temp estimated!\n")
        print("Total time till Temp= %s seconds ---" % (ostime.time() - start_time))

        fig, ax = plt.subplots(1, 3, gridspec_kw = {'wspace':0.02, 'hspace':0.02},figsize=(12, 32))
    
        cbarx = 0.13
        cbheight = 0.04
        cbary = 0.89
        cblen = 0.2
        dx1 = 0.25
        cbtitlex = 0.1
        cbtitley = 16.5
        plane = (int)(ncells[1]/2)
       

        plot = ax[0].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(rho_slice/mp),\
                            norm=mcolors.LogNorm(vmin=1.e-6, vmax=4.e2),
                            cmap='Blues')
        cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e-6,  1.e-2, 1., 1.e2))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $\rho$" + "\n" + "[g cm$^{-3}$]")
       

        # plot = ax[1].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(eint[:,plane,:]),\
        #                     norm=mcolors.LogNorm(vmin=9.e-21, vmax=5.e-10),
        #                     cmap='afmhot')
        # cax = fig.add_axes([cbarx + 3.*dx1, cbary, cblen, cbheight])
        # fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e-20, 1.e-16, 1.e-10, 1.e-8))
        # cax.xaxis.set_ticks_position('top')
        # cax.set_title(r" eint" + "\n" + " [cgs]")
        # ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=True)
        # ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=True)

        plot = ax[1].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(temp),\
                    norm=mcolors.LogNorm(vmin=2.e2, vmax=5.e7),
                    cmap='afmhot')
        cax = fig.add_axes([cbarx + dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e2, 1.e4,1.e5, 1.e6,1.e8))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" T" + "\n" + " [K]")
        ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=True)

        ax[0].text(0.5, 0.5, '%.2f'%(ds.current_time.to('Myr')))

        plot = ax[2].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(vz[:,plane,:])/kmps,\
                            vmin=-500., vmax=500.,
                            cmap='seismic')
        cax = fig.add_axes([cbarx + 2.*dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-500.,0.0, 500.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" vz" + "\n" + " [kmps]")
        ax[2].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=True)
        
        plt.setp(ax, 'ylim', (-7.9, 8.0))
        
        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
        
        outputfile_name =os.path.join(output_folder, 'density-vz-' + f.split('plt')[1] + '.jpeg')
        plt.savefig(outputfile_name, bbox_inches='tight', dpi=160)
        print('Written file----->',outputfile_name)
        
        
queue      = Queue()
start_time = ostime.time()
listfile = list_file[17:27]
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
