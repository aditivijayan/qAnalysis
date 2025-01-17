#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *


folder_list = [ 'SummitData/GasGravity/Production2pc/R8/', 'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
               'SummitData/GasGravity/Production2pc/R4/', 'SummitData/GasGravity/Production2pc/R4-0.2Zsol/',\
                  'SetonixData/R16/4pc/', 'SetonixData/R16-0.2Zsol/4pc/', \
                    'SummitData/GasGravity/Production2pc/R16-h300-Zsol/']


folder_list = ['SetonixData/R16/4pc/', 'SetonixData/R16-0.2Zsol/4pc/',\
                  'SummitData/GasGravity/Production2pc/R16-h300-Zsol/']

Zunit = Msun/1.e3
input_folder = folder_list[-1]
output_folder = os.path.join(fig_path, input_folder + '/ProjPlot/')
data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', input_folder)
os.chdir(data_path)
list_file = glob.glob("proj*")



cbarx = 0.15
cbheight = 0.02
cbary = 0.89
cblen = 0.2
dx1 = 0.25
cbtitlex = 0.1
cbtitley = 16.5


def makeSlicePlot(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        
        f = item[0]
        infile = item[1]
        
        dom_min, dom_max, ncells = getdomain(infile)
        fac = 1
        zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
        xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
        yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))

        inputfile = os.path.join(data_path, f)

        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
        print(inputfile)
        curr_timestep = ds.current_time.to('Myr')
        
        tot_mass_flux = np.array(data['mass_outflow'])[:,0,:]
        tot_rho = np.array(data['rho'])[:,0,:]

        hot_mass_flux = np.array(data['hot_mass_outflow'])[:,0,:]
        warm_mass_flux = np.array(data['warm_mass_outflow'])[:,0,:]

        hot_scal_flux = np.array(data['hot_scalar_outflow'])[:,0,:]*Zunit
        warm_scal_flux = np.array(data['warm_scalar_outflow'])[:,0,:]*Zunit
        tot_scal = np.array(data['scalar'])[:,0,:]*Zunit

        fig, ax = plt.subplots(1, 3, gridspec_kw = {'wspace':0.02, 'hspace':0.02},figsize=(18, 48))
        i=0


        plot = ax[0].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(tot_rho),\
                            norm=mcolors.LogNorm(vmin=1.e-8, vmax=1.e-2),
                            cmap='Blues')
        cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e-6,  1.e-2, 1.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $\rho$" + "\n" + "[g cm$^{-2}$]")


        plot = ax[1].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(hot_mass_flux)/1.e3,\
                            vmin=-0.5, vmax=0.5,
                            cmap='seismic')
        cax = fig.add_axes([cbarx + dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-1., 0.0, 1.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $\rho$" + "\n" + "[g cm$^{-3}$]")

        plot = ax[2].pcolormesh(yrange/kpc,zrange/kpc, np.transpose(warm_mass_flux)/1.e2,\
                            vmin=-1.15, vmax=1.15,
                            cmap='seismic')
        cax = fig.add_axes([cbarx + 2.*dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-1., 0.0, 1.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $\rho$" + "\n" + "[g cm$^{-3}$]")

        ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=True, left=True)
        ax[2].tick_params(axis='y', labelleft=False, labelright=False, right=True, left=True)
        plt.setp(ax, 'ylim', (-7.9, 8.))
        plt.setp(ax[0:-1], 'xlim', (0.0, 0.95))
        ax[0].text(0.4, 0.9, '%.2f'%(curr_timestep) + ' Myr', transform = ax[0].transAxes, color='black')

        outputfile_name =os.path.join(output_folder, 'proj_plot_' + f.split('plt')[1] + '.jpeg')
        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
        plt.savefig(outputfile_name, bbox_inches='tight', dpi=160)
            
        
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
