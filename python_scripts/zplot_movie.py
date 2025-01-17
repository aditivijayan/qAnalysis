#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns

plt.rcParams['font.size']=42
plt.rcParams['axes.linewidth']=0.2
plt.rcParams['xtick.major.size']=4
plt.rcParams['xtick.minor.size']=2
plt.rcParams['xtick.major.width']=4.1
plt.rcParams['xtick.minor.width']=2.
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.major.size']=12
plt.rcParams['ytick.minor.size']=6
plt.rcParams['ytick.major.width']=2
plt.rcParams['ytick.minor.width']=2
plt.rcParams['ytick.direction']='in'

folder_list = [ 'SummitData/GasGravity/Production2pc/R8/', 'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
               'SummitData/GasGravity/Production2pc/R4/', 'SummitData/GasGravity/Production2pc/R4-0.2Zsol/',\
               'SummitData/GasGravity/Production2pc/R4-h75-0.2Zsol',\
                  'SetonixData/R16/4pc/', 'SetonixData/R16-0.2Zsol/',\
              'SummitData/GasGravity/Production2pc/R16-h300-Zsol/',\
                'SummitData/GasGravity/Production2pc/R16-h1kpc-Zsol/',\
                    'SummitData/GasGravity/Production2pc/R16-h2kpc-Zsol/',
                    'SummitData/GasGravity/Production2pc/R16-h1.5kpc-Zsol/']


input_folder = folder_list[5]

data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', input_folder)

os.chdir(data_path)
list_file = glob.glob("proj_y_plt*/")
output_folder = os.path.join(fig_path, input_folder + '/ZplotMov/')

Zunit = Msun/1.e3
infile   = os.path.join(data_path, 'metal_uniform.in')

dom_min, dom_max, ncells = getdomain(infile)
fac = 1
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))


dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))


max_time=0.0
i=0
k=0

def makeSlicePlot(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        
        f = item[0]
        infile = item[1]

        inputfile = os.path.join(data_path, f)

        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
        print(inputfile)
        curr_timestep = ds.current_time.to('Myr')

        tot_rho = np.array(data['rho'])[:,0,:]
        tot_scal = np.array(data['scalar'])[:,0,:]*Zunit
        
    
        
        fig, ax = plt.subplots(1, 2, gridspec_kw = {'wspace':0.03, 'hspace':0.0},figsize=(24, 144))
        i=0
        plt.style.use('dark_background')

        cbarx = 0.05
        cbheight = 0.755
        cbary = 0.125
        cblen = 0.06
        dx1 = 0.86
        dx2 = 0.8
        cbtitlex = 0.15
        cbtitley = 16.5



        plot = ax[0].pcolormesh(xrange/kpc,zrange/kpc, np.transpose(tot_scal/tot_rho),\
                            norm=mcolors.LogNorm(vmin=1.e-4, vmax=1.e-1),
                            cmap=sns.color_palette("Spectral", as_cmap=True))
        cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='vertical', ticks=(1.e-4, 1.e-3, 1.e-2, 1.e-1))
        cax.xaxis.set_tick_params(rotation=45, labelleft=True, left=True, right=False)
        cax.xaxis.set_ticks_position('both')
        cax.yaxis.set_ticks_position('left')


        ax[0].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=False)
        ax[0].tick_params(axis='x', labelbottom=False, bottom=False)
        ax[0].text(0.05,0.96, '%.1f'%(curr_timestep) + ' Myr', color='purple', \
                   transform=ax[0].transAxes, fontsize=100)
        ax[0].text(-0.4,0.5, r" $\langle Z \rangle $ [$Z_{\odot}$]", color='white',\
                    transform=ax[0].transAxes, fontsize=100, rotation=90)

        plot = ax[1].pcolormesh(xrange/kpc,zrange/kpc, np.transpose(tot_rho)/mp,\
                            norm=mcolors.LogNorm(vmin=1.e16, vmax=1.e22),
                            cmap=sns.color_palette("cubehelix", as_cmap=True))
        cax = fig.add_axes([cbarx+dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='vertical', ticks=(1.e16, 1.e18, 1.e20, 1.e22))
        cax.xaxis.set_ticks_position('top')
        cax.xaxis.set_tick_params(rotation=45)
        cax.xaxis.set_ticks_position('bottom')


        ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=False, left=False)
        ax[1].tick_params(axis='x', labelbottom=False, bottom=False)
        ax[1].text(1.2,0.5, r" $\mathrm{N}_H $ [cm$^{-2}$]", color='white', \
                   transform=ax[1].transAxes, fontsize=100, rotation=270)

        # plt.setp(ax, 'ylim', (0.0, 6.0))


        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
                
        outputfile_name =os.path.join(output_folder, 'zmov-' + f.split('proj_y_plt')[-1].split('/')[0] + '.jpeg')
        if os.path.isfile(outputfile_name):
            print("File already exists-->", outputfile_name)
        else:
            plt.savefig(outputfile_name, bbox_inches='tight', dpi=160)
            print('Created file--', outputfile_name)


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