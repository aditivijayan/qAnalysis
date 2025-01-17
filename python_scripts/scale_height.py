#!/usr/bin/env python
# coding: utf-8


import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns

plt.rcParams['font.size']=28
plt.rcParams['axes.linewidth']=3.
plt.rcParams['xtick.major.size']=12
plt.rcParams['xtick.minor.size']=1
plt.rcParams['xtick.major.width']=2
plt.rcParams['xtick.minor.width']=2.
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.major.size']=12
plt.rcParams['ytick.minor.size']=6
plt.rcParams['ytick.major.width']=2
plt.rcParams['ytick.minor.width']=2
plt.rcParams['ytick.direction']='in'


folder_list = ['SetonixRuns/R16/4pc/Redo/',\
               'SummitData/GasGravity/Production2pc/R16-h300-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h1kpc-Zsol',\
              'SummitData/GasGravity/Production2pc/R16-h1.5kpc-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h2kpc-Zsol']


sigma_sfr = 1.58e-6/yr_to_sec

hSN = [150., 300., 1000., 1500., 2000.]

color =sns.color_palette("plasma",n_colors=len(folder_list)) 

sfr = np.asarray(sigma_sfr) * Msun * 1.e2

fig, ax = plt.subplots(1, 1, gridspec_kw = {'wspace':0.0, 'hspace':0.0},figsize=(12, 6))
kk=0
ls = ['-', '--', '-.', ':', '--', '-']
lw=2.0

label = [r'$150$ pc', r'$300$ pc',  r'$1$ kpc', r'$1.5$ kpc', r'$2$ kpc', ]

for folder in folder_list:
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    if(kk==0):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder)
    print(data_path)
    
    os.chdir(data_path)
    list_file1 = glob.glob("proj_y_plt*")
    infile   = os.path.join(data_path, 'metal_uniform.in')

    
    list_file = list_file1
    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))


    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))

    tot_mass = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
    tot_mass_avg = np.zeros((len(list_file), zrange.shape[0]))
    scale_hgt = np.zeros(len(list_file))
    timestep = np.zeros(len(list_file))

    i=0
    for f in list_file:
        inputfile = os.path.join(data_path, f)
        file_size = get_folder_size(inputfile) / (1024 * 1024)
        
        if(file_size<1.):
            print('Skipping ', inputfile, '\n')
            continue
            
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions)
        
        timestep[i] = ds.current_time.to('Myr')
        tot_mass[i,:] = np.array(data['rho'])[:,0,:]
        i+=1  
    
    
    tot_mass_x = np.sum(tot_mass, axis=1)*dx*dz/kpc/kpc
    
    for i in range(tot_mass_x.shape[0]):
        tot_mass_avg[i] = (tot_mass_x[i,:] + tot_mass_x[i, ::-1])/2.
    
    size = int (tot_mass_avg.shape[1]/2)
    pos_z = tot_mass_avg[:, size:]
    cumsum = np.cumsum(pos_z, axis=1)
    
    for tt in range(cumsum.shape[0]):
        index = np.where(cumsum[tt]>= (np.sum(pos_z[tt])/np.exp(1)))[0][0]
        scale_hgt[tt] = zrange[size+index]
        
    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = list(map(timestep.__getitem__, indexes))
    
    sorted_h  = list(map(scale_hgt.__getitem__, indexes))
    
    ax.plot(sorted_time , np.asarray(sorted_h)/pc, ls=ls[kk], lw=lw+0.5, color=color[kk], label=label[kk])

    print('Max scale_h = ', np.amax(scale_hgt)/pc)
    lw=lw+0.5
    kk+=1
    
    
ax.set_ylabel(r'$h_{\rm gas}$ [pc]')
ax.set_xlabel('t [Myr]')
ax.set_xlim(0.5,300.)
ax.set_ylim(0.0, 3.e3)
ax.legend(frameon=False, ncol=3)
# plt.setp(ax, ('ylim'), (0.0, 10.))
filename = fig_path + '/Paper/scale_height.jpeg'
plt.savefig(filename, bbox_inches='tight', dpi=160)
print('Created file--', filename)