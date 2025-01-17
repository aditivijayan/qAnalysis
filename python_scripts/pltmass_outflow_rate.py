
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns


source = 'Summit' 
# source = 'Gadi'

if(source == 'Summit'):
    # folder_list = ['SummitData/GasGravity/Production2pc/R8/', 
    folder_list = ['SummitData/GasGravity/Production2pc/R8/', 'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
               'SummitData/GasGravity/Production2pc/R4/', 'SummitData/GasGravity/Production2pc/R4-0.2Zsol/',\
                  'SetonixData/R16/4pc/', 'SetonixData/R16-0.2Zsol/4pc/']
    data_path0 = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder_list[0])
    
else:
    folder_list = ['GasGravity/', 'GasGravity/R8-0.2Zsol/', 'GasGravity/R4/', 'GasGravity/R4-0.2Zsol/']
    data_path0 = os.path.join(scratch, 'sims/', folder_list[0])


k=0

fig, ax = plt.subplots(1, 1, gridspec_kw = {'wspace':0.02, 'hspace':0.02},figsize=(12, 8))


label= ['R8', r'$0.2Z_{\odot}$', 'R4', '', 'R16', '']
ls = ['-', '--', '-', '--', '-', '--' ]

clors =  sns.color_palette("rocket", np.ceil(len(folder_list)/2).astype(int))
colors = [clors[0], clors[0], clors[1], clors[1], clors[2], clors[2]]

for folder in folder_list:
    if(source=='Summit'):
        data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    else:
        data_path = os.path.join(scratch, 'sims/', folder)

    os.chdir(data_path)
    list_file = glob.glob("proj_y_*")
    
    infile   = os.path.join(data_path, 'metal_uniform.in')
    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))

    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    
    dom_min, dom_max, ncells = getdomain(infile)
    
    fac = 1
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    munit = dx*yr_to_sec/Msun/1.e-2
    
    timestep = np.zeros((len(list_file)))
    tot_mass_flux = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))

    tot_scalar_flux = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))

    
    height = 2. * kpc
    if(k>=4):height = 6. * kpc
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    
    
    i = 0 
    for f in list_file:
        inputfile = os.path.join(data_path, f)
        
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
        timestep[i] = ds.current_time.to('Myr')
        tot_mass_flux[i,:] = np.array(data['mass_outflow'])[:,0,:]
        i+=1
        
    print('Max timestep = ', np.amax(timestep))
    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = list(map(timestep.__getitem__, indexes))
    
    flux_sum = np.sum(tot_mass_flux[:,:,index], axis=1)
    sorted_mdot  = list(map(flux_sum.__getitem__, indexes))
    ax.plot(sorted_time, np.asarray(sorted_mdot)*munit, ls=ls[k],\
            label=label[k],markersize=12, lw=2.5, color=colors[k])
    k+=1

ax.legend()
ax.set_ylabel('$\dot{M}$' + '\n' + '[$10^{-2} $M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$]')
ax.set_title('Total Mass Outflow')
ax.axhline(0.0, ls='--', color='darkslategray')
ax.set_yscale('Symlog', linthresh=1.e-4)
ax.set_ylim(-10., 10.)    
plt.savefig(fig_path + '/tot_mass_outflow.jpeg', bbox_inches='tight')