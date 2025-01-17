import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns

folder_list = ['SummitData/GasGravity/Production2pc/R4/', 'SummitData/GasGravity/Production2pc/R8/',\
               'SetonixData/R16/4pc/']


kk = 0

fig, ax = plt.subplots(1, 2, gridspec_kw = {'wspace':0.02, 'hspace':0.02},figsize=(22, 6))
colors =  sns.color_palette("rocket", len(folder_list))

sigma_sfr = [0.000398107/yr_to_sec, 6.e-5/yr_to_sec , \
             1.58e-6/yr_to_sec, 1.58e-6/yr_to_sec]   

label = [r'4$\Sigma_{\rm Fid}$', r'$\Sigma_{\rm Fid}$', r'$0.2\Sigma_{\rm Fid}$']


Zunit = Msun/1.e3

for folder in folder_list[0:3]:
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    print(folder)
    os.chdir(data_path)
    list_file1 = glob.glob("proj_y_*")
    
    
    os.chdir(data_path)
    infile   = os.path.join(data_path, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    
    height = np.amax(zrange)/2.
    list_file = list_file1
    
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    vol_index = (np.abs(zrange)<height) &  (np.abs(zrange)>0.15*kpc)
    print(index, zrange.shape)
    
    timestep = np.zeros((len(list_file)))
    
    scalar_flux_xz = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_scalar = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_rho    = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    
    i = 0 
    for f in list_file:
        
        inputfile = os.path.join(data_path, f)
        file_size = get_folder_size(inputfile) / (1024 * 1024)
        
        if(file_size<100. and folder!=folder_list[-1]):
            print('Skipping ', inputfile, '\n')
            continue
        
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
        print(i, inputfile)
        timestep[i] = ds.current_time.to('Myr')
        
        scalar_flux_xz[i,:] = np.abs(np.array(data['scalar_outflow'])[:,0,:])
        tot_scalar[i,:]     = np.array(data['scalar'])[:,0,:]
        tot_rho[i,:]        = np.array(data['rho'])[:,0,:]
        i+=1
        
    
    mask1 = (np.abs(zrange)) < (height + 2*dz)
    mask2 = (np.abs(zrange)) > ( height - 2.*dz)
    mask = mask1*mask2
    
    scalar_flux_z = np.sum(scalar_flux_xz , axis=1)*Zunit*dx
    scalar_flux   = np.average(scalar_flux_z[:,mask], axis=1)  
    
    
    met_mass    = np.sum(tot_scalar[:,:,vol_index], axis=(1,2)) *dx * dz * Zunit
    total_mass  = np.sum(tot_rho[:,:,vol_index], axis=(1,2)) * dx * dz
    scalar_mass = 8.6e-3*total_mass + met_mass
    
    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = list(map(timestep.__getitem__, indexes))
    
    
    sorted_scalar  = list(map(scalar_mass.__getitem__, indexes))
    inj_scalar = sigma_sfr[kk]  * Msun * np.asarray(sorted_time)
    
    phi = 1. - np.asarray(sorted_scalar)/(inj_scalar*Myr)
    
    ax[0].plot(sorted_time, phi, ls='-',\
            label=label[kk],markersize=12, lw=2.5, color=colors[kk])
    
    
    sorted_scalar_flux  = list(map(scalar_flux.__getitem__, indexes))
    inj_scalar = sigma_sfr[kk]  * Msun * np.asarray(sorted_time) * Myr
    
    tt=1
    dtimestep = np.zeros(len(sorted_time))
    dtimestep[0] = sorted_time[0]
    while(tt<timestep.shape[0]):
        dtimestep[tt] = sorted_time[tt] - dtimestep[tt-1]
        tt+=1
    
    Mzdot_dt = np.asarray(sorted_scalar_flux) / (sigma_sfr[kk]  * Msun) #* dtimestep * Myr
    
    ax[1].plot(sorted_time,  Mzdot_dt, ls='-',\
            markersize=12, lw=2.5, color=colors[kk])
    
    kk+=1
    
ax[0].legend()
ax[0].set_xlabel('Myr')
ax[1].set_xlabel('Myr')
ax[0].tick_params(axis='y', right=True, left=True)
ax[1].tick_params(axis='y', labelleft=False, right=True, left=True, labelright=True)
ax[0].set_title('$1-\int M_Z dV/\Gamma_{SN}M_{\odot}t$')
ax[1].set_title('$ \dot{M}_Z /\Gamma_{SN}M_{\odot}$')
plt.setp(ax, 'ylim', (-0.01, 1.1))


filename = fig_path + '/mass_loading_wbg.jpeg'
plt.savefig(filename, bbox_inches='tight')
print('Created file--', filename)
