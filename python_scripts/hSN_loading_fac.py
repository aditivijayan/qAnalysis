
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns

plt.rcParams['font.size']=28
plt.rcParams['axes.linewidth']=1.
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

def get_folder_size(directory):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            # Skip if it is symbolic link
            if not os.path.islink(file_path):
                total_size += os.path.getsize(file_path)
    return total_size


folder_list = ['SetonixData/R16/4pc/',\
               'SummitData/GasGravity/Production2pc/R16-h300-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h1kpc-Zsol',\
              'SummitData/GasGravity/Production2pc/R16-h1.5kpc-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h2kpc-Zsol']

sigma_sfr = 1.58e-6/yr_to_sec
kk= 0
color =sns.color_palette("husl",n_colors=len(folder_list)) 

sfr = np.asarray(sigma_sfr) * Msun * 1.e2

fig, ax = plt.subplots(2, 1, gridspec_kw = {'wspace':0.01, 'hspace':0.01},figsize=(12, 12))

hSN = [0.15, 0.3, 1., 1.5, 2.]

eta_time_avg = [] 
eta_16_avg = []
eta_84_avg = []
phi_tavg = []
phi_16 = []
phi_84 = []

Zunit = Msun/1.e3

for folder in folder_list:
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    if(kk==0):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder)
    print(data_path)
    
    os.chdir(data_path)
    list_file1 = glob.glob("proj_y_plt*")
    list_file = list_file1
    infile   = os.path.join(data_path, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))


    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))

    tot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
    tot_scalar = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    timestep = np.zeros(len(list_file))
    i=0
    for f in list_file:
        inputfile = os.path.join(data_path, f)
    
        file_size = get_folder_size(inputfile) / (1024 * 1024)
        
        if(file_size<1.):
            print('Skipping ', inputfile, '\n')
            continue
            
        if(os.path.isfile(inputfile)):
            print('Skipping ', inputfile, '\n')
            continue
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions)
        
        timestep[i] = ds.current_time.to('Myr')

        tot_mass_flux[i,:] = np.array(data['mass_outflow'])[:,0,:]*dx
        tot_scalar[i,:]     = np.array(data['scalar'])[:,0,:]
        i+=1  
        
    sign = zrange/np.abs(zrange)
    
    eta = np.sum(tot_mass_flux, axis=1)*sign/sfr
    eta_tot = (eta + eta[::-1])
    eta_tavg = np.average(eta_tot,axis=0)
    
    eta_16 = np.percentile(eta_tot, 16, axis=0)
    eta_84 = np.percentile(eta_tot, 84, axis=0)
    
    height = np.amax(zrange)/2.
    
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    
    eta_time_avg.append(eta_tavg[index])
    eta_16_avg.append(eta_16[index])
    eta_84_avg.append(eta_84[index])
    
    vol_index = (np.abs(zrange)<height) &  (np.abs(zrange)>0.15*kpc)
    
    mask1 = (np.abs(zrange)) < (height + 2*dz)
    mask2 = (np.abs(zrange)) > ( height - 2.*dz)
    mask = mask1*mask2

    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    
    sorted_time = np.asarray(list(map(timestep.__getitem__, indexes)))
    scalar_mass     = np.sum(tot_scalar[:,:,vol_index], axis=(1,2)) *dx * dz * Zunit

    sorted_scalar  = list(map(scalar_mass.__getitem__, indexes))
    sorted_time[sorted_time<1.] = 1.e-5
    inj_scalar = sigma_sfr  * Msun * sorted_time

    phi    =  1. - np.asarray(sorted_scalar)/(inj_scalar*Myr)

    phi = np.asarray(phi)

    phi_tavg.append(np.average(phi[1:]))
    phi_16.append(np.percentile(phi[1:], 16, axis=0))
    phi_84.append((np.percentile(phi[1:], 84, axis=0)))
    
    kk+=1


yerr_low = np.asarray(eta_time_avg) - np.asarray(eta_16_avg)
yerr_high = np.asarray(eta_84_avg) - np.asarray(eta_time_avg)
ax[0].errorbar(np.asarray(hSN), eta_time_avg, yerr= [yerr_low, yerr_high],\
            marker='*', markersize=12, ls=':', lw=2.0, capsize=5, color='navy')

yerr_low = np.asarray(phi_tavg) - np.asarray(phi_16)
yerr_high = np.asarray(phi_84) - np.asarray(phi_tavg)
ax[1].errorbar(np.asarray(hSN), phi_tavg,yerr=[yerr_low, yerr_high], \
            marker='v', markersize=12, ls=':', lw=2.0, color='steelblue')

ax[0].tick_params('x', bottom=True, top=True, labelbottom=False)

ax[1].tick_params('x', bottom=True, top=True, labelbottom=True)
ax[1].set_ylim(0.0,1.2)

plt.setp(ax[1], 'xlabel', r'$h_{\rm SN}$ [kpc]')
plt.setp(ax[0], 'ylabel', r'$\eta_M$')
plt.setp(ax[1], 'ylabel', r'$\phi$')

filename = fig_path + '/Paper/hSN_loading_fac.jpeg'
plt.savefig(filename, bbox_inches='tight')
print('Created file--', filename)