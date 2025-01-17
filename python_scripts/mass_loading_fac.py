
import _init_
from constants import *
from set_path import *
from config import *
from functions import *




folder_list = [ 'SummitData/GasGravity/Production2pc/R4-h75-0.2Zsol',\
              'SummitData/GasGravity/Production2pc/R16-h300-Zsol/']

solar_met = ['SummitData/GasGravity/Production2pc/R4/', 'SummitData/GasGravity/Production2pc/R8/',\
             'SetonixData/R16/4pc/']

low_met = ['SummitData/GasGravity/Production2pc/R4-0.2Zsol/', 'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
           'SetonixData/R16-0.2Zsol/4pc/']
sigma_sfr = [ 0.000398107/yr_to_sec, 0.000398107/yr_to_sec, \
               6.e-5/yr_to_sec , 6.e-5/yr_to_sec , 1.58e-6/yr_to_sec, 1.58e-6/yr_to_sec] 
name = ['R4', 'R8', 'R16']

# sigma_sfr = [ 0.000398107/yr_to_sec, 1.58e-6/yr_to_sec] 
# name = ['R4-h75', 'R16-h300']


height = 0.5 * kpc
color =sns.color_palette("rocket",n_colors=3) 

sfr = np.asarray(sigma_sfr) * Msun * 1.e2


for folder, sigsfr, Zbg in zip(folder_list[4:], sigma_sfr[4:], Z_background):
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    
    os.chdir(data_path)
    list_file = glob.glob("proj_y_*")
    
    
    os.chdir(data_path)
    infile   = os.path.join(data_path, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    
    
    
    timestep = np.zeros((len(list_file)))
    tot_mass_flux = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_scalar_flux = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_scalar = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_rho = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    
    i = 0 
    for f in list_file:
        inputfile = os.path.join(data_path, f)
        
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)

        timestep[i] = ds.current_time.to('Myr')
        tot_mass_flux[i,:] = np.array(data['mass_outflow'])[:,0,:]
        
        i+=1
    
    eta = (np.sum(tot_mass_flux, axis=1))

    eta_tavg = np.average(eta,axis=0)*sign
    
    eta_16 = np.percentile(eta, 16, axis=0)
    eta_84 = np.percentile(eta, 84, axis=0)
    burst = (eta_84 - eta_16)/eta_tavg
    burst_avg = (burst + burst[::-1])/2.
    
    ax.plot(zrange/kpc, eta_tavg, color=color[kk], ls =ls[kk], label=label[kk])
    ax.fill_between(zrange/kpc, eta_16, eta_84, color='gray', alpha=0.2)
    kk+=1
    
    
ax.set_ylabel(r'$\eta$')
ax.set_xlim(0.5, np.amax(zrange)/kpc)
ax.legend()
plt.setp(ax, ('ylim'), (0.0, 10.))


fig, ax = plt.subplots(1, 2, gridspec_kw = {'wspace':0.0, 'hspace':0.0},figsize=(14, 6))

ax[0].plot(np.asarray(harray)/kpc, phi_avg[0], 'o', ls='-', color='mediumvioletred', label=r'$Z_{\odot}$')
ax[0].fill_between(np.asarray(harray)/kpc, phi_16[0], phi_84[0], ls='-', color='gray', alpha=0.2)

ax[1].plot(np.asarray(harray)/kpc, phi_avg[1], 'o', ls='-', color='forestgreen', label=r'$0.2Z_{\odot}$')
ax[1].fill_between(np.asarray(harray)/kpc, phi_16[1], phi_84[1], ls='-', color='gray', alpha=0.2)

plt.setp(ax, 'ylim', (-0.05, 0.5))

ax[0].tick_params(axis='y',  right=True, left=True)
ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=True, left=True)
ax[0].legend()
ax[1].legend()

filename = fig_path + '/metal_loading_R16.jpeg'
plt.savefig(filename, bbox_inches='tight')
print('Created file--', filename)


