import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns

plt.rcParams['font.size']=28
plt.rcParams['axes.linewidth']=2.
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
kk= 0
color =sns.color_palette("husl",n_colors=len(folder_list)) 

sfr = np.asarray(sigma_sfr) * Msun * 1.e2

kk= 0
lw = 5.5
Zunit = Msun/1.e3

color =sns.color_palette("husl",n_colors=3)

ETA_tavg = []
ETAW_tavg = []
ETAH_tavg = [] 

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
    
    timestep = np.zeros(len(list_file))
    tot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
    wtot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
    htot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
    tot_scalar = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    
    
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

        tot_mass_flux[i,:] = np.array(data['mass_outflow'])[:,0,:]*dx
        wtot_mass_flux[i,:] = np.array(data['warm_mass_outflow'])[:,0,:]*dx
        htot_mass_flux[i,:] = np.array(data['hot_mass_outflow'])[:,0,:]*dx
        tot_scalar[i,:]     = np.array(data['scalar'])[:,0,:]
        i+=1  
        
    sign = zrange/np.abs(zrange)
    
    tmask = (timestep>75.)
    eta = np.sum(tot_mass_flux[tmask,:], axis=1)*sign/sfr
    etaW = np.sum(wtot_mass_flux[tmask,:], axis=1)*sign/sfr
    etaH = np.sum(htot_mass_flux[tmask,:], axis=1)*sign/sfr

    eta_tavg = np.average(eta,axis=0)
    etaW_tavg = np.average(etaW,axis=0)
    etaH_tavg = np.average(etaH,axis=0)
    
    eta_16 = np.percentile(eta, 16, axis=0)
    eta_84 = np.percentile(eta, 84, axis=0)
    
    eta_16_tot = (eta_16 + eta_16[::-1])
    eta_84_tot = (eta_84 + eta_84[::-1])

    height = np.amax(zrange)
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))

    print("The mass loading for ", folder, " is -->%.3f"%eta_tavg[index])
    print("The 16th percentile is -->%.3f"%eta_16_tot[index])
    print("The 84th percentile is -->%.3f"%eta_84_tot[index])
    
    ETA_tavg.append(eta_tavg   + eta_tavg[::-1])
    ETAW_tavg.append(etaW_tavg + etaW_tavg[::-1])
    ETAH_tavg.append(etaH_tavg + etaH_tavg[::-1])


    #----------Get Phi------------#
    
    scalar_mass = np.average(tot_scalar, axis=0) *dx * dz * Zunit
    
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    vol_index = (np.abs(zrange)<height) &  (np.abs(zrange)>0.15*kpc)

    mask1 = (np.abs(zrange)) < (height + 2*dz)
    mask2 = (np.abs(zrange)) > ( height - 2.*dz)
    mask = mask1*mask2

    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = np.asarray(list(map(timestep.__getitem__, indexes)))


    scalar_mass     = np.sum(tot_scalar[:,:,vol_index], axis=(1,2)) *dx * dz * Zunit


    sorted_scalar  = list(map(scalar_mass.__getitem__, indexes))
    inj_scalar = sigma_sfr  * Msun * sorted_time

    phi    =  1. - np.asarray(sorted_scalar)/(inj_scalar*Myr)

    phi = np.asarray(phi[sorted_time>75.])

    phi_avg=np.average(phi)
    phi_16=np.percentile(phi, 16, axis=0)
    phi_84=(np.percentile(phi, 84, axis=0))
    
    print("The metal loading for ", folder, " is -->%.3f"%phi_avg)
    print("The 16th percentile is --> %.3f"%phi_16)
    print("The 84th percentile is -->%.3f"%phi_84)
    
    kk+=1
    

fig, ax = plt.subplots(3, 2, gridspec_kw = {'wspace':0.00, 'hspace':0.00},figsize=(16, 20))
eta_tavg = np.asarray(ETA_tavg)
etaw = np.asarray(ETAW_tavg)
etah = np.asarray(ETAH_tavg)
eta_int = eta_tavg - etaw - etah
lw=3.0
kk=0

label = ['Total', r'$T< 2 \times 10^4$ K', r'$2 \times 10^4$ K $<T<10^6$ K', r'$T>10^6$ K']
text = [r'$150$ pc', r'$300$ pc',  r'$1$ kpc', r'$1.5$ kpc', r'$2$ kpc', ]
for ii in range(3):
    for jj in range(2):

        if(kk==5):
            ax[ii,jj].plot(zrange/np.amax(zrange), 0.0*eta_tavg[0], lw=lw, color='black', label=label[0])
            ax[ii,jj].plot(zrange/np.amax(zrange), 0.0*etaw[0], ls='--', lw=lw, color='forestgreen', label=label[1])
            ax[ii,jj].plot(zrange/np.amax(zrange), 0.0*etaw[0], ls=((0, (3, 1, 1, 1))), lw=lw, color='magenta', label=label[2])
            ax[ii,jj].plot(zrange/np.amax(zrange), 0.0*etaH[0], ls='-.', lw=lw+2, color='darkslategrey', label=label[3])

            
        else:

            ax[ii,jj].plot(zrange/np.amax(zrange), eta_tavg[kk], lw=lw, color='black')
            ax[ii,jj].plot(zrange/np.amax(zrange), etaw[kk], ls='--', lw=lw, color='forestgreen')
            ax[ii,jj].plot(zrange/np.amax(zrange), eta_int[kk], ls=((0, (3, 1, 1, 1))), lw=lw, color='magenta')
            ax[ii,jj].plot(zrange/np.amax(zrange), etah[kk], ls='-.', lw=lw+2, color='darkslategrey')
            ax[ii,jj].text(0.7, 0.82, text[kk], transform=ax[ii,jj].transAxes) 
            kk+=1

for ii in range(3):

    ax[ii,0].tick_params(axis='y', which='both', right=True, left=True, labelleft=True)
    ax[ii,1].tick_params(axis='y', which='both', right=True, left=True, labelleft=False, labelright=True)
    ax[ii,jj].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)

ax[0,0].tick_params(axis='y',  right=True, left=True, labelleft=True)

ax[2,1].legend(loc='lower right', frameon=False, fontsize=30)

plt.setp(ax, 'xlim', (0.0,1.))
plt.setp(ax[:,0], 'xticks', [0.0, 0.5, 1.0])
plt.setp(ax[1,1], 'xticks', [0.0, 0.5, 1.0])
plt.setp(ax[0,1], 'xticks', [0.0, 0.5, 1.0])
plt.setp(ax, 'yscale', 'log')
plt.setp(ax, 'ylim', (2.e-3, 8.e1))
plt.setp(ax[:,0], 'ylabel', r'$\eta_{M}$')

ax[2,1].spines['right'].set_visible(False)
ax[2,1].spines['bottom'].set_visible(False)
ax[2,1].tick_params(axis='y', which='both',  labelright=False, labelleft=False, left=False, right=False)
ax[2,1].tick_params(axis='x',  which='both', labelbottom=False, bottom=False, top=False)

ax[2,0].set_xlabel(r'$z/L_z$')
ax[1,1].set_xlabel(r'$z/L_z$')

ax[1,1].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True)
ax[1,0].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
ax[0,0].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
ax[2,0].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True)

filename = fig_path + '/Paper/loading_fac_hSN.jpeg'
plt.savefig(filename, bbox_inches='tight')
print('Created file--', filename)
