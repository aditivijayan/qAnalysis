
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns
from scipy.ndimage import uniform_filter1d


plt.rcParams['font.size']=36
plt.rcParams['axes.linewidth']=0.5
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


# folder_list = ['SummitData/GasGravity/Production2pc/R4-0.2Zsol/',  \
            #    'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
            #   'SetonixRuns/R16-0.2Zsol/4pc/Redo/']
# label1 =[r'$\Sigma50$-Z$0.2$-H$150$', r'$\Sigma13$-Z$0.2$-H$150$',  r'$\Sigma2.5$-Z$0.2$-H$150$']

folder_list = ['SummitData/GasGravity/Production2pc/R4/',  \
               'SummitData/GasGravity/Production2pc/R8/',\
              'SetonixRuns/R16/4pc/Redo/']
label = [r'$\Sigma50$-Z$1$-H$150$', r'$\Sigma13$-Z$1$-H$150$', r'$\Sigma2.5$-Z$1$-H$150$']
# label2 =['', '', r'$\Sigma2.5$-Z$0.2$-H$150$']


sigma_sfr = [ 0.000398107/yr_to_sec,\
               6.e-5/yr_to_sec ,\
               1.58e-6/yr_to_sec]
kk= 0
color =sns.color_palette("cool",n_colors=3) 
color = ['magenta', 'darkslategrey', 'darkorange']

fig, ax = plt.subplots(2, 1, gridspec_kw = {'wspace':0.01, 'hspace':0.01},figsize=(24, 16))

Zunit = Msun/1.e3

ls = ['-', '-.', '--']
lw = 3.0

for folder in folder_list:
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    if(kk==2):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder)
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

    height = np.amax(zrange)/2.
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    
    timestep = np.zeros(len(list_file))
    scalar_flux_xz = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_mass_flux_xz  = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    tot_mass_flux  = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))

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
        tot_mass_flux_xz[i,:] = np.array(data['mass_outflow'])[:,0,:]*dx*yr_to_sec/Msun #/1.e-2
        tot_mass_flux[i,:] = np.array(data['mass_outflow'])[:,0,:]*dx
        scalar_flux_xz[i,:] = (np.array(data['scalar_outflow'])[:,0,:])*Zunit*dx*yr_to_sec/Msun #/1.e-4
        
        i+=1  
        
        
    sign = zrange/np.abs(zrange)
    
    height = np.amax(zrange)/2.
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    
    mass_flux_z = np.sum(tot_mass_flux, axis=1)*sign
    
    mask1 = (np.abs(zrange)) < (height + 2*dz)
    mask2 = (np.abs(zrange)) > ( height - 2.*dz)
    mask = mask1*mask2
    
    scalar_flux_z = np.sum(scalar_flux_xz , axis=1)*sign
    scalar_fluxz_avg   = scalar_flux_z + scalar_flux_z[::-1]
    scalar_flux = scalar_fluxz_avg[:, index]
    
    mass_flux_z = np.sum(tot_mass_flux_xz, axis=1)*sign

    mass_flux_z_avg =  (mass_flux_z +  mass_flux_z[::-1])
    mass_flux   = mass_flux_z_avg[:,index]
    
    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = np.asarray(list(map(timestep.__getitem__, indexes)))
    
    
    sorted_scalar_flux  = list(map(scalar_flux.__getitem__, indexes))
    sorted_mass_flux  = list(map(mass_flux.__getitem__, indexes))

    window_size = 5
    smoothed_array = uniform_filter1d(sorted_mass_flux, size=5)
    ax[0].plot(sorted_time, smoothed_array, color=color[kk], lw=lw, ls =ls[kk], label=label[kk])

    smoothed_array = uniform_filter1d(sorted_scalar_flux, size=5)
    ax[1].plot(sorted_time,  smoothed_array, color=color[kk], lw=lw, ls =ls[kk], label=label[kk])
    lw += 0.5
    kk+=1

ax[0].legend(frameon=True, ncol=3, fontsize=34, loc='upper left')
ax[0].set_ylabel(r'$\dot{M}$' + '\n' + r'$[M_{\odot}$ yr$^{-1}]$', fontsize =40)
ax[0].tick_params(axis='y',  right=True, left=True)
ax[0].tick_params(axis='x',  top=True, bottom=True, labelbottom=False)
ax[0].set_ylim(2.e-4, 1.)
ax[1].text(0.9, 0.85, 'Mass', transform=ax[0].transAxes, fontsize=44)

ax[1].set_ylabel(r'$\dot{M}_Z$' + '\n' + r'$[M_{\odot}$ yr$^{-1}]$', fontsize =40)
ax[1].tick_params(axis='x',  top=True, bottom=True)
ax[1].tick_params(axis='y',  left=True, right=True)
ax[1].set_ylim(8.e-8, 8.e-4)
ax[1].text(0.9, 0.85, 'Metal', transform=ax[1].transAxes, fontsize=44)

plt.setp(ax, 'yscale', 'log')
plt.setp(ax, 'xlim', (0.0, 300.))

plt.setp(ax[1], 'xlabel', 't [Myr]')
filename = fig_path + '/Paper/mass_outflow.jpeg'
plt.savefig(filename, bbox_inches='tight', dpi=160)
print('Created file--', filename)
