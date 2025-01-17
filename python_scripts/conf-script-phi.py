
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


plt.style.use('dark_background')


folder = 'SummitData/GasGravity/Production2pc/R8/'
filename = fig_path + '/Paper/loading_fac_0.2Zsol.jpeg'

sigma_sfr = [  6.e-5/yr_to_sec]
Zunit = Msun/1.e3


kk= 0
lw = 5.5

color =sns.color_palette("husl",n_colors=3)

marker = ['*', 'o', 'v']
ls = ['-', '-.', '--']

label = [['Total', '', ''],\
         ['', r'$2\times 10^4$ K$<T<10^6$ K', ''],\
         ['', '', '$T>10^6$ K']]

label1 = [['Total', '', ''],
          ['', 'Thermal', ''],
          ['', '', '']]

sfr = np.asarray(sigma_sfr) * Msun * 1.e2

fig, ax = plt.subplots(1, 1, gridspec_kw = {'wspace':0., 'hspace':0.04},figsize=(16, 10))


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


dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))

timestep = np.zeros(len(list_file))

tot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
wtot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
htot_mass_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))

tot_scal_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
wtot_scal_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
htot_scal_flux = np.zeros((len(list_file),xrange.shape[0], zrange.shape[0]))
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

    tot_scal_flux[i,:] = np.array(data['scalar_outflow'])[:,0,:]*dx
    wtot_scal_flux[i,:] = np.array(data['warm_scalar_outflow'])[:,0,:]*dx
    htot_scal_flux[i,:] = np.array(data['hot_scalar_outflow'])[:,0,:]*dx

    tot_scalar[i,:]     = np.array(data['scalar'])[:,0,:]
    i+=1  

sign = zrange/np.abs(zrange)

tmask = (timestep>0.)

height = np.amax(zrange)
index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))


# #----------Get Phi------------#
scalar_mass = np.average(tot_scalar, axis=0) *dx * dz * Zunit
harray = np.linspace(1.*kpc, np.amax(zrange), 4)    

phi_avg = []
phi_16 = []
phi_84 = []


for height in harray:

    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    vol_index = (np.abs(zrange)<height)  #&  (np.abs(zrange)>0.15*kpc)

    mask1 = (np.abs(zrange)) < (height + 2*dz)
    mask2 = (np.abs(zrange)) > ( height - 2.*dz)
    mask = mask1*mask2

    indexes = list(range(len(timestep)))
    indexes.sort(key=timestep.__getitem__)
    sorted_time = np.asarray(list(map(timestep.__getitem__, indexes)))


    scalar_mass     = np.sum(tot_scalar[:,:,vol_index], axis=(1,2)) *dx * dz * Zunit


    sorted_scalar  = list(map(scalar_mass.__getitem__, indexes))
    sorted_time[sorted_time<1.] = 1.e-5
    inj_scalar = sigma_sfr[kk]  * Msun * sorted_time

    phi    =  1. - np.asarray(sorted_scalar)/(inj_scalar*Myr)

    phi = np.asarray(phi)
    phi = phi[sorted_time>75.]

    phi_avg.append(np.average(phi))
    phi_16.append(np.percentile(phi, 16, axis=0))
    phi_84.append((np.percentile(phi, 84, axis=0)))
    print(phi_avg)

ax.plot(harray/kpc, np.asarray(phi_avg), marker='*', ls='-', \
        label=label[kk], markersize=24, lw=lw+1, color='white')

ax.fill_between(harray/kpc, phi_16, phi_84, \
                    color='cornflowerblue', alpha=0.2)

ax.tick_params(axis='y',  right=True, left=True, labelleft=True)


# #-----------------------------------------#
# ax[0,kk].legend(frameon=False, fontsize=38, loc='upper right')


# kk+=1
ax.set_xlim(0.01, 4.0)

ax.set_xlabel('$z$')

ax.set_xlim(0.5, 4.1)
ax.set_xticks([1.0, 2.0, 3.0, 4.0])

ax.set_ylim(-0.1, 1.18)

ax.set_xlabel('$z$ [kpc]')

    
plt.savefig("conf-phi.jpeg", dpi=160, bbox_inches='tight')
