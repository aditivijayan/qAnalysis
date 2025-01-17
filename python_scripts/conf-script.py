
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



#-------------------Plot eta_Z------------------#

etaZ = np.sum(tot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)
etaZW = np.sum(wtot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)
etaZH = np.sum(htot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)


etaZ_tavg_tot = np.average(etaZ[tmask],axis=0)+ np.average(etaZ[tmask],axis=0)[::-1]
etaZW_tavg_tot = np.average(etaZW[tmask],axis=0)+ np.average(etaZW[tmask],axis=0)[::-1]
etaZH_tavg_tot = np.average(etaZH[tmask],axis=0)+ np.average(etaZH[tmask],axis=0)[::-1]

etaZ_16 = np.percentile(etaZ[tmask], 16, axis=0) + np.percentile(etaZ[tmask], 16, axis=0)[::-1]
etaZ_84 = np.percentile(etaZ[tmask], 84, axis=0) + np.percentile(etaZ[tmask], 84, axis=0)[::-1]


etaZ_int_tot = etaZ_tavg_tot - etaZW_tavg_tot - etaZH_tavg_tot

ax.plot(zrange/kpc, etaZ_tavg_tot, color='white', lw=lw)
ax.plot(zrange/kpc, etaZW_tavg_tot, \
                 color='forestgreen', ls='--', lw=lw)


ax.plot(zrange/kpc, etaZ_int_tot, color='magenta',\
              ls=((0, (3, 1, 1, 1))), lw=lw, label=label[kk][1])
ax.plot(zrange/kpc, etaZH_tavg_tot, color='indianred',\
              ls='-.', lw=lw)

# ax.fill_between(zrange/kpc, etaZ_16, etaZ_84, \
#                     color='grey', alpha=0.2)


ax.set_xlim(0.01, 4.0)

ax.set_ylim(-0.1, 1.18)

ax.set_xlabel('$z$ [kpc]')

    
plt.savefig("conf-etaZ.jpeg", dpi=160, bbox_inches='tight')
