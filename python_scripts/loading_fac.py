
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns
import gc 

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



# folder_list = ['SummitData/GasGravity/Production2pc/R4/',  \
#                'SummitData/GasGravity/Production2pc/R8/',\
#               'SetonixRuns/R16/4pc/Redo/']
# name = [r'$\Sigma50$-Z$1$-H$150$', r'$\Sigma13$-Z$1$-H$150$', r'$\Sigma2.5$-Z$1$-H$150$']
# filename = fig_path + '/Paper/loading_fac_Zsol.jpeg'
# tfinal = np.asarray([75., 75., 75.])

folder_list = ['SummitData/GasGravity/Production2pc/R4-0.2Zsol/',  \
               'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
              'SetonixRuns/R16-0.2Zsol/4pc/Redo']
name = [r'$\Sigma50$-Z$0.2$-H$150$', r'$\Sigma13$-Z$0.2$-H$150$', r'$\Sigma2.5$-Z$0.2$-H$150$']
filename = fig_path + '/Paper/loading_fac_0.2Zsol.jpeg'
tfinal = np.asarray([75., 75., 75.])

sigma_sfr = [ 0.000398107/yr_to_sec,\
               6.e-5/yr_to_sec ,\
               1.58e-6/yr_to_sec]
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

fig, ax = plt.subplots(4, 3, gridspec_kw = {'wspace':0., 'hspace':0.04},figsize=(32, 32))

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
    tot_rho = np.zeros((len(list_file), xrange.shape[0], zrange.shape[0]))
    
    
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
        tot_rho[i,:] = np.array(data['rho'])[:,0,:]

        i+=1  
        
    sign = zrange/np.abs(zrange)
    
    tmask = (timestep>tfinal[kk])
    
    height = np.amax(zrange)
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    
    eta = np.sum(tot_mass_flux, axis=1)*sign/sfr[kk]
    etaW = np.sum(wtot_mass_flux, axis=1)*sign/sfr[kk]
    etaH = np.sum(htot_mass_flux, axis=1)*sign/sfr[kk]
    
    eta_tavg_tot = np.average(eta[tmask],axis=0)+ np.average(eta[tmask],axis=0)[::-1]
    etaW_tavg_tot = np.average(etaW[tmask],axis=0)+ np.average(etaW[tmask],axis=0)[::-1]
    etaH_tavg_tot = np.average(etaH[tmask],axis=0)+ np.average(etaH[tmask],axis=0)[::-1]

    eta_int_tot = eta_tavg_tot - etaW_tavg_tot - etaH_tavg_tot
    
    eta_16 = np.percentile(eta[tmask], 16, axis=0) + np.percentile(eta[tmask], 16, axis=0)[::-1]
    eta_84 = np.percentile(eta[tmask], 84, axis=0) + np.percentile(eta[tmask], 84, axis=0)[::-1]
    
    ax[0,kk].plot(zrange/np.amax(zrange), eta_tavg_tot, color='black', lw=lw,  label=label[kk][0])
    if(kk==0):
        ax[0,0].plot(zrange/np.amax(zrange), etaW_tavg_tot, \
                     color='forestgreen', ls='--', lw=lw, label=r'$T<2\times 10^4$ K')
    else:
        ax[0,kk].plot(zrange/np.amax(zrange), etaW_tavg_tot, color='forestgreen', ls='--', lw=lw, label='')
        
    ax[0,kk].plot(zrange/np.amax(zrange), eta_int_tot, color='magenta',\
                  ls=((0, (3, 1, 1, 1))), lw=lw, label=label[kk][1])
    ax[0,kk].plot(zrange/np.amax(zrange), etaH_tavg_tot, color='darkslategrey',\
                  ls='-.', lw=lw+4, label=label[kk][2])

    ax[0,kk].fill_between(zrange/np.amax(zrange), eta_16, eta_84, \
                        color='cornflowerblue', alpha=0.2)
    
    
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

    ax[1,kk].plot(zrange/np.amax(zrange), etaZ_tavg_tot, color='black', lw=lw)
    if(kk==0):
        ax[1,0].plot(zrange/np.amax(zrange), etaZW_tavg_tot, \
                     color='forestgreen', ls='--', lw=lw)
    else:
        ax[1,kk].plot(zrange/np.amax(zrange), etaZW_tavg_tot, color='forestgreen', ls='--', lw=lw, label='')
        
    ax[1,kk].plot(zrange/np.amax(zrange), etaZ_int_tot, color='magenta',\
                  ls=((0, (3, 1, 1, 1))), lw=lw, label=label[kk][1])
    ax[1,kk].plot(zrange/np.amax(zrange), etaZH_tavg_tot, color='darkslategrey',\
                  ls='-.', lw=lw+4)

    ax[1,kk].fill_between(zrange/np.amax(zrange), etaZ_16, etaZ_84, \
                        color='cornflowerblue', alpha=0.2)
    
    #----------Get Phi------------#
    scalar_mass = np.average(tot_scalar, axis=0) *dx * dz * Zunit
    harray = np.linspace(1.*kpc, np.amax(zrange), 4)    
    
    phi = []
    X, T, Z = np.meshgrid(xrange, timestep, zrange)
    
    Sign = (Z/np.abs(Z))
    sign = zrange/np.abs(zrange)
    
    
    
#     scal_time_averaged = np.average(np.ma.array(tot_scalar, mask=~Tmask), axis=0).data 
#     rho_time_averaged = np.average(np.ma.array(tot_rho, mask=~Tmask), axis=0).data 
    phi = np.zeros((timestep[tmask].shape[0], harray.shape[0]))
    jj=0
    for height in harray:
       
        vol_index = (np.abs(zrange)<height) 
        Volmask = (T>tfinal[kk]) * (np.abs(Z)<height)
        Tmask = (T>tfinal[kk])
        index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
         
        scalar_mass  = np.sum(np.ma.array(tot_scalar, mask=~Volmask), axis=(1,2)).data  * dx * dz * Zunit
        rho = np.sum(np.ma.array(tot_rho, mask=~Tmask), axis=(1,2)).data * dx * dz 
        

        mass_flux_z = np.sum(np.ma.array(tot_mass_flux, mask=~Volmask), axis=1).data 
        rev_mass_flux = mass_flux_z[:,]*sign
        rev_mass_flux = rev_mass_flux[:,::-1]
        mass_flux_tot = mass_flux_z + rev_mass_flux

        scal_flux_z = np.sum(np.ma.array(tot_scal_flux, mask=~Tmask), axis=1).data  * Zunit
        rev_scal_flux =  scal_flux_z[:,]*sign
        rev_scal_flux = rev_scal_flux[:,::-1]
        scal_flux_tot = scal_flux_z + rev_scal_flux
        
        
        print(scalar_mass.shape)
        Zinj = scalar_mass[tmask]/ rho[tmask]
        
        Mz   = scal_flux_tot[:,index]
        Mdot = mass_flux_tot[:,index]
    
        phi11 = ( (Mz[tmask] - Mdot[tmask]*Zinj)/sigma_sfr[kk]/Msun)
        phi[:,jj] = phi11
        jj+=1
        
    phi_avg = np.average(phi, axis=0)
    phi_16 = np.percentile(phi, 16, axis=0)
    phi_84 = np.percentile(phi, 84, axis=0)
        
    ax[3,kk].plot(harray/np.amax(zrange), np.asarray(phi_avg), marker='*', ls='-', \
            label=label[kk], markersize=24, lw=lw+1, color='darkslategrey')
        
    ax[3,kk].fill_between(harray/np.amax(zrange), phi_16, phi_84, \
                        color='cornflowerblue', alpha=0.2)
    
    ax[3,kk].tick_params(axis='y',  right=True, left=True, labelleft=True)
    

    #-----------------------------------------#
    ax[0,kk].legend(frameon=False, fontsize=38, loc='upper right')
   
    gc.collect() 
    kk+=1
    
 
    
#---------------Plot Energy Loading------------------------------#

kk=0

norm  = 1.e51 * np.asarray(sigma_sfr)

for folder in folder_list:
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder, 'Eloading')
    if(kk==2):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder, 'Eloading')
    print(data_path)
    
    os.chdir(data_path)
    list_file1 = glob.glob('*.h5')
    list_file = list_file1
    
    infile   = os.path.join(data_path.replace('Eloading', ''), 'metal_uniform.in')
    # print(infile)
    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))


    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    dVol = dx * dy * dz

    eout_rate = np.zeros((len(list_file),zrange.shape[0]))
    ethout_rate = np.zeros((len(list_file),zrange.shape[0]))
    timestep  = np.zeros(len(list_file))
    

    i=0
    for f in list_file:
        inputfile = os.path.join(data_path, f)
        hf = h5py.File(inputfile ,'r')
        eout_rate[i,:] = np.array(hf.get("EnergyOutflowRate"))
        ethout_rate[i,:] = np.array(hf.get("ThEnergyOutflowRate"))
        timestep[i] = np.array(hf.get("Timestep")) 

        i+=1  
        
    
    sign = zrange/np.abs(zrange)
    etaE = (eout_rate)*sign/norm[kk]
    etaEth = (ethout_rate)*sign/norm[kk]

    tmask = (timestep>tfinal[kk])
    etaE = etaE[tmask]
    etaEth = etaEth[tmask]

    eta_tavg = np.average(etaE,axis=0)
    etaTh_tavg = np.average(etaEth,axis=0)
    
    eta_tavg_tot = eta_tavg + eta_tavg[::-1]
    etaTh_tavg_tot = etaTh_tavg + etaTh_tavg[::-1]
    
    eta_16 = np.percentile(etaE, 16, axis=0)
    eta_84 = np.percentile(etaE, 84, axis=0)
    
    eta_16_tot = eta_16 + eta_16[::-1]
    eta_84_tot = eta_84 + eta_84[::-1]
    
    
    ax[2,kk].plot(zrange/np.amax(zrange), eta_tavg_tot, color='seagreen', lw=lw, label=label1[0][kk])
    ax[2, kk].plot(zrange/np.amax(zrange), etaTh_tavg_tot, color='darkorange',ls='--', lw=lw, label=label1[1][kk])
    
    ax[2, kk].fill_between(zrange/np.amax(zrange), eta_16_tot, eta_84_tot,\
                           ls='-', color='cornflowerblue', alpha=0.2)
    ax[2,0].legend(loc='upper right', frameon=False, ncol=2, fontsize=40)
    ax[2,1].legend(loc='upper right', frameon=False, ncol=2, fontsize=40)
    kk+=1
    
for ii in range(3):
    for jj in range(3):
        ax[ii,jj].tick_params(axis='y', which='both',  right=True, left=True)
        ax[0,ii].set_title(name[ii], fontsize=40, pad=40)

ax[0,0].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[0,0].set_ylabel(r'$\eta_{\rm M}$', fontsize =40, labelpad=30)
ax[1,0].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
plt.setp(ax[0,:], 'yscale', 'log')
plt.setp(ax[1,:], 'yscale', 'log')

for ii in range(4):
    ax[ii,0].tick_params(axis='y',  left=True, right=True)
    ax[ii,1].tick_params(axis='y',  left=True, right=True, labelleft=False)
    ax[ii,2].tick_params(axis='y',  left=True, right=True, labelleft=False, labelright=True)


ax[1,2].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)

ax[0,1].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[0,2].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[1,1].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[1,2].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[2,0].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[2,1].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)
ax[2,2].tick_params(axis='x',  top=True, bottom=True, labelbottom = False)

ax[-1,0].tick_params(axis='x',  top=True, bottom=True)
ax[-1,1].tick_params(axis='x',  top=True, bottom=True)
ax[-1,2].tick_params(axis='x',  top=True, bottom=True)

ax[1,0].set_ylabel(r'$\eta_Z$', fontsize =40, labelpad=30)
ax[3,0].set_ylabel(r'$\phi$', fontsize =40, labelpad=30)


ax[2,0].set_ylabel(r'$\eta_{\rm E}$', fontsize =40, labelpad=30)
plt.setp(ax[2,:], 'ylim', (1.e-2, 2.5))
plt.setp(ax[2,:], 'yscale', 'log')


plt.setp(ax, ('xlim'), (0.15/np.amax(zrange/kpc), 1.05))

# ax[0,0].set_ylim(2.e-3, 2.)
# ax[0,1].set_ylim(0.05, 8.1)

ax[0,0].set_ylim(1.e-2, 1.e2)
ax[0,1].set_ylim(1.e-2, 1.e2)
ax[0,2].set_ylim(1.e-2, 1.e2)

plt.setp(ax[1,:], 'yscale', 'log')
plt.setp(ax[1,:], 'ylim', (1.e-2,2.e1))
plt.setp(ax[2,:], 'ylim', (2.e-2,2.e1))
plt.setp(ax[3,:], 'ylim', (-0.1,1.1))
plt.setp(ax, 'xlabel', '')
plt.setp(ax[-1,:], 'xlabel', r'$z/L_z$')

plt.savefig(filename, bbox_inches='tight')
print('Created file--', filename)
