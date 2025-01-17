
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import seaborn as sns


# folder_list = ['SummitData/GasGravity/Production2pc/R4/',  \
#                'SummitData/GasGravity/Production2pc/R8/',\
#               'SetonixRuns/R16/4pc/Redo/']
# name = [r'$4\Sigma_{\rm Fid}$', r'$\Sigma_{\rm Fid}$', r'0.2$\Sigma_{\rm Fid}$']
# output_file = '/Paper/solar_met.dat'
# tfinal = np.asarray([75., 75., 75.])

# folder_list = ['SummitData/GasGravity/Production2pc/R4-0.2Zsol/',  \
#                'SummitData/GasGravity/Production2pc/R8-0.2Zsol/',\
#               'SetonixRuns/R16-0.2Zsol/4pc/Redo/']
# name = [r'$4\Sigma_{\rm Fid}^*$', r'$\Sigma_{\rm Fid}^*$', r'0.2$\Sigma_{\rm Fid}^*$']
# output_file = '/Paper/sub-solar_met.dat'
# tfinal = np.asarray([75., 75., 75.])

# sigma_sfr = [ 0.000398107/yr_to_sec,\
#                6.e-5/yr_to_sec ,\
#                1.58e-6/yr_to_sec]


folder_list = ['SetonixRuns/R16/4pc/Redo/',\
               'SummitData/GasGravity/Production2pc/R16-h300-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h1kpc-Zsol',\
              'SummitData/GasGravity/Production2pc/R16-h1.5kpc-Zsol',\
               'SummitData/GasGravity/Production2pc/R16-h2kpc-Zsol']

sigma_sfr = [1.58e-6/yr_to_sec]*len(folder_list)
tfinal = [75.] * len(folder_list)
output_file = '/Paper/solar_met_hSN.dat'

Zunit = Msun/1.e3
 

############-----Mass and Metal Loading Factors----------############

kk=0

sfr = np.asarray(sigma_sfr) * Msun * 1.e2

file = open(fig_path + output_file, 'w')

for folder in folder_list:
    # data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder)
    # if(kk==2):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder)

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
    
    height = np.amax(zrange)
    index = min(range(zrange.shape[0]), key=lambda i: abs(zrange[i] - height))
    
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

    file.write("Starting analysis for" +  folder + " \n\n")
    file.write("The mass loading is -->%.3e\n"%eta_tavg_tot[index])
    file.write("The 16th percentile is -->%.3e\n"%eta_16[index])
    file.write("The 84th percentile is -->%.3e\n\n"%eta_84[index])
    file.write("The burstiness is -->%.3e\n\n"%((eta_84[index]-eta_16[index])/eta_tavg_tot[index]))
    
    file.write("The WARM mass loading is -->%.3e\n"%etaW_tavg_tot[index])
    file.write("The HOT mass loading for is -->%.3e\n"%etaH_tavg_tot[index])

    eta_int = eta_tavg_tot - etaW_tavg_tot - etaH_tavg_tot
    file.write("The INT mass loading is -->%.3e\n\n"%eta_int[index])
    
    #----------Get Eta_Z------------#
    etaZ = np.sum(tot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)
    etaZW = np.sum(wtot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)
    etaZH = np.sum(htot_scal_flux, axis=1)*sign*Zunit/(sigma_sfr[kk]  * Msun)

    
    etaZ_tavg_tot = np.average(etaZ[tmask],axis=0)+ np.average(etaZ[tmask],axis=0)[::-1]
    etaZW_tavg_tot = np.average(etaZW[tmask],axis=0)+ np.average(etaZW[tmask],axis=0)[::-1]
    etaZH_tavg_tot = np.average(etaZH[tmask],axis=0)+ np.average(etaZH[tmask],axis=0)[::-1]
    
    etaZ_16 = np.percentile(etaZ[tmask], 16, axis=0) + np.percentile(etaZ[tmask], 16, axis=0)[::-1]
    etaZ_84 = np.percentile(etaZ[tmask], 84, axis=0) + np.percentile(etaZ[tmask], 84, axis=0)[::-1]
    
    
    etaZ_int_tot = etaZ_tavg_tot - etaZW_tavg_tot - etaZH_tavg_tot
    
    file.write("The metal loading factor is -->%.3e\n"%etaZ_tavg_tot[index])
    file.write("The 16th percentile is -->%.3e\n"%etaZ_16[index])
    file.write("The 84th percentile is -->%.3e\n\n"%etaZ_84[index])
    
    file.write("The WARM metal loading factor is -->%.3e\n"%etaZW_tavg_tot[index])
    file.write("The HOT mass loading factor for is -->%.3e\n"%etaZH_tavg_tot[index])

    etaZ_int = etaZ_tavg_tot - etaZW_tavg_tot - etaZH_tavg_tot
    file.write("The INT metal loading factor is -->%.3e\n\n"%etaZ_int[index])
    


    #----------Get Phi------------#
    vol_index = (np.abs(zrange)<height) 
    X, T, Z = np.meshgrid(xrange, timestep, zrange)
    Tmask = (T>tfinal[kk])
    sign = zrange/np.abs(zrange)
    
    mass_flux_time_averaged = np.average(np.ma.array(tot_mass_flux, mask=~Tmask), axis=0).data 
    mass_flux_z   = np.sum(mass_flux_time_averaged, axis=0) 
    
    mass1_flux_z = mass_flux_z*sign
    mass_flux_tot = mass1_flux_z + mass1_flux_z[::-1]   
    
    scal_flux_time_averaged = np.average(np.ma.array(tot_scal_flux, mask=~Tmask), axis=0).data 
    scal_flux_z   = np.sum(scal_flux_time_averaged, axis=0)  * Zunit
    scal1_flux_z = scal_flux_z*sign
    scal_flux_tot = scal1_flux_z*sign + scal1_flux_z[::-1]
    
    
    scal_time_averaged = np.average(np.ma.array(tot_scalar, mask=~Tmask), axis=0).data 
    scalar_mass  = np.sum(scal_time_averaged[:,vol_index], axis=(0,1)) * dx * dz * Zunit
    
    rho_time_averaged = np.average(np.ma.array(tot_rho, mask=~Tmask), axis=0).data 
    rho = np.sum(rho_time_averaged[:,vol_index], axis=(0,1)) * dx * dz 
    
    Zinj = scalar_mass/ rho
    
    Mz   = scal_flux_tot[index]
    Mdot = mass_flux_tot[index]
   
    phi  = (Mz - Mdot*Zinj)/sigma_sfr[kk]/Msun
    print(phi)
    
    file.write("Phi is -->%.3e\n"%phi)
    # file.write("The 16th percentile is --> %.3e\n"%phi_16)
    # file.write("The 84th percentile is -->%.3e\n\n"%phi_84)
    
    kk+=1
    
##########---------------Energy Loading Factor-------------#############

kk=0
norm  = 1.e51 * np.asarray(sigma_sfr)


for folder in folder_list:
    file.write("\n\nStarting Energy analysis for " +  folder + " \n\n")
    # data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder, 'Eloading')
    # if(kk==2):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder, 'Eloading')
    
    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', folder, 'Eloading')
    if(kk==0):    data_path = os.path.join('/scratch/jh2/av5889/sims/', folder, 'Eloading')

    os.chdir(data_path)
    list_file1 = glob.glob('*.h5')
    list_file = list_file1
    
    infile   = os.path.join(data_path.replace('Eloading', ''), 'metal_uniform.in')

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
        print(inputfile)
        hf = h5py.File(inputfile ,'r')
        eout_rate[i,:] = np.array(hf.get("EnergyOutflowRate"))
        ethout_rate[i,:] = np.array(hf.get("ThEnergyOutflowRate"))
        timestep[i] = np.array(hf.get("Timestep")) 

        i+=1  
        
    
    etaE = (eout_rate+eout_rate[::-1])/norm[kk]
    etaEth = (ethout_rate+ethout_rate[::-1])/norm[kk]

    tmask = (timestep>tfinal[kk])
    etaE = etaE[tmask]
    etaEth = etaEth[tmask]

    eta_tavg = np.average(etaE,axis=0)
    etaTh_tavg = np.average(etaEth,axis=0)
    
    
    eta_16 = np.percentile(etaE, 16, axis=0)
    eta_84 = np.percentile(etaE, 84, axis=0)
    
    file.write('\n*******************************\n')
    file.write("The total energy loading is -->%.5e\n"%eta_tavg[index])
    file.write("The therma energy loading for is -->%.5e\n"%etaTh_tavg[index])
    file.write("The 16th percentile is -->%.5e\n"%eta_16[index])
    file.write("The 84th percentile is -->%.5e\n"%eta_84[index])
    kk+=1

file.close()