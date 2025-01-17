
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import sys
import matplotlib.colors as mcolors
import matplotlib.cm

def process_file(foldername, filename):
    
    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    

    ds   = yt.load(inputfile)
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
    timestep = ds.current_time.to('Myr')
    rho_gas = np.array(data['gasDensity'])
    egas    = np.array(data['gasEnergy'])
    eint    = np.array(data['gasInternalEnergy'])
    pvz = np.array(data['z-GasMomentum'])
   
    vz = np.array(data['z-GasMomentum'])/rho_gas

    eout_rate = egas *  vz * dx * dy
    eth_rate  = eint *  vz * dx * dy

    energy_outflow_rate = np.sum(eout_rate, axis=(0,1))

    eint_outflow_rate = np.sum(eth_rate, axis=(0,1))


    output_folder = os.path.join(foldername , 'Eloading/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    outputfile_name =os.path.join(output_folder, 'eloading_' + filename.split('plt')[1] + '.h5') 
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:
        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Timestep'       , data=timestep)
        hfo.create_dataset('EnergyOutflowRate'       , data=energy_outflow_rate)
        hfo.create_dataset('ThEnergyOutflowRate'       , data=eint_outflow_rate)
        
        hfo.close()
        print("--------Written file------->>", inputfile)

def phase_mass(foldername, filename):

    file = "/ccs/home/aditiv/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5"
    grackle = h5py.File(file)
    array = grackle['CoolingRates/Primordial/MMW'][()]
    table = array[:,0,:]
    table_nH   = np.logspace(-10., 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])

    i=0
    bins = 200
    egas_arr = np.logspace(-25., -5., bins)
    nH_arr   = np.logspace(-7.0, 4.0, int(bins))
    T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))

    for egas in egas_arr:
        j=0
        for nH in nH_arr:
            C = (gamma - 1.) * egas / (boltzmann_constant_cgs*nH)
            minT = C*np.amin(table)
            maxT = C*np.amax(table)
            def func(T):
                mu = interpolate.interp2d(table_temp, table_nH, table,\
                                kind='linear', copy=True, bounds_error=False, fill_value=None)
                return C*mu(T,nH)[0] - T

            T[i,j] = scipy.optimize.toms748(func, minT, maxT)
            j+=1
        i+=1


    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
    
    
    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    dV = dx * dy * dz

    ds   = yt.load(inputfile)
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
    timestep = ds.current_time.to('Myr')
    rho_gas = np.array(data['gasDensity'])
    eint    = np.array(data['gasInternalEnergy'])
    egas0=eint
    density = rho_gas
    cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
    rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs


    logrho_arr = np.log10(nH_arr[:-1])
    logrho     = np.log10(rho0)
    delta_rho  = logrho_arr[1] - logrho_arr[0]
    idxrho     = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas     = np.log10(egas0)
    delta_egas  = logEgas_arr[1] - logEgas_arr[0]

    idxegas     = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
            wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
        (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
            wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  

    warm = temp<2.e4
    hot = temp>1.e6 

    mass = rho_gas * dV
    mass_warm = rho_gas[warm] * dV
    mass_hot  = rho_gas[hot] * dV
    vol = np.full(temp.shape, dx*dy)
    vol_warm = vol[warm]
    vol_hot = vol[hot]

    warm_mass_z  = np.sum(np.ma.array(mass, mask=~warm), axis=(0,1)).data 
    hot_mass_z  = np.sum(np.ma.array(mass, mask=~hot), axis=(0,1)).data 
    warm_vol_z = np.sum(np.ma.array(vol, mask=~warm), axis=(0,1)).data 
    hot_vol_z = np.sum(np.ma.array(vol, mask=~hot), axis=(0,1)).data 

    output_folder = os.path.join(foldername , 'Phases/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    outputfile_name =os.path.join(output_folder, 'phase_' + filename.split('plt')[1] + '.h5') 
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:
        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Timestep'       , data=timestep)
        hfo.create_dataset('WarmMass'       , data=warm_mass_z)
        hfo.create_dataset('HotMass'       , data=hot_mass_z)
        hfo.create_dataset('WarmVol'       , data=warm_vol_z)
        hfo.create_dataset('HotVol'       , data=hot_vol_z)
        
        hfo.close()
        print("--------Written file------->>", inputfile)


def temp_dens_histo(foldername, filename):
    file = "/ccs/home/aditiv/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5"
    grackle = h5py.File(file)
    array = grackledata['CoolingRates/Primordial/MMW'][()]
    table = array[:,0,:]
    table_nH   = np.logspace(-10., 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])

    i=0
    bins = 200
    egas_arr = np.logspace(-25., -5., bins)
    nH_arr   = np.logspace(-7.0, 4.0, int(bins))
    T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))

    for egas in egas_arr:
        j=0
        for nH in nH_arr:
            C = (gamma - 1.) * egas / (boltzmann_constant_cgs*nH)
            minT = C*np.amin(table)
            maxT = C*np.amax(table)
            def func(T):
                mu = interpolate.interp2d(table_temp, table_nH, table,\
                                kind='linear', copy=True, bounds_error=False, fill_value=None)
                return C*mu(T,nH)[0] - T

            T[i,j] = scipy.optimize.toms748(func, minT, maxT)
            j+=1
        i+=1


    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))


    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    dV = dx * dy * dz

    ds   = yt.load(inputfile)
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions, num_ghost_zones=0)
    timestep = ds.current_time.to('Myr')
    rho_gas = np.array(data['gasDensity'])
    eint    = np.array(data['gasInternalEnergy'])
    vz = np.array(data['z-GasMomentum'])/rho_gas
    rhoZ = np.array(data['scalar_0'])
    abund = (rhoZ/rho_gas)*Msun/1.e3
    egas0=eint
    density = rho_gas
    cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
    rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs


    logrho_arr = np.log10(nH_arr[:-1])
    logrho     = np.log10(rho0)
    delta_rho  = logrho_arr[1] - logrho_arr[0]
    idxrho     = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas     = np.log10(egas0)
    delta_egas  = logEgas_arr[1] - logEgas_arr[0]

    idxegas     = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
            wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
        (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
            wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]
    bins = 201

    (dlow, dhigh)   = (1.e-6, 1.e4)
    (tlow, thigh)   = (10., 2.e10)

    zout = (np.abs(zrange)>kpc)
    temp_outflow = temp[:,:, zout]
    nH_outflow = nH[:,:, zout]
    mass_outflow = mass[:,:, zout]

    dedges = np.logspace(np.log10(dlow), np.log10(dhigh), bins)
    tedges = np.logspace(np.log10(tlow), np.log10(thigh), bins)

    htot, binx, biny = np.histogram2d(nH.reshape(-1,1)[:,0], temp.reshape(-1,1)[:,0], \
                                    bins=[dedges, tedges], weights=mass.reshape(-1,1)[:,0])

    hout, binx, biny = np.histogram2d(nH_outflow.reshape(-1,1)[:,0], temp_outflow.reshape(-1,1)[:,0], \
                                    bins=[dedges, tedges], weights=mass_outflow.reshape(-1,1)[:,0])

    output_folder = os.path.join(foldername , 'Phases/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    outputfile_name =os.path.join(output_folder, 'temp-dens_' + filename.split('plt')[1] + '.h5')
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:

        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('DensBins'       , data=binx)
        hfo.create_dataset('TempBins'       , data=biny)
        hfo.create_dataset('TotalMass'       , data=htot)
        hfo.create_dataset('OutflowMass'       , data=hout)
        hfo.create_dataset('Timestep', data=timestep)
        hfo.close()
        print("--------Written file------->>",f)


def get_slices(foldername, filename):
    file = "/ccs/home/aditiv/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5"
    grackle = h5py.File(file)
    array = grackle['CoolingRates/Primordial/MMW'][()]
    table = array[:,0,:]
    table_nH   = np.logspace(-10., 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])
    print("Inside get slice!")
    i=0
    bins = 200
    egas_arr = np.logspace(-25., -5., bins)
    nH_arr   = np.logspace(-7.0, 4.0, int(bins))
    T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))

    for egas in egas_arr:
        j=0
        for nH in nH_arr:
            C = (gamma - 1.) * egas / (boltzmann_constant_cgs*nH)
            minT = C*np.amin(table)
            maxT = C*np.amax(table)
            print(table_nH.shape, table_temp.shape, table.shape)
            def func(T):
                #mu = interpolate.RectBivariateSpline(table_nH, table_temp, table)
                mu = interpolate.RegularGridInterpolator((table_nH, table_temp), table,\
                                                  method='nearest', bounds_error=False, fill_value=None)
                return C*mu((nH,T)) - T

            T[i,j] = scipy.optimize.toms748(func, minT, maxT)
            j+=1
        i+=1

    print("I am here!\n")
    inputfile = os.path.join(foldername, filename)
    infile   = os.path.join(foldername, 'metal_uniform.in')

    dom_min, dom_max, ncells = getdomain(infile)
    fac = 1
    xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
    zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))


    dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
    dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
    dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
    dV = dx * dy * dz
    plane = int (ncells[1]/2)
    

    print("loading file...")
    ds   = yt.load(inputfile)
    timestep = ds.current_time.to('Myr')
    ncells = ds.domain_dimensions
    dom_min = ds.domain_left_edge
    dom_max = ds.domain_right_edge

    print("current time =", timestep)

    print("reading fields...")
    data = ds.r[::ncells[0]*1j, plane, ::ncells[2]*1j]
    rho_gas = data['gasDensity']
    eint    = data['gasInternalEnergy']
    vz = data['z-GasMomentum']/rho_gas
    rhoZ = data['scalar_0']
    print('done.')


    abund = (rhoZ/rho_gas)*Msun/1.e3 

    egas0=eint
    density = rho_gas
    cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
    rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs

    logrho_arr = np.log10(nH_arr[:-1])
    logrho     = np.log10(rho0)
    delta_rho  = logrho_arr[1] - logrho_arr[0]
    idxrho     = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas     = np.log10(egas0)
    delta_egas  = logEgas_arr[1] - logEgas_arr[0]

    idxegas     = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')

    print(np.amin(idxegas), np.where(idxegas==np.amin(idxegas)))
    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
               wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
          (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
               wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  
    
    output_folder = os.path.join(foldername , 'Slices/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    outputfile_name =os.path.join(output_folder, 'slice_' + filename.split('plt')[-1].split('/')[0]  + '.h5')
    if(os.path.isfile(outputfile_name)):
        print("File already exists-->", outputfile_name)
    else:

        hfo = h5py.File(outputfile_name, 'w')
        hfo.create_dataset('Xrange'       , data=xrange)
        hfo.create_dataset('Zrange'       , data=zrange)
        hfo.create_dataset('Rho'       , data=rho_gas)
        hfo.create_dataset('Temp'       , data=temp)
        hfo.create_dataset('Vz'       , data=vz)
        hfo.create_dataset('RhoZ'       , data=rhoZ)
        hfo.create_dataset('Timestep', data=timestep)
        hfo.close()
        print("--------Written file------->>",outputfile_name)


def plot_slices(foldername, filename):

    input_folder = os.path.join(foldername , 'Slices/')
    inputfile_name =os.path.join(input_folder, 'slice_' + filename.split('plt')[-1].split('/')[0]  + '.h5')

    if(os.path.isfile(inputfile_name)):

        hf = h5py.File(data_path ,'r')

        timestep = np.array(hf.get("Timestep")) 
        zrange = np.array(hf.get("Zrange"))
        xrange = np.array(hf.get("Xrange"))
        rho_gas = np.array(hf.get("Rho")) 
        rhoZ = np.array(hf.get("RhoZ")) 
        vz = np.array(hf.get("Vz")) 
        abund = (rhoZ/rho_gas)*Msun/1.e3
        temp = np.array(hf.get("Temp"))

        fig, ax = plt.subplots(1, 4, gridspec_kw = {'wspace':0.01, 'hspace':0.0},figsize=(16, 32))
        i=0

        cbarx = 0.13
        cbheight = 0.01
        cbary = 0.89
        cblen = 0.15
        dx1 = 0.2
        cbtitlex = 0.1
        cbtitley = 16.5


        print(xrange.shape, zrange.shape, np.log10(rho_gas/mp).shape)
        plot = ax[0].pcolormesh(xrange/kpc,zrange/kpc, np.transpose(np.log10(rho_gas/mp)),\
                        vmin=-5., vmax=-1.5,
                        cmap='plasma')
        cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-6, -4, -2.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $\mathrm{log} [\rho/m_H]$" + "\n" + "[cm$^{-3}$]")
        # ax[0].text(0.4, 0.9, '%.2f'%(timestep) + ' Myr', transform = ax[0].transAxes, color='black')
        
        plot = ax[1].pcolormesh(xrange/kpc,zrange/kpc, np.transpose((temp)),\
                        norm=mcolors.LogNorm(vmin=8.e2, vmax=5.e8),
                        cmap="cubehelix")
        cax = fig.add_axes([cbarx + dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e4, 1.e6, 1.e8))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" T" + "\n" + "[K]")

        plot = ax[2].pcolormesh(xrange/kpc,zrange/kpc, np.transpose(abund)/8.6e-3,\
                        norm=mcolors.LogNorm(vmin=2.e-1, vmax=10.2),
                        cmap="BuPu_r")
        cax = fig.add_axes([cbarx + 2.*dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(1.e-1, 1., 1.e1))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r" $Z$ " + "\n" + r"$[Z_{\odot}]$")

        plot = ax[3].pcolormesh(xrange/kpc,zrange/kpc, np.transpose(vz)/kmps/1.e3,\
                        norm=mcolors.SymLogNorm(vmin=-1.2, vmax=1.2, linthresh=0.05),
                        cmap='coolwarm')
        cax = fig.add_axes([cbarx + 3.*dx1, cbary, cblen, cbheight])
        fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-1., 0.0, 1.))
        cax.xaxis.set_ticks_position('top')
        cax.set_title(r"$v_z$" + "\n" + r"$[1000$ km s$^{-1}]$")

        
        ax[0].tick_params(axis='y', labelleft=True, labelright=False, right=True, left=True)
        ax[1].tick_params(axis='y', labelleft=False, labelright=False, right=True, left=True)
        ax[2].tick_params(axis='y', labelleft=False, labelright=False, right=True, left=True)
        ax[ -1].tick_params(axis='y', labelleft=False, labelright=True, right=True, left=True)

        ax[0].tick_params(axis='x', labelbottom=True, top=True, bottom=True)
        ax[1].tick_params(axis='x', labelbottom=True, top=True, bottom=True)
        ax[2].tick_params(axis='x', labelbottom=True, top=True, bottom=True)
        ax[3].tick_params(axis='x', labelbottom=True, top=True, bottom=True) 

        #ax[0].text(0.05, 0.95, r'$\Sigma50$-$Z1$-H$150$', transform=ax[0].transAxes, color='white', fontsize=30)
        ax[0].text(0.05, 0.95, '%.1f'%(timestep) + ' Myr', transform=ax[0].transAxes, color='white', fontsize=30)

        plt.setp(ax, 'ylim', (-3.9,4.))
        plt.setp(ax, 'xlim', (0.0,1.0))

        plt.setp(ax, 'yticks', (-2., 0.0, 2.0, 4.0))
        ax[-1].set_xticks((0.0, 0.5, 1.0))

        ax[0].set_ylabel(r'$z$ [kpc]', fontsize=30)
        plt.setp(ax, 'xlabel', r'$x$ [kpc]')
        plt.setp(ax, 'xticks', (0.0,0.5))
        plt.setp(ax[-1], 'xticks', (0.0,0.5, 1.0))
        
        output_folder = os.path.join(foldername , 'Slices/')

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        outputfile_name =os.path.join(output_folder, 'slices_' + filename.split('plt')[-1].split('/')[0] + '.jpeg')
        if(os.path.isfile(outputfile_name)):
            print("File already exists-->", outputfile_name)
        else:
            plt.savefig(outputfile_name, bbox_inches='tight', dpi=160)
            print("--------Written file------->>",outputfile_name)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No file path provided.")
    else:
        foldername = sys.argv[1]
        filename = sys.argv[2]
        print(sys.argv)
        get_slices(foldername, filename)
#        plot_slices(foldername, filename)
