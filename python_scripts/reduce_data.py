
import _init_
from constants import *
from set_path import *
from config import *
from functions import *


file = '/g/data/jh2/av5889/freshquokka/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
grackle = h5py.File(file)
array = grackle['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-10, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])



i=0
bins = 100
egas_arr = np.logspace(-21., -5., bins)
nH_arr   = np.logspace(-6.0, 4.0, int(bins))
temp_arr   = np.logspace(1.0, 9.0, int(bins))
T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))
Mu_mu = np.zeros((temp_arr.shape[0],nH_arr.shape[0]))

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


mu_interp = interpolate.interp2d(table_temp, table_nH, table,\
                              kind='linear', copy=True, bounds_error=False, fill_value=None)
i=0
for ttemp in temp_arr:
    j=0
    for nH in nH_arr:
        Mu_mu[i,j] = mu_interp(ttemp, nH)
        j+=1
    i+=1    


data_path = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/N4Gpu/2pcNoAMR/PltFiles/'
infile   = os.path.join(data_path, 'metal_uniform_512.in')
dom_min, dom_max, ncells = getdomain(infile)
fac = 1
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))


f = 'plt3005000/'
f = 'plt2200000/'
inputfile = os.path.join(data_path, f)
ds   = yt.load(inputfile)
ds.current_time.to('Myr')

lev = 0
data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac, num_ghost_zones=0)
timestep = ds.current_time.to('Myr')

rho_gas = np.array(data['gasDensity'])
eint    = np.array(data['gasInternalEnergy'])
rhoZ    = np.array(data['scalar_1'])
MO = 1. * Msun
rhoOxy_inj = MO * rhoZ/1.e3
Zabund = rhoOxy_inj/rho_gas

egas0=eint
density = rho_gas
cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
X = cloudy_H_mass_fraction
Z  =  0.02
Y = 1. - X - Z
mean_metals_A = 16.

rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs
nH_val = rho0/(hydrogen_mass_cgs + m_e)

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


logTemp_arr = np.log10(temp_arr[:-1])
logTgas     = np.log10(temp)
delta_tgas  = logTemp_arr[1] -logTgas[0]
idxTgas     = (np.floor((logTgas-np.amin(logTemp_arr))/delta_tgas)).astype('int')

wgt_Tgas = (logTgas - (np.amin(logTemp_arr) + delta_tgas*idxTgas))/delta_tgas

mu_val = (1.-wgt_rho)*(1.-wgt_Tgas)* Mu_mu[tuple(idxTgas)  , tuple(idxrho)]   +\
           wgt_rho *    wgt_Tgas * Mu_mu[tuple(idxTgas+1), tuple(idxrho+1)] +\
      (1. -wgt_rho)*    wgt_Tgas * Mu_mu[tuple(idxTgas+1), tuple(idxrho)]   +\
           wgt_rho *(1.-wgt_Tgas)* Mu_mu[tuple(idxTgas)  , tuple(idxrho+1)]  

n_e = (rho_gas/(hydrogen_mass_cgs + m_e)) * (1.0 - mu_val * (X + Y / 4. + Z / mean_metals_A)) / (mu_val - (m_e / (hydrogen_mass_cgs + m_e)))

disk = (np.abs(zrange)/kpc<1.)

output_folder = os.path.join(qhome, 'quokka/sims/N4Gpu/2pcNoAMR/Setonix/Data_for_Yifei/')
if not os.path.exists(output_folder):
    print(output_folder)
    os.makedirs(output_folder)

outputfile_name =os.path.join(output_folder, 'reduce_data_20Myr.h5')

hfo = h5py.File(outputfile_name, 'w')
hfo.create_dataset('GasDensity' , data=rho_gas[:,:,disk])
hfo.create_dataset('ne' , data=n_e[:,:,disk])
hfo.create_dataset('Temp'    , data=temp[:,:,disk])
hfo.create_dataset('Metallicity'    , data=Zabund[:,:,disk])

hfo.create_dataset('X'  , data=xrange)
hfo.create_dataset('Y'  , data=yrange)
hfo.create_dataset('Z'  , data=zrange[disk])
hfo.create_dataset('Timestep', data=timestep)
hfo.close()

print("Average abundance in the disk=", np.mean(Zabund[:,:,disk]))