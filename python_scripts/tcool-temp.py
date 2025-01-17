
import _init_
from constants import *
from set_path import *
from config import *
from functions import *


cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)


file = '/g/data/jh2/av5889/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=HM2012.h5'
# file = '/g/data/jh2/av5889/quokka_myrepo/quokka/extern/grackle_data_files/input/CloudyData_UVB=FG2011.h5'
#Dim2 is redshift. Dim1 is density. Dim3 is temperature
grackle = h5py.File(file)
array = grackle['CoolingRates/Primordial/MMW'][()]
#density(1.e-6, 1.e4), redshift(0,15), temperature(10., 1.e9)
table = array[:,0,:]
table_nH   = np.logspace(-10, 4, array.shape[0])
table_temp = np.logspace(1,  9, array.shape[2])

prim_heating = grackle['CoolingRates/Primordial/Heating'][()][:,0,:]
prim_cooling = grackle['CoolingRates/Primordial/Cooling'][()][:,0,:]

met_heating  = grackle['CoolingRates/Metals/Heating'][()][:,0,:]
met_cooling  = grackle['CoolingRates/Metals/Cooling'][()][:,0,:]

netLambda_prim =  - prim_cooling 

netLambda_met  = - met_cooling
netLambda_Zsol = prim_cooling + met_cooling


LambdaZsol = interpolate.interp2d(table_temp, table_nH, netLambda_Zsol,\
                              kind='linear', copy=True, bounds_error=False, fill_value=None)

##GET MU---####


i=0
bins = 100
egas_arr = np.logspace(-21., -5., bins)
nH_arr   = np.logspace(-6.0, 4.0, int(bins))
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

    
mu_interp = interpolate.interp2d(table_temp, table_nH, table,\
                              kind='linear', copy=True, bounds_error=False, fill_value=None)


Zsim = 0.2
# Zsim = 1.0
temp = []
'''
if (Zsim == 0.2):

    infile   = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8-0.2Zsol/slice-3610000-xz.h5'
    hf = h5py.File(infile ,'r')
    rho_gas = hf['Rho'][:]
    rhoZ = hf['RhoZ'][:]
    eint = hf['IntEnergy'][:]
    xrange = hf['Xrange'][:]
    zrange = hf['Zrange'][:]
    temp = hf['Temperature'][:]

if (Zsim == 1.0):
#     infile   = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8/slice-6000000.h5'
    infile   = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8/chunk-outflow-6000000.h5'
    hf = h5py.File(infile ,'r')
    rho_gas = hf['Rho'][:]
    rhoZ = hf['RhoZ'][:]
    eint = hf['IntEnergy'][:]
    xrange = hf['Xrange'][:]
    zrange = hf['Zrange'][:]
    temp = hf['Temperature'][:]

'''

data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R4-0.2Zsol/')

infile   = os.path.join(data_path, 'metal_uniform.in')
# infile   = os.path.join(data_path, 'metal_uniform_512.in')
dom_min, dom_max, ncells = getdomain(infile)
fac = 1
zrange = np.linspace(dom_min[2], dom_max[2], (fac*int(ncells[2])))
xrange = np.linspace(dom_min[0], dom_max[0], (fac*int(ncells[0])))
yrange = np.linspace(dom_min[1], dom_max[1], (fac*int(ncells[1])))

dx = (dom_max[0]- dom_min[0])/(fac*int(ncells[0]))
dy = (dom_max[1]- dom_min[1])/(fac*int(ncells[1]))
dz = (dom_max[2]- dom_min[2])/(fac*int(ncells[2]))
dV = dx * dy * dz

f = 'plt6690000/'
inputfile = os.path.join(data_path, f)
ds   = yt.load(inputfile)
lev = 0
data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac, num_ghost_zones=0)
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




rhoH = rho_gas * cloudy_H_mass_fraction
nH =  rhoH / (mp + m_e)
mu = np.zeros(temp.shape)
for i in range(temp.shape[0]):
    for j in range(temp.shape[1]):
        mu[i,j]  = np.vectorize(mu_interp)(temp[i,j].ravel(), nH[i,j].ravel())

rhoH = rho_gas * cloudy_H_mass_fraction
nH =  rhoH / (mp + m_e)
lambdaZsol = np.zeros(temp.shape)
for i in range(temp.shape[0]):
    for j in range(temp.shape[1]):
        lambdaZsol[i,j]  = np.vectorize(LambdaZsol)(temp[i,j].ravel(), nH[i,j].ravel())


edot_tot = Zsim * (rhoH*rhoH*lambdaZsol/mp/mp)
tcool_Zsim = eint/np.abs(edot_tot)

Zmet = Zsim + rhoZ*Msun/rho_gas /1.e3 /8.6e-3
edot_tot_true = Zmet*(rhoH*rhoH*lambdaZsol/mp/mp)
tcool_true = eint/np.abs(edot_tot_true)


cs = np.sqrt(1.4*boltzmann_constant_cgs * temp/mu/rho_gas)
tdyn = 4.*kpc/cs

tcool_ratio = (tcool_Zsim/tcool_true)
tcool_tdyn = tcool_true/tdyn
weights = tcool_Zsim/tcool_true

# data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8/')
outputfile_name = os.path.join(data_path, 'tcool-temp.h5')
hfo = h5py.File(outputfile_name, 'w')
hfo.create_dataset('Temperature'       , data=temp)
hfo.create_dataset('TcoolRatio'       , data=tcool_ratio)
hfo.create_dataset('TcoolTdyn'       , data=tcool_tdyn)
hfo.create_dataset('Zrange'       , data=zrange)
hfo.create_dataset('Xrange'       , data=xrange)
hfo.close()
print("--------Written file------->>", outputfile_name)


