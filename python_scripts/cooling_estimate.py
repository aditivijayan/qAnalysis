
import _init_
from constants import *
from set_path import *
from config import *
from functions import *


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


cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
X = cloudy_H_mass_fraction
Zbg = 1.0
Z  =  0.02*Zbg
Y = 1. - X - Z
mean_metals_A = 16.
sigma_T = 6.6524e-25
T_cmb = 2.725
E_cmb = 7.5e-15 * (T_cmb * T_cmb * T_cmb * T_cmb)
electron_mass_cgs = m_e = 9.1e-28


infile   = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8/chunk-disc-6000000.h5'
hf = h5py.File(infile ,'r')   


rho_gas = hf['Rho'][:]
rhoZ = hf['RhoZ'][:]
temp = hf['Temperature'][:]
xrange = hf['Xrange'][:]
zrange = hf['Zrange'][:]

Zmet = 1. + rhoZ*Msun/1.e3/rho_gas/8.6e-3

rhoH = rho_gas * cloudy_H_mass_fraction
nH =  rhoH / (mp + m_e)
# lambdaZsol = np.zeros(temp.shape)
lambdaZsol  = LambdaZsol(temp.ravel(), nH.ravel()).reshape(temp.shape)

edot_tot = (rhoH*rhoH*lambdaZsol/mp/mp)
tcoolZsol = eint/np.abs(edot_tot)

edot_tot_true = Zmet*(rhoH*rhoH*lambdaZsol/mp/mp)
tcool_true = eint/np.abs(edot_tot_true)

cs = np.sqrt(1.4*boltzmann_constant_cgs * temp/mu/rho_gas)
tdyn = 4.*kpc/cs

bins = 100

tedges = np.linspace(2, 9, bins)  # Log-spaced bins for temperature (e.g., 10^4 to 10^8 K)
ratio_edges = np.linspace(-5, 2, bins )

tcool_ratio = tcoolZsol/tcool_true


outputfile_name = os.path.join(data_path, 'tcool-outflow-6000000.h5')
hfo = h5py.File(outputfile_name, 'w')
hfo.create_dataset('tcoolZsol'       , data=tcoolZsol)
hfo.create_dataset('tcool_true'       , data=tcool_true)
hfo.create_dataset('tcool_ratio'       , data=tcool_ratio)
hfo.create_dataset('tdyn'       , data=tdyn)
hfo.create_dataset('Timestep'       , data=timestep)

hfo.close()
print("--------Written file------->>", outputfile_name)


'''

fig, ax = plt.subplots(1, 2, gridspec_kw = {'wspace':0.00, 'hspace':0.00},figsize=(16,8))

cbarx = 0.151
cbheight = 0.04
cbary = 0.89
cblen = 0.7

dx1 = 0.4
cbtitlex = 0.1

tcool_tdyn = tcool_true/tdyn

tmask = (temp>1.e4)
outflow = (np.abs(zrange)>kpc)*tmask
disc = (np.abs(zrange)<kpc)*tmask

plot = ax[0].scatter(temp[disc].flatten(), tcool_ratio[disc].flatten(), \
                        c=np.log10(tcool_tdyn[disc]).flatten(), 
                        cmap='inferno', vmin=-8, vmax=8.)
cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-8, -6, -4, -2.,  0.0,  2., 4.0, 6.0, 8.0 ))
cax.xaxis.set_ticks_position('top')
cax.set_title(r" $\log t_{\rm cool}/ t_{\rm dyn}$")


plot = ax[1].scatter(temp[outflow].flatten(), tcool_ratio[outflow].flatten(), \
                        c=np.log10(tcool_tdyn[outflow]).flatten(), 
                        cmap='inferno', vmin=-8, vmax=8.)


ax[0].tick_params(axis='y', which='both', right=True, left=True)
ax[1].tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=True)

ax[0].tick_params(axis='x', which='both', top=True)
ax[1].tick_params(axis='x', which='both', top=True)
ax[0].set_ylabel(r"$t_{{\rm cool},{Z_{\odot}}}/t_{{\rm cool},Z}$", fontsize=30)
plt.setp(ax, 'xscale',('log'))
plt.setp(ax, 'yscale',('log'))
plt.setp(ax, 'xlim', (1.e4, 8.e8))
plt.setp(ax[1], 'xlim', (1.e4, 4.e8))
plt.setp(ax, 'ylim', (0.9, 60.))
# ax[0].set_title('Outflow')
# ax[1].set_title('Disc')
plt.setp(ax, 'xlabel', 'T [K]')

ax[0].text(0.08, 0.9, '$|z|<1$ kpc', transform=ax[0].transAxes, fontsize=30)
ax[1].text(0.08, 0.9, '$|z|>1$ kpc', transform=ax[1].transAxes, fontsize=30)
plt.savefig(os.path.join(fig_path, "Paper", "tcool_met1.jpeg"), bbox_inches='tight', dpi=160)
'''