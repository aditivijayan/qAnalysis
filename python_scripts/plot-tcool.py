
import _init_
from constants import *
from set_path import *
from config import *
from functions import *

# Zsim = 0.2
Zsim = 1.0
temp = []

if (Zsim == 1.0):
    infile   = '/g/data/jh2/av5889/quokka_myrepo/quokka/sims/SummitData/GasGravity/Production2pc/R8/tcool-temp-disc.h5'
    hf = h5py.File(infile ,'r')
    temp = hf['Temperature'][:]
    tcool_ratio = hf['TcoolRatio'][:]
    tcool_tdyn = hf['TcoolTdyn'][:]

fig, ax = plt.subplots(1, 2, gridspec_kw = {'wspace':0.00, 'hspace':0.00},figsize=(16,8))

cbarx = 0.155
cbheight = 0.04
cbary = 0.89
cblen = 0.7

dx1 = 0.4
cbtitlex = 0.1


plot = ax[0].scatter(temp.flatten(), tcool_ratio.flatten(), \
                        c=np.log10(tcool_tdyn).flatten(), 
                        cmap='RdYlBu', vmin=-8, vmax=8.)
cax = fig.add_axes([cbarx, cbary, cblen, cbheight])
fig.colorbar(plot, cax=cax, orientation='horizontal', ticks=(-8, -6, -4, -2.,  0.0,  2., 4.0, 6.0, 8.0 ))
cax.xaxis.set_ticks_position('top')
cax.set_title(r" $\log t_{\rm cool,Z}/ t_{\rm dyn}$")



plot = ax[1].scatter(temp.flatten(), tcool_ratio.flatten(), \
                        c=np.log10(tcool_tdyn).flatten(), 
                        cmap='RdYlBu', vmin=-8, vmax=8.)


ax[0].tick_params(axis='y', which='both', right=True, left=True)
ax[1].tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=True)

ax[0].tick_params(axis='x', which='both', top=True)
ax[1].tick_params(axis='x', which='both', top=True)
ax[0].set_ylabel(r"$t_{{\rm cool},{Z_{\rm bg}}}/t_{{\rm cool},Z}$", fontsize=30)
plt.setp(ax, 'xscale',('log'))
plt.setp(ax, 'yscale',('log'))
plt.setp(ax, 'xlim', (1.e4, 8.e8))
plt.setp(ax[1], 'xlim', (1.e4, 4.e8))
plt.setp(ax, 'ylim', (1.e-8, 80.))
# ax[0].set_title('Outflow')
# ax[1].set_title('Disc')
plt.setp(ax, 'xlabel', 'T [K]')

ax[0].text(0.08, 0.9, '$|z|<1$ kpc', transform=ax[0].transAxes, fontsize=30)
ax[1].text(0.08, 0.9, '$|z|>1$ kpc', transform=ax[1].transAxes, fontsize=30)
plt.savefig(os.path.join(fig_path, "Paper", "tcool_disc-full.jpeg"), bbox_inches='tight', dpi=160)

