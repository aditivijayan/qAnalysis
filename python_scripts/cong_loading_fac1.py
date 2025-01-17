import _init_
from constants import *
from set_path import *
from config import *
from functions import *

sigma_sfr = 6.e-5 / yr_to_sec
sfr = sigma_sfr * Msun * 1.e2
Zunit = Msun / 1.e3

folder = 'N4Gpu/2pcNoAMR/PltFiles/'
data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims', folder)
os.chdir(data_path)
list_file = glob.glob("plt*/")
output_folder = os.path.join(h5_path, folder, 'Eta/')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def process_file(file_name, infile_path):
    inputfile = os.path.join(data_path, file_name)

    dom_min, dom_max, ncells = getdomain(infile_path)

    zrange = np.linspace(dom_min[2], dom_max[2], int(ncells[2]))
    xrange = np.linspace(dom_min[0], dom_max[0], int(ncells[0]))
    yrange = np.linspace(dom_min[1], dom_max[1], int(ncells[1]))

    dx = (dom_max[0] - dom_min[0]) / int(ncells[0])
    dy = (dom_max[1] - dom_min[1]) / int(ncells[1])
    dz = (dom_max[2] - dom_min[2]) / int(ncells[2])

    ds = yt.load(inputfile)
    timestep = ds.current_time.to('Myr')
    data = ds.covering_grid(level=0, left_edge=dom_min, dims=ds.domain_dimensions)

    rho_gas = np.array(data['gasDensity'])
    vz = np.array(data['z-GasMomentum']) / rho_gas
    rhoZ0 = np.array(data['scalar_0'])
    sign = zrange / np.abs(zrange)

    mass_flux = np.sum(rho_gas * vz, axis=(0, 1)) * dx * dy * sign
    total_mass_flux = mass_flux + mass_flux[::-1]
    etaM = total_mass_flux / sfr

    met_sfr = sigma_sfr * Msun
    metal_flux = np.sum(rhoZ0 * vz, axis=(0, 1)) * dx * dy * Zunit * sign
    total_metal_flux = metal_flux + metal_flux[::-1]
    etaZ = total_metal_flux / met_sfr

    outputfile_name = os.path.join(output_folder, f'lfac-{file_name.split("plt")[-1].split("/")[0]}.npz')

    if os.path.isfile(outputfile_name):
        print("File already exists -->", outputfile_name)
    else:
        np.savez_compressed(outputfile_name, etaM=etaM, etaZ=etaZ, zrange=zrange, timestep=timestep)
        print("File saved -->", outputfile_name)
    
    # Clear memory explicitly
    del ds, data, rho_gas, vz, rhoZ0, mass_flux, total_mass_flux, etaM, metal_flux, total_metal_flux, etaZ
    import gc
    gc.collect()

# Sequential processing of files
infile_path = os.path.join(data_path, 'metal_uniform.in')
start_time = ostime.time()

for file_name in list_file:
    print(f"Processing file: {file_name}")
    process_file(file_name, infile_path)

print("Total time = %s seconds ---" % (ostime.time() - start_time))
