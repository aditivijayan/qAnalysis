
import _init_
from constants import *
from set_path import *
from config import *
from functions import *


sigma_sfr = 6.e-5/yr_to_sec
sfr = sigma_sfr * Msun * 1.e2
Zunit = Msun/1.e3


# folder = 'N4Gpu/2pcNoAMR/PltFiles/'
# folder = 'N4Gpu/16pcNoAMR/NewTry/'
folder = 'N4Gpu/4pcnoAMR-ScalarFix'
data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims', folder)

os.chdir(data_path)
list_file = glob.glob("plt*/")
output_folder = os.path.join(h5_path, folder, 'Eta/')

# parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
# parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
# args = parser.parse_args()


def getLoadingFac(queue):
    while True:
        item = queue.get()
        if item is None:
            break
        
        f = item[0]
        infile = item[1]
        inputfile = os.path.join(data_path, f)
        infile   = os.path.join(data_path, 'metal_uniform.in')

        dom_min, dom_max, ncells = getdomain(infile)

        zrange = np.linspace(dom_min[2], dom_max[2], (int(ncells[2])))
        xrange = np.linspace(dom_min[0], dom_max[0], (int(ncells[0])))
        yrange = np.linspace(dom_min[1], dom_max[1], (int(ncells[1])))

        dx = (dom_max[0]- dom_min[0])/(int(ncells[0]))
        dy = (dom_max[1]- dom_min[1])/(int(ncells[1]))
        dz = (dom_max[2]- dom_min[2])/(int(ncells[2]))

        Lx = (dom_max[0]- dom_min[0])
        Ly = (dom_max[1]- dom_min[1])
        Lz = (dom_max[2]- dom_min[2])

        plane = int(ncells[1]/2)
        lev = 0

        ds   = yt.load(inputfile)
        timestep = ds.current_time.to('Myr')
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions)
        rho_gas = np.array(data['gasDensity'])
        vz = np.array(data['z-GasMomentum'])/rho_gas
        rhoZ0 = np.array(data['scalar_0'])
        sign = zrange/np.abs(zrange)
        
        mass_flux = np.sum(rho_gas*vz, axis=(0,1)) * dx * dy * sign
        total_mass_flux = mass_flux + mass_flux[::-1]
        etaM = total_mass_flux/sfr

        met_sfr = sigma_sfr * Msun
        metal_flux = np.sum(rhoZ0*vz, axis=(0,1)) * dx * dy * Zunit *sign
        total_metal_flux = metal_flux + metal_flux[::-1]
        etaZ = total_metal_flux/met_sfr

        if not os.path.exists(output_folder):
            print(output_folder)
            os.makedirs(output_folder)
                    
        outputfile_name =os.path.join(output_folder, 'lfac-' + f.split('plt')[-1].split('/')[0] + '.npz')
        if os.path.isfile(outputfile_name):
            print("File already exists-->", outputfile_name)
        else:
            np.savez_compressed(outputfile_name, etaM=etaM, etaZ=etaZ, zrange=zrange, timestep=timestep)
            print("File saved-->", outputfile_name)

queue      = Queue()
start_time = ostime.time()
listfile = list_file
num = len(listfile)
infile   = os.path.join(data_path, 'metal_uniform.in')
infile_list = [infile]*num
num_workers = os.cpu_count()


the_pool = [Process(target=getLoadingFac, args=(queue,)) for i in range(num_workers)]
for p in the_pool:
    p.start()

for i in range(num):
    queue.put((listfile[i],infile_list[i]))

for i in range(num_workers):
    queue.put(None)

for p in the_pool:
    p.join()

print("Total time= %s seconds ---" % (ostime.time() - start_time))