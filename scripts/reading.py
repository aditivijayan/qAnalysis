class Reading:

    def __init__(filename, lev):
        self.file = filename
        self.lev = lev

    def getdomain(file):
        infile = open(file)
        lines = infile.readlines()
        dom_range = np.zeros((2,3))
        ncell = np.zeros(3)
        dom_min = [0.0,0.0,0.0]
        dom_min[0] = float(lines[3].split()[2])
        dom_min[1] = float(lines[3].split()[3])
        dom_min[2] = float(lines[3].split()[4])
        
        dom_max = [0.0,0.0,0.0]
        dom_max[0] = float(lines[4].split()[2])
        dom_max[1] = float(lines[4].split()[3])
        dom_max[2] = float(lines[4].split()[4])

        ncell[0]=int(lines[15].split()[2])
        ncell[1]=int(lines[15].split()[3])
        ncell[2]=int(lines[15].split()[4])
        
        return dom_min, dom_max, ncell
        
    def read(self):
        fac = 1 if self.lev == 0 else 2.*self.lev
    
        dom_min, dom_max, ncells = getdomain(infile)
        ds   = yt.load(inputfile)
        data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)

        timestep = ds.current_time.to('Myr')
        density = np.array(data['gasDensity'])
        egas0  = np.array(data['gasInternalEnergy'])
        rho0 = density/hydrogen_mass_cgs
        pz = np.array(data['z-GasMomentum'])
        vz = pz/density

        return timestep, density, egas, rho, pz, vz




read_ojb = Reading()
timestep, density, _, pz, vz = read_ojb(filename, lev)
timestep, density, _ = read_ojb(filename, lev)
