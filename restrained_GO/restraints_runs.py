import numpy as np 
from scm.plams import *

mol = Molecule('clu_lin_opt_geo.xyz')
force_field = '/data/gent/vo/000/gvo00003/vsc41842/plams_test/nilton_published_reduced.ff'
print ("molecule:", mol)


class BondScan:
    def __init__(self, name, min, max):
           self.name = name
           self.min = min
           self.max = max

class TorsScan:
    def __init__(self, name, min, max):
           self.name = name
           self.min = min
           self.max = max           

scans = []
scans.append(BondScan("18 61", 1.5, 5.0))
scans.append(BondScan("20 54", 0.5, 3.0))
scans.append(BondScan("31 54", 1.6, 3.9))
scans.append(BondScan("18 58", 1.6, 4.1))
scans.append(BondScan("20 18", 1.8, 4.7))
scans.append(BondScan("23 31", 3.9, 5.9))
scans.append(BondScan("18 31", 2.2, 5.2))
scans.append(BondScan("18 61n31 4", 1.5, 2.7))

scans.append(TorsScan("18 61 4 25", 70, 180))
scans.append(TorsScan("61 18 31 30",100, 180))
scans.append(TorsScan("72 6 58 18", 70, 180))

step = 0.1
tors_step = 10
torsion_flag = False
two_bonds_flag = False
scan_jobs = []

for scan in scans:

    # CASES: 2 bonds and dihedral;
    if 'n' in scan.name:
        two_bonds_flag = True
    elif isinstance(scan, TorsScan):  #Bond or Torsion
        torsion_flag = True
        step = tors_step

    scan_name = scan.name.replace(" ", "-")
    distances = np.arange(scan.min, scan.max + step, step)
    jobs = []

    for dist in distances:

        s = Settings()
        s.input.ams.task = 'geometryoptimization'
        s.input.ams.geometryoptimization.Convergence.Step = 0.01

        # RESTRAINTS cases: 1 & 2 bonds, and dihedral;
        if torsion_flag:
            s.input.ams.restraints.dihedral = f"{scan.name} {dist:.2f} 0.4 Harmonic"
        elif two_bonds_flag:
            bond_1 = scan.name.split('n')[0]
            bond_2 = scan.name.split('n')[1]
            s.input.ams.restraints.distance = [f"{bond_1} {dist:.2f} 0.4 Harmonic", 
                                               f"{bond_2} {dist:.2f} 0.4 Harmonic"] 
        else:
            s.input.ams.restraints.distance = f"{scan.name} {dist:.2f} 0.4 Harmonic"
        
        # REAXFF SETTINGS:
        s.input.reaxff.ForceField = force_field
        s.input.reaxff.Torsions = '2013'
       
        # DFT SETTINGS: 
        #s.input.adf.basis.type = 'TZ2P'
        #s.input.adf.xc.gga = 'PBE'
        #s.input.adf.NumericalQuality = 'Good'
        #s.input.adf.Relativity.Level = 'Scalar'

        jobs.append( AMSJob(molecule=mol, settings=s, name=f"GO_{dist:.2f}") )

    multi = MultiJob(name = scan_name, children = jobs)
    scan_jobs.append(multi)

gr = GridRunner(parallel=True, grid='slurm')

# Running JOBS;
init()
#[multi.run(jobrunner = gr, cores=8, walltime='08:00:00') for multi in scan_jobs]
[multi.run() for multi in scan_jobs]
finish()
