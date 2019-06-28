from ase import *
from ase.io import read, write
from ase.parallel import paropen, parprint
from gpaw import GPAW, setup_paths, FermiDirac, PW
from gpaw import Mixer, MixerSum, MixerDif
#from gpaw.eigensolvers import Davidson
from gpaw.poisson import PoissonSolver
from ase.units import Bohr
from gpaw import GPAW, restart
import os
from gpaw.xas import XAS

setup_paths.insert(0, './')

sys = read('EMImCCN3.xyz')

sys.center(vacuum=8)

calc1 = GPAW(h=0.16,
             charge=0,
             mode="fd",
	     xc='PBE',
             txt='EMImCCN3_gs.txt',
             eigensolver='rmm-diis',
             )
sys.set_calculator(calc1)
e1 = sys.get_potential_energy() + calc1.get_reference_energy()
rho = sys.calc.get_all_electron_density(gridrefinement=2)
write('EMImCCN3_density2.cube', sys, data=rho * Bohr**3)
calc.write('EMImCCN3.gpw')
