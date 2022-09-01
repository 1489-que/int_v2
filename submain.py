from ase import Atoms
from ase.visualize import view
from ase.build import bulk
import os
from ase.io import read, write
from ase.io import lammpsdata
from ase.lattice.cubic import BodyCenteredCubic
import numpy as np
from subprocess import Popen
import ase
from matplotlib import pyplot as plt
from ase.build import attach
from ase import io
import math
from ase.build.tools import sort
from ase import neighborlist
from ase.io import xyz
from ase.optimize.sciopt import SciPyFminBFGS as spo
from ase.optimize import FIRE
from ase.calculators.lammpsrun import LAMMPS
def distance(x0, x1, dimensions):
    delta = np.abs(x1 - x0)
    if delta[0] > ( 0.5 * dimensions ):
        delta[0] = delta[0] - dimensions
    else:
        delta[0] = delta[0]
    if delta[1] > ( 0.5 * dimensions ):
        delta[1] = delta[1] - dimensions
    else:
        delta[1] = delta[1]
    if delta[2] > ( 0.5 * dimensions ):
        delta[2] = delta[2] - dimensions
    else:
        delta[2] = delta[2]
    return np.sqrt((delta ** 2).sum())

# сам кристалл
a = 1
nsize = 6
lattice_const = 3.52
atoms = BodyCenteredCubic(directions=[[a,0,0], [0,a,0], [0,0,a]],
                          size=(nsize, nsize, nsize), symbol='Li',
                          latticeconstant=lattice_const)
c = 21.12
########################

pos = atoms.positions
atoms.append('Li')

# диагональ
d = []
for atm in atoms:
    if np.std(atm.position) <= lattice_const / 100:
        d.append(atm)
pos_int = d[6].position / 3
atoms.positions[-1] = pos_int
d = []
m = []
for atm in reversed(atoms):
    if np.std(atm.position) <= lattice_const / 100:
        d.append(atm)
        m.append(atm.index)
print(d)
print(m)

########################

lammpsdata.write_lammps_data('bulk.lammps-data', atoms)
#os.system(r"C:\Users\default.DESKTOP-J9F7KTV\OneDrive\desktop\GB_structures\CODE\run.bat")
atoms_read = io.read('E_int_1_diag2.lammps-dump-text')
pos_a_read = atoms_read.positions
print(pos_a_read)
#########################

arr = []
for i in range(len(m)):
    arr.append(pos_a_read[m[i]])
arr_mod = np.array(list(np.float_(arr)))
pos_sort = sort(arr_mod, tags=arr_mod[:, 2])
p = []
print(pos_sort)
for i in range(1, len(m)):
    p.append(distance(pos_sort[i], pos_sort[i - 1], c))
p.append(distance(pos_sort[len(m) - 1], pos_sort[0], c))
print(p)

# plot


x = np.arange(0, len(m))
plt.plot(x, p, 'o-')
plt.show()