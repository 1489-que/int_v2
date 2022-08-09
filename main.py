from ase import Atoms
from ase.visualize import view
from ase.build import bulk
from ase.io import read, write
from ase.lattice.cubic import BodyCenteredCubic
import numpy as np
import ase
from ase.build import attach
import math
from ase import neighborlist
from ase.io import xyz
from ase.optimize.sciopt import SciPyFminBFGS as spo
from ase.calculators.lammpsrun import LAMMPS

# parameters = {'pair_style': 'eam/alloy',
#               'pair_coeff': ['* * Li_Belash2013.eam.alloy Li']}
#
# files = ['Li_Belash2013.eam.alloy']
# lammps = LAMMPS(parameters=parameters, files=files)

a = 1
atoms = BodyCenteredCubic(directions=[[a,0,0], [0,a,0], [0,0,a]],
                          size=(6,6,6), symbol='Li',
                          latticeconstant=3.52)
# print(atoms.positions)
# view(atoms)

del atoms[[atom.index for atom in atoms if atom.index == 173]]

# print(pos[173], '\n')

# print(atoms.cell, '\n')
c = 21.12
r = 0
# print(len(atoms))
rel = (r, r, r)
absol = np.dot(rel, atoms.cell) + (8.9, 8.9, 8.9)
atoms.append('Li')
atoms.positions[-1] = absol

rel2 = (r, r, r)
absol2 = np.dot(rel2, atoms.cell) + (8.3, 8.3, 8.3)
atoms.append('Li')
atoms.positions[-1] = absol2

# v = 19.36
# rel3 = (r, r, r)

# absol3 = np.dot(rel2, atoms.cell) + (0, 0, v)
# atoms.append('Li')
# atoms.positions[-1] = absol3

pos = atoms.positions
m = []
for i in range(len(pos)):
    j = 0
    if pos[i][j] == pos[i][j + 1] == pos[i][j + 2]:
        m.append(i)
# print(m)

# print(len(atoms))

# view(atoms)
n_arr = []
for i in range(len(m)):
    n_arr.append(pos[m[i]])
# print(n_arr)

# atoms.calc = lammps
# spo(atoms)


# расстояние между атомами
dist = atoms.cell
print(len(n_arr))
def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

print(distance(n_arr[11], n_arr[1], dist))




