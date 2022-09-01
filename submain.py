from ase import Atoms
from ase.visualize import view
from ase.build import bulk
from ase.io import read, write
from ase.io import lammpsdata
from ase.lattice.cubic import BodyCenteredCubic
import numpy as np
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
lattice_const = 3.52
atoms = BodyCenteredCubic(directions=[[a,0,0], [0,a,0], [0,0,a]],
                          size=(6,6,6), symbol='Li',
                          latticeconstant=lattice_const)
# print(atoms.positions)
# view(atoms)
dist = 21.12
# удаления атома на диагонали
del atoms[[atom.index for atom in atoms if atom.index == 173]]


# создание атома на месте удаленного с небольшим сдвигом
c = 21.12
r = 0
# print(len(atoms))
rel = (r, r, r)
absol = np.dot(rel, atoms.cell) + (8.9, 8.9, 8.9)
atoms.append('Li')
atoms.positions[-1] = absol

# интерстишил
rel2 = (r, r, r)
absol2 = np.dot(rel2, atoms.cell) + (8.3, 8.3, 8.3)
atoms.append('Li')
atoms.positions[-1] = absol2
ff = atoms.cell
epsilon = lattice_const / 100
# вычисление диагонали  --- с точностью эпсилон
pos = atoms.positions
m = []
for i in range(len(pos)):
    if (pos[i][0] == pos[i][1] == pos[i][2]):
        m.append(i)

#   координаты атомов на диагонале, где m их индексы (включая интерстишил)
n_arr = []
for i in range(len(m)):
    n_arr.append(pos[m[i]])
# print(n_arr)

lammpsdata.write_lammps_data('bulk.lammps-data', atoms)

atoms_read = io.read('E_int_1_diag2.lammps-dump-text')
pos_a_read = atoms_read.positions

print(pos_a_read)
# print(len(pos_a_read))
print('m= ', m)
arr = []
for i in range(len(m)):
    arr.append(pos_a_read[m[i]])
print('\n')
print(arr)

arr_mod = np.array(list(np.float_(arr)))
print(arr_mod)
pos_sort = sort(arr_mod, tags=arr_mod[:, 2])
for i in range(len(m)):
    print(pos_sort[i])
p = []
c = dist
for i in range(1, len(m)):
    p.append(distance(pos_sort[i], pos_sort[i - 1], c))

print(p)

#   ГРАФИК

print(len(p))
x = np.arange(0, 12)
print(x)
plt.plot(x, p, 'o-')
plt.show()


# есть ли id атомов в асе