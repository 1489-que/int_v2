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
import math
from ase import neighborlist
from ase.io import xyz
from ase.optimize.sciopt import SciPyFminBFGS as spo
from ase.optimize import FIRE
from ase.calculators.lammpsrun import LAMMPS


def distance(x0, x1, dimensions):
    delta = np.abs(x1 - x0)
    # delta = np.where(delta[0] > 0.5 * dimensions[0], delta[0] - dimensions[0], delta[0])
    if delta[0] > ( 0.5 * dimensions ):
        delta = delta[0] - dimensions
    else:
        delta = delta[0]
    return np.sqrt((delta ** 2).sum())


# сам кристалл
a = 1
atoms = BodyCenteredCubic(directions=[[a,0,0], [0,a,0], [0,0,a]],
                          size=(6,6,6), symbol='Li',
                          latticeconstant=3.52)
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
epsilon = np.finfo(float).eps
# вычисление диагонали  --- с точностью эпсилон
pos = atoms.positions
m = []
for i in range(len(pos)):
    if (pos[i][0] == pos[i][1] == pos[i][2]) or ( (pos[i][0] == pos[i][1] == pos[i][2]) == 0 and pos[i][0] <= epsilon and pos[i][1] <= epsilon and pos[i][2] <= epsilon ):
        m.append(i)

#   координаты атомов на диагонале, где m их индексы (включая интерстишил)
n_arr = []
for i in range(len(m)):
    n_arr.append(pos[m[i]])
print(n_arr)

#   попытка сделать вывод системы в файл и в ламмпсе оптимизировать ее

lammpsdata.write_lammps_data('bulk.lammps-data', atoms)

#xyz.read_xyz('E_int_1_diag2.xyz', index=None)

path = "coords_int_diag2.txt"

# считывание данных с файла
with open(path) as textFile:
    arr = np.array([line.split() for line in textFile])
arr = np.delete(arr, 0, axis=1)
arr = arr.tolist()
arr_mod = np.array(list(np.float_(arr)))

# оптимизированный кристалл
pos_optimized = arr_mod
dim = dist

p1 = distance(pos_optimized[m[1]], pos_optimized[m[0]], dim)
p2 = distance(pos_optimized[m[2]], pos_optimized[m[1]], dim)
p3 = distance(pos_optimized[m[3]], pos_optimized[m[2]], dim)
p4 = distance(pos_optimized[m[4]], pos_optimized[m[3]], dim)
p5 = distance(pos_optimized[m[12]], pos_optimized[m[4]], dim)
p6 = distance(pos_optimized[m[11]], pos_optimized[m[12]], dim)
p7 = distance(pos_optimized[m[11]], pos_optimized[m[5]], dim)
p8 = distance(pos_optimized[m[5]], pos_optimized[m[6]], dim)
p9 = distance(pos_optimized[m[6]], pos_optimized[m[7]], dim)
p10 = distance(pos_optimized[m[8]], pos_optimized[m[7]], dim)
p11 = distance(pos_optimized[m[9]], pos_optimized[m[8]], dim)
p12 = distance(pos_optimized[m[10]], pos_optimized[m[9]], dim)
p13 = distance(pos_optimized[m[10]], pos_optimized[m[0]], dim)

diag_d = []
# diag_d.extend([p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0], p9[0], p10[0], p11[0], p12[0], p13[0]])
diag_d.extend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12])
print(diag_d)

#   ГРАФИК

print(len(diag_d))
x = np.arange(0, 12)
print(x)
plt.plot(x, diag_d)
plt.show()


# есть ли id атомов в асе