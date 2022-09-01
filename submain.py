import os
from ase.io import lammpsdata
from ase.lattice.cubic import BodyCenteredCubic
import numpy as np
from matplotlib import pyplot as plt
from ase import io
import shutil
from ase.build.tools import sort


def distance(x0, x1, dimensions):
    delta = np.abs(x1 - x0)
    if delta[0] > (0.5 * dimensions):
        delta[0] = delta[0] - dimensions
    else:
        delta[0] = delta[0]
    if delta[1] > (0.5 * dimensions):
        delta[1] = delta[1] - dimensions
    else:
        delta[1] = delta[1]
    if delta[2] > (0.5 * dimensions):
        delta[2] = delta[2] - dimensions
    else:
        delta[2] = delta[2]
    return np.sqrt((delta ** 2).sum())


# сам кристалл
a = 1
nsize = 6
lattice_const = 3.52
atoms = BodyCenteredCubic(directions=[[a, 0, 0], [0, a, 0], [0, 0, a]],
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
pos_int = d[6].position / 10 + d[6].position
atoms.positions[-1] = pos_int
d = []
m = []
for atm in reversed(atoms):
    if np.std(atm.position) <= lattice_const / 100:
        d.append(atm)
        m.append(atm.index)

########################

DESTPATH = r"C:\Users\default.DESKTOP-J9F7KTV\OneDrive\Рабочий стол\GB_structures\CODE"

lammpsdata.write_lammps_data('bulk.lammps-data', atoms)
shutil.copy2(r"bulk.lammps-data", DESTPATH)
os.system(DESTPATH + "\run.bat")
# os.system("pause")
atoms_read = io.read(DESTPATH + "\E_int_1_diag2.lammps-dump-text")
pos_a_read = atoms_read.positions

#########################

arr = []
for i in range(len(m)):
    arr.append(pos_a_read[m[i]])
arr_mod = np.array(list(np.float_(arr)))
pos_sort = sort(arr_mod, tags=arr_mod[:, 2])
p = []
for i in range(1, len(m)):
    p.append(distance(pos_sort[i], pos_sort[i - 1], c))
p.append(distance(pos_sort[len(m) - 1], pos_sort[0], c))
print(p)

# plot
x = np.arange(0, len(m))
plt.plot(x, p, 'o-')
plt.show()
