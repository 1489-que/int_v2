import os
from ase.io import lammpsdata
from ase.lattice.cubic import BodyCenteredCubic
import numpy as np
from matplotlib import pyplot as plt
from ase import io
from matplotlib.backends.backend_pdf import PdfPages
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
nsize = 21
lattice_const = 3.52
atoms = BodyCenteredCubic(directions=[[a, 0, 0], [0, a, 0], [0, 0, a]],
                          size=(nsize, nsize, nsize), symbol='Li',
                          latticeconstant=lattice_const)

c = atoms.cell[0][0]

########################

pos = atoms.positions
atoms.append('Li')

del_index = nsize

temp = lattice_const / 6

# диагональ
d = []
m = []
for atm in atoms:
    if np.std(atm.position) <= lattice_const / 100:
        d.append(atm)
        m.append(atm.index)

# print(len(m))
pos_int = temp + d[del_index].position
pos_b_int = d[del_index].position - temp
# print(pos_b_int, ' ', pos_int, ' ', d[del_index - 1], d[del_index + 1])
del (d[del_index])
atoms.positions[-1] = pos_int
atoms.positions[-1] = pos_b_int

d = []
m = []
for atm in reversed(atoms):
    if np.std(atm.position) <= lattice_const / 100:
        d.append(atm)
        m.append(atm.index)

########################

DESTPATH = r"C:\Users\default.DESKTOP-J9F7KTV\OneDrive\PROJ\CODE\\"
lammpsdata.write_lammps_data(DESTPATH + "\\bulk.lammps-data", atoms)
os.chdir(DESTPATH)
os.system("run.bat")
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

# отклонение от идеального значения

dev_p = np.array(p) - lattice_const / 2 * 3**0.5
print(dev_p)

# plot
x = np.arange(0, len(m))
fig, ax = plt.subplots(figsize=(8, 6))
fig.subplots_adjust(bottom=0.15, left=0.2)
pdf = PdfPages('bulks' + str(nsize) + '.pdf')
plt.plot(x, dev_p, 'o-')

ax.set_xlabel('Atom id')
ax.set_ylabel('MD [Å]')
pdf.savefig()
pdf.close()

plt.title(str(nsize) + '×' + str(nsize) + '×' + str(nsize))

plt.show()
