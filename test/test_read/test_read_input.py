from src.read_input import Atom
import numpy as np


atom1 = Atom.from_incar()
print('atom1:\n',atom1)
print('atom1.norb:\n',atom1.norb)
print('atom1.nmin:\n',atom1.nmin)
print('atom1.nmax:\n',atom1.nmax)
print('atom1.int_type:\n',atom1.int_type)
print('atom1.int_val:\n',atom1.int_val)
print('atom1.soc_type:\n',atom1.soc_type)
print('atom1.soc_val:\n',atom1.soc_val)
print('atom1.cfd_mat:\n',atom1.cfd_mat)
print('atom1.point_group:\n',atom1.point_group)
print('atom1.vpm_type:\n',atom1.vpm_type)
atom1.check_incar()
