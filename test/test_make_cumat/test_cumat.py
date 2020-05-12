#! /usr/bin/env python
from module.atomic_stream import atomic_make_cumat
import numpy as np
from module.mod_dump import dump_4dc
from module.mod_read import read_4dc

norb = 10
U  = 5.0
Up = 3.2
J  = 0.9

cumat = atomic_make_cumat(norb,U,Up,J,J,J)
cumat_read = read_4dc(norb,norb,norb,norb,'/share/home/pengsy/test/test_rtgw2020/test/LaNiO2/gutz_j17_3/atom_1/solver.cumat.out')
hd_dump = dump_4dc(norb,norb,norb,norb,cumat,'test_LaNiO2_cumat.dat')

print("abs(sum(cumat-cumat_read)) = ",np.abs(cumat-cumat_read).sum(axis=(0,1,2,3)))
