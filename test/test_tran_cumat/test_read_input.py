from src.read_input import Atom
import numpy as np
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc
from module.mod_read import read_4dc


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
cumat=atomic_make_cumat(atom1.norb,atom1.int_val['U'],atom1.int_val['Up'],atom1.int_val['Jz'],atom1.int_val['Js'],atom1.int_val['Jp'])
#cumat_t = np.zeros((atom1.norb,atom1.norb),dtype=np.complex128)
cumat_t = atomic_tran_cumat(atom1.norb,atom1.amat,cumat)
dump_cumat_t = dump_4dc(atom1.norb,atom1.norb,atom1.norb,atom1.norb,cumat_t,'test-tumat.dat')

cumat_old = read_4dc(atom1.norb,atom1.norb,atom1.norb,atom1.norb,'/share/home/pengsy/test/test_rtgw2020/test/LaNiO2/gutz_j17_3/atom_1/solver.tumat.out')
print("abs(sum(cumat-cumat_old)) = ",np.abs(cumat_t-cumat_old).sum(axis=(0,1,2,3)))

