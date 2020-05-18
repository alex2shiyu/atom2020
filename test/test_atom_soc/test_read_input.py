from src.read_input import Atom
import numpy as np
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc
from module.mod_read import read_4dc
from module.atomic_hmat import atomic_make_hmtrx#(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
from src.atomic_subs   import atomic_state  
from scipy.special import comb
from module.atomic_basis import atomic_make_basis#(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)

# construct atom object
atom1 = Atom.from_incar()
atom1.make_ncfgs()
atom1.check_incar()
atom1.make_soc()
atom1.show_incar()
#
totncfgs = atom1.make_totncfgs()

# construct cumat and transform it to natural basis
cumat   = atomic_make_cumat(atom1.norb,atom1.int_val['U'],atom1.int_val['Up'],atom1.int_val['Jz'],atom1.int_val['Js'],atom1.int_val['Jp'])
cumat_t = atomic_tran_cumat(atom1.norb,atom1.amat,cumat)

# make basis
nstat = atomic_state(atom1.norb)
basis,invcd,invsn = atomic_make_basis(atom1.norb,totncfgs,atom1.ncfgs,atom1.norb,atom1.nmin,atom1.nmax,nstat)

# make umat

# make hmat
# (norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
hmat  =  atomic_make_hmtrx(atom1.norb,totncfgs,atom1.ncfgs,basis,invcd,invsn,atom1.cfd_mat,cumat_t)


