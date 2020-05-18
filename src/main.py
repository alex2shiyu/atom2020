import numpy as np
from src.read_input import Atom
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc
from module.mod_read import read_4dc
from module.atomic_hmat import atomic_make_hmtrx#(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
from src.atomic_subs   import atomic_state  
from scipy.special import comb
from module.atomic_basis import atomic_make_basis#(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)
from pg.PG import DPG, TranOrb, Mb_Pg
from pg.PG_util import Point_ID,pmutt,decompose_vec,get_TranOrb_param
from pg.read_IR_DSG import *
from wanntb.tran import tran_op


# debug
test = True
test_jy = False
test_jy_type = False

# construct atom object
atom1 = Atom.from_incar()
atom1.make_instance()

# construct cumat and transform it to natural basis
cumat   = atomic_make_cumat(atom1.norb,atom1.int_val['U'],atom1.int_val['Up'],atom1.int_val['Jz'],atom1.int_val['Js'],atom1.int_val['Jp'])
cumat_t = atomic_tran_cumat(atom1.norb,atom1.amat,cumat)

# make basis
nstat = atomic_state(atom1.norb)
print('nstat=\n',nstat) if test else 0
basis,invcd,invsn = atomic_make_basis(atom1.norb,atom1.totncfgs,atom1.ncfgs,atom1.norb,atom1.nmin,atom1.nmax,nstat)

# make hmat
# (norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
hmat  =  atomic_make_hmtrx(atom1.norb,atom1.totncfgs,atom1.ncfgs,basis,invcd,invsn,atom1.cfd_mat,cumat_t)



# make point group
# <just for the matrix representation(MRC: MRCar) in cartisian basis of generators of corresponding point group as well as irreps>
pgid  = Point_ID(atom1.point_group)
dpg71 = DPG(atom1.point_group,pgid)
sgid  = atom1.space_group 
dpg71.get_data(sgid)
#dpg71.groupclass_irrep()
dpg71.show_attribute()
 
# <make matrix representation(MRO: MROrb) in d f orbitals space> 
npoly1, dim1, npower1, nfunc1 = get_TranOrb_param(atom1.soc_type)
Oprt_PG = TranOrb(dpg71.rep_vec,dpg71.rep_spin,npoly1,dim=dim1,npower=npower1,nfunc=nfunc1,shell=atom1.soc_type)
Oprt_PG.show_attribute() if test else 0

# <transform MRO into natural basis(MRN)>
umat_so_natural = tran_op(Oprt_PG.umat_so, atom1.amat) 

# construct instance of Mb_Pg
pg_manybody = Mb_Pg(len(umat_so_natural),umat_so_natural,dpg71.irreps)
pg_manybody.show_attribute()

# construct projectors 

