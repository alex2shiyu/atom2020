import numpy as np
from src.read_input import Atom
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc, dump_2dc
from module.mod_read import read_4dc
from module.atomic_hmat import atomic_make_hmtrx#(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
from src.atomic_subs   import atomic_state  
from scipy.special import comb
from module.atomic_basis import atomic_make_basis#(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)
from pg.PG import DPG, TranOrb, MBPG, MBPGsubs
from pg.PG_util import Point_ID,pmutt,decompose_vec,get_TranOrb_param
from pg.read_IR_DSG import *
from wanntb.tran import tran_op
import copy
#from module.atomic_angular import atomic_make_sp2np
from module.gw_make_new import gw_make_newui 


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
##dpg71.show_attribute()
 
# <make matrix representation(MRO: MROrb) in d f orbitals space> 
npoly1, dim1, npower1, nfunc1 = get_TranOrb_param(atom1.soc_type)
Oprt_PG = TranOrb(dpg71.rep_vec,dpg71.rep_spin,npoly1,dim=dim1,npower=npower1,nfunc=nfunc1,shell=atom1.soc_type)
##Oprt_PG.show_attribute() if test else 0

# <transform MRO into natural basis(MRN)>
print('umat_so_original:',Oprt_PG.umat_so)
umat_so_natural = tran_op(Oprt_PG.umat_so, atom1.amat) 
manybody_umat   = []
for imat in umat_so_natural:
#for imat in Oprt_PG.umat_so:
#   umat_mb_tmp = atomic_make_sp2np(atom1.norb, atom1.totncfgs, atom1.ncfgs, basis, invsn, invcd, imat)
    umat_mb_tmp = gw_make_newui(atom1.norb,atom1.ncfgs,atom1.totncfgs,atom1.nmin,atom1.nmax,imat,basis,invcd,invsn,1.0E-8)
    manybody_umat.append(umat_mb_tmp)
# construct instance of Mb_Pg
pg_manybody = MBPG(len(manybody_umat),manybody_umat)
##pg_manybody.show_attribute()

# construct projectors 
sta = int(0)
for inn in range(atom1.nmin,atom1.nmax+1):
    # _sp means <sub space>
    len_sp   = int(comb(atom1.norb,inn))
    print('sta: ',sta,'len: ',len_sp)
    umat_sp  = []
    for iop in range(pg_manybody.nop):
        umat_tmp = pg_manybody.op[iop][sta:sta+len_sp,sta:sta+len_sp]
#       filename = 'test_opmatrx_'+str(iop+1)+'.dat'
#       dump_2dc(umat_tmp.shape[0],umat_tmp.shape[1],umat_tmp,filename)
#       print('umat.shape  :',umat_tmp.shape)
#       print('trace <',iop+1,'>  :',np.trace(umat_tmp))
        umat_sp.append(umat_tmp)
    pg_mb_sp = MBPGsubs(pg_manybody.nop,umat_sp)
    pg_mb_sp.Cal_ReductRep(dpg71.irreps)
    for irr in pg_mb_sp.irrep:
        print('characters:\n',irr.characters.real)
        print('multi for ',irr.label,'is ',irr.multi)
    # renew the start point of the next sub space  
    sta += len_sp 



