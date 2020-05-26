import numpy as np
from src.read_input import Atom
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc, dump_2dc, dump_1dr, dump_1dc
from module.mod_read import read_4dc
from module.atomic_hmat import atomic_make_hmtrx#(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
from src.atomic_subs   import atomic_state, decide_unitary_len 
from scipy.special import comb,perm
from module.atomic_basis import atomic_make_basis#(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)
from pg.PG import DPG, TranOrb, MBPG, MBPGsubs
from pg.PG_util import Point_ID,pmutt,decompose_vec,get_TranOrb_param
from pg.read_IR_DSG import *
from wanntb.tran import tran_op, tran_unitary
import copy
from module.atomic_angular import atomic_make_sp2np
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
hmat  =  atomic_make_hmtrx(atom1.norb,atom1.totncfgs,atom1.ncfgs,basis,invcd,invsn,atom1.onsite,cumat_t)


# make point group
# <just for the matrix representation(MRC: MRCar) in cartisian basis of generators of corresponding point group as well as irreps>
pgid  = Point_ID(atom1.point_group)
dpg71 = DPG(atom1.point_group,pgid)
sgid  = atom1.space_group 
dpg71.get_data(sgid)
#dpg71.groupclass_irrep()
#dpg71.show_attribute()
 
# <make matrix representation(MRO: MROrb) in d f orbitals space> 
npoly1, dim1, npower1, nfunc1 = get_TranOrb_param(atom1.soc_type)
Oprt_PG = TranOrb(dpg71.rep_vec,dpg71.rep_spin,npoly1,dim=dim1,npower=npower1,nfunc=nfunc1,shell=atom1.soc_type)
#Oprt_PG.show_attribute() if test else 0

# <transform MRO into natural basis(MRN)>
##print('umat_so_original:',Oprt_PG.umat_so)
umat_so_natural = tran_op(Oprt_PG.umat_so, atom1.amat) #umat_so_natrual  transforms in rows
#umat_so_natural  = Oprt_PG.umat_so #umat_so_natrual  transforms in rows


#umat_so_natural = tran_unitary(Oprt_PG.umat_so,atom1.amat,transpose=True) 

# construct projectors 
Basis_list = []
sta = int(0)
pg_manybody = MBPG(len(umat_so_natural),[])
pg_manybody.dim = atom1.ncfgs
for inn in range(atom1.nmin,atom1.nmax+1):
    print(40*'*')
    print('* nocc = ',inn)
    len_sp   = int(comb(atom1.norb,inn))

    basis1,invcd1,invsn1 = atomic_make_basis(atom1.norb,atom1.totncfgs,len_sp,atom1.norb,inn,inn,nstat)
    manybody_umat   = []
    unitary_len = perm(atom1.norb, inn, exact=True)
#   unitary_len =  decide_unitary
    cnt_op = int(-1)
    for imat in umat_so_natural:
        print('     ',30*'*')
        cnt_op += 1
        print('     * nop = ',cnt_op)
# the input unitary matrix of gw_make_newui should be transformed in columns, and output is also transformed in columns.
#       umat_mb_tmp1 = atomic_make_sp2np(atom1.norb, atom1.totncfgs, atom1.ncfgs, basis, invsn, invcd, imat)
        umat_mb_tmp = gw_make_newui(atom1.norb,len_sp,atom1.totncfgs,unitary_len,inn,inn,np.transpose(imat),basis1,invcd1,invsn1,1.0E-8)
# the manybody's unitary transform matrix should transfrom in columns like that in group theory
        manybody_umat.append(umat_mb_tmp) # 
# construct instance of Mb_Pg
#   pg_manybody = MBPG(len(manybody_umat),manybody_umat)
##pg_manybody.show_attribute()
#   pg_mb_sp = MBPGsubs(pg_manybody.nop,umat_sp)
    print(5*' ',20*'* * ')
    print('     * I start to check the unitarity of operatros *')
    for iop in range(len(manybody_umat)) :
        print(10*' ',5*'- ')
        check_u_tmp = np.dot(np.transpose(np.conjugate(manybody_umat[iop])),manybody_umat[iop])
        print(10*' ','diagonal elements for iop =',iop)
        if np.abs(np.sum(np.abs(check_u_tmp)) - check_u_tmp.shape[0]) < 1.0e-4:
            print(10*' ','   it is a identity matrix : Pass')
#       print(check_u_tmp.diagonal())
        print('')
    print(5*' ',20*' * * ')

    pg_mb_sp = MBPGsubs(len(manybody_umat),manybody_umat)
    pg_mb_sp.dim = len_sp
    pg_mb_sp.ham = hmat[sta:sta+len_sp,sta:sta+len_sp]
    pg_mb_sp.Cal_ReductRep(dpg71.irreps)
    pg_mb_sp.check_projectors()
    pg_mb_sp.collect_basis()#the eigenwaves arranges in rows
    pg_mb_sp.trans_ham()
    pg_mb_sp.diag_ham() 
    pg_mb_sp.check_basis_irreps1()# should be executed before pg_mb_sp.cal_degeneracy()
    pg_mb_sp.cal_degeneracy()
    dump_1dr(len_sp,pg_mb_sp.ham_eig,path='eig_n_'+str(inn)+'.dat',prec=1.0e-6)
    dump_2dc(len_sp,len_sp,pg_mb_sp.ham_evc,path='evc_n_'+str(inn)+'.dat',prec=1.0e-6)
    print(10*'* * ')
#   print('|*| final basis of irreps :\n',pg_mb_sp.allbasis['matrix'])
    print('|*| the label of final basis of irreps :\n',pg_mb_sp.allbasis['irreplabel'])
    print(10*'* * ')
#   pg_mb_sp.decompose_degenerate() 
    pg_mb_sp.check_basis_irreps()
#   print(5*' ','check for degeneracy of irreps :')
#   for j in range(len_sp):
#       cnt_t = 0
#       label_t = None
#       logic_t = True
#       for i in range(len_sp):
#           if np.abs(pg_mb_sp.ham_evc[i,j]) > 1.0E-6:
#               cnt_t += 1
#               if label_t == None:
#                   label_t = pg_mb_sp.allbasis['irreplabel'][i]
#               else:
#                   if label_t != pg_mb_sp.allbasis['irreplabel'][i]:
#                       logic_t = False
#       print(5*' ','j =',j,'deg=',cnt_t,'logic:',logic_t,'label = ',label_t)

    for irr in pg_mb_sp.irrep:
        print(5*' ','characters:\n',irr.characters.real)
        print(5*' ','multi for ',irr.label,'is ',irr.multi)
    # renew the start point of the next sub space  
    Basis_list.append(pg_mb_sp)
    sta += len_sp 

# transform hamiltonian into basis of irreps



