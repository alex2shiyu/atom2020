#! /usr/bin/env python
import numpy as np
import copy
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
from module.atomic_angular import atomic_make_sp2np
from module.gw_make_new import gw_make_newui 
from timeit import default_timer as timer
from datetime import timedelta
from src.atomic_subs import atom_timer, show_header, show_subheader, show_error, show_success
from src.atomic_subs import atom_timer_Nsubs,show_subsubheader, show_sub3header

# time 
timedata = atom_timer()

# show header
show_header()

# construct atom object
timedata.reset()
show_subheader('Read Incar(atom instance)')
try:
    atom1 = Atom.from_incar()
    atom1.make_instance()
except:
    show_error('Read Incar')
else:
    show_success()
    print('\n\n')
timedata.atom = timedata.count()

# construct cumat and transform it to natural basis
timedata.reset()
show_subheader('Cumat')
try :
    cumat   = atomic_make_cumat(atom1.norb,atom1.int_val['U'],atom1.int_val['Up'],atom1.int_val['Jz'],atom1.int_val['Js'],atom1.int_val['Jp'],atom1.iprint)
    cumat_t = atomic_tran_cumat(atom1.norb,atom1.amat,cumat)
except :
    show_error('Cumat')
else:
    show_success()
    print('\n\n')
timedata.cumat = timedata.count()

# make basis
timedata.reset()
show_subheader('Basis')
try :
    nstat = atomic_state(atom1.norb)
    print('nstat=\n',nstat) if atom1.iprint >= 2 else 0
    basis,invcd,invsn = atomic_make_basis(atom1.norb,atom1.totncfgs,atom1.ncfgs,atom1.norb,atom1.nmin,atom1.nmax,nstat,2)
except :
    show_error('Basis')
else:
    print('basis:\n',basis) if atom1.iprint >= 2 else 0
    show_success()
    print('\n\n')
timedata.basis = timedata.count()

# make hmat
timedata.reset()
show_subheader('Hamiltonian')
try :
    hmat  =  atomic_make_hmtrx(atom1.norb,atom1.totncfgs,atom1.ncfgs,basis,invcd,invsn,atom1.onsite,cumat_t,atom1.iprint)
except :
    show_error('Hamiltonian')
else:
    show_success()
    print('\n\n')
timedata.ham = timedata.count()


# make point group
# <just for the matrix representation(MRC: MRCar) in cartisian basis of generators of corresponding point group as well as irreps>
timedata.reset()
show_subheader('PG(instance)')
try :
    pgid  = Point_ID(atom1.point_group)
    dpg = DPG(atom1.point_group,pgid)
    sgid  = atom1.space_group 
    dpg.get_data(sgid)
    dpg.show_attribute() if atom1.iprint >= 2 else 0
except :
    show_error('PG')
else:
    show_success()
    print('\n\n')
timedata.dpg = timedata.count()
 
# <make matrix representation(MRO: MROrb) in d f orbitals space> 
timedata.reset()
show_subheader('Trans matrix in atomic orbitals')
try :
    npoly1, dim1, npower1, nfunc1 = get_TranOrb_param(atom1.soc_type)
    Oprt_PG = TranOrb(dpg.rep_vec,dpg.rep_spin,npoly1,dim=dim1,npower=npower1,nfunc=nfunc1,shell=atom1.soc_type,iprint=atom1.iprint)
    Oprt_PG.show_attribute() if atom1.iprint == 3 else 0
    Oprt_PG.check_symm_crystal(atom1.cfd_mat)
    Oprt_PG.check_symm_soc(atom1.soc_mat)
except :
    show_error('TranOrb')
else:
    show_success()
    print('\n\n')
timedata.TranOrb = timedata.count()

# <transform MRO into natural basis(MRN)>
# should notice that this may be wrong
#umat_so_natural = tran_op(Oprt_PG.umat_so, atom1.amat) #umat_so_natrual  transforms in rows
#umat_so_natural  = Oprt_PG.umat_so #umat_so_natrual  transforms in rows


#umat_so_natural = tran_unitary(Oprt_PG.umat_so,atom1.amat,transpose=True) 

# construct projectors 
timedata.reset()
show_subheader('DPG reduction in every N subspace')
Basis_list = []
sta = int(0)
pg_manybody = MBPG(len(Oprt_PG.umat_so),[],atom1.iprint)
pg_manybody.dim = atom1.ncfgs
for inn in range(atom1.nmin,atom1.nmax+1):
    timeNsubs = atom_timer_Nsubs() 
    timeNsubs.noc = inn
    show_subsubheader('Nocc = ' + str(inn))
    len_sp   = int(comb(atom1.norb,inn))
#
    timeNsubs.reset()
    basis1,invcd1,invsn1 = atomic_make_basis(atom1.norb,atom1.totncfgs,len_sp,atom1.norb,inn,inn,nstat,1)
    timeNsubs.basis = timeNsubs.count()

    print('* basis:\n',basis1) if atom1.iprint >= 2 else 0
    manybody_umat   = []
    unitary_len = perm(atom1.norb, inn, exact=True)
#   unitary_len =  decide_unitary
    cnt_op = int(-1)
    timeNsubs.reset()
    for imat in Oprt_PG.umat_so:
        cnt_op += 1
        if atom1.iprint >= 3 :
            print('     ',30*'*')
            print('     * nop = ',cnt_op)
# the input unitary matrix of gw_make_newui should be transformed in columns, and output is also transformed in columns.
#       umat_mb_tmp1 = atomic_make_sp2np(atom1.norb, atom1.totncfgs, atom1.ncfgs, basis, invsn, invcd, imat)
        umat_mb_tmp = gw_make_newui(atom1.norb,len_sp,atom1.totncfgs,unitary_len,inn,inn,np.transpose(imat),basis1,invcd1,invsn1,atom1.iprint,1.0E-8)
# the manybody's unitary transform matrix should transfrom in columns like that in group theory
        manybody_umat.append(umat_mb_tmp) # 
    timeNsubs.operators = timeNsubs.count()

    pg_mb_sp = MBPGsubs(len(manybody_umat),manybody_umat,atom1.iprint)
    pg_mb_sp.dim = len_sp
    pg_mb_sp.ham = hmat[sta:sta+len_sp,sta:sta+len_sp]
#
    timeNsubs.reset()
    show_sub3header('Check Ham(symm)')
    pg_mb_sp.check_symm_ham()
    timeNsubs.check_ham = timeNsubs.count()
#
    timeNsubs.reset()
    show_sub3header('Reduction')
    pg_mb_sp.Cal_ReductRep(dpg.irreps)
    timeNsubs.reduct = timeNsubs.count()
#
    timeNsubs.reset()
    show_sub3header('Check projectors')
    pg_mb_sp.check_projectors()
    timeNsubs.check_pro = timeNsubs.count()
#
    show_sub3header('Collect basis of irreps')
    pg_mb_sp.collect_basis()#the eigenwaves arranges in rows
#
    timeNsubs.reset()
    show_sub3header('Check validity of finnal basis of irreps')
    pg_mb_sp.check_irrepbasis_final()
    timeNsubs.check_irrepbasis = timeNsubs.count()

#test
    if atom1.iprint == 3 :
        pg_mb_sp.trans_ham(trans=False)
        pg_mb_sp.diag_ham(trans=False) 
        dump_1dr(len_sp,pg_mb_sp.ham_eig,path='eig_n_before_tran_'+str(inn)+'.dat',prec=1.0e-6)
        dump_2dc(len_sp,len_sp,pg_mb_sp.ham_evc,path='evc_n_before_tran_'+str(inn)+'.dat',prec=1.0e-6)
#test
#
    timeNsubs.reset()
    show_sub3header('trans ham to irreps basis')
    pg_mb_sp.trans_ham(trans=True)
    show_sub3header('diagonalize ham')
    pg_mb_sp.diag_ham(trans=True) 
    timeNsubs.ham_trandiag = timeNsubs.count()

#
    timeNsubs.reset()
    show_sub3header('check (degenerate space) <----> (irrep space)')
    pg_mb_sp.check_basis_irreps1()# should be executed before pg_mb_sp.cal_degeneracy()
    show_sub3header('check (eigenwave) <----> (column of irrep)')
    pg_mb_sp.check_basis_irreps()
    show_sub3header('cal degeneracy of ham')
    timeNsubs.check_eigenwave_irrep = timeNsubs.count()
#
    timeNsubs.reset()
    pg_mb_sp.cal_degeneracy()
    if atom1.iprint == 3:
        print(10*'* * ')
        print('|*| the label of final basis of irreps :\n',pg_mb_sp.allbasis['irreplabel'])
        print(10*'* * ')
    show_sub3header('decompose eigenwaves in every degenerate space')
    pg_mb_sp.decompose_degenerate() 
    dump_1dr(len_sp,pg_mb_sp.ham_eig,path='eig_n_'+str(inn)+'_after.dat',prec=1.0e-6)
    dump_2dc(len_sp,len_sp,pg_mb_sp.ham_evc,path='evc_n_'+str(inn)+'_after.dat',prec=1.0e-6)
    show_sub3header('check (eigenwave) <----> (column of irrep) [after decomposition]')
    pg_mb_sp.check_basis_irreps()
    timeNsubs.decompose = timeNsubs.count()
#
    print('\n\n') 
    for irr in pg_mb_sp.irrep:
        print(5*' ','characters:\n',irr.characters.real)
        print(5*' ','multi for ',irr.label,'is ',irr.multi)
    # renew the start point of the next sub space  
    Basis_list.append(pg_mb_sp)
    sta += len_sp 
#
    timedata.Nsubs.append(timeNsubs)

timedata.Nsubstot = timedata.count()

# transform hamiltonian into basis of irreps

# show info of elapsed time
timedata.timeover()
timedata.show()
