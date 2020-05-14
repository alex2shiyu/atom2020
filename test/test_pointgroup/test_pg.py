from src.read_input import Atom
import numpy as np
from module.atomic_stream import atomic_make_cumat,atomic_tran_cumat#(norbs, amtrx, cumat, cumat_t)
from module.mod_dump import dump_4dc
from module.mod_read import read_4dc
from module.atomic_hmat import atomic_make_hmtrx#(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, hmat)
from src.atomic_subs   import atomic_state  
from scipy.special import comb
from module.atomic_basis import atomic_make_basis#(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)
from pg.PG import DPG
from pg.read_IR_DSG import *

# construct atom object
atom1 = Atom.from_incar()
atom1.make_ncfgs()
atom1.check_incar()
atom1.show_incar()

# test or not
test = True
test_jy = False
test_jy_type = False

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

# test reading pg data from yijiang's database
dpg71 = DPG('D2h',pgid=8)
print("just test the DPG:\n",dpg71) if test else 0
print("dpg71.name_sc",dpg71.name_sc) if test else 0
print("dpg71.pgid",dpg71.pgid) if test else 0

# pull data from yijiang's code
for gid in range(71,72):
    lgrps = loadIR(gid)
    print(gid)
    for grp in lgrps:
        print(grp.klabel) if test_jy else 0
        print(grp.rotC) if test_jy else 0
        print(grp.tauC) if test_jy else 0
        print('generators:') if test_jy else 0
        print(grp.generators) if test_jy else 0
        
        print("grp.klabel",type(grp.klabel))  if test_jy_type else 0
        print("grp.rotC",type(grp.rotC))  if test_jy_type else 0
        print("grp.tauC",type(grp.tauC))  if test_jy_type else 0
        print("grp.generators",type(grp.generators))  if test_jy_type else 0
        for ir in grp.irreps:
            print(ir.label) if test_jy else 0
            print(ir.matrices) if test_jy else 0
            print(ir.characters) if test_jy else 0
            print("ir.label",type(ir.label))  if test_jy_type else 0
            print("ir.matrices",type(ir.matrices))  if test_jy_type else 0
            print("ir.characters",type(ir.characters))   if test_jy_type else 0
for grp in lgrps:
    if grp.klabel == 'GM' :
        dpg71.nrank    = len(grp.rotC)
        dpg71.rep_vec  = grp.rotC
        dpg71.rep_spin = grp.su2s
        dpg71.irreps   = grp.irreps
        print('dpg71.nrank:\n',dpg71.nrank)
        print('dpg71.rep_vec:\n',dpg71.rep_vec)
        print('dpg71.rep_spin:\n',dpg71.rep_spin)
        for ir in dpg71.irreps:
#           print('dpg71.irreps.label:\n',ir.label)
#           print('dpg71.irreps.matrices:\n',ir.matrices)
#           print('dpg71.irreps.characters:\n',ir.characters[0])
            print(ir.label,ir.characters.real)

# define generators
dpg71.gen_pos = np.array([0,1,2,3,4,5,6,7,8,12],dtype=np.int32)
dpg71.gen_mul = np.array([1,2,2,2,1,2,2,2,1,1],dtype=np.int32)
print('dpg71.gen_pos\n',dpg71.gen_pos)
print('dpg71.gen_mul\n',dpg71.gen_mul)

# check the orthogonality of irreps
for ir in dpg71.irreps:
    character_reduce = np.array([ir.characters[dpg71.gen_pos[i]].real for i in range(len(dpg71.gen_pos))])
    print('reduced character for ',ir.label ,'th irreps:\n',character_reduce)
