import numpy as np
from pg.PG import TranOrb
from pg.PG_util import pmutt,decompose_vec

npoly=np.array([1,1,1],dtype=np.int32)
umat = [np.array([[0,1,0],[-1,0,0],[0,0,1]],dtype=np.float64)]
t2g_test = TranOrb(umat,npoly,3,2,3,shell='t2g')
t2g_test.def_func()
t2g_test.make_mapbasis()
print('dim=',t2g_test.dim)
print('npower=',t2g_test.npower)
print('npoly=',t2g_test.npoly)
print('nfunc=',t2g_test.nfunc)
print('umat=',t2g_test.umat)
print('shell=',t2g_test.shell)
print('func_o=',t2g_test.func_o)
print('umat_mb=',t2g_test.umat_mb)
print('mapbasis=',t2g_test.mapbasis)


# test pmutt
newstruct = np.array([[1,2,0],[1,2,0],[1,2,3]],dtype=np.float64)

pmutt_basis , pmutt_coeff = pmutt(newstruct)
print('newstruct : \n',newstruct)
print('pmutt_basis : \n',pmutt_basis)
print('pmutt_coeff : \n',pmutt_coeff)
basis_vec  = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype=np.float64)
target_vec = np.array([4.5,-2.1,0.1],dtype=np.float64)
rep_new = decompose_vec(basis_vec,target_vec)
print('rep_new\n',rep_new)

