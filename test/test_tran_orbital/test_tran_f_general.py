import numpy as np
from pg.PG import TranOrb
from pg.PG_util import pmutt,decompose_vec

npoly=np.array([3,3,3,2,1,2,2],dtype=np.int32)
# C4
#umat = [np.array([[0,1,0],[-1,0,0],[0,0,1]],dtype=np.float64)]
# P(inversion)
umat = [np.array([[-1,0,0],[0,-1,0],[0,0,-1]],dtype=np.float64)]
# C_{2*pi/3}
#umat = [np.array([[-0.5,np.sqrt(3)/2,0],[-np.sqrt(3)/2,-0.5,0],[0,0,1]],dtype=np.float64)]
t2g_test = TranOrb(umat,npoly,dim=3,npower=3,nfunc=7,shell='f')
t2g_test.show_attribute()

#t2g_test.def_func()
#t2g_test.make_mapbasis()
#t2g_test.make_vec_oldbasis()
##print('dim=',t2g_test.dim)
##print('npower=',t2g_test.npower)
##print('npoly=',t2g_test.npoly)
##print('nfunc=',t2g_test.nfunc)
##print('umat=',t2g_test.umat)
##print('shell=',t2g_test.shell)
##print('func_o=',t2g_test.func_o)
##print('mapbasis=',t2g_test.mapbasis)
##print('umat_new=',t2g_test.umat_new)
#t2g_test = TranOrb(umat,npoly,3,2,5,shell='d')
#t2g_test.def_func()
#t2g_test.make_mapbasis()
#t2g_test.make_vec_oldbasis()
#print('dim=',t2g_test.dim)
#print('npower=',t2g_test.npower)
#print('npoly=',t2g_test.npoly)
#print('nfunc=',t2g_test.nfunc)
#print('umat=',t2g_test.umat)
#print('shell=',t2g_test.shell)
#print('func_o=',t2g_test.func_o)
#print('mapbasis=',t2g_test.mapbasis)
#t2g_test.make_func_new()
#print('umat_new=',t2g_test.umat_new)





