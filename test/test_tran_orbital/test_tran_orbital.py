import numpy as np
from pg.PG import TranOrb


npoly=np.array([1,1,1],dtype=np.int32)
umat = np.array([[0,1,0],[-1,0,0],[0,0,1]],dtype=np.float64)
t2g_test = TranOrb(umat,npoly,3,2,3,shell='t2g')
t2g_test.def_func()
print('dim=',t2g_test.dim)
print('npower=',t2g_test.npower)
print('npoly=',t2g_test.npoly)
print('nfunc=',t2g_test.nfunc)
print('umat=',t2g_test.umat)
print('shell=',t2g_test.shell)
print('func_o=',t2g_test.func_o)
print('umat_mb=',t2g_test.umat_mb)






