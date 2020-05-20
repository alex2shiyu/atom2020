from pg.PG_util import isOrthogonal
import numpy as np

a  = np.array([1,-1,0],dtype=np.complex128)
b1 = np.array([1,1,1],dtype=np.complex128)
b2 = [b1,np.array([0,0,1],dtype=np.complex128)]
b3 = [b1,np.array([0,0,1],dtype=np.complex128),np.array([0,1,1],dtype=np.complex128)]

print('is a orthogonal to b1? ',isOrthogonal(a,b1))
print('is a orthogonal to b2? ',isOrthogonal(a,b2))
print('is a orthogonal to b3? ',isOrthogonal(a,b3))
