import numpy as np
import math
from src.atomic_subs   import atomic_state  
from scipy.special import comb
from module.atomic_basis import atomic_make_basis
#(norbs, totncfgs, ncfgs, ntots, nmin, nmax, nstat, basis, invcd, invsn)


norbs    = np.int(10)
totncfgs = math.pow(2, norbs)
ncfgs    = np.int(comb(norbs,7) + comb(norbs,8) + comb(norbs,9))
ntots    = np.int(8)
nmin     = np.int(7)
nmax     = np.int(9)
nstat    = atomic_state(norbs)
print('nstat=\n',nstat)
basis,invcd,invsn = atomic_make_basis(norbs,totncfgs,ncfgs,ntots,nmin,nmax,nstat)
print('norbs=',norbs)
print('totncfgs=',totncfgs)
print('ncfgs=',ncfgs)
print('ntots=',ntots)
print('nmin=',nmin)
print('nmax=',nmax)
try:
    with open('test_atom.basis.dat','w') as f:
        for i in range(ncfgs):
            line='{:^8d}{:^8d}{:^8d}'.format(i+1,basis[i],invsn[basis[i]])
            f.write(line)
            for j in range(norbs):
                line1='{:>2d}'.format(invcd[j,i])
                f.write(line1)
            f.write('\n')
except IOError:
    print("open the file {} failed", 'test_atom.basis.dat')
