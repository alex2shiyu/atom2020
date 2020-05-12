import numpy as np

def atomic_state(norbs):
    from scipy.special import comb
    nstat = np.zeros(norbs+1,dtype=np.int32)
    for i in range(norbs+1):
        nstat[i] = comb(norbs,i)
    return nstat

