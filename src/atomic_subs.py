import numpy as np
from atomic_constants import * 
from scipy.special import comb, perm

def atomic_state(norbs):
    from scipy.special import comb
    nstat = np.zeros(norbs+1,dtype=np.int32)
    for i in range(norbs+1):
        nstat[i] = comb(norbs,i)
    return nstat

def decide_unitary_len(nocc,unitary_mat):
    '''
    aim : to decide the unitary_len needed by fortran warperred function gw_make_newui
    with the precision cutoff setting
    intput : 
            nocc : [integer] to decide the subspace
            unitary : [numpy.ndarray] the unitary matrix in single particle basis
    '''
    norb = unitary_mat.shape[0]
    norb_eff = int(0)
    norb_eff_set = []
    for iorb in range(norb):
        norb_eff_set[iorb] = sum(np.abs(unitary_mat[:,iorb]) > unitary_manybody_prec)
    norb_eff = int(np.max(norb_eff_set))

    unitary_len = perm(norb_eff,nocc) + int(100)
    return unitary_len

