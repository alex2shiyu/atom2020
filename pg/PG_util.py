from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
import os
import sys
import copy


def check_ham_wave(ham,eig,wave):
    '''
    aim : check whether the three inputs are consistent with each other.
    input : 
            ham  : 2d numpy.ndarray 
            eig  : 1d numpy.ndarray
            wave : 2d numpy.ndarray with "i"th columns corrosponding to "i"th eigen 
    '''
    print('check eigenfunctions :\n')
    dim = ham.shape[0]
    istrue = [True]
    for j in range(dim):
        for i in range(dim):
            wave_new = np.einsum('ij,j->i',ham,wave[:,j])
            value    = np.dot(np.conjugate(wave[:,i]),wave_new)
            if i == j  and np.abs(eig[j] - value) < 1.0e-10:
                istrue.append(True)
            elif i != j and np.abs(value) < 1.0e-10:
                istrue.append(True)
            else:
                print(5*' ','<',i,'|',j,'> =',value)
                istrue.append(False)
    if all(istrue) :
        print(5*' ','success')
    else :
        print(5*' ','failed')

#def decompose_vec1(basis):
#    '''
#    aim : to decompose [[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]
#              ------>  [[a'11, 0,  0],[  0,a'22, 0],[  0,  0,a'33]]
#              where aij is one dimension vector
#    '''
#    if len(basis) == 1:
#        return basis
#    else:
def find_rep(basis,op):
    '''
    aim    : find the representation matrix for <op> in the set of <basis> 
    op     : [numpy.ndarray] operators in a certain space, should transform functions in columns 
    basis  : [list of 1D-numpy.ndarray] make it a close space for <op>: should transfrom in rows
    output : rep transform in rows
    '''
    norb = len(basis)
    ndim = op.shape[0]
    rep  = np.zeros((norb,norb),dtype=np.complex128)
    for irow in range(norb):
        for icol in range(norb):
            # <1|O|1>=a_11,  <2|O|1> = a_12 
            rep[irow,icol] = np.einsum('i,ij,j->',np.conjugate(basis[icol]),op,basis[irow])
    # make sure the basis is colse
    is_close = []
    for irow in range(norb):
        basis_after_trans  = np.einsum('ij,j->i',op,basis[irow])
        basis_recover = np.zeros(ndim,dtype=np.complex128)
        for iorb in range(norb):
            basis_recover += rep[irow,iorb]*basis[iorb]
        error = np.sum(np.abs(basis_recover - basis_after_trans))
        if error < 1.0e-6 :
            is_close.append(True)
        else:
            print('ERROR: <find_rep>')
            print('basis[i],i=',irow)
            print('    \n',basis[irow])
            print('op\n',op)
            print('basis_after_trans:\n',basis_after_trans)
            print('basis_recover:\n',basis_recover)
            print('error = ',error)
            is_close.append(False)
    if not all(is_close):
        raise ValueError('ERROR: basis is not complete : ', error)
    else:
        print('<find_rep> success in <check_irrepbasis_final>')
        return rep



def check_unitary(manybody_umat):
    '''
    manybody_umat : list of numpy.ndarray
    '''
    print('>>>> check manybody unitary matrix :')
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


def RepBasisNorm(basis):
    '''
    normalize a basis for a certain rep which is actually a normalization of a vector
    '''
    if isinstance(basis,np.ndarray):
        return basis/np.linalg.norm(basis)
    elif isinstance(basis,list):
        tmp_list = []
        for iba in basis:
            tmp_list.append(iba/np.linalg.norm(iba))
        return tmp_list
    else:
        raise ValueError("RepBasisNorm just support <list of numpy.ndarray> and <numpy.ndarray> not yet this type: ",type(basis))


def isindependent(basis,basis_set):
    '''
    to judge whether the function <basis> is independent to the other basis set
    input : 
        basis : 1D numpy-ndarray
        basis_set : list of 1D numpy-ndarray or just numpy-ndarray
    '''
    print('  ')
    print('   entering isindependent ....')
    if isinstance(basis_set, list):
        isOrtho = True
        param_tmp = np.zeros(len(basis_set),dtype=np.complex128)
        basis_t   = copy.deepcopy(basis)
        print('      I will show1 the components ...')
        for ir in range(len(basis_set)):
            param_tmp[ir] = np.dot(np.conjugate(basis_set[ir]),basis)/np.dot(np.conjugate(basis_set[ir]),basis_set[ir])
            print('       i=',ir,'components = ',param_tmp[ir])
            basis_t = basis_t - param_tmp[ir] * basis_set[ir]
            print('basis    :\n',basis)
            print('basis_set:\n',basis_set)
            print('para     :\n',param_tmp)
            print('results  :\n',basis_t)
        if np.sum(np.abs(basis_t)) < 1.0 : 
            isOrtho = False
            print('         Sad (isindependent value) : ',np.sum(np.abs(basis_t)))
        else :
            print('     Not Sad (isindependent value) : ',np.sum(np.abs(basis_t)))
        print("")
        return [np.sum(np.abs(basis_t))], [isOrtho]
    elif isinstance(basis_set,np.ndarray):
        print('      I will show2 the components ...')
        isOrtho = True
        basis_t   = copy.deepcopy(basis)
        param_tmp = np.dot(np.conjugate(basis_set),basis)
        basis_t   = basis_t - param_tmp * basis_set
        if np.sum(np.abs(basis_t)) < 0.1 : 
            isOrtho = False
            print('    Sad (isindependent value) : ',np.sum(np.abs(basis_t)))
        else :
            print('Not Sad (isindependent value) : ',np.sum(np.abs(basis_t)))
        return [np.sum(np.abs(basis_t))],[isOrtho]
    else:
        raise ValueError("Just support list of numpy-1Darray or numpy-1Darray not yet :",type(basis_set))

def isOrthogonal(basis,basis_set):
    '''
    to judge whether the function <basis> is orthogonal to the other basis set
    input : 
        basis : 1D numpy-ndarray
        basis_set : list of 1D numpy-ndarray or just numpy-ndarray
    '''
    if isinstance(basis, np.ndarray) and isinstance(basis_set, list):
        isOrtho = [True for i in range(len(basis_set))]
        print("")
        for ir in range(len(basis_set)):
            print(5*'-')
            print('max_value:',np.max(np.abs(basis)),np.max(np.abs(basis_set[ir])))
            if np.abs(np.dot(np.conjugate(basis),basis_set[ir])) > 1.0E-6:
                print('Sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(basis),basis_set[ir])))
                isOrtho[ir] = False
            else :
#           print('wave1:\n',basis)
#           print('wave2:\n',basis_set)
                print('Not sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(basis),basis_set[ir])))
        return isOrtho
    elif isinstance(basis, np.ndarray) and isinstance(basis_set,np.ndarray):
        isOrtho = True
        print(5*'-')
        print('max_value:',np.max(np.abs(basis)),np.max(np.abs(basis_set)))
        if np.abs(np.dot(np.conjugate(basis),basis_set)) > 1.0E-6:
            print('Sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(basis),basis_set)))
            isOrtho = False
        else:
            print('Not sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(basis),basis_set)))
        print('wave1:\n',basis)
        print('wave2:\n',basis_set)
        print('isOrtho value:',np.abs(np.dot(np.conjugate(basis),basis_set)))
        return isOrtho
    elif isinstance(basis,list) and isinstance(basis_set,list):
        print('basis:\n',basis)
        print('basis_set:\n',basis_set)
        isOrtho = True
        cnt1 = 0
        cnt2 = 0
        for ir1 in basis:
            cnt1 += 1
            cnt2 = 0
            for ir2 in basis_set:
                cnt2 += 1
                print(5*'-')
                print('ir1=',cnt1,' ir2=',cnt2)
                print('max_value:',np.max(np.abs(ir1)),np.max(np.abs(ir2)))
                print('ir1=\n',ir1)
                print('ir2=\n',ir2)
                if np.array_equal(ir1,ir2):
                    print('-----> equal')
                    continue 
                elif np.abs(np.dot(np.conjugate(ir1),ir2)) > 1.0E-6:
                    print('Sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(ir1),ir2)))
                    isOrtho = False
                else:
                    print('Not sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(ir1),ir2)))
                print('')
        return isOrtho
    elif isinstance(basis,list) and isinstance(basis_set,np.ndarray):
        isOrtho = True
        for ir1 in basis:
            print(5*'-')
            print('max_value:',np.max(np.abs(ir1)),np.max(np.abs(basis_set)))
            if np.array_equal(ir1,basis_set):
                print('-----> equal')
                continue 
            elif np.abs(np.dot(np.conjugate(ir1),basis_set)) > 1.0E-6:
                print('Sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(ir1),basis_set)))
                isOrtho = False
            else:
                print('Not sad (isOrtho value) : ',np.abs(np.dot(np.conjugate(ir1),basis_set)))
        return isOrtho
    else:
        raise ValueError("Just support list of numpy-1Darray or numpy-1Darray not yet :",type(basis_set))

def get_TranOrb_param(orb_case):
    if orb_case == 'f' :
        npoly1=np.array([3,3,3,2,1,2,2],dtype=np.int32)
        dim1    = 3
        npower1 = 3
        nfunc1  = 7
    elif orb_case == 'd' :
         npoly1=np.array([3,1,1,2,1],dtype=np.int32)
         dim1    = 3
         npower1 = 2
         nfunc1  = 5
    elif orb_case == 't2g' :
         npoly1=np.array([1,1,1],dtype=np.int32)
         dim1    = 3
         npower1 = 2
         nfunc1  = 3
    return npoly1, dim1, npower1, nfunc1

def Point_ID(sch):
    '''
    convert schoenflies to serial number ID
    '''
    if sch == 'C1' or '1':
        return int(1)
    elif sch == 'Ci' or '-1':
        return int(2)
    elif sch == 'C2' or '2':
        return int(3)
    elif sch == 'Cs' or 'm':
        return int(4)
    elif sch == 'C2h' or '2/m':
        return int(5)
    elif sch == 'D2' or '222':
        return int(6)
    elif sch == 'C2v' or 'mm2':
        return int(7)
    elif sch == 'D2h' or 'mmm':
        return int(8)
    elif sch == 'C4' or '4':
        return int(9)
    elif sch == 'S4' or '-4':
        return int(10)
    elif sch == 'C4h' or '4/m':
        return int(11)
    elif sch == 'D4' or '422':
        return int(12)
    elif sch == 'C4v' or '4mm':
        return int(13)
    elif sch == 'D2d' or '-42m':
        return int(14)
    elif sch == 'D4h' or '4/mmm':
        return int(15)
    elif sch == 'C3' or '3':
        return int(16)
    elif sch == 'C3i' or '-3':
        return int(17)
    elif sch == 'D3' or '32':
        return int(18)
    elif sch == 'C3v' or '3m':
        return int(19)
    elif sch == 'D3d' or '-3m':
        return int(20)
    elif sch == 'C6' or '6':
        return int(21)
    elif sch == 'C3h' or '-6':
        return int(22)
    elif sch == 'C6h' or '6/m':
        return int(23)
    elif sch == 'D6' or '622':
        return int(24)
    elif sch == 'C6v' or '6mm':
        return int(25)
    elif sch == 'D3h' or '-6m2':
        return int(26)
    elif sch == 'D6h' or '6/mmm':
        return int(27)
    elif sch == 'T' or '23':
        return int(28)
    elif sch == 'Th' or 'm-3':
        return int(29)
    elif sch == 'O' or '432':
        return int(30)
    elif sch == 'Td' or '-43m':
        return int(31)
    elif sch == 'Oh' or 'm-3m':
        return int(32)


def decompose_vec(basis_vec, target_vec):
    '''
    aim   : I will decompose the target target_vec into a linear combination of basis_vec
    input : basis_vec : ndarray[drow*dcol] 
                        and a column means a basis,the number of column(dcol) means the number of basis
            target_vec: one dimension ndarray
    output: one dimension ndarray with length of <dcol> of basis_vec
    formula : V = \sum_i a_i A_i, A_i = \sum_j b_{ji} B_j, 
              thus, V = \sum_{ij} a_i b_{ji} B_j = \sum_j v_j B_j
              here A_i means funciton basis, like seven f orbitals, while B_j is the polynomia basis like {x3,x2y,...}
              ten in total of f orbitals. Thus b is a retangular matrix not square matrix.
              problem: is I know v_j, how to get a_i in the case that {A_i} is not orthogonal and SO DOES {B_i}?
              b a = v  --->  b^T b a = b^T v  --->  b_new = (b^T b) is a square matrix which can be inversed 
              ---> a = (b^T b)^{-1} b^T v 

    '''
    drow,dcol = basis_vec.shape
    rep_new   = np.zeros(dcol,dtype=np.float64)
    basis_vec_tmp = np.dot(np.transpose(basis_vec),basis_vec)
    basis_vec_inv = np.linalg.inv(basis_vec_tmp)
    basis_util    = np.dot(basis_vec_inv,np.transpose(basis_vec))
    rep_new   = np.einsum('ij,j->i', basis_util, target_vec)
#   for i in range(dcol):
#       basis_tmp = basis_vec[:,i]
#       vec_inner = np.dot(basis_tmp,target_vec)
#       basis_mod = np.dot(basis_tmp,basis_tmp)
#       rep_new[i] = vec_inner / basis_mod
    # check the completeness of decomposition
    check_o = np.einsum('ij,j->i', basis_vec, rep_new)
    error = np.abs(target_vec - check_o)
    if error.sum() < 1.0E-6:
        print('>>> decomposition succeed ...')
        return rep_new
    else:
        print("error=\n",error)
        print("basis_vec \n",basis_vec)
        print("target_vec \n",target_vec)
        print("check_o \n",check_o)
        print("rep_new \n",rep_new)
        raise Exception("decomposition not successfull in PG_util.decompose_vec")


# permutation of new_struct to construct new polynomias
def pmutt(newstruct):
    '''
    permutation and combination of an array to a serials vector which will be also arranged as an array
    for example:
        ("x2yz")       (tran_func)                                      x      y      z      z      z
        [1,1,3]  -----------------------> newstruct=coeff_t*np.array([u[0,:],u[1,:],u[2,:],u[2,:],u[2,:]])
        if umat = [[0.1,0.2,0],[0,1,2],[0,0,1]],(here I drop off the coeff_t which doesn't change the result and
            also explictly the umat here is not orthogonal and normalized , sorry for that because it's just an example)
        then newstruct is : newstruct = [[0.1,0.2,0],[0,1,2],[0,0,1],[0,0,1],[0,0,1]]
        and the output result of pmutt is :
        [0.1,1,1,1,1] -> xyz3 -> [113] -> "113"
        [0.1,2,1,1,1] -> xz4  -> [104] -> "104"
        [0.2,1,1,1,1] -> y2z3 -> [023] -> "023"
        [0.2,2,1,1,1] -> yz4  -> [014] -> "014"
    '''
    drow, dcol = newstruct.shape
#
    pmutt_basis = []
    pmutt_coeff = {}
    cnt_xyz     = np.zeros(dcol,np.int32)
    if drow == 2 :
        for i1 in range(dcol):
            if abs(newstruct[0,i1]) > 1.0E-10 :
                for i2 in range(dcol):
                    if abs(newstruct[1,i2]) > 1.0E-10:
                        cnt_xyz[i1] += 1
                        cnt_xyz[i2] += 1
                        name_tmp = str(cnt_xyz[0]) + str(cnt_xyz[1]) + str(cnt_xyz[2])
                        value_tmp = newstruct[0,i1]*newstruct[1,i2]
                        if name_tmp in pmutt_basis :
                            pmutt_coeff[name_tmp] +=  value_tmp
                        else :
                            pmutt_basis.append(name_tmp)
                            pmutt_coeff[name_tmp] = value_tmp
                        cnt_xyz     = np.zeros(dcol,np.int32)
    elif drow == 3 :
        for i1 in range(dcol):
            if abs(newstruct[0,i1]) > 1.0E-10 :
                for i2 in range(dcol):
                    if abs(newstruct[1,i2]) > 1.0E-10:
                        for i3 in range(dcol):
                            if abs(newstruct[2,i3]) > 1.0E-10:
                                cnt_xyz[i1] += 1
                                cnt_xyz[i2] += 1
                                cnt_xyz[i3] += 1
                                name_tmp  = str(cnt_xyz[0]) + str(cnt_xyz[1]) + str(cnt_xyz[2])
                                value_tmp = newstruct[0,i1]*newstruct[1,i2]*newstruct[2,i3]
                                if name_tmp in pmutt_basis :
                                    pmutt_coeff[name_tmp] +=  value_tmp
                                else :
                                    pmutt_basis.append(name_tmp)
                                    pmutt_coeff[name_tmp] = value_tmp
                                cnt_xyz     = np.zeros(dcol,np.int32)
    # very rare 
    elif drow == 4 :
        for i1 in range(dcol):
            if abs(newstruct[0,i1]) > 1.0E-10 :
                for i2 in range(dcol):
                    if abs(newstruct[1,i2]) > 1.0E-10:
                        for i3 in range(dcol):
                            if abs(newstruct[2,i3]) > 1.0E-10:
                                for i4 in range(dcol):
                                    if abs(newstruct[3,i4]) > 1.0E-10:
                                        cnt_xyz[i1] += 1
                                        cnt_xyz[i2] += 1
                                        cnt_xyz[i3] += 1
                                        cnt_xyz[i4] += 1
                                        name_tmp = str(cnt_xyz[0]) + str(cnt_xyz[1]) + str(cnt_xyz[2])
                                        value_tmp = newstruct[0,i1]*newstruct[1,i2]*newstruct[2,i3]\
                                                *newstruct[3,i4]
                                        if name_tmp in pmutt_basis :
                                            pmutt_coeff[name_tmp] +=  value_tmp
                                        else :
                                            pmutt_basis.append(name_tmp)
                                            pmutt_coeff[name_tmp] = value_tmp
                                        cnt_xyz     = np.zeros(dcol,np.int32)
    return pmutt_basis, pmutt_coeff
