from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
import os
import sys

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

def isOrthogonal(basis,basis_set):
    '''
    to judge whether the function <basis> is orthogonal to the other basis set
    input : 
        basis : 1D numpy-ndarray
        basis_set : list of 1D numpy-ndarray
    '''
    isOrtho = True
    for ir in basis_set:
        if np.dot(basis,ir) > 1.0E-6:
            isOrtho = False
    return isOrtho

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
