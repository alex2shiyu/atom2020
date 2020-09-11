from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
from wanntb.tran import tran_op
import os
import sys
from pg.PG_util import pmutt, decompose_vec, RepBasisNorm, isOrthogonal, isindependent,find_rep,check_ham_wave
from pg.read_IR_DSG import loadIR 
import copy
from pg.gramSchmidt import GramSchmidt
from module.mod_dump import dump_4dc, dump_2dc, dump_1dr, dump_1dc
from src.atomic_subs import show_sub4header
import itertools

class Irrep:
    '''
    Irreps class 
    '''
    def __init__(self,label='',dim=1):
        self.label = label  # label for the representation. 'd' as postfix indicates double-valued irreps
        self.dim   = dim  # dimension of the representation
        self.matrices   = []  # representation matrices, ordered as the operations of the belonging group
        self.characters = []  # characters as a list, ordered the same as matrices
        self.rank  = 0

    def cal_multi(self,rep):
        '''
        calculate the multiplicity of the irrep in the given rep.
        rep:  should be a list of numpy.adarray which is the matrix representation of operators of PG 
        '''
        self.nrank = len(self.matrices) # the rank of the PG
        assert isinstance(rep,list)
        assert isinstance(rep[0],np.ndarray)
        character_t = complex(0.0)
        for iop in range(self.nrank):
            character_t += np.conjugate(self.characters[iop])*np.trace(rep[iop]) 
        multi_tmp = np.round(character_t.real/float(self.nrank))
        if np.abs(multi_tmp - character_t.real/float(self.nrank)) < 1.0E-6 :
            print("multi success!") if self.iprint == 3 else 0
            return int(multi_tmp)
        else:
            raise ValueError('multi is not integer : ',multi_tmp,character_t.real/float(self.nrank),np.abs(multi_tmp - character_t.real/float(self.nrank)))

class ReductRep(Irrep):
    def __init__(self,label='', dim=1, iprint=1):
        super().__init__(label='',dim=1)
        self.multi      = 0 # dictionary to contain various projectors for a certain Irrep
        self.projector  = {}  # dictionary to contain various projectors for a certain Irrep
#       self.reductstate= ''
        self.basis      = {}  # a dict of different set of basis for different <multi>, each set contain
                              # a list of basis whose number is the dimension of the irrep. For example
                              # self.multi = 2, self.dim = 3 , then the self.basis = dict{'multi0'=[np.array([.1D.the
                              # first basis.]),[..the second one.],[...the third one ...]],
                              # 'multi1'=[np.array([...1D..the first..]),np.array([...the second...]),np.array([...the
                              # third...])]} 
        self.basis_normortho   = {}  # it is orthogonal and normalized. Structure of this attribute is same with that of
                                     # self.basis
        self.iprint = iprint


    def irrep_normortho2(self):
        '''
        aim : we assume that the basis space of different type of irrep constructed with the help
              of projectors is already orthogonal to each other. It should be right and I make sure that 
              with the method MBPGsubs.check_projectors(). The results for D2h DPG in N=8 subspace of d orbitals
              shows that the value for inner product of basis of different irreps is less than 1.0E-4, there are
              about 10 inner products is ~ 1.0E-5 and most of them is less than 1.0E-9. I don't know why there are 
              still about 10 of them is as large as 1.0E-5. I guess because of the way we construct unitary matrix 
              in many-body space. However this is already a good approximation.
        note : should use this method after the reduction has been used
        difference : the difference with self.irrep_normortho2() is that I will orthogonalize the basis traversing 
                     every dimension instead of orthogonalize only one and then othogonalize others with the unitary 
                     matrix obtained from the one before.
        '''
        vectors_new = []
        for idim in range(self.dim):
            vectors_tmp = []
            for i in range(self.multi):
                name_t = 'multi' + str(i)
                vectors_tmp.append(self.basis[name_t][idim])
            vectors_n1 = GramSchmidt(vectors_tmp)
            vectors_new.append(vectors_n1)

        # normalized and orthogonalized
        vectors_normortho = []
        vectors_new_all = []
        for i in range(self.dim):
            print('basis before norm :\n',vectors_new[i]) if self.iprint ==3 else 0
            vectors_normortho_t = RepBasisNorm(list(vectors_new[i]))
            print('basis after norm :\n',vectors_normortho_t) if self.iprint ==3 else 0
            vectors_new_all += list(vectors_normortho_t)
            vectors_normortho.append(vectors_normortho_t)
       
        # check for the orthogonalize of all the produced basis of a certain irreps
        if self.iprint == 3:
            print(' ')
            print(10*'* * ')
            print(5*' ','>>> orthogonality[all] of the whole newly produced basis of the irrep ...')
        vectors_new_all_copy = copy.deepcopy(vectors_new_all)
        if self.iprint == 3:
            print('check before entering isOrthogonal :\n')
            print('                   vectors_new_all :',vectors_new_all)
            print('                   vectors_new_all_copy :',vectors_new_all_copy)
        is_allortho = isOrthogonal(vectors_new_all,vectors_new_all_copy,self.iprint) 
        if self.iprint == 3:
            print(5*' ','All is Orthogonal ?   ',is_allortho)

        
        # assign self.basis_normortho
        for i in range(self.multi):
            name_t = 'multi' + str(i)
            basis_or_t = []
            for j in range(self.dim):
                basis_or_t.append(vectors_normortho[j][i])
            self.basis_normortho[name_t] = copy.deepcopy(basis_or_t)
        if self.iprint == 3:
            print(' ')
            print(10*'* * ')
            print(5*' ','>>> orthogonality[multi] between different multi ...')
        for imul1 in range(self.multi):
            name1 = 'multi' + str(imul1)
            for imul2 in range(self.multi):
                name2 = 'multi' + str(imul2)
                if self.iprint == 3:
                    print(5*'--')
                is_ortho_multi = isOrthogonal(self.basis_normortho[name1],self.basis_normortho[name2],self.iprint)
                if self.iprint == 3:
                    print('Ortho: multi1 = ',imul1,' multi2 = ',imul2,'  -> ',is_ortho_multi)
                    print('')



    def irrep_normortho(self):
        '''
        aim : we assume that the basis space of different type of irrep constructed with the help
              of projectors is already orthogonal to each other. It should be right and I make sure that 
              with the method MBPGsubs.check_projectors(). The results for D2h DPG in N=8 subspace of d orbitals
              shows that the value for inner product of basis of different irreps is less than 1.0E-4, there are
              about 10 inner products is ~ 1.0E-5 and most of them is less than 1.0E-9. I don't know why there are 
              still about 10 of them is as large as 1.0E-5. I guess because of the way we construct unitary matrix 
              in many-body space. However this is already a good approximation.
        note : should use this method after the reduction has been used
        '''
        vectors_new = []
        vectors_org = []
        for idim in range(self.dim):
            vectors_tmp = []
            for i in range(self.multi):
                name_t = 'multi' + str(i)
                vectors_tmp.append(self.basis[name_t][idim])
            vectors_org.append(vectors_tmp)
        vectors_n1 = GramSchmidt(vectors_org[0])
        vectors_new.append(vectors_n1)
        # get transform matrix 
        v_o = np.array(vectors_org[0])
        v_o_d = np.dot(v_o, np.transpose(np.conjugate(v_o)))
        Schmidt_umat = np.dot(np.dot(vectors_new[0],np.transpose(np.conjugate(v_o))),np.linalg.inv(v_o_d))
        # check the self-consistency of Schmidt_umat 
        check_umat = np.sum(np.abs(np.dot(Schmidt_umat,np.array(vectors_org[0]))-vectors_new[0]))
        check_max  = np.max(np.abs(np.dot(Schmidt_umat,np.array(vectors_org[0]))-vectors_new[0]))
        if check_umat < 1.0E-8:
            print(5*" ","check the self-consistency of Schmidt_umat of ",self.label," : OK")
        else:
            print(5*" ","check the self-consistency of Schmidt_umat of ",self.label," : Fail --> ",check_umat)
            print(5*" ","                                              ",self.label," : max  --> ",check_max)
        if self.dim > 0 :
            print(5*' ',10*'- . ')
            for idim in range(1,self.dim):
                vectors_n_t = np.dot(Schmidt_umat,np.array(vectors_org[idim]))
                vectors_new.append(vectors_n_t)
                isOrtho_no = isOrthogonal(list(vectors_n_t),list(vectors_n_t),self.iprint)
                print(" * check for orthogonality of new produced basis")
                print(5*' ',10*'- .. ')
                print(5*' ',' *[New Ortho] dim =',idim)
                print(5*' ','                   ',isOrtho_no)
                print('')
                # check orthogonal


        # which has been normalized and orthogonalized
        vectors_normortho = []
        for i in range(self.dim):
            vectors_normortho_t = RepBasisNorm(list(vectors_new[i]))
            vectors_normortho.append(vectors_normortho_t)
        
        # assign self.basis_normortho
        for i in range(self.multi):
            name_t = 'multi' + str(i)
            basis_or_t = []
            for j in range(self.dim):
                basis_or_t.append(vectors_normortho[j][i])
            self.basis_normortho[name_t] = copy.deepcopy(basis_or_t)


    # multi should be assigned with the help of method of super class 
    def reduction(self,op,irrep_prev):
        '''
        will reduce the reducible representation op
        input : 
            op [list]  : is a list of numpy.ndarray which is the representation matrix of every operators of the PG
            irrep_prev : other irrep whose multiplicity is greater than 0 and has been reduced successfully. 
            <self.projector and self.multi refers to information about "this irrep" and "this input op">
        '''
        for imul in range(self.multi):
            name_key = 'multi'+str(imul)
            if self.iprint == 3:
                print('')
                print('')
                print(20*'=')
                print('* * imul :',imul)
            phi_t1 = self.make_phi1(imul,op[0].shape[0],irrep_prev)
            if self.iprint == 3 :
                print(20*'- -') 
                print('* WAVE(phi1) Done :  of multi<',imul,'> of <',self.label,'> is :',np.linalg.norm(phi_t1))
            basis_list = []
            basis_list += phi_t1 # phi_t1 is list of numpy.ndarray
            if self.dim > 1 :
                phi_other = self.make_phi_other1(phi_t1[0],imul)
                basis_list += phi_other
#           if self.iprint == 2 :
            if self.iprint >= 2 :
                print('* basis set for multi',imul,'of ',self.label,'is\n',basis_list)
            is_multi_ortho = isOrthogonal(basis_list,basis_list,self.iprint)
#           if self.iprint == 2 :
            if self.iprint >= 2 :
                print('* is this set basis orthogonal to each other :',is_multi_ortho)
                print(20*'- -') 
            self.basis[name_key] = basis_list
                



    def make_phi1(self,multi_now,ndim,irrep_prev):
        '''
        according to intro. in Alttman's point group tables,
        W^i_11 \phi = \phi^i_1 
        however, if self.multi > 1, when make_phi1 for multi_now > 1, you should make sure the 
        \phi^i_1 is orthogonal to the previous set of bases
        input : 
                multi_now : Integer <counts from 0>
                ndim :      dimension of the subspace which is also the dimension of op matrix
                            I will traverse the basis of the space like 
                            [1,0,0,0, ...]
                            [0,1,0,0, ...]
                            [0,0,1,0, ...]
                            until the useful \phi^i_1 has been found
                irrep_prev: other irrep whose multiplicity is greater than 0 and has been reduced successfully. 
                            this input is used to make the new produced \phi^i_1 is orthogonal to the basis set of 
                            other irreps of this PG
        output: func_1 [list of a <1D-numpy.ndarray>] 
        '''
        label_list = []
        name_ip = 'P00'
        isindependt_all = [False]
        for i in range(ndim):
            if self.iprint == 3 :
                print("\n",10*'*','\n',2*'* ','i :',i)
            func_all    = np.zeros(ndim,dtype=np.complex128)
            func_all[i] = complex(1.0)
#           for ip in range(1):
            if self.iprint == 3 :
                print("\n",'. . . . . . . . . ','\n','i(P_ii)=',0,'|','\n','. . . . . . . . . ')
            Proj    = self.projector[name_ip]
            func_1  = np.einsum('ji,i->j',Proj,func_all)
            # to make sure phi2 = P10 phi1 has nonzero element
            if self.dim > 1 :
                phi2 = self.make_phi_other1(func_1,multi_now)
#               if np.sum(np.abs(phi2)) > 1e-2 :
                max_flag = 100.0
                for i in range(len(phi2)):
                    max_current = np.abs(np.array(phi2[i])).max()
                    if max_current < max_flag:
                        max_flag = copy.deepcopy(max_current) 
                if max_flag > 1e-2 :
                    phi2_flag = True
                else:
                    phi2_flag = False
            else :
                phi2_flag = True 
            if np.sum(np.abs(func_1)) > 1.0e-4 and phi2_flag:
                func_1  = RepBasisNorm(func_1) 
                basis_set = []
                if self.iprint == 3:
                    print(5*'- ','\n','   check independence ...')
                if int(multi_now) > 0 :
                    for imul in range(multi_now):
                        mul_name  = 'multi' + str(imul)
                        basis_set += self.basis[mul_name]
                if len(irrep_prev) > 0:
                    for irr in irrep_prev :
                        if irr.multi > 0 :
                            for imul in range(irr.multi):
                                mul_name = 'multi' + str(imul)
                                basis_set += irr.basis[mul_name]
                if basis_set == [] :
                    if self.iprint == 3:
                        print(5*'  ','NO NEED : because imul=',multi_now,'len(irrep_prev)=',len(irrep_prev))
                    label_list.append(np.sum(np.abs(func_1)))
                    isindependt_all += [True]
                else :
                    label_t, independt   = isindependent(func_1,basis_set,self.iprint) 
                    if self.iprint == 3:
                        print('is independent? : ',independt)
                    isindependt_all += independt
                    if not independt :
                        label_list.append(float(0.0))
                    else :
                        label_list.append(np.sum(np.abs(np.array(label_t))))
            else:
                isindependt_all += [False]
                label_list.append(float(0.0))
        
        if any(isindependt_all) :
            pos_target = label_list.index(max(label_list))
            if self.iprint >= 2 :
                print('    ---> the proper phi for phi1 of mulit=',multi_now,'of irreps ',self.label,'is',pos_target)
            func_all    = np.zeros(ndim,dtype=np.complex128)
            func_all[pos_target] = complex(1.0)
            Proj    = self.projector[name_ip]
            func_1  = np.einsum('ji,i->j',Proj,func_all)
            func_1  = RepBasisNorm(func_1) 
            return [func_1]
        else:
            print('ERROR-info1:isindependt:\n',isindependt_all)
            print('ERROR-info1:label_list:\n',label_list)
            raise ValueError('ERROR in make_phi1')
#####               print('   ??? is independent to pervious set: \n','   ',isindependt_all)
#####               print('\n')
#####               if all(isindependt_all) == True :
#####                   print(10*'*-')
#####                   print('* the proper phi for phi1 of mulit=',multi_now,'of irreps ',self.label,'is',i+1)
#####                   print(10*'*-')
#####                   return [func_1]
#####               elif i == ndim-1  and ip == self.dim-1 :
#####                   raise ValueError("Can not find proper phi_1 orthogonal to previous basis \
#####                           set:",multi_now,isindependt_all)
#####           elif i == ndim-1 and ip == self.dim-1 :
#####               raise ValueError("Can not find phi_1 with finite elements  in <make_phi1>")
#####           else :
#####               print('Warning : norm of phi_1 is too small ...')
        
    def make_phi_other(self,phi_1,multi_now):
        '''
        aim : once there is phi_1, produce the other basis using projectors W^i_10, W^i_21, ... when the dimension
              of the irrepresentation is larger than 1
        input : 
                phi_1 : 1D-numpy.ndarray
                multi_now : Integer <counts from 0>
        output
                phi_other : list of     
        '''
        basis_other = []
        for i in range(self.dim-1):
            name_proj = 'P' + str(i+1) + str(i)
#           name_proj = 'P' + str(i) + str(i+1)
            proj      = self.projector[name_proj]
#           basis_t   = np.dot(proj,phi_1)
            basis_t   = np.einsum('ji,i->j',proj,phi_1)
            if self.iprint >= 2 :
                print('* WAVE(phi'+str(i+2)+') Done :  of multi<',multi_now,'> of <',self.label,'> is :',np.linalg.norm(basis_t))
            basis_other.append(basis_t)
        return basis_other

    def make_phi_other1(self,phi_1,multi_now):
        '''
        version : advanced version for dimension of irreps > 2
        aim : once there is phi_1, produce the other basis using projectors W^i_10, W^i_21, ... when the dimension
              of the irrepresentation is larger than 1
        input : 
                phi_1 : 1D-numpy.ndarray
                multi_now : Integer <counts from 0>
        output
                phi_other : list of     
        '''
        basis_other = []
        basis_other.append(phi_1)
        for i in range(self.dim-1):
            name_proj = 'P' + str(i+1) + str(i)
#           name_proj = 'P' + str(i) + str(i+1)
            proj      = self.projector[name_proj]
#           basis_t   = np.dot(proj,phi_1)
            basis_t   = np.einsum('ji,i->j',proj,basis_other[len(basis_other)-1])
            if self.iprint >= 2 :
                print('* WAVE(phi'+str(i+2)+') Done :  of multi<',multi_now,'> of <',self.label,'> is :',np.linalg.norm(basis_t))
            basis_other.append(basis_t)
        basis_other.pop(0) # delete the first one which is phi_1
        return basis_other
        

    def make_projector(self,op):
        '''
        aim : make projectors of every irreps for subspace of N operators 
        op  : [list of numpy.ndarray] the matrix rep. of operators of PG
        '''
        assert isinstance(op,list)
        assert isinstance(op[0],np.ndarray)
        dim_op =op[0].shape[0] 
# firstly make character projector
        P_character = np.zeros((dim_op,dim_op),dtype=np.complex128)
        for i in range(self.nrank):
            P_character += self.dim/self.nrank * np.conjugate(self.characters[i])*op[i]
        self.projector['Pcharacter'] = P_character
# secondly make quasi projector
        P_quasi = np.zeros((self.dim,dim_op,dim_op),dtype=np.complex128)
        for i in range(self.dim):
            for j in range(self.nrank):
                P_quasi[i,:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][i,i])*op[j]
        self.projector['Pquasi'] = P_quasi
# thirdly make normal projector
        # first W_1p (p = 1, ..., self.dim)
        for i in range(self.dim):
            P_normal = np.zeros((dim_op,dim_op),dtype=np.complex128)
            name = 'P' + str(0) + str(i)
            for j in range(self.nrank):
#               P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][0,i])*op[j]
                P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][0,i])*op[j]
            self.projector[name] = P_normal
            if i > 0:
                P_normal = np.zeros((dim_op,dim_op),dtype=np.complex128)
                name1 = 'P' + str(i) + str(i)
                for j in range(self.nrank):
                    P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][i,i])*op[j]
                self.projector[name1] = P_normal

        # then W10, W21, ...
        if self.dim > 1 :
            for i in range(self.dim-1):
                P_normal = np.zeros((dim_op,dim_op),dtype=np.complex128)
                name = 'P' + str(i+1) + str(i)
                if self.iprint == 3 :
                    print(5*'** ','check for projectors matrix value : ')
                for j in range(self.nrank):
#                   print(10*' ')
#                   print('Projectors : j=',j)
#                   print(5*' ','characters : \n',self.matrices[j])
#                   print(5*' ','op matrix  : \n',op[j])
#                   print('')
#                   P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][i+1,i])*op[j]
                    P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][i+1,i])*op[j]
                self.projector[name] = P_normal
#       print('Pro name   :',name)
        if self.iprint == 3 :
            print('>>> Projectors  :',self.projector)




# double Point group
class DPG():
    '''
    class for double point group

    attributes:
        pgid[I]          : serial number of point group
   [NY] name_hm[str]     : Hermann-Mauguin symbol for point group
        name_sc[str]     : schonflies notation for point group
        nrank[I]         : rank of the double point group
        rep_vec[list]    : vec rep in (x,y,z) space
        rep_spin[list]   : vec rep in (up,dn) space
        mb_shell[str]    : shell of atom -> s p d f
        mb_norb[I]       : number of orbitals of the corresponding shell
        mb_rep[ndarray]  : the operator matrix in joint space of s or p or d or f
        nclass[I]        : number of class
        irrep[list]      : irreducible representation of double point group
        gen_pos[ndarray] : position of generator in rep  
        gen_mul[ndarray] : number of elements of irrep in rep  

    comment: [NY] not yet implemented
    '''
    def __init__(self, name_sc='', pgid=0):
        self.name_sc   = name_sc
        self.pgid      = pgid
        self.nrank     = None
        self.rep_vec   = []    # the rep_vec has been transposed compared with data from Bilbao in the method get_data
        self.rep_spin  = None
        self.mb_shell  = None
        self.mb_norb   = None
        self.mb_rep    = None
        self.nclass    = None
        self.irreps    = None
        self.gen_pos   = None
        self.gen_mul   = None
        self.irreps_class = []
  
    def show_attribute(self):
        print(10*'-')
        print('name_sc:',self.name_sc)
        print('pgid:',self.pgid)
        print('nrank:',self.nrank)
        print('nclass:',self.nclass)
        print('mb_norb:',self.mb_norb)
        print('rep_vec:\n',self.rep_vec)
        print('rep_spin:\n',self.rep_spin)
        print('mb_shell:\n',self.mb_shell)
        print('irreps:\n',self.irreps)
        print('gen_pos:\n',self.gen_pos)
        print('gen_mul:\n',self.gen_mul)
        print('irreps_class:\n',self.irreps_class)
    # should first assign self.irreps
    def groupclass_irrep(self):
        self.gen_pos = np.array([0,1,2,3,4,5,6,7,8,12],dtype=np.int32)
        self.gen_mul = np.array([1,2,2,2,1,2,2,2,1,1],dtype=np.int32)
        self.nclass  = len(self.gen_pos)
        for ir in self.irreps:
            irr = Irrep()
            irr.label = ir.label
            irr.dim   = ir.dim
#           character_reduce = np.array([ir.characters[self.gen_pos[i]].real for i in range(self.nclass)])
            character_reduce = []
            matrix_reduce = []
            for igen in range(self.nclass):
                character_reduce.append(ir.characters[self.gen_pos[igen]])
                matrix_reduce.append(ir.matrices[self.gen_pos[igen]])
            irr.matrices = matrix_reduce
            irr.characters =np.array(character_reduce)
            self.irreps_class.append(irr)

        
        

    def def_mb_byhand(self,case='f'):
        if case == 'f' :
            num_oprt = 10 # just for D2h double point group
            a = np.zeros((7,7),dtype = np.float64)# the orbitals number is same with that of wannier90

        elif case == 'd' :
            a == np.zeros((5,5),dtype = np.float64)
        elif case == 't2g' :
            a == np.zeros((3,3),dtype = np.float64)
        elif case == 'p' :
            a == np.zeros((3,3),dtype = np.float64)
        else :
            raise ValueError('PG can\'t recognize the value: ',case)

    def get_data(self,sgid):
        for gid in range(sgid,sgid+1):
            lgrps = loadIR(gid)
        for grp in lgrps:
            if grp.klabel == 'GM' :
                self.nrank    = len(grp.rotC)
                self.rep_vec  =[np.transpose(np.conjugate(i)) for i in grp.rotcar] # I find the bilbao's data is also transfor in
                                                         # rows, for example dsp#139, what ever grp.rotcar should be 
                                                         # made sure transform in rows. Then why I need to transform it
                                                         # using transpose, becasue we need R^(-1) not R, that is mean
                                                         # that P_R f[r] = f[R^(-1) r]
                                                         #----------------below is written down before
                                                         # bilbao's data is transfromed in columns and we transpose it
                                                         # here but you should know when we use the manybody version of
                                                         # it to construct the projectors, we will use the form which
                                                         # transform in columns as same as the convention of group
                                                         # theory
#               self.rep_spin = [a[::-1,::-1] for a in grp.su2s]
#               self.rep_spin = [np.transpose(np.conjugate(a)) for a in grp.su2s]
                self.rep_spin = grp.su2s # it transform in rows for vectors
                self.irreps   = grp.irreps
#               print('dpg71.nrank:\n',self.nrank)
#               print('rotC:\n',grp.rotC)
#               print('dpg71.rep_vec:\n',self.rep_vec)
#               print('dpg71.rep_spin:\n',self.rep_spin)
#               for ir in self.irreps:
#                   print(ir.label,ir.characters.real)

# many body point group
class MBPG():
    '''
    after get all the prerequisite, a class of essential information of point group need by atomic's problem is desired
    '''
    def __init__(self, nop, op,iprint):
        '''
        attributes:
            nop [I]       : rank of the point group
            op  [list]    : a list of numpy.ndarray which is the matrix rep. of operators in natural basis, transfrom in
                            # columns
            Nsubs [list]   : to contain the instance of MBPGsubs of different subspace of \hat{N}
            basis [list]  : to contain basis 
#           irrep [list]  : a list of Irrep 
        '''
        self.nop      = nop
        self.op       = op
        self.iprint   = iprint
        self.ham      = None 
        self.ham_eig  = None
        self.ham_evc  = None
        self.ham_natural = None # which has been transfromed into natural basis
        self.Nsubs    = []
        self.dim      = 100  
        self.op_irrep = []
        self.vpmnum   = 0
        self.vpmsy    = None

    def show_results(self):
        print('\n\n')
        print(' *'+75*'='+'*')
        print(' |'+29*' ','SymmAtom Results',28*' '+'|')
        print(' *'+75*'='+'*')
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  DPG',' ',' ','|')) 
        print(' |'+75*'-'+'|')
        print('{:>2}{}{:<14}{:>35}{:>23}{:>2}'.format('|',2*' ','TotNumCfgs',':',self.dim,'|')) 
        print('{:>2}{}{:<14}{:>35}{:>23}{:>2}'.format('|',2*' ','TotNumVpms',':',self.vpmnum,'|')) 
        for nsubs in self.Nsubs :
            print('{:>2}{}{:<20}{:>29}{:>23}{:>2}'.format('|','  ',11*'-',' ',' ','|')) 
            print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  N='+str(nsubs.nocc),' ',' ','|')) 
            if nsubs.vpmtype[0] == 'g' :
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','NumCfgs',':',nsubs.dim,'|')) 
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','NumVpms',':',nsubs.vpmnum,'|')) 
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','-------',' ',' ','|')) 
                for irr in nsubs.irrep :
                    if irr.multi > 0 :
                        print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ',irr.label,':',irr.multi,'|')) 
            elif nsubs.vpmtype[0] == 'd' :
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','NumCfgs',':',nsubs.dim,'|')) 
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','NumVpms',':',nsubs.vpmnum,'|')) 
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ','-------',' ',' ','|')) 
                print('{:>2}{}{:<14}{:>54}{:>2}'.format('|',6*' ','diagonal vpms',' ','|')) 
#           if nsubs != self.Nsubs[len(self.Nsubs)-1]:
        # print number of total configurations and vpms

        print(' *'+75*'-'+'*')
# hamiltonian
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  Hamiltonain',' ',' ','|')) 
        print(' |'+75*'-'+'|')
        for nsubs in self.Nsubs :
            print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  N='+str(nsubs.nocc),' ',' ','|')) 
            if nsubs.vpmtype[0] == 'g' :
                list_t = [nsubs.deginfo[ii]['degeneracy'] for ii in range(nsubs.deg)]
                list_tt= [nsubs.irrep_accidental_deg[ii] for ii in range(nsubs.deg)]
                dim_t = len(list_t)
                remainder_t = dim_t % 10
                if remainder_t == 0:
                    times_t = int(dim_t/10)
                else :
                    times_t = int(dim_t /10) + 1
                if remainder_t > 0:
                    for i in range(10 - remainder_t):
                        list_t.append(' ')
                        list_tt.append(' ')
                ibase = int(0)
                for iline in range(times_t):
                    if iline == 0 :
                        print('{:>2}{}{:<13}{:<4}'.format('|',6*' ','Deg',':'),end=' ')
                    else :
                        print(' |'+23*' ',end=' ')
                    for icnt in range(10):
                        print('{0:<{width}}'.format(list_t[icnt + ibase], width=4), end=' ')
                    print('{:>2}'.format('|'))
#
                    if iline == 0 :
                        print('{:>2}{}{:<13}{:<4}'.format('|',6*' ','Accidental?',':'),end=' ')
                    else :
                        print(' |'+23*' ',end=' ')
                    for icnt in range(10):
                        print('{0:<{width}}'.format(list_tt[icnt + ibase], width=4), end=' ')
                    print('{:>2}'.format('|'))
                    ibase = (iline+1)*10
                print('{:>2}{}{:<14}{:>31}{:>23}{:>2}'.format('|',6*' ',14*'-',' ',' ','|')) 
                for keys in nsubs.ham_evc_irrepindex['group'] :
#                   print('{:>2}{}{:<8}{:<5}{:>36}{:>23}{:>2}'.format('|',6*' ',keys,':',nsubs.ham_evc_irrepindex['group'][keys],'|')) 
#                   print('{:>2}{}{:<8}{:<4}'.format('|',6*' ',keys,':'),\
#                           *(map('{0[1]}'.format, enumerate(nsubs.ham_evc_irrepindex['group'][keys]))),'{:>2}'.format('|')) 
                    list_t = nsubs.ham_evc_irrepindex['group'][keys]
                    dim_t = len(list_t)
                    remainder_t = dim_t % 11
                    if remainder_t == 0:
                        times_t = int(dim_t/11)
                    else :
                        times_t = int(dim_t /11) + 1
                    if remainder_t > 0:
                        for i in range(11 - remainder_t):
                            list_t.append(' ')
                    ibase = int(0)
                    for iline in range(times_t):
                        if iline == 0 :
                            print('{:>2}{}{:<8}{:<4}'.format('|',6*' ',keys,':'),end=' ')
                        else :
                            print(' |'+18*' ',end=' ')
                        for icnt in range(11):
                            print('{0:<{width}}'.format(list_t[icnt + ibase], width=4), end=' ')
                        print('{:>2}'.format('|'))
                        ibase = (iline+1)*11
#                   print('{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}'.format(*(nsubs.ham_evc_irrepindex['group'][keys][0:10])),'{:>2}'.format('|')) 
            elif nsubs.vpmtype[0] == 'd' :
                print('{:>2}{}{:<14}{:>54}{:>2}'.format('|',6*' ','diagonal vpms',' ','|')) 
            print('{:>2}{}{:<68}{:>2}'.format('|',6*' ',68*'-','|')) 
        print(' *'+75*'='+'*')


    def show_attribute(self):
        print('nop:',self.nop)
        print('operators:\n',self.op)
#       print('irreps:\n',self.irrep)

    def collect_eig(self,nstat):
        '''
        aim : stack bases from self.Nsubs to an array
        input : nstat ：length of every block(to check)
        '''
        
        if len(nstat) == len(self.Nsubs):
            ndim = sum(nstat)

            self.ham_eig = np.zeros(ndim,dtype=np.float64)
            ibase = int(0)
            for i in range(len(nstat)):
                self.ham_eig[ibase:ibase+int(nstat[i])] = self.Nsubs[i].ham_eig
                ibase += int(nstat[i])
        else :
            raise ValueError('Error in collect_eig of MBPG :','len(nstat)=',len(nstat),'len(self.Nsubs)=',len(self.Nsubs))
    
    def collect_evc(self,nstat):
        '''
        aim : stack bases from self.Nsubs to an array
        input : nstat ：length of every block(to check)
        '''
        
        if len(nstat) == len(self.Nsubs):
            ndim = sum(nstat)

            self.ham_evc = np.zeros((ndim,ndim),dtype=np.complex128)
            ibase = int(0)
            for i in range(len(nstat)):
                self.ham_evc[ibase:ibase+int(nstat[i]),ibase:ibase+int(nstat[i])] = self.Nsubs[i].ham_evc_natural
                ibase += int(nstat[i])
        else :
            raise ValueError('Error in collect_evc of MBPG :','len(nstat)=',len(nstat),'len(self.Nsubs)=',len(self.Nsubs))
    
    def collect_ham(self,nstat):
        '''
        aim : stack bases from self.Nsubs to an array
        input : nstat ：length of every block(to check)
        '''
        
        if len(nstat) == len(self.Nsubs):
            ndim = sum(nstat)

            self.ham_natural = np.zeros((ndim,ndim),dtype=np.complex128)
            ibase = int(0)
            for i in range(len(nstat)):
                self.ham_natural[ibase:ibase+int(nstat[i]),ibase:ibase+int(nstat[i])] = self.Nsubs[i].ham_natural
                ibase += int(nstat[i])
        else :
            raise ValueError('Error in collect_ham of MBPG :','len(nstat)=',len(nstat),'len(self.Nsubs)=',len(self.Nsubs))
    
    def collect_vpm(self,nstat):
        '''
        aim : stack bases from self.Nsubs to an array
        input : nstat ：length of every block(to check)
        '''
        self.vpmnum = int(0)
        if len(nstat) == len(self.Nsubs):
            ndim = sum(nstat)
            if ndim != self.dim :
                raise ValueError('Error in collect_vpm : dim =',ndim,'ncfgs =',self.dim)
            self.vpmsy = np.zeros((ndim,ndim),dtype=np.int64)
            ibase = int(0)
            threshhold = int(0)
            for i in range(len(nstat)):
                vpmsy_copy = copy.deepcopy(self.Nsubs[i].vpmsy)
                mask_array = vpmsy_copy > 0
                vpmsy_copy[mask_array] += threshhold
                self.vpmsy[ibase:ibase+int(nstat[i]),ibase:ibase+int(nstat[i])] = vpmsy_copy
                ibase += int(nstat[i])
                threshhold = np.max(self.vpmsy)
                self.vpmnum += self.Nsubs[i].vpmnum
        else :
            raise ValueError('Error in collect_vpm of MBPG :','len(nstat)=',len(nstat),'len(self.Nsubs)=',len(self.Nsubs))
class MBPGsubs(MBPG): 
    '''
    every subspace of N will have their own information when do reduciton
    '''
    def __init__(self, nop, op, nocc, iprint, vpmtype):
        super().__init__(nop, op, iprint)
        self.nocc = nocc
        self.irrep = []
        self.allbasis = {}
        self.allbasis['matrix'] = []# due to it's a list, thus the eigenwaves arranged in rows
        self.allbasis['irreplabel'] = []
        self.ham_irrep = np.zeros(op[0].shape,dtype=np.complex128) if len(op) > 0 else None
        self.ham_evc_natural = None # which has been transfromed into natural basis
        self.ham_evc_irrepindex = {}
        self.irrepindex = {}
        self.ham_evc_decompose = []
        self.Focknatural = None
        self.vpmtype = vpmtype  # str
        self.deg     = 100
        self.deginfo = []
        self.irrep_accidental_deg = [] # record whether there is accidential degeneracy in every degenerate space
 
    def group_irrep_evc(self):
        '''
        aim : group evc with same irrep as well as same columns together which is convenient to make atom.vpmsy.in

        note : here the keys is named with the index starting from 1 not 0 just for unambiguous output
               GM5d_1 GM5d_2 ...
        '''
        self.ham_evc_irrepindex['group'] = {}
        for i in range(len(self.ham_evc_irrepindex['sequential'])):
            keys = self.ham_evc_irrepindex['sequential'][i][0] + '_' + str(self.ham_evc_irrepindex['sequential'][i][2]+1)
            if keys in self.ham_evc_irrepindex['group'] :
                self.ham_evc_irrepindex['group'][keys].append(i+1) 
            else :
                self.ham_evc_irrepindex['group'][keys] = [i+1]
 
    def make_vpmsy(self):
        self.vpmsy = np.zeros((self.dim,self.dim),dtype=np.int32)
        if self.vpmtype[0] == 'g':
            num_vpm = int(0)
            for idict in dict.items(self.ham_evc_irrepindex['group']):
                for idim in idict[1]:
                    for jdim in idict[1]:
                        num_vpm += 1
                        self.vpmsy[idim-1,jdim-1] =  num_vpm
        elif self.vpmtype[0] == 'd' :
            num_vpm = int(0)
            for idim in range(self.dim):
                num_vpm += 1
                self.vpmsy[idim,idim] = num_vpm
        self.vpmnum = num_vpm


    def getirrepindex(self):
        '''
        aim   : to get the index in self.irrep
        usage : self.irrepindex['GM1+'] = 0
        '''
        cnt_t = int(-1)
        for ir in self.irrep :
            cnt_t += 1
            self.irrepindex[ir.label] = int(cnt_t)


    # check the invariant of hamiltonian 
    def check_symm_ham(self):
        ham_new = np.zeros((self.dim,self.dim),dtype=np.complex128)
        for iop in range(self.nop):
            ham_new += np.dot(np.conjugate(np.transpose(self.op[iop])),np.dot(self.ham,self.op[iop]))
        ham_new = ham_new / float(self.nop)

        error = np.sum(np.abs(self.ham - ham_new))
        if error > 1.0e-6 :
            dump_2dc(self.dim,self.dim,ham_new,path='ham_new.dat',prec=1.0e-6)
            dump_2dc(self.dim,self.dim,self.ham,path='ham.dat',prec=1.0e-6)
            dump_2dc(self.dim,self.dim,self.ham - ham_new,path='ham_diff.dat',prec=1.0e-6)
            raise ValueError('error in check_symm_ham and error value : ',error)
        elif self.iprint >= 2 :
            print('    ---> success : Check for the symmetry of hamiltonian :',error)

    # transfrom operators into basis of irreps after collecting bases
#   def trans_operators(self,trans):
#       if trans :
#           print('Check : the transformation matrix(basis matrix) in rows :\n',self.allbasis['matrix'])
#           for iop in range(self.nop):
#               op = tran_op(self.ham,np.transpose(np.array(self.allbasis['matrix'])))
#       else :
#           self.ham_irrep = copy.deepcopy(self.ham)



    # transfrom ham into basis of irreps after collecting bases
    def trans_ham(self,trans=True):
        if trans :
            if self.iprint == 3 :
                print('Check : the transformation matrix(basis matrix) in rows :\n',self.allbasis['matrix'])
            self.ham_irrep = tran_op(self.ham,np.transpose(np.array(self.allbasis['matrix'])))
            print(4*' '+'---> transformation done.')
        else :
            self.ham_irrep = copy.deepcopy(self.ham)
            print(4*' '+'---> not transformation.')



    # diagonalize ham after transforming the hamiltonian into irreps bases
    def diag_ham2(self,trans=True):
        self.ham_eig, self.ham_evc = np.linalg.eigh(self.ham_irrep)
        print(20*'** ')
        print('')
        check_ham_wave(self.ham_irrep,self.ham_eig,self.ham_evc,self.iprint)
        print('Check for the symmetry-kept eigenfunctions :')
        print('[index from 0]')
        Bmat_col = np.transpose(np.array(self.allbasis['matrix']))
        print('op : \n',self.op)
        if trans :
            ham_evc_original = np.dot(Bmat_col,self.ham_evc)
            print('Check ham :\n',np.dot(np.conjugate(np.transpose(ham_evc_original)),np.dot(self.ham,ham_evc_original)))
        for idim in range(self.dim):
            for iop in range(self.nop):
                if trans :
# transform functions basis in columns
                    op = np.dot(np.dot(np.transpose(np.conjugate(Bmat_col)),self.op[iop]),Bmat_col)
                    print('operators in irreps basis :\n',iop)
                    print(op)
                    op1 = copy.deepcopy(self.op[iop])
#                   op = np.dot(np.dot(np.transpose(Bmat),self.op[iop]),np.conjugate(Bmat))
                else :
                    op = copy.deepcopy(self.op[iop])
#               if trans :
#                   error_ham_irr = np.sum(np.abs(np.dot(self.ham_irrep,op)-np.dot(op,self.ham_irrep)))
#                   wave_new = np.einsum('ij,j->i',op,ham_evc[:,idim])
                error_ham_irr = np.sum(np.abs(np.dot(self.ham_irrep,op)-np.dot(op,self.ham_irrep)))
                wave_new = np.einsum('ij,j->i',op,self.ham_evc[:,idim])
#               else :
#                   error_ham_irr = np.sum(np.abs(np.dot(self.ham_irrep,op)-np.dot(op,self.ham_irrep)))
#                   wave_new = np.einsum('ij,j->i',op,self.ham_evc[:,idim])
                for jdim in range(self.dim):
#                   if trans :
#                       inner_product = np.dot(np.conjugate(self.ham_evc[:,jdim]),wave_new)
                    inner_product = np.dot(np.conjugate(self.ham_evc[:,jdim]),wave_new)
#                   else :
#                       inner_product = np.dot(np.conjugate(self.ham_evc[:,jdim]),wave_new)
                    if np.abs(inner_product) > 1.0e-6  :
                        print('WARNING:  m=',jdim,'O_i',iop,'n=',idim,'<m|O_i|n> =',inner_product)
                        print(5*' ','error_ham_irr =',error_ham_irr)
                        print(5*' ')
        print('')
        print(20*'** ')

    # diagonalize ham after transforming the hamiltonian into irreps bases
    def diag_ham(self,trans=True):
        self.ham_eig, self.ham_evc = np.linalg.eigh(self.ham_irrep)

    # to collect all the basis after self.Cal_ReductRep()
    def collect_basis(self):
        for irr in self.irrep:
            if irr.multi > 0 :
                for imul in range(irr.multi):
                    name = 'multi' + str(imul)
                    for idim in range(irr.dim):
                        self.allbasis['matrix'].append(irr.basis[name][idim])
                        self.allbasis['irreplabel'].append([irr.label,irr.dim,idim,irr.multi])

    def Cal_ReductRep(self,irrep_input):
        '''
        traverse every possible irreps of the PG to calculate multiplicities and corresponding projectors
        '''
        true_pre_irrep = [] # to contain the irreps which has self.multi greater than one for sake of reduction for the
                            # following irreps
        for ir in irrep_input:
            irr = ReductRep(iprint=self.iprint)
            irr.label = copy.deepcopy(ir.label)
            irr.dim   = copy.deepcopy(ir.dim)
            irr.matrices   = copy.deepcopy(ir.matrices) #[np.transpose(imat) for imat in ir.matrices]
            irr.characters = copy.deepcopy(ir.characters)
            print('>>>>> irrep :',irr.label) if self.iprint >= 2 else 0
            irr.multi      = irr.cal_multi(self.op)
            print('      multi :',irr.multi) if self.iprint >= 2 else 0
            irr.make_projector(self.op)
            print('      projectors :\n',irr.projector.keys()) if self.iprint >= 2 else 0
            if irr.multi > 0 :
#               irr.reduction(self.op,self.irrep)
                irr.reduction(self.op,true_pre_irrep)
                irr.irrep_normortho2()
                true_pre_irrep.append(irr)
#           irr.reductstate= 'ok'
            self.irrep.append(irr)
        #
        self.getirrepindex()
        #
        # test for the orthogonal of different irreps
        show_sub4header('orthogonality[irreps] of all bases')
        isOrtho_tmp = [True]
        for i in range(1,len(self.irrep)):
            if self.irrep[i].multi == 0 :
                continue
            else:
                for imul in range(self.irrep[i].multi):
                    name_i = 'multi' + str(imul)
                    basis_i = self.irrep[i].basis[name_i]
                    for j in range(i):
                        if self.irrep[j].multi == 0:
                            continue
                        else:
                            for jmul in range(self.irrep[j].multi):
                                name_j  = 'multi' + str(jmul)
                                basis_j = self.irrep[j].basis[name_j]
                                isOrtho = isOrthogonal(basis_i,basis_j,self.iprint)
                                isOrtho_tmp.append(isOrtho)
                                if self.iprint == 3:
                                    print('')
                                    print('multi:',imul,'of',self.irrep[i].label,'<-->','multi:',jmul,'of',self.irrep[j].label,\
                                            ' : ', isOrtho)
                                    print('')
        if all(isOrtho_tmp) :
            print(4*' '+'---> '+'success')
    # should execute after the self.diag_ham has been executed
    def cal_degeneracy(self):
        '''
        aim : extract the information of degeneracy
        '''
        eng_flag  = 10000.0
        deg_cnt   = 0 
        deginfo_t = {} # 
        for i in range(self.dim):
            if np.abs(self.ham_eig[i]-eng_flag) < 1.0e-6:
                degeneracy_cnt += 1
                if i == self.dim -1 :
                    deginfo_t['degeneracy'] = degeneracy_cnt
                    self.deginfo.append(deginfo_t)
            else:
                deg_cnt += 1
                if i > 0 :
                    deginfo_t['degeneracy'] = int(degeneracy_cnt)
                    self.deginfo.append(deginfo_t)
                deginfo_t = {} # 
                eng_flag = self.ham_eig[i]
                deginfo_t['start'] = int(i)
                deginfo_t['energy']   = self.ham_eig[i] 
                degeneracy_cnt = 1
                if i == self.dim -1 :
                    deginfo_t['degeneracy'] = degeneracy_cnt
                    self.deginfo.append(deginfo_t)
            if self.iprint == 3:
                print('deginfo -> i(band index,from 0)=',i)
                print(5*' ',deginfo_t)
        self.deg = deg_cnt 
        if self.iprint == 3: 
            print('deg = ',self.deg)
        # check the self-consistence:
        if self.deg != len(self.deginfo):
            print('self.deg = ',self.deg,'len(self.deginfo) = ',len(self.deginfo))
            raise ValueError("ERROR in PG.MBPGsubs.cal_degeneracy")
        for i in range(self.deg):
            deg_irreplabel = None
            self.deginfo[i]['irrep'] = [] 
            for j in range(self.deginfo[i]['start'],self.deginfo[i]['start']+self.deginfo[i]['degeneracy']):
                for jj in range(self.dim):
                    if np.abs(self.ham_evc[jj,j]) > 1.0e-3 : 
                        deg_irreplabel = self.allbasis['irreplabel'][jj][0]
                        if deg_irreplabel not in self.deginfo[i]['irrep']:
                            self.deginfo[i]['irrep'].append(copy.deepcopy(deg_irreplabel))
            print('self.deginfo[i][irrep] for i=',i,' is \n',self.deginfo[i]['irrep']) if self.iprint >=2 else 0
        # print the irreps adn accidential degeneracy of every degenerate space 
        # and judge whether there are accidential degeneracy of one irreps
        print('\n')
        print('Irreps of degenerate space and information of accidential degeneracy:\n')
        print('{:15}{:10}{:10}{:20}'.format('Deg Space','Deg','Taste','Irreps'))
        for i in range(self.deg):
            deg_from_irreps = 0
            taste1   = False # record whether there is accidential degeneracy, True is yes False means One degenerate
                             # space just contains only one irrep 
            taste    = True  # record whether there is accidential degeneracte space which contain a certain degenerate 
                             # space more than once which we currently can not handle 
            break_dump = False
            if len(self.deginfo[i]['irrep']) > 1 :
                taste1 = True
            for ii in range(len(self.deginfo[i]['irrep'])):
                deg_from_irreps += self.irrep[self.irrepindex[self.deginfo[i]['irrep'][ii]]].dim
            if deg_from_irreps != self.deginfo[i]['degeneracy']:
                break_dump = True
                taste      = False
            line_t = '{:<15}{:<10}{}{}{}'.format(i+1,self.deginfo[i]['degeneracy'],taste,'   ',\
                    self.deginfo[i]['irrep'])
            print(line_t)
            self.irrep_accidental_deg.append(copy.deepcopy(taste1))
        print('\n')
        if self.iprint == 3 :
            print('Degeneracy:  \n',self.deg)
            print('Degeneracy(info):  \n',self.deginfo)
        if break_dump :
            raise ValueError('SymmAtom can not handle a accidential degenerate space due to contain a same Irreps more than once')
       
    def check_irrepbasis_final(self):
        '''
        aim : check whether it's the basis of irreps which we want after all 
              bases have been found.
        '''
        isTrue = [True]
        for ir in self.irrep :
            if ir.multi > 0:
                for imul in range(ir.multi):
                    name      = 'multi' + str(imul)
                    rep_ir_imul = []
                    basis_set = ir.basis[name]
                    error = float(0.0)
                    for iop in range(self.nop):
#                       op = np.transpose(self.op[iop])
                        op = self.op[iop]
                        rep_matrix  = find_rep(basis_set,op,self.iprint) 
                        rep_ir_imul.append(rep_matrix)
                        error = np.sum(np.abs(rep_matrix - ir.matrices[iop]))
                        if np.abs(error) > 1.0e-6:
                            isTrue.append(False)
                            print('Error : ir=',ir.label,'imul=',imul,'iop',iop,'character:',ir.matrices[iop],'rep',rep_matrix,' error =',error)
                            print('   op:\n',op)
                            print('   basis_set:\n',basis_set)
                            raise ValueError('ERROR: basis is not complete : ', error)
                        else:
                            isTrue.append(True)
                            if self.iprint >= 2:
                                print('Success : ir=',ir.label,'imul=',imul,'iop',iop,'character:',ir.matrices[iop],'rep',rep_matrix,' error =',error)
        if all(isTrue) :
            print('\n',4*' '+'--->'+'PASS')
        else :
            print(4*' '+'--->'+'Failed')

    # decompose eigen wavefunction mixing irreps belongs to different columns of a irrep because of energy degeneracy
    def decompose_degenerate(self):
        '''
        aim : decompose eigen wavefunction mixing irreps belongs to different columns of a irrep because of energy degeneracy
        '''
        logic_all = []
        self.ham_evc_decompose = copy.deepcopy(self.ham_evc)
        for i in range(self.deg):
            basis_list = []
            cnt_dict = int(-1)
            label_dict = {}
            # collect basis
            for j in range(self.deginfo[i]['start'],self.deginfo[i]['start']+self.deginfo[i]['degeneracy']):
                cnt_dict += 1
                basis_list.append(self.ham_evc[:,j])
                label_dict[str(cnt_dict)] = j
            work_array  = np.transpose(np.array(basis_list))

            if self.iprint == 3:
                print('work array :',work_array)
                print('Deg(from 0)             : ', i)
#           print('the shape of work_array : ',work_array.shape)

            len_sp = self.ham_evc.shape[0]
            logic_list = []
            cnt_list = []
            for jj in range(len(basis_list)): 
#               print('IMPORTANT CHECK : ',basis_list[jj][0:5])
                if self.iprint == 3 :
                    print('     band index : ',label_dict[str(jj)])
                label_t = None
                logic_t = True
                cnt_t   = 0
                for ii in range(len_sp):
                    if np.abs(basis_list[jj][ii]) > 1.0E-7:
                        cnt_t += 1
                        if label_t == None:
                            label_t = self.allbasis['irreplabel'][ii]
                        else:
                            if label_t != self.allbasis['irreplabel'][ii]:
                                logic_t = False
                logic_list.append(logic_t)
                cnt_list.append(cnt_t)
            print('logical                 : ',logic_list,'  Deg: ',cnt_list) if self.iprint >= 2 else 0
            # decompose : 
            # check whether its nonzero elements only lies in a single irreps space
            # if yes : pass
            # if no  : error
            # now the irrep_t may contain many irreps
            irrep_t = [self.irrep[self.irrepindex[self.deginfo[i]['irrep'][irrr]]] for irrr in range(len(self.deginfo[i]['irrep']))]

            if self.iprint == 3:
                print(' irrep now :',irrep_t)
            if not all(logic_list) :
                multi_list = []
                for iir in range(len(irrep_t)):
                    multi_list.append(irrep_t[iir].multi)
                # make permutation and combination
                permu_list = []
                permu_t = []
                for ilen in range(len(multi_list)):
                    permu_t1 = [imm for imm in range(multi_list[ilen])]
                    permu_t.append(permu_t1)
                permu_t2 = itertools.product(*permu_t)
                for iele in permu_t2:
                    permu_list.append(list(iele))

                # base 
                ibase = np.zeros(len(irrep_t),dtype=np.int32)
                for iirr in range(len(irrep_t)):
                    for kk in range(self.irrepindex[irrep_t[iirr].label]):
                        if self.irrep[kk].multi > 0:
                            ibase[iirr] += self.irrep[kk].multi * self.irrep[kk].dim
                    if self.iprint == 3:
                        print('ibase for ',irrep_t[iirr].label, 'is : ',ibase[iirr])

                # make the dimension for work space
                len_eff_basis = 0
                for iirr in range(len(irrep_t)):
                    len_eff_basis += irrep_t[iirr].dim

                for ilist in permu_list :
                    # initial
                    find_label = False

#                   for k in range(irrep_t[iirr].multi):
                   
                    basis_array = np.zeros((len_eff_basis,len_eff_basis),dtype=np.complex128)
                    base_t = 0
                    for icol in range(len(irrep_t)):
                        basis_array[base_t:base_t+irrep_t[icol].dim,:] = \
                                work_array[(ibase[icol]+ilist[icol]*irrep_t[icol].dim):(ibase[icol]+(ilist[icol]+1)*irrep_t[icol].dim),:]
                        base_t += irrep_t[icol].dim
                    # judge whether the basis_array is reversible?
                    if self.iprint == 3:
                        print('>>> matrix start to make decomposition :\n',basis_array)
                    if np.abs(np.linalg.det(basis_array)) > 1e-3 :
                        find_label = True
                        print(">>> index for multi of basis decomposition(from 0) : ",ilist) if self.iprint == 3 else 0
                        break
                    elif ilist == permu_list[-1] :
                        print(">>> decomposition(from 0) : Failed in multi=",ilist,'determinant =',np.linalg.det(basis_array))
                        print(">>> index for multi of basis decomposition(from 0) : Failed")
                        raise ValueError(">>> Can not find proper matrix to decompose")
                    else:
                        if self.iprint >= 2:
                            print(">>> decomposition(from 0) : Failed in multi=",'determinant =',np.linalg.det(basis_array))
                if find_label :
                    if self.iprint ==3 :
                        print('>>> decomposition matrix original :\n',basis_array)
                        print('>>> decomposition matrix :\n',np.linalg.inv(basis_array))
                    work_array_new = np.dot(work_array,np.linalg.inv(basis_array))
                    self.ham_evc_decompose[:,self.deginfo[i]['start']:self.deginfo[i]['start']+self.deginfo[i]['degeneracy']] =\
                            copy.deepcopy(work_array_new[:,:])
            else:
                if self.iprint >= 2 :
                    print(">>> basis decomposition(from 0) : no Need")
                self.ham_evc_decompose[:,self.deginfo[i]['start']:self.deginfo[i]['start']+self.deginfo[i]['degeneracy']] =\
                        copy.deepcopy(work_array[:,:])
        if self.ham_evc.shape == self.ham_evc_decompose.shape :
            for icnt in range(self.ham_evc_decompose.shape[1]):
                self.ham_evc[:,icnt] =RepBasisNorm(self.ham_evc_decompose[:,icnt]) 
        else :
            raise ValueError('Error in decompose_degenerate: dim(self.ham_evc)=',self.ham_evc.shape,\
                    'dim(self.ham_evc_decompose.shape)=',self.ham_evc_decompose.shape)
#       self.ham_evc = copy.deepcopy(self.ham_evc_decompose)

    def check_basis_irreps1(self):
        '''
        aim : check whether the nonzero elements only lies in a single irreps space and do not cross between different
              irreps space
        '''
        len_sp = self.ham_evc.shape[0]
        logic_list = []
        label_list = []
        number_list= []
        for j in range(len_sp):
            cnt_t   = 0
            label_t = []
            logic_t = True
            num_t   = 0
            for i in range(len_sp):
                if np.abs(self.ham_evc[i,j]) > 1.0E-10:
                    cnt_t += 1
                    if label_t == []:
                        label_t.append(self.allbasis['irreplabel'][i][0])
                        num_t = 1
                    else:
                        if self.allbasis['irreplabel'][i][0] not in label_t :
                            label_t.append(self.allbasis['irreplabel'][i][0])
                            num_t += 1
                            logic_t = False
            logic_list.append(logic_t)
            label_list.append(label_t)
            number_list.append(num_t)
            if self.iprint >= 2 :
                print('number of nonzeros elements of ',j,'th column is : ',cnt_t)
                print('irrep label in ',j,'th column is \n',label_t)
                print('logical label in ',j,'th column is \n',logic_t)
                print('number of different irreps in ',j,'th column is \n',num_t)
        if max(number_list) == 1 or logic_list == []:
            print('\n',4*' '+'---> PASS')
        else:
            print('irrep label in is \n',label_list)
            print('logical label is \n',logic_list)
            print('number of different irreps is \n',number_list)
            print('\n',4*' '+' need to decompose ...')
#           raise ValueError("ERROR in PG.MBPGsubs.check_basis_irreps1")


    # check the structure of irreps of basis
    def check_basis_irreps(self):
        '''
        aim : check whether the nonzero elements only lies in a single column of a single irreps space and do not cross between different
              irreps space and different columns in a irreps.
        '''
        len_sp = self.ham_evc.shape[0]
        logic_list = [True]
        self.ham_evc_irrepindex['sequential'] = []
        for j in range(len_sp):
#           print('IMPORTANT CHECK : ',self.ham_evc[0:5,j])
            if self.iprint ==3 :
                print('             band index : ',j)
            cnt_t   = 0
            label_t = None
            logic_t = True
            for i in range(len_sp):
                if np.abs(self.ham_evc[i,j]) > 1.0E-10:
                    cnt_t += 1
                    if label_t == None:
                        label_t = self.allbasis['irreplabel'][i]
                        self.ham_evc_irrepindex['sequential'].append(label_t)
                    else:
                        if label_t != self.allbasis['irreplabel'][i]:
                            logic_t = False
            logic_list.append(logic_t)
            if self.iprint >= 2 :
                print(5*' ','j =',j,'deg=',cnt_t,'logic:',logic_t,'label = ',label_t)
        if self.iprint >= 2 :
            print('irreps info of eigen wavefunctions :\n',self.ham_evc_irrepindex['sequential'])
        if all(logic_list) : 
            print('\n',4*' '+'---> PASS')
        else :
            print('\n',4*' '+'---> Failed')



    # this method is valid to use only in the case that the self.irrep.projectors() has been assigned
    def check_projectors(self):
        proj_is_good_all = [True]
        if self.iprint >= 2 :
            print("")
            print(10*'* *')
            print('    ---> I am checking the properties of projectors : W^i_{np}W^j_{qr} = W^i_{nr}\delta_{ij}\delta_{pq}')
            print(10*'* *')
            print("")
        for i in range(len(self.irrep)):
            Proi = []
            for idim in range(self.irrep[i].dim):
                Proi.append(self.irrep[i].projector['P'+str(0)+str(idim)])
                if idim < self.irrep[i].dim-1:
                    Proi.append(self.irrep[i].projector['P'+str(idim+1)+str(idim)])
                
            for j in range(len(self.irrep)):
                proj_is_good = []
                Proj = []
                if j != i :
                    for jdim in range(self.irrep[j].dim):
                        Proj.append(self.irrep[j].projector['P'+str(0)+str(jdim)])
                        if jdim < self.irrep[j].dim-1:
                            Proj.append(self.irrep[j].projector['P'+str(jdim+1)+str(jdim)])

                    for ipro in Proi:
                        for jpro in Proj:
                            if np.sum(np.abs(np.dot(ipro,jpro))) > 1.0E-6:
                                proj_is_good.append(False)
                            else : 
                                proj_is_good.append(True)
                if j == i :
                    for ipro in range(self.irrep[i].dim):
                        pro_t = self.irrep[i].projector['P'+str(ipro)+str(ipro)]
                        if np.sum(np.abs(np.dot(pro_t,pro_t) - pro_t)) > 1.0E-6:
                            proj_is_good.append(False)
                        else : 
                            proj_is_good.append(True)
                proj_is_good_all += proj_is_good
                print('i = ',self.irrep[i].label,'j = ',self.irrep[j].label) if self.iprint == 3 else 0
                if all(proj_is_good) == True or proj_is_good == [] :
                    print('     : beautiful~~~', proj_is_good) if self.iprint >= 2 else 0
                else:
                    print('     : sad ~~~',proj_is_good) if self.iprint >= 2 else 0
        if all(proj_is_good_all) == True:
            print(4*' '+'---> Pass')
        else:
            print('\n',4*' '+'---> Failed')


#   def Fock2irreps(self):
#       '''
#       aim : after MBPGsubs.Cal_ReductRep() has been well executed,
#             we need to construct the unitary matrix which transfroming the Fock basis to 
#             the irreps basis in each N subspace
#       '''
#       for i in range(len(self.irrep)):
#           if self.irrep[i].multi == 0 :
#               continue
#           else:
#               for imul in range(self.irrep[i].multi):
#                   name_i  = 'multi' + str(imul)
#                   basis_i = self.irrep[i].basis[name_i]

class TranOrb():
    '''
    aim:  once we know the representation of an operator in the space of (x,y,z),
          how to get the representation of that operator in the space of p or d or f shell
    input: 
          dim     =  3  -> decide the dimension of coordinates
          npower  =  2  -> decide the maximum power of basis function: for example npower=2 for example 
                           dxy= sqrt(15/4*pi)*{xy \over r^2}
          npoly   =  np.array([1,1,1])  -> the number of polynomials of each basis function for example
                                           np.array([1,1,1]) for {dyz, dxz, dxy}
          nfunc   =  3  -> the number of basis function, for example n=3 for t2g of d shell, namely "dxy, dyz, dxz"
                        respectively
          umat    =  [ndarray,ndarray] -> the unitary transformation matrxi of the target operator of corresponding point group
                                3*3 numpy array for 3-D real space and 2*2 for 2-D real space
                     note that the umat is not R(\hat{n},\phi) but R^{-1}(\hat{n},\phi) = R(\hat{n},-\phi)
                     the former one evaluate the new corrdinates of a point in the old cartisian system(r) after
                     transformation, but the latter one means the new corrdinates of a point in the new cartisian
                     system(r'=Rr) after the same transformation
                     the formula of R(\hat{n},\phi) operate on a vec is :
                     r' = Rr or r = R^{-1}r' not r'^T = r^T R like that of representation theory of group theory
      [O] shell   = "d" -> decide the power of basis function
      [O] func_o  = dict{'struct'=[[[]],[[]],[[]]],'coeff'=[[],[],[]]} first : nfunc * npoly * dim, the structure, 
                                                      second: nfunc * npoly      , the coefficients(you can drop the
                                                      over coefficients)
    output: umat_new [the most important output]
    note: 
          1. at least one of shell and func_o should be give
          2. orbital order should better be arranged as that in userguide of wannier90 if you just give "shell" instead
             of both "shell" and "fun_o" because we will automatically produce the fun_o in the order of wannier90's
             convention
    '''
    def __init__(self,umat,su2,npoly,dim=3,npower=2,nfunc=3,shell=None,func_o={},iprint=1):
        '''
        input: 
               mapbasis : map index (the index of "x2" is "200", and index of "xyz" is "111") to their position in the
                       basis function space 
        '''
        self.dim       = int(dim)
        self.npower    = int(npower)
        self.npoly     = np.array(npoly,dtype=np.int32)
        self.nfunc     = int(nfunc)
        self.umat      = umat
        self.su2       = su2
        self.shell     = shell
        self.func_o    = func_o
        self.iprint    = iprint
        self.mapbasis  = {}
        self.nop       = len(self.umat)
        self.vec_oldbasis = None
        self.umat_new     = [] # the resulting matrix transforms in rows, should be noticed when using it
        self.umat_so      = [] # also in rows
        
        assert self.shell in ['s','p','t2g','d','f'],   "I can not recognize shell"
        assert self.dim > 0, "dim should great than zero as input of tran_orbital"
        assert self.func_o != {} or self.shell != None, "you should give at least one of <shell> and <func_o>"
#
        self.show_attribute() if self.iprint == 3 else 0
        self.def_func()
        self.make_mapbasis()
        if self.iprint == 3:
            print('mapbasis=',self.mapbasis)
        self.make_vec_oldbasis()
        if self.iprint == 3:
            print('vec_oldbasis=',self.vec_oldbasis)
        self.make_func_new()
        self.make_umat_so()

    # the most last method to be use in this class:
    def make_umat_so(self):
        for i in range(self.nop):
            mat = np.kron(self.umat_new[i],np.transpose(self.su2[i]))
            self.umat_so.append(mat)

    def check_symm_soc(self,soc_matrix):
        ndim = int(self.nfunc) * 2
        soc_new = np.zeros((ndim,ndim),dtype=np.complex128)
        # below is proper for double group but not for single group
        for i in range(self.nop):
            # note rep_vec is R^{-1}
            soc_new += np.dot(np.conjugate(self.umat_so[i]),np.dot(soc_matrix,np.transpose(self.umat_so[i])))
        soc_new = soc_new / float(self.nop)

        error = np.sum(np.abs(soc_matrix - soc_new))
        if error > 1.0e-6 :
            dump_2dc(ndim,ndim,soc_new,path='soc_new.dat',prec=1.0e-6)
            dump_2dc(ndim,ndim,soc_matrix,path='soc.dat',prec=1.0e-6)
            dump_2dc(ndim,ndim,soc_matrix - soc_new,path='soc_diff.dat',prec=1.0e-6)
            raise ValueError('error in check_symm_soc and error value : ',error)
        elif self.iprint >= 2 :
            print('    ---> success : Check for the symmetry of soc :',error)

    def check_symm_crystal(self,crystal_matrix):
        ndim = int(self.nfunc)
        crystal_new = np.zeros((ndim,ndim),dtype=np.complex128)
        if self.iprint == 3:
            print('ndim:\n',ndim)
            print('nop:\n',self.nop)
            print('self.umat_new:\n',self.umat_new[0].shape)
            print('crystal_matrix:\n',crystal_matrix.shape)
        # below is proper for double group but not for single group
        for i in range(int(self.nop/2)):
            # note rep_vec is R^{-1}
            crystal_new += np.dot(np.conjugate(self.umat_new[i]),np.dot(crystal_matrix[0:2*ndim:2,0:2*ndim:2],np.transpose(self.umat_new[i])))
        crystal_new = crystal_new*2 / float(self.nop)

        error = np.sum(np.abs(crystal_matrix[0:2*ndim:2,0:2*ndim:2] - crystal_new))
        if error > 1.0e-6 :
            dump_2dc(ndim,ndim,crystal_new,path='cfd_new.dat',prec=1.0e-6)
            dump_2dc(ndim,ndim,crystal_matrix[0:2*ndim:2,0:2*ndim:2],path='cfd.dat',prec=1.0e-6)
            dump_2dc(ndim,ndim,crystal_matrix[0:2*ndim:2,0:2*ndim:2] - crystal_new,path='cfd_diff.dat',prec=1.0e-6)
            raise ValueError('error in check_symm_crystal and error value : ',error)
        elif self.iprint >= 2 :
            print('    ---> success : Check for the symmetry of crystal :',error)

    def show_attribute(self):
        print('dim=',self.dim)
        print('npower=',self.npower)
        print('npoly=',self.npoly)
        print('nfunc=',self.nfunc)
        print('umat=',self.umat)
        print('nop=',self.nop)
        print('shell=',self.shell)
        print('func_o=',self.func_o)
        print('mapbasis=',self.mapbasis)
        print('vec_oldbasis=',self.vec_oldbasis)
        print('umat_new=',self.umat_new)
        print('su2=',self.su2)
        print('umat_so=',self.umat_so)

    def def_func(self):
        '''
        define func_o respect to the given value of shell
        '''
        ph1=np.sqrt(5/2);ph2=2*np.sqrt(15);ph3=np.sqrt(3/2)
        ph5=np.sqrt(3/2);ph6=np.sqrt(15);ph7=np.sqrt(5/2)
        struct = np.zeros((self.nfunc,np.max(self.npoly),self.dim),dtype=np.int32)
        coeff  = np.zeros((self.nfunc,np.max(self.npoly)),dtype=np.float64)
        if self.shell == None:
            print("you should not use this method if you do not give <shell>!")
        elif self.shell == 'f':
            # orbital order : fz3,fxz2,fyz2,fz(x2-y2),fxyz,fx(x2-3y2),fy(3x2-y2)
            # the polynomial order : as hamonic sphere function table in wikipedia
            # for example : fz3 = 2z3 -3zx2 -3zy2 (I always like drop the overall coefficients, for example
            # 1/4*sqrt(7/pi)/r3 here.
            # fz3
            struct[0,0,0:3] = np.array([0,0,3])
            struct[0,1,0:3] = np.array([2,0,1])
            struct[0,2,0:3] = np.array([0,2,1])
            coeff[0,0] = 2
            coeff[0,1] =-3
            coeff[0,2] =-3
            # fxz2
            struct[1,0,0:3] = np.array([1,0,2])
            struct[1,1,0:3] = np.array([3,0,0])
            struct[1,2,0:3] = np.array([1,2,0])
            coeff[1,0] = ph5 * 4
            coeff[1,1] =-ph5 
            coeff[1,2] =-ph5
            # fyz2
            struct[2,0,0:3] = np.array([0,1,2])
            struct[2,1,0:3] = np.array([2,1,0])
            struct[2,2,0:3] = np.array([0,3,0])
            coeff[2,0] = ph3 * 4
            coeff[2,1] =-ph3 
            coeff[2,2] =-ph3
            # fz(x2-y2)
            struct[3,0,0:3] = np.array([2,0,1])
            struct[3,1,0:3] = np.array([0,2,1])
            coeff[3,0] = ph6
            coeff[3,1] =-ph6 
            # fxyz
            struct[4,0,0:3] = np.array([1,1,1])
            coeff[4,0] = ph2
            # fx(x2-3y2)
            struct[5,0,0:3] = np.array([3,0,0])
            struct[5,1,0:3] = np.array([1,2,0])
            coeff[5,0] = ph7
            coeff[5,1] =-ph7 * 3 
            # fy(3x2-y2)
            struct[6,0,0:3] = np.array([2,1,0])
            struct[6,1,0:3] = np.array([0,3,0])
            coeff[6,0] = ph1 * 3
            coeff[6,1] =-ph1 
        elif self.shell == 't2g':
            # orbital order : wannier90 m=2,3,5 : dxz, dyz, dxy
            struct[0,0,0:3] = np.array([1,0,1])
            struct[1,0,0:3] = np.array([0,1,1])
            struct[2,0,0:3] = np.array([1,1,0])
            coeff[0,0] = 1
            coeff[1,0] = 1
            coeff[2,0] = 1
        elif self.shell == 'd':
            ph1 = float(1.0);ph2 = 2.0*np.sqrt(3);ph3=ph2
            ph4 = np.sqrt(3);ph5=ph2
            # orbital order : wannier90  : dz2, dxz, dyz, dx2-y2, dxy
            # z2
            struct[0,0,0:3] = np.array([2,0,0])
            struct[0,1,0:3] = np.array([0,2,0])
            struct[0,2,0:3] = np.array([0,0,2])
            # xz
            struct[1,0,0:3] = np.array([1,0,1])
            # yz
            struct[2,0,0:3] = np.array([0,1,1])
            # x2-y2
            struct[3,0,0:3] = np.array([2,0,0])
            struct[3,1,0:3] = np.array([0,2,0])
            # xy
            struct[4,0,0:3] = np.array([1,1,0])
            # z2
            coeff[0,0] = -1
            coeff[0,1] = -1
            coeff[0,2] =  2
            # xz
            coeff[1,0] = ph2
            # yz
            coeff[2,0] = ph3
            # x2-y2
            coeff[3,0] = ph4
            coeff[3,1] =-ph4
            # xy
            coeff[4,0] = ph5
        else:
            raise ValueError("I can not recognize the value of <shell>:",shell)
        #
        self.func_o['struct'] = struct
        self.func_o['coeff']  = coeff
            
        
    def make_mapbasis(self):
        '''
        aim    : I map every basis to a vector and then comes a problem:
                 if I know a polynomia(x2 in dz2), how do I know the position of x2 in the set of basis of the vection rep.
        output : self.label(x2) = 2 (2 means the third position in vection)
        '''
        cnt = int(-1)
        for k in range(self.npower+1):
            for j in range(self.npower+1):
                for i in range(self.npower+1):
                    if abs(k+j+i - self.npower) < 1.0E-6 :
                        cnt = cnt + 1
                        lab_tmp = str(i) + str(j) + str(k)
                        self.mapbasis[lab_tmp] = cnt


    # should use this func after self.func_o has been assigned
    def make_vec_oldbasis(self):
        '''
        I will construct the vec representation of the old function
        '''
        vec_old = np.zeros((len(list(self.mapbasis.values())),self.nfunc),dtype=np.float64)
        for ifunc in range(self.nfunc):
            for ipoly in range(self.npoly[ifunc]):
                name_tmp = str(self.func_o['struct'][ifunc,ipoly,0]) \
                        + str(self.func_o['struct'][ifunc,ipoly,1]) \
                        + str(self.func_o['struct'][ifunc,ipoly,2])
                vec_old[self.mapbasis[name_tmp],ifunc] = self.func_o['coeff'][ifunc,ipoly]
        self.vec_oldbasis = vec_old


    def make_func_new(self):
        '''
        <the main part of this class>
        aim: make transfromation matrix of operators of point group in the space of s p d f orbitals
        '''
        umat_new_list = []
        for iop  in range(self.nop):
#           the components arrange as rows
            umat_tmp = self.umat[iop]
            umat_new = np.zeros((self.nfunc,self.nfunc),dtype=np.float64)
            for ifunc in range(self.nfunc):
                vec_ifunc = np.zeros(len(list(self.mapbasis.values())),np.float64)
                for ipoly in range(self.npoly[ifunc]):
                    vec_ifunc_ipoly = np.zeros(len(list(self.mapbasis.values())),np.float64)
                    poly_struct = self.func_o['struct'][ifunc,ipoly,:]
                    poly_coeff  = self.func_o['coeff'][ifunc,ipoly]
#                   print('ifunc=',ifunc,'ipoly=',ipoly)
#                   print('poly_struct:',poly_struct)
#                   print('poly_coeff:',poly_coeff)
#                   print('umat_tmp:\n',umat_tmp)
                    vec_ifunc_ipoly = self.tran_func(poly_struct,poly_coeff,umat_tmp)
#                   print('vec_ifunc_ipoly\n',vec_ifunc_ipoly)
                    vec_ifunc += vec_ifunc_ipoly
#                   print(2*'-')
#               print(5*'-')
                # I will write a new function to decompose 
#               print('ifunc=',ifunc,'ipoly=',ipoly)
#               print('vec_ifunc\n',vec_ifunc)
                rep_new_tmp = decompose_vec(self.vec_oldbasis,vec_ifunc,iprint = self.iprint)
#               the components arrange as rows
                umat_new[ifunc,:] = rep_new_tmp
            umat_new_list.append(umat_new)
        self.umat_new = umat_new_list


    def tran_func(self,struct_t,coeff_t,umat):
        '''
        aim : I will tranfrom (for example <xy>) to a linear combination of the basis (for example <dxy> <dyz> <dxz>)
        input :
               struct_t = [1,1,0] (for xy function)
               coeff_t  = real number which is the coefficient of the polynomia (xy for example) in the basis function
                          (for example dxy here, which in this case coeff_t is 1.0) 
               umat     = the representation matrix of operators of point group in cartesian space
                          (this function cope with just one operator once, for loop over every umat
                          of self.umat will be done in make_func_new function)
        note : * should first assign self.mapbasis
               * the umat is transfrom in rows(which means r' = umat * r) not like that in group theory in columns
        '''
        # construct new_struct,for example struct_t = np.array([2,1,1]), umat=np.array([[1,1,0],[0,1,1],[1,1,1]])
        # new_struct = np.array([[1,1,0],[1,1,0],[0,1,1],[1,1,1]])
        basis_vec = np.zeros(len(list(self.mapbasis.values())),np.float64)
        new_struct = np.zeros((self.npower,self.dim),dtype=np.float64)
        cnt_t = int(-1)
        for idim in range(self.dim):
            if struct_t[idim] > 0 :
                for imulti in range(struct_t[idim]):
                    cnt_t = cnt_t + 1
                    new_struct[cnt_t,:]  = umat[idim,:]
        new_struct = new_struct 
        # pmutt
        pmutt_basis,pmutt_coeff = pmutt(new_struct) 
        for i in range(len(pmutt_basis)):
            corrd_i = self.mapbasis[pmutt_basis[i]]
            basis_vec[corrd_i] = pmutt_coeff[pmutt_basis[i]]

        return basis_vec * coeff_t

        
#   @staticmethod
#   def pmutt_index_list(drow,dcol=int(2)):
#       '''
#       it's an auxiliary function of pmutt
#       aim : to produce all the index of a permutation which means all the possibility of taking only one element
#             from every row of drow * dcol array
#       input : (the shape of the array)
#               drow: (integer)
#           [O] dcol: (integer) if not given, it will be set as 3(x,y,z) automatically 
#       output  : ndarray of all the possible permutation
#                 list of all the possible permutation
#       example : input  : drow = 2, dcol = 2
#                 output : np.array([[0,0],[0,1],[1,0],[1,1]])
#                          ['00','01','10','11']
#       '''
#       cnt = int(0)
#       if  int(drow) == 1 : 
#           for i in range(dcol):
