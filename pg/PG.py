from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
import os
import sys
from pg.PG_util import pmutt, decompose_vec, RepBasisNorm, isOrthogonal, isindependent
from pg.read_IR_DSG import loadIR 
import copy
from pg.gramSchmidt import GramSchmidt

class Irrep:
    '''
    Irreps class 
    '''
    def __init__(self,label='',dim=1):
        self.label = label  # label for the representation. 'd' as postfix indicates double-valued irreps
        self.dim   = dim  # dimension of the representation
        self.matrices   = []  # representation matrices, ordered as the operations of the belonging group
        self.characters = []  # characters as a list, ordered the same as matrices

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
            print("multi success!")
            return int(multi_tmp)
        else:
            print('multi is not integer : ',multi_tmp,character_t.real/float(self.nrank),np.abs(multi_tmp - character_t.real/float(self.nrank)))

class ReductRep(Irrep):
    def __init__(self,label='',dim=1):
        super().__init__(label='',dim=1)
        self.multi      = 100 # dictionary to contain various projectors for a certain Irrep
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
            vectors_normortho_t = RepBasisNorm(list(vectors_new[i]))
            vectors_new_all += list(vectors_normortho_t)
            vectors_normortho.append(vectors_normortho_t)
       
        # check for the orthogonalize of all the produced basis of a certain irreps
        print(' ')
        print(10*'* * ')
        print(5*' ','>>> orthogonality[all] of the whole newly produced basis of the irrep ...')
        is_allortho = isOrthogonal(vectors_new_all,vectors_new_all) 
        print(5*' ','All is Orthogonal ?   ',is_allortho)

        
        # assign self.basis_normortho
        for i in range(self.multi):
            name_t = 'multi' + str(i)
            basis_or_t = []
            for j in range(self.dim):
                basis_or_t.append(vectors_normortho[j][i])
            self.basis_normortho[name_t] = copy.deepcopy(basis_or_t)
        print(' ')
        print(10*'* * ')
        print(5*' ','>>> orthogonality[multi] between different multi ...')
        for imul1 in range(self.multi):
            name1 = 'multi' + str(imul1)
            for imul2 in range(self.multi):
                name2 = 'multi' + str(imul2)
                print(5*'--')
                is_ortho_multi = isOrthogonal(self.basis_normortho[name1],self.basis_normortho[name2])
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
                isOrtho_no = isOrthogonal(list(vectors_n_t),list(vectors_n_t))
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
            print('')
            print('')
            print(20*'=')
            print('* * imul :',imul+1)
            phi_t1 = self.make_phi1(imul,op[0].shape[0],irrep_prev)
            print(20*'- -') 
            print('* WAVE(phi1) Done :  of multi<',imul,'> of <',self.label,'> is :',np.linalg.norm(phi_t1))
            basis_list = []
            basis_list += phi_t1
            if self.dim > 1 :
                phi_other = self.make_phi_other(phi_t1[0],imul)
                basis_list += phi_other
            print('* basis set for multi',imul,'of ',self.label,'is\n',basis_list)
            is_multi_ortho = isOrthogonal(basis_list,basis_list)
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
        for i in range(ndim):
            label_list_tmp = []
            print("\n",10*'*','\n',2*'* ','i :',i)
            func_all    = np.zeros(ndim,dtype=np.complex128)
            func_all[i] = complex(1.0)
#           for ip in range(1):
            print("\n",'. . . . . . . . . ','\n','i(P_i+1_i)=',0,'|','\n','. . . . . . . . . ')
            Proj    = self.projector[name_ip]
            func_1  = np.einsum('ij,i->j',Proj,func_all)
            func_1  = RepBasisNorm(func_1) 
            isindependt_all = []
            print(5*'- ','\n','   check independence ...')
            if int(multi_now) > 0 or len(irrep_prev) > 0 :
                print('multi_now(from 0) =',multi_now,'/',self.multi,'  ','num of irrep_prev :',len(irrep_prev))
                print(5*'- ')
                if int(multi_now) > 0 :
                    print('')
                    for imul in range(multi_now):
                        print('   independent:  multi(from 0)=',imul,'multi(now)=',multi_now)
                        mul_name  = 'multi' + str(imul)
                        basis_set = self.basis[mul_name]
                        label_t, independt   = isindependent(func_1,basis_set) 
                        isindependt_all += independt
                        label_list_tmp += label_t
                if len(irrep_prev) > 0 :
                    for irr in irrep_prev :
                        if irr.multi > 0:
                            for imul in range(irr.multi):
                                print('   independent:  multi(from 0)=',imul,'/',irr.multi)
                                mul_name = 'multi' + str(imul)
                                basis_set = irr.basis[mul_name]
                                label_t, independt   = isindependent(func_1,basis_set) 
                                isindependt_all += independt
                                label_list_tmp += label_t
                label_list.append(np.sum(np.array(label_list_tmp)))
            else:
#               return [func_1]
                label_list.append(np.sum(np.abs(func_1)))
        
        pos_target = label_list.index(max(label_list))
        print(10*'*-')
        print('* the proper phi for phi1 of mulit=',multi_now,'of irreps ',self.label,'is',pos_target+1)
        print(10*'*-')
        func_all    = np.zeros(ndim,dtype=np.complex128)
        func_all[pos_target] = complex(1.0)
        Proj    = self.projector[name_ip]
        func_1  = np.einsum('ij,i->j',Proj,func_all)
        func_1  = RepBasisNorm(func_1) 
        return [func_1]
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
            proj      = self.projector[name_proj]
#           basis_t   = np.dot(proj,phi_1)
            basis_t   = np.einsum('ij,i->j',proj,phi_1)
            print('* WAVE(phi'+str(i+2)+') Done :  of multi<',multi_now,'> of <',self.label,'> is :',np.linalg.norm(basis_t))
            basis_other.append(basis_t)
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
#               print(5*'** ','check for projectors matrix value : ')
                for j in range(self.nrank):
#                   print(10*' ')
#                   print('Projectors : j=',j)
#                   print(5*' ','characters : \n',self.matrices[j])
#                   print(5*' ','op matrix  : \n',op[j])
#                   print('')
                    P_normal[:,:] += self.dim/self.nrank * np.conjugate(self.matrices[j][i+1,i])*op[j]
                self.projector[name] = P_normal




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
        self.nclass = len(self.gen_pos)
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
                self.rep_vec  = [np.transpose(i) for i in grp.rotcar] # bilbao's data is transfromed in columns and we transpose it
                                                         # here but you should know when we use the manybody version of
                                                         # it to construct the projectors, we will use the form which
                                                         # transform in columns as same as the convention of group
                                                         # theory
                self.rep_spin = grp.su2s
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
    def __init__(self, nop, op):
        '''
        attributes:
            nop [I]       : rank of the point group
            op  [list]    : a list of numpy.ndarray which is the matrix rep. of operators in natural basis, transfrom in
                            # columns
#           irrep [list]  : a list of Irrep 
        '''
        self.nop = nop
        self.op  = op
#       self.irrep = []
#       for ir in irrep:
#           irr = Irrep()
#           irr.label = ir.label
#           irr.dim   = ir.dim
#           irr.matrices = ir.matrices
#           irr.characters = ir.characters
#           self.irrep.append(irr) 
#       print("check->irrep:\n",type(self.irrep[0])) 
#       if len(self.irrep) > 0 :
#           assert isinstance(self.irrep[0], Irrep)

    def show_attribute(self):
        print('nop:',self.nop)
        print('operators:\n',self.op)
#       print('irreps:\n',self.irrep)

class MBPGsubs(MBPG): 
    '''
    every subspace of N will have their own information when do reduciton
    '''
    def __init__(self, nop, op):
        super().__init__(nop, op)
        self.irrep = []
    
    def Cal_ReductRep(self,irrep_input):
        '''
        traverse every possible irreps of the PG to calculate multiplicities and corresponding projectors
        '''
        for ir in irrep_input:
            irr = ReductRep()
            irr.label = copy.deepcopy(ir.label)
            irr.dim   = copy.deepcopy(ir.dim)
            irr.matrices   = copy.deepcopy(ir.matrices)
            irr.characters = copy.deepcopy(ir.characters)
            print('>>>>> irrep :',irr.label)
            irr.multi      = irr.cal_multi(self.op)
            print('>>>>> multi :',irr.multi)
            irr.make_projector(self.op)
            print('>>>>> projectors :\n',irr.projector)
            if irr.multi > 0 :
                irr.reduction(self.op,self.irrep)
                irr.irrep_normortho2()
#           irr.reductstate= 'ok'
            self.irrep.append(irr)
        # test for the orthogonal of different irreps
        print(5*"* ",' orthogonality[irreps] between different irreps')
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
                                isOrtho = isOrthogonal(basis_i,basis_j)
                                print('')
                                print('multi:',imul,'of',self.irrep[i].label,'<-->','multi:',jmul,'of',self.irrep[j].label,\
                                ' : ', isOrtho)
                                print('')

    # this method is valid to use only in the case that the self.irrep.projectors() has been assigned
    def check_projectors(self):
        print("")
        print(10*'* *')
        print(' * I am checking the properties of projectors : W^i_{np}W^j_{qr} = W^i_{nr}\delta_{ij}\delta_{pq}')
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
                print('i = ',self.irrep[i].label,'j = ',self.irrep[j].label)
                if all(proj_is_good) == True or proj_is_good == [] :
                    print('     : beautiful~~~', proj_is_good)
                else:
                    print('     : sad ~~~',proj_is_good)


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
    def __init__(self,umat,su2,npoly,dim=3,npower=2,nfunc=3,shell=None,func_o={}):
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
        self.mapbasis  = {}
        self.nop       = len(self.umat)
        self.vec_oldbasis = None
        self.umat_new   = [] # the resulting matrix transforms in rows, should be noticed when using it
        self.umat_so    = []
        
        assert self.shell in ['s','p','t2g','d','f'],   "I can not recognize shell"
        assert self.dim > 0, "dim should great than zero as input of tran_orbital"
        assert self.func_o != {} or self.shell != None, "you should give at least one of <shell> and <func_o>"
#
        self.show_attribute()
        self.def_func()
        self.make_mapbasis()
        print('mapbasis=',self.mapbasis)
        self.make_vec_oldbasis()
        print('vec_oldbasis=',self.vec_oldbasis)
        self.make_func_new()
        self.make_umat_so()

    # the most last method to be use in this class:
    def make_umat_so(self):
        for i in range(self.nop):
            mat = np.kron(self.umat_new[i],self.su2[i])
            self.umat_so.append(mat)


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
                rep_new_tmp = decompose_vec(self.vec_oldbasis,vec_ifunc)
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
