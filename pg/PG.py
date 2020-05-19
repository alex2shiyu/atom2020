from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
import os
import sys
from pg.PG_util import pmutt, decompose_vec
from pg.read_IR_DSG import loadIR 

class Irrep:
    '''
    class of Irrep copied from yijiang's scripts
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
            print('multi is not integer : ',multi_tmp)

class ReductRep(Irrep):
    def __init__(self,label='',dim=1):
        super().__init__(label='',dim=1)
        self.multi      = 0   # dictionary to contain various projectors for a certain Irrep
        self.projector  = {}  # dictionary to contain various projectors for a certain Irrep
        self.basis      = {}  # a dict of different set of basis for different <multi>, each set contain
                              # a list of basis whose number is the dimension of the irrep. For example
                              # self.multi = 2, self.dim = 3 , then the self.basis = dict{'multi1'=[np.array([.1D.the
                              # first basis.]),[..the second one.],[...the third one ...]],
                              # 'multi2'=[np.array([...1D..the first..]),np.array([...the second...]),np.array([...the
                              # third...])]} 
    # multi should be assigned with the help of method of super class 
    def reduction(self,op):
        '''
        will reduce the reducible representation op
        input : 
            op [list]: is a list of numpy.ndarray which is the representation matrix of every operators of the PG
            <self.projector and self.multi refers to information about "this irrep" and "this input op">
        '''
        for imul in range(self.multi):
            pass

    def make_phi1(self,multi_now):
        '''
        according to intro. in Alttman's point group tables,
        W^i_11 \phi = \phi^i_1 
        however, if self.multi > 1,
        '''
        

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
        # then W10, W21, ...
        if self.dim > 1 :
            for i in range(self.dim-1):
                P_normal = np.zeros((dim_op,dim_op),dtype=np.complex128)
                name = 'P' + str(i+1) + str(i)
                for j in range(self.nrank):
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
        self.rep_vec   = []
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
                self.rep_vec  = grp.rotcar
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
            op  [list]    : a list of numpy.ndarray which is the matrix rep. of operators in natural basis
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
            irr.label = ir.label
            irr.dim   = ir.dim
            irr.matrices   = ir.matrices
            irr.characters = ir.characters
            irr.multi      = irr.cal_multi(self.op)
            irr.projector  = irr.make_projector(self.op)
            self.irrep.append(irr)


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
        self.umat_new   = []
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
        note : should first assign self.mapbasis
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
