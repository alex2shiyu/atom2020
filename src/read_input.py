#! /usr/bin/env python

import numpy as np
from module.mod_read import read_2dc
from wanntb.soc import atom_hsoc
from wanntb.tran import tran_op, tmat_c2r
from scipy.special import perm

class Atom():
    """
    Atom class which need to do the atomic problem
    property:
             norb [I] : number of orbitals(with spin, for example, norb=10 for d robitals )
             nmin [I] : for cut-off : use for caution
             nmax [I] : 
             int_type [char] : support "kanamori" now 
             int_val  [dict] : kanamori interaction parameters
                U    : kanamori interaction parameters
                Up   
                Jz
                Js
                Jp
             soc_type   [char]       : construct soc for which type of orbitals : eg? t2g? d ? p ?
             soc_val    [float64]    : soc strength
             soc_mat    [complex128] : soc matrix(single particles)
             cfd_mat    [complex128] : cfd matrix(single particles)
    """
    def __init__(self,norb,nmin,nmax,int_type,int_val,soc_type,soc_val,cfd,cfd_mat,gqn,point_group,space_group,basis_tran,amat,vpm_type):
        self.norb = norb
        self.nmin = nmin
        self.nmax = nmax
        self.int_type = int_type
        self.int_val  = int_val
        self.soc_type = soc_type
        self.soc_val  = soc_val
        self.soc_mat  = None
        self.cfd      = cfd    
        self.cfd_mat  = cfd_mat
        self.gqn      = gqn
        self.basis_tran  = basis_tran
        self.amat     = amat
        self.point_group = point_group
        self.space_group = int(space_group)
        self.vpm_type = vpm_type
        self.ncfgs    = None
        self.totncfgs = None # 2^norb, set this value just because f2py can not coye with 2**norb as the dimension of
                             # input parameters
        self.onsite   = None

    @staticmethod
    def from_incar(incar_file = 'atom2020.incar',cemat_file='atom2020.cemat.in',amat_file='atom2020.amat.in'):
        """
        Generate an Atom instance from input file : atom2020.incar
        Args: atom2020.incar
            -----
            norb  10
            nmin  7
            nmax  9
            int_type  kanamori
            int_val   3.0 2.4 0.3 0.3 0.3 # must be in the order of "U, Up, Jz, Js, Jp"
            soc_type  p
            soc_val   0.28  (unit in ev)
            cfd       yes
            gqn       3 
            point_group C4v
            space_group 71
            vpm_type    1
            -----
            #gqn:
            #1[pure]: n ;  
            #2[no soc, no cfd]: n, ps, s, sz, l, lz; 
            #3[soc, no cfd]: n, j, jz (double group of SO(3))
            #4[no soc, cdf]: n, s, sz, irrep of point group
            #5[soc, cfd]: n, irrep of double of point group
        """
        norb     = None
        nmin     = None
        nmax     = None
        int_type = None
        int_val  = {}
        soc_type = None
        soc_val  = None
        cfd      = None
        gqn      = None
        point_group = None
        vpm_type    = None
        basis_tran  = None
        amat        = None
        try:
            with open(incar_file, 'r') as f:
                fileread = f.readlines()
                for cnt, line in enumerate(fileread):
                    line1 = line.strip()
                    if line1[0:4] == "norb" :
                        norb = np.int32(line1.split()[1])
                    elif line1[0:4] == "nmin" :
                        nmin = np.int32(line1.split()[1])
                    elif line1[0:4] == "nmax" :
                        nmax = np.int32(line1.split()[1])
                    elif line1[0:8] == "int_type" :
                        int_type = line1.split()[1]
                    elif line1[0:7] == "int_val" :
                        int_val['U'] = np.float64(line1.split()[1])
                        int_val['Up'] = np.float64(line1.split()[2])
                        int_val['Jz'] = np.float64(line1.split()[3])
                        int_val['Js'] = np.float64(line1.split()[4])
                        int_val['Jp'] = np.float64(line1.split()[5])
                    elif line1[0:8] == "soc_type" :
                        soc_type = line1.split()[1]
                    elif line1[0:7] == "soc_val" :
                        soc_val = np.float64(line1.split()[1])
                    elif line1[0:3] == "cfd" :
                        cfd = line1.split()[1]
                    elif line1[0:3] == 'gqn' :
                        gqn = np.int32(line1.split()[1])
                    elif line1[0:11] == "point_group" :
                        point_group = line1.split()[1]
                    elif line1[0:11] == "space_group" :
                        space_group = int(line1.split()[1])
                    elif line1[0:10] == "basis_tran" :
                        basis_tran = line1.split()[1]
                    elif line1[0:8] == "vpm_type" :
                        vpm_type = np.int32(line1.split()[1])
                    elif line1 == '':
                        print("omit a blank line")
                    elif line1[0] in ['#','!','%']:
                        print("omit a comment line")
                    else :
                        print("I can not recognize this input option:",line1,' and it will be omitted')
#               return Atom(norb,nmin,nmax,int_type,int_val,soc_type,soc_val,cfd_mat)
#           return Atom(norb,nmin,nmax,int_type,int_val,soc_type,soc_val)
        except IOError:
            print("File:" + "\"" + incar_file + "\"" +  " doesn't exist!")
        #
        if basis_tran == 'yes' :
            try : 
                amat = read_2dc(norb,norb,amat_file)
            except IOError:
                print("File:" + "\"" +amat_file+ "\"" + " doesn't exist!")
        else : 
            amat = np.identity(norb,dtype=np.complex128) 

        if cfd == 'yes' :
            try : 
                cfd_mat_t = read_2dc(int(norb/2),int(norb/2),cemat_file)
                spin_t = np.identity(2,dtype=np.complex128)
                cfd_mat = np.kron(cfd_mat_t,spin_t)
#               cfd_mat = tran_op(cfd_mat,amat)
                return Atom(norb,nmin,nmax,int_type,int_val,soc_type,soc_val,cfd,cfd_mat,gqn,point_group,space_group,basis_tran,amat,vpm_type)
            except IOError:
                print("File:" + "\"" +cemat_file+ "\"" + " doesn't exist!")
        else : 
#           cfd_mat = None
            cfd_mat = np.zeros((norb,norb),dtype=np.complex128)
            return Atom(norb,nmin,nmax,int_type,int_val,soc_type,soc_val,cfd,cfd_mat,gqn,point_group,space_group,basis_tran,amat,vpm_type)
   
    def make_instance(self):
        self.make_ncfgs()
        self.check_incar()
        self.make_soc()
        self.make_onsite()
        self.make_totncfgs()
        self.show_incar()

#   def make_unitary_lenth(self):
#       '''
#       aim: there is function(gw_make_newui in gw_make_new.f90 file) warpered using fortran to get the unitary transformation matrix in given many body space in
#       the case we know their transformation matrix in single particle orbitals, the algorithm is very simple thus one
#       type of the array (1D or 2D) is huge and the dimension depends on the cutoff space of Fock space. This function
#       is to calculate the upper bound of the dimension using self.nmin, self.nmax, self.norb
#       '''
#       return perm(10, 3, exact=True)

    def make_onsite(self):
        '''
        I will sum atomic soc and cfd to onsite, both are on original basis
        '''
        self.onsite = self.soc_mat + self.cfd_mat
    
    # this method will use the library of pylibwyl/wanntb
    def make_soc(self):
# soc in complex basis 
        soc_c = atom_hsoc(self.soc_type,self.soc_val)
# transformation matrix from complex basis to real basis(wannier90 orbitals order)
        tmat  = tmat_c2r(self.soc_type,True)
# soc in real orbitals(wannier90 orbitals order)
        soc_r = tran_op(soc_c,tmat)
# soc in real natural basis(gutzwiller), need input transfrom matrix 
#       soc_r_natural = tran_op(soc_r,self.amat)
# assign attribute of instance 
#       self.soc_mat  = soc_r_natural
        self.soc_mat  = soc_r

    def make_ncfgs(self):
        '''
        just calculate the ncfgs of working space
        '''
        from scipy.special import comb
        ncfgs = np.int(0)
        for i in range(self.nmin,self.nmax+1):
            ncfgs = comb(self.norb, i) + ncfgs
        self.ncfgs = np.int(ncfgs)

    def make_totncfgs(self):
        '''
        calculate the total number of configurations of the whole shell
        for example : p: 64
                      d: 1024
        '''
        import math
        totncfgs = math.pow(2, self.norb)
        self.totncfgs = int(totncfgs)
#       return totncfgs

    def check_incar(self):
        """
        destination : check the reasonablity of the parameters from incar
        """
        if self.norb < 0 or self.norb > 14:
            raise IOError("<norb> has been set wrong!")
        elif self.nmin < 0 or self.nmin > self.norb or self.nmax < 0 or self.nmax < self.nmin or self.nmax > self.norb :
            raise IOError("<nmin> or <nmax> has been set wrong!")
        if self.int_type == None or self.int_val == None :
            raise IOError("<int_type> or <int_val> can not be found")
        if self.soc_type == None :
            self.soc_val = None
            print("soc will be absent in atomic hamiltonian, but may still exist in gqn!")
        elif self.soc_type not in ['s','p','d','f']:
            raise IOError("<soc_type> is not recognize and should be one of <s, p, d, f> !")
        elif self.soc_val == None :
            self.soc_val = np.float64(0.0)
            print("<soc_val> is None and has beed set to <0.0>")
        if self.gqn == None :
            self.gqn = np.int32(1)
            print("<gqn> has been set to default value: gqn= 1")
        elif self.gqn == 4 or self.gqn == 5 and self.point_group == None :
            raise IOError("if you set gqn to <4> or <5>, the <point_group> should be given")
        elif self.gqn == 2 or self.gqn ==3 and self.point_group != None :
            print("if you set gqn to <2> or <3>, the <point_group> should not be given")
        if self.cfd == 'yes' :
            print('<cfd> will be read from atom2020.cemat.in')
        elif self.cfd == 'no' :
            self.cfd = None
            print("cfd will be absent in atomic hamiltonian, but may still exist in gqn!")
        if self.basis_tran == 'no' :
            print("Are you sure the basis need not to transform to natural basis from atomic basis!")

    def show_incar(self):
        print('self.norb:\n',self.norb)
        print('self.nmin:\n',self.nmin)
        print('self.nmax:\n',self.nmax)
        print('self.ncfgs:\n',      self.ncfgs)
        print('self.int_type:\n',   self.int_type)
        print('self.int_val:\n',    self.int_val)
        print('self.soc_type:\n',   self.soc_type)
        print('self.soc_val:\n',    self.soc_val)
        print('self.soc_mat:\n',    self.soc_mat)
        print('self.cfd_mat:\n',    self.cfd_mat)
        print('self.onsite:\n',     self.onsite)
        print('self.point_group:\n',self.point_group)
        print('self.vpm_type:\n',   self.vpm_type)
