#! /usr/bin/env python

import numpy as np
from module.mod_read import read_2dc

class Atom():
    """
    Atom class which need to do the atomic problem
    property:
             norb [I] : number of orbitals(with spin, for example, norb=10 for d robitals )
             nmin [I] : for cut-off : use for caution
             nmax [I] : 
             int_type [char] : support "kanamori" now 
             int_para [dict] : kanamori interaction parameters
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
    def __init__(self,norb,nmin,nmax,int_type,int_para,soc_type,soc_val,cfd_mat):
        self.norb = norb
        self.nmin = nmin
        self.nmax = nmax
        self.int_type = int_type
        self.int_para = int_para
        self.soc_type = soc_type
        self.soc_val  = soc_val
        self.soc_mat  = None
        self.cfd_mat  = cfd_mat

    @staticmethod
    def from_incar(incar_file = 'atom2020.incar',cemat_file='atom2020.cemat.in'):
        """
        Generate an Atom instance from input file : atom2020.incar
        Args: atom2020.incar
            -----
            norb  10
            nmin  7
            nmax  9
            int_type  kanamori
            int_para  3.0 2.4 0.3 0.3 0.3 # must be in the order of "U, Up, Jz, Js, Jp"
            soc_type  p
            soc_val   0.28 # (unit in ev)
            -----
        """
        int_para = {}
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
                    elif line1[0:8] == "int_para" :
                        int_para['U'] = np.float64(line1.split()[1])
                        int_para['Up'] = np.float64(line1.split()[2])
                        int_para['Jz'] = np.float64(line1.split()[3])
                        int_para['Js'] = np.float64(line1.split()[4])
                        int_para['Jp'] = np.float64(line1.split()[5])
                    elif line1[0:8] == "soc_type" :
                        soc_type = line1.split()[1]
                    elif line1[0:7] == "soc_val" :
                        soc_val = np.float64(line1.split()[1])
#               return Atom(norb,nmin,nmax,int_type,int_para,soc_type,soc_val,cfd_mat)
#           return Atom(norb,nmin,nmax,int_type,int_para,soc_type,soc_val)

        except IOError:
            print("File:" + "\"" + incar_file + "\"" +  " doesn't exist!")
        #
        try : 
            cfd_mat = read_2dc(norb,norb,cemat_file)
            return Atom(norb,nmin,nmax,int_type,int_para,soc_type,soc_val,cfd_mat)
        except IOError :
            print("File:" + "\"" +cemat_file+ "\"" + " doesn't exist!")
