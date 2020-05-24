#! /usr/bin/env python

import numpy as np

OPR1 = np.diag(np.array([1,1,1,1,1],dtype=np.complex128))
OPS1 = np.array([[1,0],[0,1]],dtype=np.complex128)
IRC1 = np.array([[1,0],[0,1]],dtype=np.complex128)

OPR2 = np.diag(np.array([1,-1,-1,1,1],dtype=np.complex128))
OPS2 = np.array([[-1j,0],[0,1j]],dtype=np.complex128)
IRC2 = np.array([[0,-1],[1,0]],dtype=np.complex128)

OPR3 = np.diag(np.array([1,1,-1,1,-1],dtype=np.complex128))
OPS3 = np.array([[0,-1],[1,0]],dtype=np.complex128)
IRC3 = np.array([[0,-1j],[-1j,0]],dtype=np.complex128)

OPR4 = np.diag(np.array([1,-1,1,1,-1],dtype=np.complex128))
OPS4 = np.array([[0,-1j],[-1j,0]],dtype=np.complex128)
IRC4 = np.array([[-1j,0],[0,1j]],dtype=np.complex128)

OPR5 = np.diag(np.array([1,1,1,1,1],dtype=np.complex128))
OPS5 = np.array([[1,0],[0,1]],dtype=np.complex128)
IRC5 = np.array([[1,0],[0,1]],dtype=np.complex128)

OPR6 = np.diag(np.array([1,-1,-1,1,1],dtype=np.complex128))
OPS6 = np.array([[-1j,0],[0,1j]],dtype=np.complex128)
IRC6 = np.array([[0,-1],[1,0]],dtype=np.complex128)

OPR7 = np.diag(np.array([1,1,-1,1,-1],dtype=np.complex128))
OPS7 = np.array([[0,-1],[1,0]],dtype=np.complex128)
IRC7 = np.array([[0,-1j],[-1j,0]],dtype=np.complex128)

OPR8 = np.diag(np.array([1,-1,1,1,-1],dtype=np.complex128))
OPS8 = np.array([[0,-1j],[-1j,0]],dtype=np.complex128)
IRC8 = np.array([[-1j,0],[0,1j]],dtype=np.complex128)

OPR9 = np.diag(np.array([1,1,1,1,1],dtype=np.complex128))
OPS9 = np.array([[-1,0],[0,-1]],dtype=np.complex128)
IRC9 = np.array([[-1,0],[0,-1]],dtype=np.complex128)

OPR10= np.diag(np.array([1,-1,-1,1,1],dtype=np.complex128))
OPS10= np.array([[1j,0],[0,-1j]],dtype=np.complex128)
IRC10= np.array([[0,1],[-1,0]],dtype=np.complex128)

OPR11= np.diag(np.array([1,1,-1,1,-1],dtype=np.complex128))
OPS11= np.array([[0,1],[-1,0]],dtype=np.complex128)
IRC11= np.array([[0,1j],[1j,0]],dtype=np.complex128)

OPR12= np.diag(np.array([1,-1,1,1,-1],dtype=np.complex128))
OPS12= np.array([[0,1j],[1j,0]],dtype=np.complex128)
IRC12= np.array([[1j,0],[0,-1j]],dtype=np.complex128)

OPR13= np.diag(np.array([1,1,1,1,1],dtype=np.complex128))
OPS13= np.array([[-1,0],[0,-1]],dtype=np.complex128)
IRC13= np.array([[-1,0],[0,-1]],dtype=np.complex128)

OPR14= np.diag(np.array([1,-1,-1,1,1],dtype=np.complex128))
OPS14= np.array([[1j,0],[0,-1j]],dtype=np.complex128)
IRC14= np.array([[0,1],[-1,0]],dtype=np.complex128)

OPR15= np.diag(np.array([1,1,-1,1,-1],dtype=np.complex128))
OPS15= np.array([[0,1],[-1,0]],dtype=np.complex128)
IRC15= np.array([[0,1j],[1j,0]],dtype=np.complex128)

OPR16= np.diag(np.array([1,-1,1,1,-1],dtype=np.complex128))
OPS16= np.array([[0,1j],[1j,0]],dtype=np.complex128)
IRC16= np.array([[1j,0],[0,-1j]],dtype=np.complex128)

OPR_list = [OPR1,OPR2,OPR3,OPR4,OPR5,OPR6,OPR7,OPR8,OPR9,OPR10,OPR11,OPR12,OPR13,OPR14,OPR15,OPR16]
OPS_list = [OPS1,OPS2,OPS3,OPS4,OPS5,OPS6,OPS7,OPS8,OPS9,OPS10,OPS11,OPS12,OPS13,OPS14,OPS15,OPS16]
IRC_list = [IRC1,IRC2,IRC3,IRC4,IRC5,IRC6,IRC7,IRC8,IRC9,IRC10,IRC11,IRC12,IRC13,IRC14,IRC15,IRC16]

# constants
ndimir  = IRC1.shape[0]
nrank = len(OPR_list)
ndim = OPR1.shape[0]

OPSO_list = []
for i in range(nrank):
    OPSO_list.append(np.kron(OPR_list[i],OPS_list[i]))


P00 = np.zeros((ndim*2,ndim*2),dtype=np.complex128)
P10 = np.zeros((ndim*2,ndim*2),dtype=np.complex128)
for i in range(nrank):
    print(10*' ')
    print('Projectors : j=',i)
    print(5*' ','characters : \n',IRC_list[i])
    print(5*' ','op, matrix : \n',OPSO_list[i])
    P00 += float(ndimir)/float(nrank) * np.conjugate(IRC_list[i][0,0]) * OPSO_list[i]
    P10 += float(ndimir)/float(nrank) * np.conjugate(IRC_list[i][1,0]) * OPSO_list[i]
print('P00:\n',P00)
print('P10:\n',P10)



