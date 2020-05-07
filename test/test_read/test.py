#!/usr/bin/env python 

import numpy as np
#from src.incar import *
#from src.func  import *
from  module.MakeFock import *
import requests,os

# test fermi distribution
#print("the fermi distribution is : ",  FdisA(0.1,0.0))
#print("the fermi distribution is : ", dFdEnA(0.1,0.0))
#print("the fermi distribution is : ",  dFduA(0.1,0.0))
#print("{:.8f}".format(FdisA(0.02,0.0)))
#print("{:.8f}".format(dFdEnA(0.02,0.0)))
#print("{:.8f}".format(dFduA(0.02,0.0)))

# test m0fock and dMdN0
NumOrb = 3
NumMin = 1
NumMax = 2
NumCfg = 6
#print(dir(MakeFock))#it lets us to determine all the attributes of a module

cg = np.zeros(NumCfg,dtype=np.int16)
cg = makecg(NumOrb, NumMin, NumMax, NumCfg)
print("Cg",cg)
n0=np.array([0.5,0.2,0.3])
print("N0 = ",n0)
#print(dir(gs))
#NumCfg = MakeFock.verif(10,5)
#print("results is : ", NumCfg)
