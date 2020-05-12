#! /usr/bin/env python

from module.mod_read import read_2dr, dump_2dr
import numpy as np
from  MakeFock import *

a = np.zeros((2,3),dtype=np.float64)
a[0,0] = 1.0
a[1,2] = 3.0
dimen = np.array([2,3],dtype=np.int16)
nline = 0
ioflag =2
prec = 0.0001
NumOrb = 3
NumMin = 1
NumMax = 2
NumCfg = 6
#print(dir(MakeFock))#it lets us to determine all the attributes of a module

cg = np.zeros(NumCfg,dtype=np.int16)
cg = makecg(NumOrb, NumMin, NumMax, NumCfg)
print('cg\n',cg)
a=dump_2dr(dimen[0],dimen[1],a,'a.dat')
b = read_2dr(dimen[0],dimen[1],'a.dat')
print('size of b',b.shape)
print('b\n',b)
#a=dump_2dr(dimen[0],dimen[1],a,'a.dat',nline,ioflag,prec)

