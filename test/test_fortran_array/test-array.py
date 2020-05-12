#! /usr/bin/env python

from test import test_array
import numpy as np

a = np.array([0,1,2,3,4],dtype=np.int32)
norb = np.int(4)
b = test_array(norb,a)
print('b=',b)
