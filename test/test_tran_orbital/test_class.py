import numpy as np

class ADD():
    '''
    in order to test calling one function of the class from another function of the same class
    '''
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.c = None


    def method1(self):
        a1 = 7
        x = self.method2(a1)
        self.c = float(self.b) + x

    def method2(self,xx):
        tmp = float(xx) + float(self.a)
        return tmp
