import os
import re
import h5py
from functools import reduce
import numpy as np
from math import cos, sin, acos, asin, pi, sqrt

Library=os.path.abspath(os.path.dirname(__file__))

# =====================================================================================
# Lattice
# =====================================================================================
SGTricP = [1, 2]
SGMonoP = [3, 4, 6, 7, 10, 11, 13, 14]
SGMonoB = [5, 8, 9, 12, 15]
SGOrthP = [16, 17, 18, 19, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 47, \
           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]
SGOrthB1 = [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]
SGOrthB2 = [38, 39, 40, 41]
SGOrthF = [22, 42, 43, 69, 70]
SGOrthI = [23, 24, 44, 45, 46, 71, 72, 73, 74]
SGTetrP = [75, 76, 77, 78, 81, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, \
           95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 111, 112, \
           113, 114, 115, 116, 117, 118, 123, 124, 125, 126, 127, \
           128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138]
SGTetrI = [79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120, \
           121, 122, 139, 140, 141, 142]
SGTrigP = [146, 148, 155, 160, 161, 166, 167]
SGHexaP = [143, 144, 145, 147, 149, 150, 151, 152, 153, 154, 156, \
           157, 158, 159, 162, 163, 164, 165, 168, 169, 170, 171, \
           172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, \
           183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194]
SGCubcP = [195, 198, 200, 201, 205, 207, 208, 212, 213, 215, 218, \
           221, 222, 223, 224]
SGCubcF = [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]
SGCubcI = [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]


def PrimInConv(gid):
    if gid in SGTricP + SGMonoP + SGOrthP + SGTetrP + SGHexaP + SGCubcP:
        return np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    elif gid in SGMonoB + SGOrthB1:
        return np.array([[1. / 2, -1. / 2, 0.], [1. / 2, 1. / 2, 0.], [0., 0., 1.]])
    elif gid in SGOrthB2:
        return np.array([[1., 0., 0.], [0., 1. / 2, 1. / 2], [0., -1. / 2, 1. / 2]])
    elif gid in SGOrthI + SGTetrI + SGCubcI:
        return np.array([[-1. / 2, 1. / 2, 1. / 2], [1. / 2, -1. / 2, 1. / 2], [1. / 2, 1. / 2, -1. / 2]])
    elif gid in SGOrthF + SGCubcF:
        return np.array([[0., 1. / 2, 1. / 2], [1. / 2, 0., 1. / 2], [1. / 2, 1. / 2, 0.]])
    elif gid in SGTrigP:
        return np.array([[2. / 3, 1. / 3, 1. / 3], [-1. / 3, 1. / 3, 1. / 3], [-1. / 3, -2. / 3, 1. / 3]])
    else:
        raise ValueError('wrong gid !!!')


# =============================================================================
# Little group and irreps
# =============================================================================
class LittleGroup:  # little group at a special k point
    def __init__(self):
        self.gid = ''  # space group number as a string
        self.klabel = ''  # label for the k point
        self.kvecC = np.zeros(3)  # k vector
        self.rotC = []  # the rotation part of all operations as a list
        self.tauC = []  # the translation vector of all operations as a list
        self.su2s = []  # the SU2 matrices of all operations as a list, read from Bilbao
        self.irreps = []  # the irreducible representations as a list, each being an 'Irrep' object
        #
        self.kvecP = []
        self.rotP = []
        self.tauP = []
        #
        self.shiftC = np.array([0, 0, 0])
        self.shiftP = np.array([0, 0, 0])
        #
        self.basisP = []  # primitive basis, read from POSCAR
        self.su2c = []  # SU2 matrix calculated from SO3 matrix
        self.index = []  # index of SU2 matrix need to time -1, in order to align with Bilbao

        self.degenerate_pair=[] # list to collect degenerate irrep pairs
        #
        self.generators = [] # list of index of group generators
        self.antiunitary_op = antiunitary_op()
        self.coirreps = []

    def prim(self):
        # C : direct,     prim. lattices in conv.
        #   : reciprocal, conv. lattices in prim.
        # D : reciprocal, prim. lattices in conv.
        #   : direct,     conv. lattices in prim.
        # Cij*Di'j=delta_i'i
        #
        # ai = Cij*Aj,  bi = Dij*Bj
        # Aj = ai*Dij,  Bj = bi*Cij
        C = PrimInConv(self.gid)
        D = np.linalg.inv(C).T
        #
        self.kvecP = np.dot(C, self.kvecC)
        self.tauP = [np.dot(D, T) for T in self.tauC]
        rotP = [np.dot(np.dot(D, R), C.T) for R in self.rotC]
        self.rotP = [np.array(np.round(R), dtype=int) for R in rotP]
        if any([np.max(abs(self.rotP[i] - rotP[i])) > 1e-4 for i in range(len(rotP))]):
            raise TypeError('Non-integer rotation matrix !!!')
        #

    def shift(self, shift=np.zeros(3), coord='conv'):
        # g*(r-r0) + t = r'-r0
        # => g*r + (r0-g*r0+t) = r'

        coord = 'conv'
        sgdir = 'sg'+str(self.gid)
        if (sgdir == 'sg48' or sgdir == 'sg86' or sgdir == 'sg126' or sgdir == 'sg210'):
            shift = np.array([-0.25, -0.25, -0.25])
        elif (sgdir == 'sg70'):
            shift = np.array([0.125, 0.125, 0.125])
            coord = 'prim'
        elif (sgdir == 'sg201'):
            shift = np.array([-0.375, -0.375, -0.375])
        elif (sgdir == 'sg85' or sgdir == 'sg129' or sgdir == 'sg130'):
            shift = np.array([0.25, -0.25, 0.0])
        elif (sgdir == 'sg50' or sgdir == 'sg125' or sgdir == 'sg59'):
            shift = np.array([-0.25, -0.25, 0.0])
        elif (sgdir == 'sg133' or sgdir == 'sg137' or sgdir == 'sg134'):
            shift = np.array([0.25, -0.25, 0.25])
        elif sgdir == 'sg68':
            shift = np.array([0.25, -0.25, 0.25])
            coord = 'p'
        elif (sgdir == 'sg88'):
            shift = np.array([0.375, 0.125, 0.25])
            coord = 'p'
        elif (sgdir == 'sg141' or sgdir == 'sg142'):
            shift = np.array([0.5, 0.25, 0.125])
        elif (sgdir == 'sg138'):
            shift = np.array([0.25, -0.25, -0.25])
        elif (sgdir == 'sg222' or sgdir == 'sg224'):
            shift = np.array([0.25, 0.25, 0.25])
        elif (sgdir == 'sg227'):
            shift = np.array([0.125, 0.125, 0.125])
        else:
            shift = np.array([0, 0, 0])

        if np.max(np.abs(shift)) > 1e-4:
            C = PrimInConv(self.gid)
            D = np.linalg.inv(C).T
            #
            if coord[0] == 'c':
                self.shiftC = np.array(shift)
                self.shiftP = np.dot(D, self.shiftC)
            else:
                self.shiftP = np.array(shift)
                self.shiftC = np.dot(C.T, self.shiftP)
            #
            # The tau
            dtauC = [self.shiftC - np.dot(self.rotC[i], self.shiftC) for i, t in enumerate(self.tauC)]
            self.tauC = [t + dtauC[i] for i, t in enumerate(self.tauC)]
            self.tauP = [np.dot(D, T) for T in self.tauC]
            #
            # The irreps
            # for irrep in self.irreps:
            #    irrep.matrices  = [ irrep.matrices[i]*np.exp(2j*np.pi*np.dot(dtauC[i],self.kvecC)) \
            #                        for i in range(len(dtauC)) ]
            #    irrep.characters= [ irrep.characters[i]*np.exp(2j*np.pi*np.dot(dtauC[i],self.kvecC)) \
            #                        for i in range(len(dtauC)) ]

    def SU2(self):  # SU2 matrix calculated from SO3 matrix, may be different from su2s read from Bilbao
        sigma0 = np.array([[1, 0], [0, 1]], dtype=complex)  # Pauli matrix
        sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
        sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)

        A = np.array([self.basisP[0], self.basisP[1], self.basisP[2]], dtype=float).T
        B = np.linalg.inv(A)
        # print('A',A)
        # print('B',B)
        for iop in range(len(self.rotP)):
            rotCart = np.dot(np.dot(A, self.rotP[iop]), B)  # SO3 matrix in Cartesian coordinates
            angle, axis = get_rotation(rotCart)
            su2 = cos(angle / 2) * sigma0 - 1j * sin(angle / 2) \
                  * (axis[0] * sigma1 + axis[1] * sigma2 + axis[2] * sigma3)
            self.su2c.append(su2)
        # print('klabel', self.klabel)
        # print('SO3',self.rotP[iop])
        # print('SO3c',rotCart)
        # print('SU2',su2)
        # print('angle',angle*360/2/pi)
        # print('axis',axis,'\n')

    def align_su2_GM(self):
        n = round(len(self.rotP) / 2)
        su2 = [su2s.copy() for su2s in self.su2s]
        su2c = self.su2c
        multi_table1 = 100 * np.ones((n, n))  # multiplication table of Bilbao su2s, times 100 to distinguish error
        multi_table2 = 100 * np.ones((n, n))  # multiplication table of calculated su2c

        for i in range(n):
            for j in range(n):
                tmp1 = np.dot(su2[i], su2[j])
                tmp2 = np.dot(su2c[i], su2c[j])

                for k, op3 in enumerate(su2):
                    if np.max(np.array(abs(tmp1 - op3))) < 1e-4:
                        multi_table1[i, j] = k
                        break
                    elif np.max(np.array(abs(tmp1 + op3))) < 1e-4:
                        multi_table1[i, j] = -k
                        break
                assert multi_table1[i, j] != 100, ('can\'t find multiplication table for Bilbao SU2!', \
                                                   i, su2[i], j, su2[j], 'i*j', tmp1)

                for k, op3 in enumerate(su2c):
                    if np.max(np.array(abs(tmp2 - op3))) < 1e-4:
                        multi_table2[i, j] = k
                        break
                    elif np.max(np.array(abs(tmp2 + op3))) < 1e-4:
                        multi_table2[i, j] = -k
                        break

        minus = np.zeros(n)
        for i in range(n):
            for j in range(n):
                if multi_table2[i, j] == -multi_table1[i, j] and multi_table1[i, j] != 0:
                    minus[i] += 1

        max = np.max(minus)
        if max > 0:
            for i in range(n):
                if minus[i] == max:
                    self.index.append(i)

        # check whether or not  the multiplication tables match
        if len(self.index) > 0:
            for i in self.index:
                su2[i] *= -1
            for i in range(n):
                for j in range(n):
                    tmp1 = np.dot(su2[i], su2[j])
                    for k, op3 in enumerate(su2):
                        if np.max(np.array(abs(tmp1 - op3))) < 1e-4:
                            multi_table1[i, j] = k
                            break
                        elif np.max(np.array(abs(tmp1 + op3))) < 1e-4:
                            multi_table1[i, j] = -k
                            break
                    assert multi_table1[i, j] != 100, ('can\'t find multiplication table for Bilbao SU2!', \
                                                       i, su2[i], j, su2[j], 'i*j', tmp1)
        assert np.max(np.abs(multi_table1 - multi_table2)) < 1e-4, 'can\'t find the right SU2 matrix to time -1!'


    def __str__(self):
        strg = 'gid:'+str(self.gid)+'  k label: ' + self.klabel + '\n'
        strg += '  k = %6.3f %6.3f %6.3f (conv)  %6.3f %6.3f %6.3f (prim) \n' % (tuple(self.kvecC) + tuple(self.kvecP))
        for ii in range(round(len(self.rotP) / 2)):
            strg += 'op %3d \n' % ii
            strg += '    %4d %4d %4d\n' % tuple(self.rotC[ii][0])
            strg += '    %4d %4d %4d\n' % tuple(self.rotC[ii][1])
            strg += '    %4d %4d %4d\n' % tuple(self.rotC[ii][2])
            strg += '     %10.7f %10.7f %10.7f\n' % tuple(self.tauC[ii])
            strg += '    (%10.7f %10.7f)  (%10.7f %10.7f) \n' % (self.su2s[ii][0][0].real, self.su2s[ii][0][0].imag, \
                                                                 self.su2s[ii][0][1].real, self.su2s[ii][0][1].imag)
            strg += '    (%10.7f %10.7f)  (%10.7f %10.7f) \n' % (self.su2s[ii][1][0].real, self.su2s[ii][1][0].imag, \
                                                                 self.su2s[ii][1][1].real, self.su2s[ii][1][1].imag)
        return strg

    def find_generators(self):
        # given a group G of SO(3) matrices, find generators of it, and return the index list.
        if len(self.rotC)==2: #only identity
            self.generators = [0]
        else:
            element_set = [np.eye(3)]
            generator_set = []
            for ig in range(round(len(self.rotC)/2)):
                g = self.rotC[ig]
                if sum([np.linalg.norm(g-tmpg)<1e-8 for tmpg in element_set]) == 0:
                    generator_set.append(ig)

                for nfold in [1,2,3,4,6]:
                    gn = np.linalg.matrix_power(g,nfold) 
                    if np.linalg.norm(gn - np.eye(3)) > 1e-8 and \
                            sum([np.linalg.norm(gn - tmpg)<1e-8 for tmpg in element_set]) == 0:
                        element_set.append(gn)
                    else:
                        break                   

                for ie1, e1 in enumerate(element_set):
                    for ie2, e2 in enumerate(element_set):
                        tmp = np.dot(e1, e2)
                        if np.linalg.norm(tmp - np.eye(3)) > 1e-8 and \
                                sum([np.linalg.norm(tmp - tmpg)<1e-8 for tmpg in element_set]) == 0:
                            element_set.append(tmp)
                for ie1, e1 in enumerate(element_set):
                    for ie2, e2 in enumerate(element_set):
                        for ie3, e3 in enumerate(element_set):
                            tmp = np.dot(np.dot(e1, e2),e3)
                            if np.linalg.norm(tmp - np.eye(3)) > 1e-8 and \
                                    sum([np.linalg.norm(tmp - tmpg)<1e-8 for tmpg in element_set]) == 0:
                                element_set.append(tmp)
                        
            self.generators = generator_set
       #print(self.generators)



class Irrep:
    def __init__(self):
        self.label = ''  # label for the representation. 'd' as postfix indicates double-valued irreps
        self.dim = 1  # dimension of the representation
        self.matrices = []  # representation matrices, ordered as the operations of the belonging group
        self.characters = []  # characters as a list, ordered the same as matrices

class coIrrep:
    def __init__(self):
        self.label = ''  # label for the representation. 'd' as postfix indicates double-valued irreps
        self.dim = 1  # dimension of the representation
        self.matrices = []  # representation matrices, ordered as the operations of the belonging group
        self.characters = []  # characters as a list, ordered the same as matrices
        self.TRmatrix = []

class antiunitary_op:
    def __init__(self):
        self.exist = False
        self.label = None
        self.rotC = np.eye(3)
        self.tauC = np.zeros(3)
        self.rotP = np.eye(3)
        self.tauP = np.zeros(3)
        self.su2  = np.eye(2)
        self.square_index = 0 # the index of (Tg)^2 in the little group
        self.square_rotC = np.eye(3)
        self.square_tauC = np.zeros(3)
        self.square_su2  = np.eye(2)

def get_rotation(R):
    det = np.linalg.det(R)
    tmpR = det * R
    arg = (np.trace(tmpR) - 1) / 2
    if arg > 1:
        arg = 1
    elif arg < -1:
        arg = -1
    angle = acos(arg)
    axis = np.zeros((3, 1))
    if abs(abs(angle) - pi) < 1e-4:
        for i in range(3):
            axis[i] = 1
            axis = axis + np.dot(tmpR, axis)
            if max(abs(axis)) > 1e-1:
                break
        assert max(abs(axis)) > 1e-1, 'can\'t find axis'
        axis = axis / np.linalg.norm(axis)
    elif abs(angle) > 1e-3:
        # standard case, see Altmann's book
        axis[0] = tmpR[2, 1] - tmpR[1, 2]
        axis[1] = tmpR[0, 2] - tmpR[2, 0]
        axis[2] = tmpR[1, 0] - tmpR[0, 1]
        axis = axis / sin(angle) / 2
    elif abs(angle) < 1e-4:
        axis[0] = 1

    return angle, axis


def loadIR(gid, fname=None, test=False, shift=[0, 0, 0], coord='c'):
    if not fname:
#       fname = Library + '/IR_DSG.dat'
        fname = '/share/home/pengsy/work/atom2020/pg' + '/IR_DSG.dat'

   #if dfttype=='abinit':
   #    from . import readabinit
   #    pra, prb, prc = readabinit.prvec()
   #else:
   #    from . import readvasp
   #    pra, prb, prc = readvasp.prvec()
    pra, prb, prc = [1,0,0], [0,1,0], [0,0,1]

    lgrps = []
    with open(fname) as file:
        # find the space group
        while True:
            line = file.readline()
            print(line, end='') if test else 0
            if line.startswith('gid ' + str(gid)):
                print('****') if test else 0
                break
            elif not line:
                raise ValueError('I can\'t find ``gid ' + str(gid) + '\" in IR_DSG.dat')
        # Read data
        #
        # k points ========================================================
        #
        line = file.readline()
        print(line, end='') if test else 0
        assert line.startswith('nk')
        nk = int(line.split()[1])
        for i in range(nk):
            grp = LittleGroup()  ####
            grp.basisP = [pra, prb, prc]  # primitive basis
            line = file.readline()
            print(line, end='') if test else 0
            assert line.startswith('k_label')
            grp.klabel = line.split()[1]
            grp.gid = gid
            line = file.readline()
            print(line, end='') if test else 0
            assert line.startswith('k_vec')
            s = re.search('\[(.*?)\]', line).groups()[0]
            grp.kvecC = np.array([float(x) for x in s.split()])
            #
            # operations ==================================================
            #
            line = file.readline()
            print(line, end='') if test else 0
            assert line.startswith('nop')
            nop = int(line.split()[1])
            grp.rotC = []
            grp.tauC = []
            for i in range(nop):
                line = file.readline()
                print(line, end='') if test else 0
                assert line.startswith('operation')
                rot, trans, su2 = re.findall('\[(.*?)\]', line)
                rot = np.array([float(x) for x in rot.split()])
                rot.shape = (3, 3)
                grp.rotC.append(rot)
                trans = np.array([float(x) for x in trans.split()])
                grp.tauC.append(trans)
                su2 = np.array([complex(x) for x in su2.split()])
                su2.shape = (2, 2)
                grp.su2s.append(su2)
            # the below line is added by sypeng
#           grp.find_generators()
            #
            # irreps ======================================================
            #
            line = file.readline()
            print(line, end='') if test else 0
            assert line.startswith('nir')
            nir = int(line.split()[1])
            irreps = []
            for j in range(nir):
                ir = Irrep()
                line = file.readline()
                print(line, end='') if test else 0
                assert line.startswith('ir_label')
                ir.label = line.split()[1]
                line = file.readline()
                print(line, end='') if test else 0
                assert line.startswith('ir_dim')
                ir.dim = int(line.split()[1])
                line = file.readline()
                print(line, end='') if test else 0
                assert line.startswith('ir_matrices')
                mats = re.findall('\[(.*?)\]', line)
                mats = [np.array([complex(x) for x in m.split()]) for m in mats]
                for m in mats:
                    m.shape = (ir.dim, ir.dim)
                ir.matrices = mats
                line = file.readline()
                print(line, end='') if test else 0
                assert line.startswith('ir_characters')
                chrcts = re.search('\[(.*?)\]', line).groups()[0]
                ir.characters = np.array([complex(x) for x in chrcts.split()])
                irreps.append(ir)
            #
            # collect
            #
            grp.irreps = irreps
            lgrps.append(grp)
        #
        # Post-processing
        #
 
        for grp in lgrps:
            grp.prim()
            grp.shift(shift, coord)
            grp.SU2()
 
#       igamma = [i for i in range(len(lgrps)) if max(abs(lgrps[i].kvecC)) < 1e-4][0]
#       lgrps[igamma].align_su2_GM()  # find operations which need to time -1 on its SU2 matrix
#
#       GM_index = lgrps[igamma].index  # list of operations need to time -1
#       GM_operation = [rotP for irot, rotP in enumerate(lgrps[igamma].rotP) if irot in GM_index]
#
#       if len(GM_index) > 0:
#           for grp in lgrps:
#               nop_nsoc = round(len(grp.rotP) / 2)  # number of operations in nsoc case
#               for i in range(nop_nsoc):
#                   for R in GM_operation:
#                       if np.max(np.abs(R - grp.rotP[i])) < 1e-6:
#                           grp.index.append(i)
#                           grp.su2c[i] *= -1  # align calculated su2c to bilbao su2
 
    lgrps = read_TR_degeneracy(gid, lgrps)
    return lgrps


def read_TR_degeneracy(gid, lgrps):
    file = h5py.File(Library + '/TR_degeneracy.h5','r')
    for ik, grp in enumerate(lgrps):
        TR_pair_nsoc = np.array(file[str(gid) + 'nsoc' + grp.klabel]) # like [[1,2],[3,3]]
        TR_pair_soc = np.array(file[str(gid) + 'soc' + grp.klabel])   # like [[4,4],[5,6]]

        TR_pair = []
        for ip, pair in enumerate(TR_pair_nsoc):
            TR_pair.append(pair)
        for ip, pair in enumerate(TR_pair_soc):
            TR_pair.append(pair)
            
        grp.degenerate_pair = TR_pair  #like  [[1,2],[3,3],[4,4],[5,6]]
       #print('TRpair',TR_pair)
    
    file.close()
    return lgrps

def latt_home(vec,tol=1e-4):
    vec0=[]
    for x in vec:
        if abs(x-round(x))<tol:
            vec0.append(round(x))
        else:
            vec0.append(x-np.floor(x))
    return np.array(vec0)

def latt_eq(v1,v2,tol=1e-6):
    return np.linalg.norm(v1-v2)<tol



#if __name__ == '__main__':
#    for gid in range(71,72):
#        lgrps = loadIR(gid)
#        print(gid)
#        for grp in lgrps:
#            print(grp.klabel)
#            print(grp.rotC)
#            print(grp.tauC)
#            print('generators:')
#            print(grp.generators)
#            for ir in grp.irreps:
#                print(ir.label)
#                print(ir.matrices)
#                print(ir.characters)




