from __future__ import division, print_function, absolute_import

import numpy as np
import pickle as pk
import os
import sys


# lattvec ================================================================================
class lattvec(np.ndarray) : 
    '''
    The class for lattice position.
        lattvec is a 1D or 2D ndarray, where each line represents a vector.
    Class attributes:
        tol, the tolerance for two lattvec to be equal with each other.
    '''
    tol=1e-4
    def __new__(cls, pos):
        array = np.asarray(pos,dtype=np.float)
        if array.shape[-1] != 3 or len(array.shape)>2 :
            raise TypeError('wrong shape of lattvec !!!')
        return array.view(cls)
        #
    def norm(self):
        'Calculate the norm of the vector'
        return np.linalg.norm(self)
        #
    def __str__(self):
        return '[%8.4f,%8.4f,%8.4f]'%(self[0],self[1],self[2])
    def __repr__(self):
        return self.__str__()
        #
    def zero(self) :
        'Whether self is a zero vector'
        return float(np.max(np.abs(self)))<lattvec.tol
        #
    def __eq__(self, other) :
        'Overload =='
        if isinstance(other, np.ndarray) :
            return (self-other).zero()
        else :
            return False
        #
    def __ne__(self, other) :
        return not(self==other)
        #
    '''def __contains__(self, vv) :
        '
        Overload "in" and "not in"
        '
        # self shoud be lattvec with shape (n,3), i.e. a set of vectors
        if len(self.shape)==2 and len(vv.shape)==1 :
            return bool( np.min(np.max(np.abs(self-vv),1))<lattvec.tol )
        else :
            return False
    '''
        #
    def home(self, AA=np.array([[1,0,0],[0,1,0],[0,0,1]]), \
                   BB=np.array([[1,0,0],[0,1,0],[0,0,1]]) ) :
        'Move the vector to home cell'
        #
        direc = np.dot(self,BB.T)
        if len(self.shape)==1 :
            for ii, dd in enumerate(direc) :
                if abs(dd-round(dd))<lattvec.tol :
                    direc[ii] = 0.0
                else :
                    direc[ii] = dd-np.floor(dd)
        elif len(self.shape)==2 :
            print(direc)
            for ii, dd in enumerate(direc) :
                for jj, dd_ in enumerate(dd) :
                    if abs(dd_-round(dd_))<lattvec.tol :
                        direc[ii][jj] = 0.0
                    else :
                        direc[ii][jj] = dd_-np.floor(dd_)
        return lattvec(np.dot(direc,AA))
        #
    def order(self, AA=np.array([[1,0,0],[0,1,0],[0,0,1]]), \
                    BB=np.array([[1,0,0],[0,1,0],[0,0,1]]), maxorder=10) :
        'Order of a fractional vector'
        find = False
        direc = np.dot(self,BB.T)
        for ii in range(1,maxorder+1):
            if (ii*direc).home().zero():
                find = True
                break
        if find:
            return ii
# END of lattvec ========================================================================

# symop =================================================================================
class symop(object):
    '''
    The class for symmetry operation
    Attributes:
        mat, the operation matrix (in cart)
        tau, the translation part (in cart)
        det, the determination of mat
        alph, the rotation angle of det*mat
        nfold, the order of rotation
        axis, the rotation axis of det*mat
        AA,  AA[i] is the i-th lattice vector
        BB,  BB[i] is the i-th reciprocal lattice vector
    Class attributes:
        tol, the tolerance
    '''
    tol = 1.0e-4
    # __init__ ==========================================================================
    def __init__(self, mat, tau, AA=lattvec([[1,0,0],[0,1,0],[0,0,1]]), \
                                 BB=lattvec([[1,0,0],[0,1,0],[0,0,1]]), AAC=None, BBC=None ):
        '''
        Notice: tau will be moved to the home cell'''
        #
        # primitive cell
        self.AA = AA
        self.BB = BB
        #
        AAC = AA if AAC is None else AAC
        BBC = BB if BBC is None else BBC
        self.AAC = AAC
        self.BBC = BBC
        #
        # check cell
        if np.max(np.abs(np.dot(BB.T,AA)-np.eye(3)))>1e-4  or \
           np.max(np.abs(np.dot(BBC.T,AAC)-np.eye(3)))>1e-4 :
            raise ValueError('wrong reciprocal lattice !!!')
        #
        # transform matrix and translation vector
        self.mat = mat
        self.tau = tau.home(AA=AA, BB=BB)
        self.tauC= np.dot(self.BBC,self.tau)
        self.matP= np.dot(self.BB,np.dot(self.mat,self.AA.T))
        self.tauP= np.dot(self.BB,self.tau)
        assert np.linalg.norm(self.matP-np.round(self.matP))<1e-4, \
        'rotation matrix in primitive coordinate are not integers \n'+str(self.matP)
        self.matP=np.round(self.matP)
        #
        # determinant and rotational part
        self.det = np.linalg.det( self.mat )
        if abs(abs(self.det) - 1.0) > symop.tol :
            raise ValueError('wrong determinant !!!')
        self.det = np.sign(self.det)
        rot = self.det*self.mat
        #
        # rotation angle
        self.alph = np.arccos( max(min(0.5*(rot.trace()-1),1),-1) )
        if abs(self.alph)<symop.tol :
            self.nfold = 1
        else :
            self.nfold = round(2*np.pi/self.alph)
        #
        # rotation axis
        if abs(self.alph - np.pi)<symop.tol :
            ww, vv = np.linalg.eig(rot)
            for ii, wwi in enumerate(ww):
                if abs(wwi-1.0)<symop.tol:
                    break
            if abs(wwi-1.0)>=symop.tol:
                raise ValueError('can not find axis !!!')
            self.axis = lattvec( vv[:,ii].real )
        elif abs(self.alph)>symop.tol:
            vv = np.array( [ rot[2,1]-rot[1,2], rot[0,2]-rot[2,0], rot[1,0]-rot[0,1] ] )
            self.axis = lattvec( vv/np.sin(self.alph)/2.0 )
        else:
            self.axis = lattvec( [0.0, 0.0, 1.0] )
        #
        # Indices in conventional cell  ------------------------------------------
        #
        #  Mirror & glide
        if self.nfold==2 and self.det<0 :
            #
            # the minimal reciprocal lattice vector along axis
            found = False
            direc_= np.dot(self.AAC, self.axis)
            direc_= direc_/max(abs(direc_))
            for ii in range(1,5) :
                direc = ii*direc_
                if all([ abs(dd-round(dd))<1e-4 for dd in direc ]) :
                    found = True
                    break
            if not found :
                raise ValueError('Can not find lattice in the axis direction !!!')
            #
            self.indice = [ round(a) for a in direc ]
        #
        #   Rotation, screw, rotoreflection
        elif self.nfold>1 :
            #
            # the minimal lattice vector along axis
            found = False
            direc_= np.dot(self.BBC, self.axis)
            direc_= direc_/max(abs(direc_))
            for ii in range(1,5) :
                direc = ii*direc_
                if all([ abs(dd-round(dd))<1e-4 for dd in direc ]) :
                    found = True
                    break
            if not found :
                raise ValueError('Can not find lattice in the axis direction !!!')
            #
            self.indice = [ round(a) for a in direc ]
            #
            if self.nfold==2 and sum(self.indice)<0 :
                self.axis   = -self.axis
                self.indice = [ -a for a in self.indice ]
            #
        else :
            self.indice = [ 0,0,0 ]

    # END of __init__ ===============================================================
    #
    def __mul__(self,other) :
        '''
        Overload *
            if other is a symop, then return the group multiplication
            if other is a lattvec, then return the lattvec after rotation
        '''
        if ( isinstance(other,symop) ) :
            errA = np.max(np.abs(self.AA-other.AA))
            errB = np.max(np.abs(self.BB-other.BB))
            if (errA<symop.tol and errB<symop.tol ) :
                mat = np.dot(self.mat, other.mat)
                tau = self.tau + np.dot(self.mat, other.tau)
                return symop(mat, tau, self.AA, self.BB, self.AAC, self.BBC)
            else :
                raise ValueError('The lattices of two entrances of group multiplication should be identical !!!')
        elif (isinstance(other,lattvec)) :
            return np.dot(self.mat,other) + self.tau 
        #
    def __eq__(self,other) :
        'Overload =='
        errA = np.max(np.abs(self.AA-other.AA))
        errB = np.max(np.abs(self.BB-other.BB))
        errm = np.max(np.abs(self.mat-other.mat))
        errt = np.max(np.abs(self.tau-other.tau))
        return( errA<1.0e-4 and errB<1.0e-4 and errm<1.0e-4 and errt<1.0e-4 )
        #return( errm<1.0e-4 and errt<1.0e-4 )
        #
    def __ne__(self,other) :
        'Overload !='
        return( not self==other )
        
# END of symop ==========================================================================

        
# Spacegroup ============================================================================
class SpaceGroup(object):
    '''
    Class for space group
    Attributes:
        name, the name of the space group
        nop, the number of symops
        op[nop], the symops
        AA, AA[i] is the i-th lattice vector
        BB, BB[i] is the i-th reciprocal lattice vector
        id0, the index of identity
        mtb[48][48], the multiplication table
        inv[48], the inversion table
        ncls, the number of classes
        cls[ncls], the operation indexes in each class
    '''
    # __init__ ==========================================================================
    def __init__(self, gen, name='', spgid=0) :
        #
        self.name = name
        self.spgid= int(spgid)
        self.nop = len(gen)
        assert self.nop>0
        assert isinstance(gen[0],symop)
        #
        # get all the symmetry operations -------------------------------------
        #
        self.op = gen
        self.AA = self.op[0].AA
        self.BB = self.op[0].BB
        self.AAC= self.op[0].AAC
        self.BBC= self.op[0].BBC
        #
        for op1 in self.op :
            for op2 in self.op :
                op3 = op1*op2
                if all([op3 != op0 for op0 in self.op]) :
                    self.op.append(op3)
        self.nop = len(self.op)
        #
        # get the multiplication table ----------------------------------------
        #
        self.mtb =[ [None]*self.nop for ii in range(self.nop) ]
        for ii in range(self.nop):
            for jj in range(self.nop):
                self.mtb[ii][jj] = list(filter(lambda k: self.op[k]==self.op[ii]*self.op[jj],\
                                               range(self.nop)))[0]
        #
        # get the identity and inversion ------------------------------------------------
        #
        self.inv = [None]*self.nop
        for kk, op0 in enumerate(self.op) :
            if op0.det>0 and abs(op0.alph)<1.0e-3 :
                self.id0 = kk
                #print(self.id0)
                break
        for ii in range(self.nop) :
            #self.inv[ii] = list(filter(lambda k: self.mtb[ii][k]==self.id0, range(self.nop)))[0]
            self.inv[ii] = [k for k in range(self.nop) if self.mtb[ii][k]==self.id0][0]
        #
        # get the classes ---------------------------------------------------------------
        #
        self.cls   = []
        self.ncls  = 0
        unassign = [True]*self.nop
        for ii in range(self.nop) :
            clstmp   = []
            for jj in range(self.nop) :
                kk = self.mtb[jj][ii]
                kk = self.mtb[kk][self.inv[jj]]
                if unassign[kk] :
                    clstmp.append( kk )
                    unassign[kk] = False
            if len(clstmp)>0 :
                self.cls.append( clstmp )
                self.ncls += 1
    # END of __init__ ===================================================================
    #
    # get_PhysClass =====================================================================
    def get_PhysClass(self):
        '''
        Get the physical classes. An operation is physically equivalent with its inversion
        Output:
            pclassrep, the representative of each physical class. We choose the operation 
                with max indice[2]
            pclassidx, the indexes of physical classes
        '''
        pclassrep = []
        pclassidx = []
        phys = [True]*self.ncls
        #
        # Loop for class
        for icls, cls in enumerate(self.cls) :
            #
            # choose i0 as the one with maximal indice[2]
            idz= [ self.op[i].indice[2] for i in cls ]
            i0 = cls[idz.index(max(idz))] 
            #
            if i0==self.id0 :
                phys[icls] = False
            #
            # Enumerate all the other classes
            for jcls, cls_ in enumerate(self.cls) :
                if jcls != icls and phys[jcls] :
                    #
                    # Check whether icls is physical
                    for kk in cls_ :
                        
                        if i0==self.inv[kk] :
                            if self.op[i0].indice[2]<=self.op[kk].indice[2] :
                                phys[icls] = False
                            break
                if not phys[icls] :
                    break
            #
            # Basic class
            if phys[icls] :
                pclassrep.append(self.op[i0])
                pclassidx.append(icls)
        return pclassrep, pclassidx
    # END of get_PhysClass ==============================================================
    #
    # subgroups =========================================================================
    def get_subgroup(self,gen):
        'Get the subgroup from a set of generators'
        sub = gen
        for ii in sub:
            for jj in sub:
                kk = self.mtb[ii][jj]
                if kk not in sub:
                    sub.append(kk)
        sub.sort()
        return sub
        #
    def get_allsubgroups(self):
        'Get all the subgroups'
        #
        # firstly get the cyclic subgroups
        subset = []
        for iop in range(self.nop) :
            sub = self.get_subgroup([iop])
            #
            if sub not in subset :
                subset.append(sub)
        #
        # secondly get all the subgroups by combining existed subgroups
        for sub1 in subset:
            for sub2 in subset:
                sub3 = self.get_subgroup( list(set(sub1 + sub2)) )
                if sub3 not in subset :
                    subset.append(sub3)
        return subset
        #
    def get_coset(self,sub,cotype='L'):
        #if set(sub) != set(self.get_subgroup(sub)) :
        #    raise ValueError('The input sub is not a sub group !!!')
        set0 = sorted(sub, key=lambda x : -1 if x==self.id0 else x)
        coset= [set0]
        for iop in range(self.nop) :
            if cotype[0]=='L' :
                set_ = [ self.mtb[iop][jop] for jop in sub ]
            else :
                set_ = [ self.mtb[jop][iop] for jop in sub ]
            set_ = sorted(set_, key=lambda x : -1 if x==self.id0 else x)
            if set_ not in coset:
                coset.append(set_)
        return coset

    # END of subgroups ==================================================================
    #
    def str_lattice(self) :
        string = 'Lattice :\n'
        string+= '  conv :\n'
        string+= '  A1 = %6.3f  %6.3f  %6.3f ' % tuple(self.AAC[0])
        string+= '  B1 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BBC[0])
        string+= '  A2 = %6.3f  %6.3f  %6.3f ' % tuple(self.AAC[1])
        string+= '  B2 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BBC[1])
        string+= '  A3 = %6.3f  %6.3f  %6.3f ' % tuple(self.AAC[2])
        string+= '  B3 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BBC[2])
        string+= '  prim :\n'
        string+= '  a1 = %6.3f  %6.3f  %6.3f ' % tuple(self.AA[0])
        string+= '  b1 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BB[0])
        string+= '  a2 = %6.3f  %6.3f  %6.3f ' % tuple(self.AA[1])
        string+= '  b2 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BB[1])
        string+= '  a3 = %6.3f  %6.3f  %6.3f ' % tuple(self.AA[2])
        string+= '  b3 = %6.3f  %6.3f  %6.3f \n' % tuple(self.BB[2])
        return string
        #
    def str_symop(self):
        string = 'Operations :\n'
        string+= '    op det   angle  indice (conv)     tau (conv)\n'
        for ii, op in enumerate(self.op):
            string+= '  %4d  %2d  %6.2f  ( %3d %3d %3d )   ( %7.4f %7.4f %7.4f )\n'  \
             %(ii, op.det, op.alph/np.pi*180, op.indice[0], op.indice[1], op.indice[2], \
               op.tauC[0], op.tauC[1], op.tauC[2] )
        string+= '\n'
        string+= '    * Indices of identity and inversion are set to zero \n'
        string+= '      Indices of mirror & glide operations are defined in conventional reciprocal lattice \n'
        string+= '      Indices of other operations are defined in conventional lattice \n'
        string+= '      All tau (include the glides\') are defined in conventional lattice \n'
        string+= '\n'
        return string
        #
    def str_class(self):
        string = 'Classes :\n'
        string+= '  class  op\n'
        for ii, cls in enumerate(self.cls) :
            string+= ('     %2d  ' % ii )+cls.__str__()+'\n'
        return string
        #
    def __str__(self) :
        return ('SPG: %3d' % self.spgid) + '\n' + 'Name: '+self.name + '\n' \
            + self.str_lattice() + self.str_symop()
# End of Space group ====================================================================

#
# Get the space group from space group number ===========================================
#
def get_SPGfromID(spgid):
    if (spgid==1):
        aa = symop( np.array([ [1,0,0], [0,1,0], [0,0,1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa], name='P1', spgid=spgid)
    if (spgid==2):
        aa = symop(-np.array([ [1,0,0], [0,1,0], [0,0,1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa], name='P-1', spgid=spgid)
    if (spgid==3):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa], name='P2', spgid=spgid)
    if (spgid==4):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]) )
        return SpaceGroup([aa], name='P2_1', spgid=spgid)
    if (spgid==5):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T 
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC )
        return SpaceGroup([aa], name='C2', spgid=spgid)
    if (spgid==6):
        aa = symop( np.array([ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa], name='Pm', spgid=spgid)
    if (spgid==7):
        aa = symop( np.array([ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]) )
        return SpaceGroup([aa], name='Pc', spgid=spgid)
    if (spgid==8):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T 
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa], name='Cm', spgid=spgid)
    if (spgid==9):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T 
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC )
        return SpaceGroup([aa], name='Cc', spgid=spgid)
    if (spgid==10):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='P2/m', spgid=spgid)
    if (spgid==11):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]) )
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='P2_1/m', spgid=spgid)
    if (spgid==12):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T 
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='C2/m', spgid=spgid)
    if (spgid==13):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]) )
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='P2/c', spgid=spgid)
    if (spgid==14):
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]) )
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='P2_1/c', spgid=spgid)
    if (spgid==15):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='C2/c', spgid=spgid)
    if (spgid==16):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='P222', spgid=spgid)
    if (spgid==17):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]) )
        return SpaceGroup([aa,bb], name='P222_1', spgid=spgid)
    if (spgid==18):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]) )
        return SpaceGroup([aa,bb], name='P2_1 2_1 2', spgid=spgid)
    if (spgid==19):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]) )
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]) )
        return SpaceGroup([aa,bb], name='P2_1 2_1 2_1', spgid=spgid)
    if (spgid==20):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='C222_1', spgid=spgid)
    if (spgid==21):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        return SpaceGroup([aa,bb], name='C222', spgid=spgid)
    if (spgid==22):
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB,AAC=AAC, BBC=BBC )
        return SpaceGroup([aa,bb], name='F222', spgid=spgid)
    if (spgid==23):
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='I222', spgid=spgid)
    if (spgid==24):
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ]).T
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='I2_1 2_1 2_1', spgid=spgid)
    if (spgid==25):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa,bb], name='Pmm2', spgid=spgid)
    if (spgid==26):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        return SpaceGroup([aa,bb], name='Pmc2_1', spgid=spgid)
    if (spgid==27):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        return SpaceGroup([aa,bb], name='Pcc2', spgid=spgid) 
    if (spgid==28):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]) )
        return SpaceGroup([aa,bb], name='Pma2', spgid=spgid) 
    if (spgid==29):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]) )
        return SpaceGroup([aa,bb], name='Pca2_1', spgid=spgid)
    if (spgid==30):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0.5, 0.5]) )
        return SpaceGroup([aa,bb], name='Pnc2', spgid=spgid)
    if (spgid==31):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]) )
        return SpaceGroup([aa,bb], name='Pmn2_1', spgid=spgid)
    if (spgid==32):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]) )
        return SpaceGroup([aa,bb], name='Pba2', spgid=spgid)
    if (spgid==33):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]) )
        return SpaceGroup([aa,bb], name='Pna2_1', spgid=spgid)
    if (spgid==34):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]) )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0.5]) )
        return SpaceGroup([aa,bb], name='Pnn2', spgid=spgid)
    if (spgid==35):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB,AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Cmm2', spgid=spgid)
    if (spgid==36):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB,AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Cmc2_1', spgid=spgid)
    if (spgid==37):
        AA = lattvec([ [0.5,0.5,0], [-0.5,0.5,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Ccc2', spgid=spgid)
    if (spgid==38):
        AA = lattvec([ [1,0,0], [0,0.5,0.5], [0,-0.5,0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        return SpaceGroup([aa,bb], name='Amm2', spgid=spgid)
    if (spgid==39):
        AA = lattvec([ [1,0,0], [0,0.5,0.5], [0,-0.5,0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Aem2', spgid=spgid)
    if (spgid==40):
        AA = lattvec([ [1,0,0], [0,0.5,0.5], [0,-0.5,0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Ama2', spgid=spgid)
    if (spgid==41) :
        AA = lattvec([ [1,0,0], [0,0.5,0.5], [0,-0.5,0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC )
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Aea2', spgid=spgid)
    if (spgid==42) :
        AA = lattvec([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Fmm2', spgid=spgid)
    if (spgid==43) :
        AA = lattvec([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.25, 0.25, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb], name='Fdd2', spgid=spgid)
    # 44: Imm2
    if (spgid==44) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='Imm2', spgid=spgid)
    # 45: Iba2
    if (spgid==45) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='Iba2', spgid=spgid)
    # 46: Ima2
    if (spgid==46) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [ 1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='Ima2', spgid=spgid)
    # 47: Pmmm
    if (spgid==47) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pmmm', spgid=spgid)
    # 48: Pnnn
    if (spgid==48) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pnnn', spgid=spgid)
    # 49: Pccm
    if (spgid==49) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pccm', spgid=spgid)
    # 50: Pban
    if (spgid==50) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pban', spgid=spgid)
    # 51: Pmma
    if (spgid==51) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pmma', spgid=spgid)
    # 52: Pnna
    if (spgid==52) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pnna', spgid=spgid)
    # 53: Pmna
    if (spgid==53) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pmna', spgid=spgid)
    # 54: Pcca
    if (spgid==54) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pcca', spgid=spgid)
    # 55: Pbam
    if (spgid==55) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pbam', spgid=spgid)
    # 56: Pccn
    if (spgid==56) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pccn', spgid=spgid)
    # 57: Pbcm
    if (spgid==57) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pbcm', spgid=spgid)
    # 58: Pnnm
    if (spgid==58) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pnnm', spgid=spgid)
    # 59: Pmmn
    if (spgid==59) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pmmn', spgid=spgid)
    # 60: Pbcn
    if (spgid==60) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pbcn', spgid=spgid)
    # 61: Pbca
    if (spgid==61) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pbca', spgid=spgid)
    # 62: Pnma
    if (spgid==62) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='Pnma', spgid=spgid)
    # 63: Cmcm
    if (spgid==63) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Cmcm', spgid=spgid)
    # 64: Cmce
    if (spgid==64) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Cmce', spgid=spgid)
    # 65: Cmmm
    if (spgid==65) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Cmmm', spgid=spgid)
    # 66: Cccm
    if (spgid==66) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Cccm', spgid=spgid)
    # 67: Cmme
    if (spgid==67) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Cmme', spgid=spgid)
    # 68: Ccce
    if (spgid==68) :
        AA = lattvec([ [0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Ccce', spgid=spgid)
    # 69: Fmmm
    if (spgid==69) :
        AA = lattvec([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fmmm', spgid=spgid)
    # 70: Fddd
    if (spgid==70) :
        AA = lattvec([ [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.75, 0.75, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.75, 0, 0.75]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fddd', spgid=spgid)
    # 71: Immm
    if (spgid==71) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5,-0.5, 0.5], [0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Immm', spgid=spgid)
    # 72: Ibam
    if (spgid==72) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5,-0.5, 0.5], [0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Ibam', spgid=spgid)
    # 73: Ibca
    if (spgid==73) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5,-0.5, 0.5], [0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Ibca', spgid=spgid)
    # 74: Imma
    if (spgid==74) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [0.5,-0.5, 0.5], [0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0, 1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array([ [-1,0,0], [0, 1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Imma', spgid=spgid)
    # 75: P4
    if (spgid==75) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P4', spgid=spgid)
    # 76: P4_1
    if (spgid==76) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0, 0.25]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P4_1', spgid=spgid)
    # 77: P4_2
    if (spgid==77) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P4_2', spgid=spgid)
    # 78: P4_3
    if (spgid==78) :
        AA = lattvec([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0, 0.75]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P4_3', spgid=spgid)
    # 79: I4
    if (spgid==79) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB,AAC=AAC, BBC=BBC)
        return SpaceGroup([aa],name='I4', spgid=spgid)
    # 80: I4_1
    if (spgid==80) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array([ [0,-1,0], [1,0,0], [0,0, 1] ]), lattvec([0, 0.5, 0.25]), \
                    AA=AA, BB=BB,AAC=AAC,BBC=BBC)
        return SpaceGroup([aa],name='I4_1', spgid=spgid)
    # 81: P-4
    if (spgid==81) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P-4', spgid=spgid)
    # 82: I-4
    if (spgid==82) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB,AAC=AAC,BBC=BBC)
        return SpaceGroup([aa],name='I-4', spgid=spgid)
    # 83: P4/m
    if (spgid==83) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4/m', spgid=spgid)
    # 84: P4_2/m
    if (spgid==84) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2/m', spgid=spgid)
    # 85: P4/n
    if (spgid==85) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4/n', spgid=spgid)
    # 86: P4_2/n
    if (spgid==86) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2/n', spgid=spgid)
    # 87: I4/m
    if (spgid==87) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        return SpaceGroup([aa,bb],name='I4/m', spgid=spgid)
    # 88: I4_1/a
    if (spgid==88) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.75, 0.25, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC,BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4_1/a', spgid=spgid)
    # 89: P422
    if (spgid==89) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P422', spgid=spgid)
    # 90: P4 2_1 2
    if (spgid==90) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4 2_1 2', spgid=spgid)
    # 91: P4_1 22
    if (spgid==91) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.25]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_1 22', spgid=spgid)
    # 92: P4_1 2_1 2
    if (spgid==92) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.25]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.25]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_1 2_1 2', spgid=spgid)
    # 93: P4_2 22
    if (spgid==93) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 22', spgid=spgid)
    # 94: P4_2 2_1 2
    if (spgid==94) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 2_1 2', spgid=spgid)
    # 95: P4_3 2 2
    if (spgid==95) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.75]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_3 2 2', spgid=spgid)
    # 96: P4_3 2_1 2
    if (spgid==96) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.75]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.75]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_3 2_1 2', spgid=spgid)
    # 97: I422 
    if (spgid==97) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB,AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I422', spgid=spgid)
    # 98: I4_1 22 
    if (spgid==98) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.5, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.75]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4_1 22', spgid=spgid)
    # 99: P4mm
    if (spgid==99) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4mm', spgid=spgid)
    # 100: P4bm
    if (spgid==100) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4bm', spgid=spgid)
    # 101: P4_2 cm
    if (spgid==101) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 cm', spgid=spgid)
    # 102: P4_2 nm
    if (spgid==102) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 cm', spgid=spgid)
    # 103: P4cc
    if (spgid==103) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4cc', spgid=spgid)
    # 104: P4nc
    if (spgid==104) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4nc', spgid=spgid)
    # 105: P4_2 mc
    if (spgid==105) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 mc', spgid=spgid)
    # 106: P4_2 bc
    if (spgid==106) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P4_2 bc', spgid=spgid)
    # 107: I4mm
    if (spgid==107) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4mm', spgid=spgid)
    # 108: I4cm
    if (spgid==108) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4cm', spgid=spgid)
    # 109: I4_1 md
    if (spgid==109) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.5, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4_1 md', spgid=spgid)
    # 110: I4_1 cd
    if (spgid==110) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.5, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4_1 cd', spgid=spgid)
    # 111: P-42m
    if (spgid==111) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-42m', spgid=spgid)
    # 112: P-42c
    if (spgid==112) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-42c', spgid=spgid)
    # 113: P-4 2_1 m
    if (spgid==113) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4 2_1 m', spgid=spgid)
    # 114: P-4 2_1 c
    if (spgid==114) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4 2_1 c', spgid=spgid)
    # 115: P-4m2
    if (spgid==115) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4m2', spgid=spgid)
    # 116: P-4c2
    if (spgid==116) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4c2', spgid=spgid)
    # 117: P-4b2
    if (spgid==117) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4b2', spgid=spgid)
    # 118: P-4n2
    if (spgid==118) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-4n2', spgid=spgid)
    # 119: I-4m2
    if (spgid==119) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-4m2', spgid=spgid)
    # 120: I-4c2
    if (spgid==120) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-4c2', spgid=spgid)
    # 121: I-42m
    if (spgid==121) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-42m', spgid=spgid)
    # 122: I-42d
    if (spgid==122) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,1,0], [-1,0,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.75]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-42d', spgid=spgid)
    # 123: P4/mmm
    if (spgid==123) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/mmm', spgid=spgid)
    # 124: P4/mcc
    if (spgid==124) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/mcc', spgid=spgid)
    # 125: P4/nbm
    if (spgid==125) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/nbm', spgid=spgid)
    # 126: P4/nnc
    if (spgid==126) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/nnc', spgid=spgid)
    # 127: P4/mbm
    if (spgid==127) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/mbm', spgid=spgid)
    # 128: P4/mnc
    if (spgid==128) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/mnc', spgid=spgid)
    # 129: P4/nmm
    if (spgid==129) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/nmm', spgid=spgid)
    # 130: P4/ncc
    if (spgid==130) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4/ncc', spgid=spgid)
    # 131: P4_2/mmc
    if (spgid==131) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/mmc', spgid=spgid)
    # 132: P4_2/mcm
    if (spgid==132) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/mcm', spgid=spgid)
    # 133: P4_2/nbc
    if (spgid==133) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/nbc', spgid=spgid)
    # 134: P4_2/nnm
    if (spgid==134) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/nnm', spgid=spgid)
    # 135: P4_2/mbc
    if (spgid==135) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/mbc', spgid=spgid)
    # 136: P4_2/mnm
    if (spgid==136) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/mnm', spgid=spgid)
    # 137: P4_2/nmc
    if (spgid==137) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/nmc', spgid=spgid)
    # 138: P4_2/ncm
    if (spgid==138) :
        AA = lattvec([ [1, 0, 0], [ 0, 1, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P4_2/ncm', spgid=spgid)
    # 139: I4/mmm
    if (spgid==139) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='I4/mmm', spgid=spgid)
    # 140: I4/mcm
    if (spgid==140) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='I4/mcm', spgid=spgid)
    # 141: I4_1/amd
    if (spgid==141) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.25, 0.75, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='I4_1/amd', spgid=spgid)
    # 142: I4_1/acd
    if (spgid==142) :
        AA = lattvec([ [-0.5, 0.5, 0.5], [ 0.5,-0.5, 0.5], [ 0.5, 0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.25, 0.75, 0.25]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,1,0], [0,0,-1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='I4_1/acd', spgid=spgid)
    # 143: P3
    if (spgid==143) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P3', spgid=spgid)
    # 144: P3_1
    if (spgid==144) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P3_1', spgid=spgid)
    # 145: P3_2
    if (spgid==145) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P3_2', spgid=spgid)
    # 146: R3
    if (spgid==146) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa],name='R3', spgid=spgid)
    # 147: P-3
    if (spgid==147) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-3', spgid=spgid)
    # 148: R-3
    if (spgid==148) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb],name='R-3', spgid=spgid)
    # 149: P312
    if (spgid==149) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P312', spgid=spgid)
    # 150: P321
    if (spgid==150) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P321', spgid=spgid)
    # 151: P3_1 12 
    if (spgid==151) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3_1 12', spgid=spgid)
    # 152: P3_1 21 
    if (spgid==152) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3_1 21', spgid=spgid)
    # 153: P3_2 12 
    if (spgid==153) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3_2 12', spgid=spgid)
    # 154: P3_2 21 
    if (spgid==154) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3_2 21', spgid=spgid)
    # 155: R32
    if (spgid==155) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb],name='R32', spgid=spgid)
    # 156: P3m1 
    if (spgid==156) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3m1', spgid=spgid)
    # 157: P31m
    if (spgid==157) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P31m', spgid=spgid)
    # 158: P3c1
    if (spgid==158) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/2]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P3c1', spgid=spgid)
    # 159: P31c
    if (spgid==159) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/2]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P31c', spgid=spgid)
    # 160: R3m
    if (spgid==160) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb],name='R3m', spgid=spgid)
    # 161: R3c
    if (spgid==161) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 1/2]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb],name='R3c', spgid=spgid)
    # 162: P-31m
    if (spgid==162) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.dot(np.dot(AA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P-31m', spgid=spgid)
    # 163: P-31c
    if (spgid==163) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [-1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.dot(np.dot(AA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P-31c', spgid=spgid)
    # 164: P-3m1
    if (spgid==164) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.dot(np.dot(AA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P-3m1', spgid=spgid)
    # 165: P-3c1
    if (spgid==165) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.dot(np.dot(AA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P-3c1', spgid=spgid)
    # 166: R-3m
    if (spgid==166) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        cc = symop( np.dot(np.dot(CAA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb,cc],name='R-3m', spgid=spgid)
    # 167: R-3c
    if (spgid==167) :
        CAA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        CBB = np.linalg.inv(CAA).T
        PAA = np.dot(lattvec([ [2/3,1/3,1/3], [-1/3,1/3,1/3], [-1/3,-2/3,1/3] ]),CAA)
        PBB = np.linalg.inv(PAA).T
        aa = symop( np.dot(np.dot(CAA.T,np.array( [ [0,-1,0], [1,-1,0], [0,0,1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        bb = symop( np.dot(np.dot(CAA.T,np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0.5]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        cc = symop( np.dot(np.dot(CAA.T,np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ])),CBB), \
                    np.dot(lattvec([0, 0, 0]),CAA), \
                    AA=PAA, BB=PBB, AAC=CAA, BBC=CBB)
        return SpaceGroup([aa,bb,cc],name='R-3c', spgid=spgid)
    # 168: P6
    if (spgid==168) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6', spgid=spgid)
    # 169: P6_1
    if (spgid==169) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/6]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6_1', spgid=spgid)
    # 170: P6_5
    if (spgid==170) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 5/6]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6_5', spgid=spgid)
    # 171: P6_2
    if (spgid==171) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/6]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6_2', spgid=spgid)
    # 172: P6_4
    if (spgid==172) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 4/6]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6_4', spgid=spgid)
    # 173: P6_3
    if (spgid==173) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T,np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 3/6]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P6_3', spgid=spgid)
    # 174: P-6
    if (spgid==174) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, -np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa],name='P-6', spgid=spgid)
    # 175: P6/m
    if (spgid==175) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6/m', spgid=spgid)
    # 176: P6_3/m
    if (spgid==176) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_3/m', spgid=spgid)
    # 177: P622
    if (spgid==177) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P622', spgid=spgid)
    # 178: P6_1 22
    if (spgid==178) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/6]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_1 22', spgid=spgid)
    # 179: P6_5 22
    if (spgid==179) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 5/6]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_5 22', spgid=spgid)
    # 180: P6_2 22
    if (spgid==180) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/6]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 2/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_2 22', spgid=spgid)
    # 181: P6_4 22
    if (spgid==181) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 4/6]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 1/3]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_4 22', spgid=spgid)
    # 182: P6_3 22
    if (spgid==182) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_3 22', spgid=spgid)
    # 183: P6mm
    if (spgid==183) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6mm', spgid=spgid)
    # 184: P6cc
    if (spgid==184) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6cc', spgid=spgid)
    # 185: P6_3 cm
    if (spgid==185) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_3 cm', spgid=spgid)
    # 186: P6_3 mc
    if (spgid==186) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P6_3 mc', spgid=spgid)
    # 187: P-6m2
    if (spgid==187) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, -np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-6m2', spgid=spgid)
    # 188: P-6c2
    if (spgid==188) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, -np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,-1,0], [-1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-6c2', spgid=spgid)
    # 189: P-62m
    if (spgid==189) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, -np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-62m', spgid=spgid)
    # 190: P-62c
    if (spgid==190) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, -np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb],name='P-62c', spgid=spgid)
    # 191: P6/mmm
    if (spgid==191) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P6/mmm', spgid=spgid)
    # 192: P6/mcc
    if (spgid==192) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB)
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P6/mcc', spgid=spgid)
    # 193: P6_3/mcm
    if (spgid==193) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB )
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P6_3/mcm', spgid=spgid)
    # 194: P6_3/mmc
    if (spgid==194) :
        AA = lattvec([ [0,-1, 0], [np.sqrt(3)/2, 1/2, 0], [ 0, 0, 1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.dot(np.dot(AA.T, np.array( [ [1,-1,0], [1,0,0], [0,0,1] ])),BB), \
                    np.dot(lattvec([0, 0, 0.5]),AA), \
                    AA=AA, BB=BB )
        bb = symop( np.dot(np.dot(AA.T, np.array( [ [0,1,0], [1,0,0], [0,0,-1] ])),BB), \
                    np.dot(lattvec([0, 0, 0]),AA), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB)
        return SpaceGroup([aa,bb,cc],name='P6_3/mmc', spgid=spgid)
    # 195: P23
    if (spgid==195) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P23', spgid=spgid)
    # 196: F23
    if (spgid==196) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='F23', spgid=spgid)
    # 197: I23
    if (spgid==197) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB, AAC=AAC, BBC=BBC )
        return SpaceGroup([aa,bb],name='I23', spgid=spgid)
    # 198: P2_1 3
    if (spgid==198) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P2_1 3', spgid=spgid)
    # 199: I2_1 3
    if (spgid==199) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I2_1 3', spgid=spgid)
    # 200: Pm-3
    if (spgid==200) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pm-3', spgid=spgid)
    # 201: Pn-3
    if (spgid==201) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0.5, 0]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pn-3', spgid=spgid)
    # 202: Fm-3
    if (spgid==202) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fm-3', spgid=spgid)
    # 203: Fd-3
    if (spgid==203) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.75, 0.75, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fd-3', spgid=spgid)
    # 204: Im-3
    if (spgid==204) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Im-3', spgid=spgid)
    # 205: Pa-3
    if (spgid==205) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pa-3', spgid=spgid)
    # 206: Ia-3
    if (spgid==206) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0, 0.5]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Ia-3', spgid=spgid)
    # 207: P432
    if (spgid==207) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P432', spgid=spgid)
    # 208: P4_2 32
    if (spgid==208) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P4_2 32', spgid=spgid)
    # 209: F432
    if (spgid==209) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='F432', spgid=spgid)
    # 210: F4_1 32
    if (spgid==210) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.75, 0.75, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='F4_1 32', spgid=spgid)
    # 211: I432
    if (spgid==211) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I432', spgid=spgid)
    # 212: P4_3 32
    if (spgid==212) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.75, 0.25, 0.75]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P4_3 32', spgid=spgid)
    # 213: P4_1 32
    if (spgid==213) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.25, 0.75, 0.25]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P4_1 32', spgid=spgid)
    # 214: I4_1 32
    if (spgid==214) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.25, 0.75, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I4_1 32', spgid=spgid)
    # 215: P-43m
    if (spgid==215) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P-43m', spgid=spgid)
    # 216: F-43m
    if (spgid==216) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='F-43m', spgid=spgid)
    # 217: I-43m
    if (spgid==217) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-43m', spgid=spgid)
    # 218: P-43n
    if (spgid==218) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb],name='P-43n', spgid=spgid)

    # 219: F-43c
    if (spgid==219) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='F-43c', spgid=spgid)
    # 220: I-43d
    if (spgid==220) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop(-np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.75, 0.25, 0.75]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb],name='I-43d', spgid=spgid)
    # 221: Pm-3m
    if (spgid==221) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pm-3m', spgid=spgid)
    # 222: Pn-3n
    if (spgid==222) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0, 0]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pn-3n', spgid=spgid)
    # 223: Pm-3n
    if (spgid==223) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pm-3n', spgid=spgid)
    # 224: Pn-3m
    if (spgid==224) :
        AA = lattvec([ [1,0,0], [0,1,0], [0,0,1] ])
        BB = np.linalg.inv(AA).T
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.5, 0.5]), \
                    AA=AA, BB=BB )
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB )
        return SpaceGroup([aa,bb,cc],name='Pn-3m', spgid=spgid)
    # 225: Fm-3m
    if (spgid==225) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fm-3m', spgid=spgid)
    # 226: Fm-3c
    if (spgid==226) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.5, 0.5]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fm-3c', spgid=spgid)
    # 227: Fd-3m
    if (spgid==227) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.5, 0.75, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fd-3m', spgid=spgid)
    # 228: Fd-3c
    if (spgid==228) :
        AA = lattvec([ [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0.75, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Fd-3c', spgid=spgid)
    # 229: Im-3m
    if (spgid==229) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Im-3m', spgid=spgid)
    # 230: Ia-3d
    if (spgid==230) :
        AA = lattvec([ [-0.5,0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,-0.5] ])
        BB = np.linalg.inv(AA).T
        AAC= lattvec([ [1.0,0.0,0], [ 0.0,1.0,0], [0,0,1] ])
        BBC= np.linalg.inv(AAC).T 
        aa = symop( np.array( [ [0,0,1], [1,0,0], [0,1,0] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        bb = symop( np.array( [ [0,-1,0], [1,0,0], [0,0,1] ]), lattvec([0.25, 0.75, 0.25]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        cc = symop( np.array( [ [-1,0,0], [0,-1,0], [0,0,-1] ]), lattvec([0, 0, 0]), \
                    AA=AA, BB=BB , AAC=AAC, BBC=BBC)
        return SpaceGroup([aa,bb,cc],name='Ia-3d', spgid=spgid)
    # -2: p2
    if (spgid==-2):
        aa = symop( np.array([ [-1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0, 0, 0]) )
        return SpaceGroup([aa],name='p2',spgid=spgid)
    # -4: p1g1
    if (spgid==-4):
        aa = symop( np.array([ [1,0,0], [0,-1,0], [0,0,1] ]), lattvec([0.5, 0, 0]) )
        return SpaceGroup([aa],name='pg',spgid=spgid)
    #
    wg2sg = {-1:1, -2:3, -3:6, -5:8, -6:25, -7:28, -8:32, -9:35, -10:75, -11:99, -12:100,\
             -13:143, -14:156, -15:157, -16:168, -17:183}
    wg2sg_nm = {-1:'p1', -2:'p2', -3:'pm', -5:'cm', -6:'p2mm', -7:'p2mg', -8:'p2gg', -9:'c2mm',\
                -10:'p4', -11:'p4mm', -12:'p4gm', -13:'p3', -14:'p3m1', -15:'p31m', -16:'p6',\
                -17:'p6mm'}
    #
    if spgid in wg2sg:
        sg = get_SPGfromID(wg2sg[spgid])
        sg.name = wg2sg_nm[spgid]
        sg.spgid= spgid
        return sg
        
    return None


def getSPG(ispg,new):
    #
    if not os.path.exists('SPG') :
        os.mkdir('SPG')
    name = './SPG/'+str(ispg)+'.dat'
    new_= new or (not os.path.isfile(name))
    if new_:
        spg = get_SPGfromID(ispg)
        fd = open(name, 'wb')
        pk.dump(spg,fd)
        fd.close()
    else :
        fd = open(name, 'rb')
        spg = pk.load(fd)
        fd.close()
        if not isinstance(spg,SpaceGroup) :
            return getSPG(ispg,True)
    return spg


if __name__ == '__main__' :
    #
    if len(sys.argv)==1 :
        print('Please give the space group number')
    elif len(sys.argv)==2 :
        argv1= sys.argv[1]
        new0=False
    elif len(sys.argv)>=3 :
        argv1= sys.argv[1]
        new0=True if sys.argv[2][0]=='n' else False
    #
    # One SG
    if argv1[0] not in 'sw':
        ispg = int(argv1)
        spg  = getSPG(ispg,new0)
        print(spg)
    #
### # All SGs
### elif argv1[0] == 's':
###     for ispg in range(1,231) :
###         spg  = getSPG(ispg,new0)
###         print(spg)
### elif argv1[0] == 'w':
###     for ispg in range(1,18) :
###         spg  = getSPG(-ispg,new0)
###         print(spg)
