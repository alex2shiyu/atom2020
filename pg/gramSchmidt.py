#!/usr/bin/env python
import numpy as np

"""
Example:
    input :
            vspace =  [[1,2,3],[6,3,0],[5,1,6]]
            onorm = gramSchmidt.gramSchmidt(vspace)
                >> array([[ 0.26726124,  0.53452248,  0.80178373],
                        [ 0.87287156,  0.21821789, -0.43643578],
                        [ 0.40824829, -0.81649658,  0.40824829]])
"""


def proj(x,u):
    ## Can't hurt
    u = unit_vec(u)
    return np.dot(np.conjugate(x),u)/np.dot(np.conjugate(u),u) * u

def unit_vec(x):
    """Get unit vector of x. Same direction, norm==1"""
#   return x/np.linalg.norm(x)
    return x

def GramSchmidt(vectors):
    """ _correct_ recursive implementation of Gram Schmidt algo that is not subject to 
    rounding erros that the original formulation is. 

    Function signature and usage is the same as gramSchmidt()
    """
    ###### Ensure the input is a 2d array (or can be treated like one)
    vectors = np.atleast_2d(vectors)

    ###### Handle End Conditions
    if len(vectors) == 0:
        return [[]]

    ## Always just take unit vector of first vector for the start of the basis
    u1 = unit_vec(vectors[0])

    if len(vectors) == 1:
        return np.array([u1])

    ###### Orthonormalize the rest of the vectors
    #           | easy row stacking
    #           |                                            | Get the orthagonal projection of each subsequent vector onto u1 (ensures whole space is now orthagonal to u1)                  
    #                       | Recurse on the projections     |
    basis = np.vstack( (u1, GramSchmidt( list(map(lambda v: v - proj(v,u1), vectors[1:])))) ) # not explicit list(map) conversion, need for python3+

    return np.array(basis)

def _is_orthag(vectors):
    """Simple check, sees if all of the vectors in v are orthagonal to eachother.

    Takes the dot product of each vector pair, sees if the result is close to zero
    """
    orthag = True
    vectors = np.atleast_2d(vectors)
    for vector in vectors:
        for vector2 in vectors:
            ## Don't dot itself
            if np.array_equal(vector,vector2):
                continue
            ## Dot product alwys has some numerical precision remainder
            if abs(np.dot(vector,vector2)) > 1e-5:
                orthag = False
    return orthag

