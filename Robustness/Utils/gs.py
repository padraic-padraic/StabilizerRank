"""Functions to perform the GramSchmidt Orthogonalisation procedure for a collection
of Stabiliser states, and return the Projector on to the space spanned by them.
"""
import numpy as np

from scipy.linalg import lu

__all__ = ['Projector', 'OrthoProjector','GramSchmidt']

def check_lin_independence(vectors):
    """Gram-Schmidt only applies if vectors are linearly independent.
    We expect this to be tha case given Stab States for a mutually unbiased 
    basis but it's worth checking anyway."""
    if len(vectors) == 1:
        return True
    M = np.zeros([len(vectors), len(vectors[0])], dtype=np.complex_)
    for i in range(len(vectors)):
        M[i] = vectors[i].T
    pl, u = lu(M, permute_l=True)
    if any([np.count_nonzero(M[i]) == 0 for i in range(len(vectors))]):
        return False #M must be full rank for linear independence
    return True

def gs_prj(base, target):
    return np.sum((base.H*target)) / np.sum((base.H*base)) * base

def GramSchmidt(vectors):
    dim = vectors[0].size
    if not check_lin_independence(vectors):
        raise ValueError('Oh god what even')
    V = np.matrix(np.zeros([dim,dim], dtype=np.complex_))
    U = np.matrix(np.zeros([dim,dim], dtype=np.complex_))
    for i in range(len(vectors)):
        V[:,i] = vectors[i]
    for i in range(len(vectors)):
        U[:,i] = V[:,i]
        for j in range(i):
            U[:,i] -= gs_prj(U[:,j], V[:,i])
    for i in range(len(vectors)):
        norm = np.linalg.norm(np.matrix(U[:,i]), 2)
        U[:,i] /= norm
    return [np.matrix(U[:,i]) for i in range(len(vectors))]

def Projector(vectors):
    dim = len(vectors[0])
    A = np.matrix(np.zeros([dim, len(vectors)], dtype=np.complex_))
    for i in range(len(vectors)):
        A[:,i] = vectors[i]
    return A*(A.H*A).I*A.H

def OrthoProjector(vectors):
    dim  = len(vectors[0])
    ortho_vecs = GramSchmidt(vectors)
    A = np.matrix(np.zeros([dim, len(ortho_vecs)], dtype=np.complex_))
    for i in range(len(ortho_vecs)):
        A[:,i] = ortho_vecs[i]
    P = A*A.H
    return P
