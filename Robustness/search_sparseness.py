from Utils import n_stab, stab_states, OrthoProjector, SL0

import numpy as np
import qutip as qt

from cmath import exp as cexp #Differentiate from real-valued exp
from itertools import combinations
from math import acos, cos, sin, sqrt
from math import pow as fpow #Differentiate from stdlib pow
BETA = acos(1/sqrt(3)) /2
PI = acos(0.)*2
H = (cos(PI/8)*qt.basis(2,0) + sin(PI/8)*qt.basis(2,1)).unit()
F = (cos(BETA)*qt.basis(2,0) + cexp(1j*PI/4)*sin(BETA)*qt.basis(2,1)).unit()

def brute_force_sparseness(target, stabs):
    for i in range(len(stabs)):
        for basis in combinations(stabs, i+1): #Try all combinations of i stabiliser states
            print(i+1)
            proj = OrthoProjector([b.full() for b in basis])
            projection = np.trace(proj*np.matrix(target.full()))
            print(projection)
            if np.allclose(projection, 1.):
                return i+1
    return 'This probably shouldn\'t have happened you dolt.'

def SL0_estimate(target, stabs, n_qubits):
    h_dim = pow(2,n_qubits)
    basis = [qt.basis(h_dim, i) for i in range(h_dim)]
    A = np.matrix(np.zeros([h_dim, len(stabs)], dtype=np.complex_))
    for i in range(h_dim): ##Build the h_dim x n_stab dimensional A matric
        for j in range(len(stabs)):
            A[i,j] = basis[i].overlap(stabs[j])
    b = np.matrix(np.zeros(h_dim, dtype=np.complex_)).T
    for i in range(h_dim): #Build the state vector as a numpy matrix
        b[i] = basis[i].overlap(target)
    x = SL0(A, b, 1e-15) #Run the SL0 routine defined in Utils.SL0
    return np.count_nonzero(x)

if __name__ == '__main__':
    out_str = """For the {0} state, the SL0 algorithm gives a sparseness of {1},
                 and the brute force search gives {2} for {3} qubits"""
    stabs = stab_states(2)
    target = qt.tensor(H,H)
    print(out_str.format('H', SL0_estimate(target, stabs, 2), 
                         brute_force_sparseness(target, stabs), 2))
    target = qt.tensor(F,F)
    print(out_str.format('F', SL0_estimate(target, stabs, 2), 
                         brute_force_sparseness(target, stabs), 2))
    target = qt.rand_ket(4)
    print(out_str.format('random', SL0_estimate(target, stabs, 2), 
                         brute_force_sparseness(target, stabs), 2))

