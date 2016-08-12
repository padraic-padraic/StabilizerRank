from Utils import n_stab, Projector, OrthoProjector, SL0, stab_states
from Utils.dispatcher import star_execution

import datetime
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
            proj = OrthoProjector([b.full() for b in basis])
            projection = np.linalg.norm(proj*target.full(), 2)
            if np.allclose(projection, 1.):
                return basis, i+1
    return 'This probably shouldn\'t have happened you dolt.'

def SL0_estimate(target, stabs, n_qubits):
    h_dim = pow(2,n_qubits)
    basis = [qt.basis(h_dim, i) for i in range(h_dim)]
    A = np.matrix(np.zeros([h_dim, len(stabs)], dtype=np.complex_))
    for i in range(h_dim): ##Build the h_dim x n_stab dimensional A matric
        for j in range(len(stabs)):
            try:
                A[i,j] = basis[i].overlap(stabs[j])
            except TypeError:
                print(basis[i])
                print(stabs[j])
                raise TypeError('Quitting for real now')
    b = np.matrix(np.zeros(h_dim, dtype=np.complex_)).T
    for i in range(h_dim): #Build the state vector as a numpy matrix
        b[i] = basis[i].overlap(target)
    x = SL0(A, b, 1e-15) #Run the SL0 routine defined in Utils.SL0
    for i in range(x.size):
        x[i] = np.abs(x[i])
    # print(x)
    return np.count_nonzero(x)

def do_for_n_qubits(n, **kwargs):
    print('Doing for ' + str(n) + ' qubits')
    out_str = """For the {0} state, the SL0 algorithm gives a sparseness of {1},
                 and the brute force search gives {2} for {3} qubits"""
    out_str2 = """Minimal stabiliser decomposition found was:\n"""
    stabs = stab_states(n)
    target = qt.tensor([H]*n).unit()
    l_norm = SL0_estimate(target, stabs, n)
    basis, sparsity = brute_force_sparseness(target, stabs)
    res = [out_str.format('H', l_norm, sparsity, n), 
           out_str2 + "\n".join([str(b) for b in basis])]
    target = qt.tensor([F]*n).unit()
    l_norm = SL0_estimate(target, stabs, n)
    basis, sparsity = brute_force_sparseness(target, stabs)
    res.append(out_str.format('F', l_norm, sparsity, n)) 
    res.append(out_str2 + "\n".join([str(b) for b in basis]))
    target = qt.rand_ket(pow(2,n))
    basis, sparsity = brute_force_sparseness(target, stabs)
    res.append(out_str.format('random', l_norm, sparsity, n)) 
    res.append(out_str2 + "\n".join([str(b) for b in basis]))
    print('Done for ' + str(n) + ' qubits')
    return res

if __name__ == '__main__':
    # ns = [[1], [2], [3], [4], [5], [6], [7]]
    # kwargs_list = [{}]*len(ns)
    # outputs = star_execution(do_for_n_qubits, ns, kwargs_list)
    ostring = datetime.datetime.now().strftime('%d%m%Y_%H%M%S')+".txt"
    for i in range(1,4):
        with open(ostring, 'w') as f:
            out = do_for_n_qubits(i)
            for res in out:
                f.write(res)
                f.write("\n")
    