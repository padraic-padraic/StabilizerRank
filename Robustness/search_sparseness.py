from Utils import gen_stabiliser_groups, n_stab, Projector, OrthoProjector, SL0, stab_states
# from Utils.dispatcher import star_execution

import datetime
import multiprocessing
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
plus = (qt.basis(2,0)+qt.basis(2,1)).unit()
S = qt.sigmaz().sqrtm()
T = S.sqrtm()
RT = T.sqrtm() * plus
RRT = T.sqrtm().sqrtm() * plus

def brute_force_sparseness(target, stabs):
    print('Start brute force')
    for i in range(len(stabs)):
        print('Test with ' + str(i+1) + ' states')
        for basis in combinations(stabs, i+1): #Try all combinations of i stabiliser states
            proj = OrthoProjector([b.full() for b in basis])
            projection = np.linalg.norm(proj*target.full(), 2)
            if np.allclose(projection, 1.):
                print('Done brute force for target')
                return basis, i+1
    return 'This probably shouldn\'t have happened you dolt.'

def SL0_estimate(targets, stabs, n_qubits):
    res = []
    h_dim = pow(2,n_qubits)
    basis = [qt.basis(h_dim, i) for i in range(h_dim)]
    A = np.matrix(np.zeros([h_dim, len(stabs)], dtype=np.complex_))
    for i in range(h_dim): ##Build the h_dim x n_stab dimensional A matrix once and only once
        for j in range(len(stabs)):
            try:
                A[i,j] = basis[i].overlap(stabs[j])
            except TypeError:
                print(basis[i])
                print(stabs[j])
                raise TypeError('Quitting for real now')
    b = np.matrix(np.zeros(h_dim, dtype=np.complex_)).T
    for target in targets:
        for i in range(h_dim): #Build the state vector as a numpy matrix
            b[i] = basis[i].overlap(target)
        x = SL0(A, b, 1e-15) #Run the SL0 routine defined in Utils.SL0
        for i in range(x.size):
            x[i] = np.abs(x[i])
        # print(x)
        res.append(np.count_nonzero(x))
        print('SL0 done for target')
    return res

FILE_NAME = datetime.datetime.now().strftime('%d%m%Y_%H%M%S')+".txt"
OUT_STRING_1 = "For the {0} state on {1} qubits,"
OUT_STRING_2 = "the SL0 estimate gives a sparsenes of {}"
OUT_STRING_3 = "the brute force search returns a sparseness of {}."
OUT_STRING_4 = "The resulting basis states found were:"

if __name__ == '__main__':
    FILE_NAME = 'RandSparse.txt'
    # ns = [1,2,3,4]
    # strs = ['H', 'F']
    # strs = ['Root T', 'Root Root T']
    r = qt.rand_ket(2)
    r2 = qt.rand_ket(4)
    targets = [r, qt.tensor(r,r), r2, qt.tensor(r,r,r)]
    ns = [1, 2, 2, 3]
    strs = ['Random', '2fold product', '2 qubit random', '3 fold product']
    for i, n in enumerate(ns):
        # targets = [qt.tensor([H]*n), qt.tensor([F]*n)]
        # targets = [qt.tensor([RT]*n), qt.tensor([RRT]*n)]
        stabs = stab_states(n)
        SL0s = SL0_estimate([targets[i]], stabs, n)
        basis, val = brute_force_sparseness(targets[i], stabs)
        with open(FILE_NAME, 'a') as f:
            out = OUT_STRING_1.format(strs[i], n)+"\n"
            out += OUT_STRING_2.format(SL0s[0])
            out += OUT_STRING_3.format(val) + "\n"
            out += "\n".join(str(b) for b in basis) +"\n"
            f.write(out)

    #for n in ns: #Running in parallel is too resource intensive
                 #Calling in a subprocess guarantees allocated memory is freed on completion
        # p = multiprocessing.Process(target=do_for_n_qubits, args=(n, ostring))
        # p.start()
        # p.join()


