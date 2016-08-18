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

def brute_force_sparseness(target, stabs):
    for i in range(len(stabs)):
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

def do_for_n_qubits(n, queue=None):
    strs = ['H', 'F', 'random']
    targets = (qt.tensor([H]*n), qt.tensor([F]*n),
               qt.rand_ket(pow(2,n), dims=[[2]*n, [1]*n]))
    stabs = stab_states(n)
    print('Do SL0 estimate')
    SL0s = SL0_estimate(targets, stabs, n)
    print('Do brute force sparseness')
    sparseness = [tuple(brute_force_sparseness(target, stabs)) for target in targets]
    out_string = """For the {0} state on {1} qubits, the SL0 estimate gives a sparseness
    of {2}, and the brute force search returns a sparseness of {3}.
    The resulting basis states found were:
    """
    res = []
    for i, s in enumerate(strs):
        res.append(out_string.format(s, n, SL0s[i], sparseness[i][1]) 
                   + "\n".join([str(b) for b in sparseness[i][0]]))
    if queue:
        queue.put(res)
        return
    return res
if __name__ == '__main__':
    ostring = datetime.datetime.now().strftime('%d%m%Y_%H%M%S')+".txt"
    queue = multiprocessing.Queue()
    ns = [1,2,3,4]
    for n in ns: #Running in parallel is too resource intensive
                 #Calling in a subprocess guarantees allocated memory is freed on completion
        p = multiprocessing.Process(target=do_for_n_qubits, args=(n, queue))
        p.start()
        p.join()
        res = queue.get()
        # res = do_for_n_qubits(n)
        with open(ostring, 'a') as f:
            for bit in res:
                f.write(bit+"\n")


    
