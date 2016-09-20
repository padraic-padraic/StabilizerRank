import multiprocessing
import numpy as np
import qutip as qt

from bitarray import bitarray
from cmath import exp as cexp #Differentiate from real-valued exp
from copy import deepcopy
from math import acos, cos, exp, sin, sqrt
from math import pow as fpow #Differentiate from stdlib pow
from random import randrange, random
from Utils import gen_random_stabilisers, OrthoProjector, string_to_pauli
from Utils.dispatcher import star_execution

BETA = acos(1/sqrt(3)) /2
PI = acos(0.)*2
H = (cos(PI/8)*qt.basis(2,0) + sin(PI/8)*qt.basis(2,1)).unit()
F = (cos(BETA)*qt.basis(2,0) + cexp(1j*PI/4)*sin(BETA)*qt.basis(2,1)).unit()
PLUS = (qt.basis(2,0)+qt.basis(2,1)).unit()
S = qt.sigmaz().sqrtm()
T = S.sqrtm()
RT = T.sqrtm() * PLUS
RRT = T.sqrtm().sqrtm() * PLUS

def do_anneal(n_qubits, target, chi, **kwargs):
    beta = kwargs.pop('beta_init', 1)
    beta_max = kwargs.pop('beta_max', 4000)
    anneal_steps = kwargs.pop('M', 1000)
    b_diff = (beta_max-beta)/anneal_steps
    walk_steps = kwargs.pop('steps', 100)
    states = gen_random_stabilisers(n_qubits, chi)
    identity = qt.qeye(pow(2,n_qubits))
    identity.dims = [[2]*n_qubits, [2]*n_qubits]
    while beta <= beta_max:
        for i in range(walk_steps):
            projector = OrthoProjector([s.full() for s in states])
            projection = np.linalg.norm(projector*target.full(), 2)
            if np.allclose(projection, 1.):
                return chi, states
            a = randrange(0,len(states))
            p_string = bin(randrange(0,pow(2,2*n_qubits)))[2:]
            p_bits = bitarray(2*n - len(p_string))
            p_bits.extend(p_string)
            pauli = string_to_pauli(p_bits)
            new_states = deepcopy(states)
            new_states[a] = ((identity + pauli).unit() * states[a])
            new_proj = OrthoProjector([s.full() for s in new_states])
            new_projection = np.linalg.norm(new_proj*target.full(), 2)
            if new_projection > projection:
                states = new_states
            else:
                p_accept = exp(-beta * (projection-new_projection))
                if random() >  p_accept:
                    states = new_states
        beta += b_diff
    return 'No decomposition found for chi={}\n'.format(chi)

if __name__ == '__main__':
    targets = [qt.tensor([RT]*2), qt.tensor([RT]*3)]
    ns = [2, 3]
    zip_list = [list(pair) for pair in zip(targets, ns)]
    args_list = []
    for state, n in zip_list:
        for i in range(2, n):
            args_list.append([n,state,i])
    kwargs_list = [{}]*2
    results = star_execution(do_anneal, args_list, kwargs_list)
    for r in results:
        if isinstance(r, tuple):
            print('Decomposition found with chi={}'.format(r[0]))
        else:
            print(r)
    # for n, target in zip(ns, targets):
    #     for chi in range(2, pow(2,n)+1):
    #         res = do_anneal(n, target, chi)
    #         if isinstance(res, tuple):
    #             print('Found decomposition with rank {}').format(res[0])
