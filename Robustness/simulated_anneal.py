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
from Utils.dispatcher import OutputToQueue
from Utils.pretty_print import AnalysisResult

BETA = acos(1/sqrt(3)) /2
PI = acos(0.)*2
H = (cos(PI/8)*qt.basis(2,0) + sin(PI/8)*qt.basis(2,1)).unit()
F = (cos(BETA)*qt.basis(2,0) + cexp(1j*PI/4)*sin(BETA)*qt.basis(2,1)).unit()
PLUS = (qt.basis(2,0)+qt.basis(2,1)).unit()
S = qt.sigmaz().sqrtm()
T = S.sqrtm()
RT = T.sqrtm() * PLUS
RRT = T.sqrtm().sqrtm() * PLUS

STATES = {'T':T,
          'RT':RT,
          'RRT':RRT,
          'F':F,
          'H':H}

result_queue = multiprocessing.Queue()

@OutputToQueue(result_queue)
def format_for_output(**kwargs):
    func = kwargs.pop('func')
    fname = kwargs.pop('fname')
    ostring = kwargs.pop('ostring')
    n = kwargs.pop('n_qubits')
    target = kwargs.pop('target_string')
    output = func(kwargs.pop('func_inputs'))
    if output is None:
        return None
    return AnalysisResult(n, target, output, fname=fname, ostring=ostring)

def do_anneal(**kwargs):
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
    states = ['RT', 'RRT', 'F']
    ns = [2, 3, 4, 5]
    jobs = []
    for state, n in zip(states, ns):
        details = {'ostring':ostring,
                   'fname':state+".txt",
                   'n_qubits':n,}
                   'target_string':state}
        func_inputs = {'target':qt.tensor([STATES[state]]*n),
                       'n_qubits':n}
        for i in range(2, n):
            job = deepcopy(details)
            job['func_inputs'] = deepcopy(func_inputs)
            job['func_inputs']['chi'] = i
            jobs.append(job)
    pool = multiprocessing.Pool(Utils.dispatcher.N_PROCESSORS)
    results = pool.map_async(format_for_output, jobs)
    pool.close()
    pool.join()
    
