#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import multiprocessing
import numpy as np

from bitarray import bitarray
from cmath import exp as cexp #Differentiate from real-valued exp
from copy import deepcopy
from math import acos, cos, exp, sin, sqrt
from math import pow as fpow #Differentiate from stdlib pow
from random import randrange, random
from qutip import basis, qeye, tensor, sigmaz
from Utils import gen_random_stabilisers, OrthoProjector, string_to_pauli
from Utils.dispatcher import N_PROCESSORS
from Utils.pretty_print import AnalysisResult

BETA = acos(1/sqrt(3)) /2
PI = acos(0.)*2
H = (cos(PI/8)*basis(2,0) + sin(PI/8)*basis(2,1)).unit()
F = (cos(BETA)*basis(2,0) + cexp(1j*PI/4)*sin(BETA)*basis(2,1)).unit()
PLUS = (basis(2,0)+basis(2,1)).unit()
S = sigmaz().sqrtm()
T = S.sqrtm()
RT = T.sqrtm() * PLUS
RRT = T.sqrtm().sqrtm() * PLUS

STATES = {'T':T,
          'RT':RT,
          'RRT':RRT,
          'F':F,
          'H':H}

success_string = """A decomposition was found with rank {chi}. The list of 
                    states found looked like:\n {states}"""

write_string = """Checking for the existence of a decompositon of the state 
                  {target} on {n_qubits} qubit(s).\n {}"""

def format_for_output(job_dict):
    func = job_dict.pop('func')
    fname = job_dict.pop('fname')
    ostring = job_dict.pop('ostring')
    n = job_dict.pop('n_qubits')
    target = job_dict.pop('target_string')
    output = func(**job_dict.pop('func_inputs'))
    return AnalysisResult(n, target, output, fname=fname, ostring=write_string)

def do_anneal(**kwargs):
    beta = kwargs.pop('beta_init', 1)
    beta_max = kwargs.pop('beta_max', 4000)
    anneal_steps = kwargs.pop('M', 1000)
    b_diff = (beta_max-beta)/anneal_steps
    walk_steps = kwargs.pop('steps', 100)
    n_qubits = kwargs.pop('n_qubits')
    target = kwargs.pop('target')
    chi = kwargs.pop('chi')
    states = gen_random_stabilisers(n_qubits, chi)
    identity = qeye(pow(2,n_qubits))
    identity.dims = [[2]*n_qubits, [2]*n_qubits]
    while beta <= beta_max:
        for i in range(walk_steps):
            projector = OrthoProjector([s.full() for s in states])
            projection = np.linalg.norm(projector*target.full(), 2)
            if np.allclose(projection, 1.):
                return success_string.format(chi=chi, states=states)
            a = randrange(0,len(states))
            p_string = bin(randrange(0,pow(2,2*n_qubits)))[2:]
            p_bits = bitarray(2*n_qubits - len(p_string))
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

def run_analysis():
    states = ['RT', 'RRT', 'F']
    ns = [3,4]*len(states)
    results = []
    print('Create Pool')
    pool = multiprocessing.Pool(8)
    print('Created Pool')
    for state, n in zip(states, ns):
        print('Prepping {} with {} qubit(s).'.format(state,n))
        details = {'ostring': write_string,
                   'fname':"/home/ucaphdc/Scratch/output/"+state+".txt",
                   'n_qubits':n,
                   'target_string':state,
                   'func':do_anneal}
        func_inputs = {'target':tensor([STATES[state]]*n),
                       'n_qubits':n}
        for i in range(2, pow(2,n)):
            job = deepcopy(details)
            job['func_inputs'] = deepcopy(func_inputs)
            job['func_inputs']['chi'] = i
            results.append(pool.apply_async(format_for_output, [job])) #Pool is so dumb
    pool.close()
    pool.join()
    while results:
        for index, res in enumerate(results):
            if res.ready():
                out = results.pop(index).get()
                out.write()
                break

if __name__ == '__main__':
    run_analysis()

