"""Function to generate the full set of n-qubit stabiliser groups
using the set of n-qubit paulis with real-valued phase.
Note this set excludes the n-fold identity, and -1*n-fold identity."""
from .binary_subspace import BinarySubspace
from .n_stab import n_stab
from bitarray import bitarray
from functools import reduce
from itertools import combinations
from qutip import qeye, sigmax, sigmay, sigmaz, tensor
from random import sample, randrange, random

import numpy as np
import operator as op
import os
import pickle

I = qeye(2)
X = sigmax()
Y = sigmay()
Z = sigmaz()

__all__ = ['gen_stabiliser_groups', 'BinarySubspace', 'string_to_pauli',
           'stab_states', 'get_proj_eigenstate', 'gen_random_stabilisers']

def ncr(n, r):
    """Efficient evaluation of ncr, taken from StackOverflow
    http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python"""
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom

# def random_combination(iterable, r):
#     """Random selection from itertools.combinations(iterable, r)
#     Taken from itertools.recipes on pydoc"""
#     pool = tuple(iterable)
#     n = len(pool)
#     indices = sorted(sample(range(n), r))
#     return tuple(pool[i] for i in indices)

def string_to_pauli(bits):
    n = len(bits)//2
    pauli_chain = []
    for x, z in zip(bits[:n], bits[n:]):
        if not x and not z:
            pauli_chain.append(I)
        elif x and z:
            pauli_chain.append(Y)
        elif x and not z:
            pauli_chain.append(X)
        else:
            pauli_chain.append(Z) 
    return tensor(pauli_chain)

def symplectic_inner_product(n, a, b):
    x_a, z_a = a[:n], a[n:]
    x_b, z_b = b[:n], b[n:]
    # count = (x_a&z_b).count() + (x_b&z_a).count()
    count = np.sum((x_a&z_b)) + np.sum((x_b&z_a))
    return count%2

def test_commutivity(n, bits1, bits2):
    return symplectic_inner_product(n, bits1, bits2) == 0 #1 if they anticommute, 0 if they commute

def gen_bitstrings(n):
    bitstrings = []
    for i in range(1, pow(2,2*n)): #We ignore the all 0 string as it corresponds to I^{n}
        bin_string = bin(i)[2:] #strip the 0b from the string
        # a = bitarray(2*n - len(bin_string))
        # a.extend(bin_string)
        # bitstrings.append(a)
        bin_string = '0'*(2*n - len(bin_string)) + bin_string
        a = np.array([b == '1' for b in bin_string])
        bitstrings.append(a)
    return bitstrings

def find_generators(bitstrings, **kwargs):
    n = len(bitstrings[0])//2
    subspaces = []
    generators = []
    target = kwargs.pop('target', n_stab(n) // pow(2,n)) #We add in the phase later
    i=0
    for group in combinations(bitstrings, n):
        i+=1
        if len(group) == 2:
            if not test_commutivity(n, group[0], group[1]):
                continue
        if len(group) > 2:
            if not all([test_commutivity(n, pair[0], pair[1]) 
                        for pair in combinations(group, 2)]): 
                continue
        candidate = BinarySubspace(*group)
        if len(candidate.generators) < n:
            continue
        if len(candidate._items) < pow(2,n):
            continue
        # res = tuple(i for i in sorted(candidate._items))
        res = tuple(i for i in sorted(candidate._items, key=np.sum))
        # if not res in subspaces:
        for space in subspaces:
            if np.all([np.array_equal(_el1, _el2) for _el1, _el2 in zip(res, space)]):
                continue
        subspaces.append(res)
        generators.append(tuple(candidate.generators))
        if len(generators) == target:
            break
    return generators

def phaseify_groups(n, generators):
    # Add phase 'by hand'
    phase_strings = []
    for i in range(1, pow(2,n)): #2^n different phase strings exist
        base = bin(i)[2:]
        a = bitarray(n-len(base))
        a.extend(base)
        phase_strings.append(a)
    for i in range(len(generators)):
        for ps in phase_strings:
            generators.append([-1*p if b else p 
                                    for p, b in zip(generators[i], ps)])
    return generators

def gen_stabiliser_groups(n):
    path = os.path.join(os.path.dirname(__file__),
                        str(n)+'_generators.pkl')    
    if os.path.isfile(path):
        return pickle.load(open(path, 'rb'))
    bitstrings = gen_bitstrings(n)
    print('Found {} binary strings'.format(len(bitstrings)))
    generators = find_generators(bitstrings)
    nophase_total = len(generators)
    print('Found {} unique generators, sans phase'.format(nophase_total))
    pauli_generators = []
    for group in generators:
        pauli_generators.append([string_to_pauli(p) for p in group])
    pauli_generators = phaseify_groups(n, pauli_generators)
    print('Added phase to give {} generatprs'.format(len(pauli_generators)))
    pickle.dump(pauli_generators, open(path, 'wb'))
    return pauli_generators

def find_projectors(generating_sets):
    n_qubits = len(generating_sets[0])
    Id= qeye(pow(2, n_qubits))
    dims = [[2]*n_qubits, [2]*n_qubits]
    Id.dims = dims
    projectors = []
    for genset in generating_sets:
        res = qeye(pow(2, n_qubits))
        res.dims = dims
        for g in genset:
            res *= (Id+g)/2
        projectors.append(res)
    return projectors

def get_proj_eigenstate(projector):
    eigs, vecs = projector.eigenstates(sort='high')
    for n, eig in enumerate(eigs):
        if np.allclose(eig, complex(1)) or np.allclose(eig, 1.):
            return vecs[n]
        else:
            print(eigs)
            return None

def stab_states(n):
    path = os.path.join(os.path.dirname(__file__),
                        str(n)+'_stabs.pkl')
    if os.path.isfile(path):
        print('Loaded')
        return pickle.load(open(path, 'rb'))
    projectors = find_projectors(gen_stabiliser_groups(n))
    if len(projectors) != n_stab(n):
        raise ValueError('Your code is bad and you should feel bad')
    vals = [get_proj_eigenstate(proj) for proj in projectors]
    if any([v is None for v in vals]):
        print('Eep')
    pickle.dump(vals, open(path, 'wb'), -1)
    return vals

def gen_random_stabilisers(n_qubits, chi):
    print('Doing for chi={0} on {1} qubit(s)'.format(chi, n_qubits))
    bitstrings = gen_bitstrings(n_qubits)
    generators = find_generators(bitstrings, target=chi)
    pauli_generators = []
    for group in generators:
        pauli_generators.append([string_to_pauli(g) for g in group])
    for i in range(chi):
        if random() > (1 / pow(2, n_qubits)): # Add a phase! Randomly...
            phase_num = bin(randrange(1,pow(2,n_qubits)))[2:]
            p_bits = bitarray(n_qubits-len(phase_num))
            p_bits.extend(phase_num)
            pauli_generators[i] = [-1*p if b else p 
                                    for p, b in zip(pauli_generators[i], p_bits)]
    return [get_proj_eigenstate(proj) for proj in 
            find_projectors(pauli_generators)]
