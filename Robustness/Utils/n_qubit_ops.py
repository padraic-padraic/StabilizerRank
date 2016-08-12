"""Function to generate the full set of n-qubit stabiliser groups
using the set of n-qubit paulis with real-valued phase.
Note this set excludes the n-fold identity, and -1*n-fold identity."""
from .n_stab import n_stab
from bitarray import bitarray
from functools import reduce
from itertools import combinations
from qutip import commutator, qeye, Qobj, sigmax, sigmay, sigmaz, tensor
from random import sample, shuffle

import numpy as np
import operator as op
import os
import pickle

I = qeye(2)
X = sigmax()
Y = sigmay()
Z = sigmaz()

__all__ = ['gen_stabiliser_groups', 'BinarySubspace', 'string_to_pauli',
           'stab_states', 'get_proj_eigenstate']

def ncr(n, r):
    """Efficient evaluation of ncr, taken from StackOverflow
    http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python"""
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom

def random_combination(iterable, r):
    """Random selection from itertools.combinations(iterable, r)
    Taken from itertools.recipes on pydoc"""
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(sample(range(n), r))
    return tuple(pool[i] for i in indices)

def xor(a,b):
    return (a|b)&~(a&b)

def xnor(a,b):
    return (a&b)^(~a&~b)

class BinarySubspace(object):
    """Set-like class for bitarray objects to generate a closed subspace."""
    def __init__(self, *data):
        self.order = 0
        self._items = []
        self.generators = []
        for val in data:
            # if not isinstance(val, np.ndarray):
            #     raise ValueError('This class works for numpy arrays only!')
            self.add(val)
    
    def __contains__(self, it):
        for _el in self._items:
            # if np.array_equal(_el, it):
            if all(xnor(_el, it)):
                return True
        return False

    def __iter__(self):
        for item in self._items:
            yield item

    def _generate(self, obj):
        for item in self._items:
            new = item^obj
            # new = xor(item, obj)
            if new in self:
                continue
            else:
                self.order +=1
                self._items.append(new)
                self._generate(new)
        return

    def __eq__(self, other):
        return all([_el in other for _el in self._items])

    def add(self, obj):
        for _el in self._items:
            # if np.array_equal(obj, _el):
            if all(xnor(obj, _el)):
                return self
        self.order +=1
        self.generators.append(obj)
        self._items.append(obj)
        self._generate(obj)
        return self

def string_to_pauli(n, bits):
    phase = bits[-1]
    bits = np.delete(bits, -1)
    pauli_chain = []
    for x, z in zip(bits[:n], bits[n:]):
        if ~x and ~z:
            pauli_chain.append(I)
        elif x and z:
            pauli_chain.append(Y)
        elif x and ~z:
            pauli_chain.append(X)
        else:
            pauli_chain.append(Z) 
    if phase:
        return -1 * tensor(pauli_chain)
    return tensor(pauli_chain)

def test_commutivity(n, bits1, bits2):
    if all(xnor(bits1[:-1], bits2[:-1])): ##Additional check, removes groups
    ## Supposedly 'stabilised' by e.g IX, -IX....
        return False
    p1, p2 = string_to_pauli(n, bits1), string_to_pauli(n, bits2)
    return np.sum(np.abs(commutator(p1,p2).full())) == 0.
    # pbits1 = bits1[:n]|bits1[n:-1]
    # pbits2 = bits2[:n]|bits2[n:-1]
    # p_count = 0
    # for b1, b2 in zip(pbits1, pbits2):
    #     if b1 and b2:
    #         p_count += 1
    # return p_count%2 == 0

def find_generators(n, bitstrings):
    path = os.path.join(os.path.dirname(__file__),
                        str(n)+'_generators.pkl')    
    if os.path.isfile(path):
        return pickle.load(open(path, 'rb'))
    subspaces = set()
    generators = []
    target = n_stab(n)
    for group in random_combination(combinations(bitstrings, n),
                                    ncr(len(bitstrings), n)):
        if len(group) == 2:
            if not test_commutivity(n, group[0], group[1]):
                continue
        if len(group) > 2:
            if not all([test_commutivity(n, pair[0], pair[1]) 
                        for pair in combinations(group, 2)]): 
                continue
        candidate = BinarySubspace(*group)
        found = len(subspaces)
        subspaces.add(tuple([i.to01() for i in sorted(candidate._items)]))
        if len(subspaces) == found+1:
            generators.append(candidate.generators)
        # print(len(subspaces))
        if len(subspaces) == target:
            break
    res =  [tuple(gen_set) for gen_set in generators]
    pickle.dump(res, open(path, 'wb'))
    return res

def gen_stabiliser_groups(n):
    bitstrings = []
    for i in range(2, pow(2, 2*n+1)): #We ignore the strings for 0 and 1 as these correspond to I^{n} and -I^{n}.
        bin_string = bin(i)[2:] #strip the 0b from the string
        a = bitarray(2*n+1 - len(bin_string))
        a.extend(bin_string)
        # bin_string = '0'*(2*n+1 - len(bin_string)) + bin_string
        # a = np.array([b == '1' for b in bin_string])
        bitstrings.append(a)
    print('Found {} binary strings'.format(len(bitstrings)))
    # shuffle(bitstrings)
    generators = find_generators(n, bitstrings)
    # print(generators)
    print('Found {} unique generators'.format(len(generators)))
    pauli_generators = []
    for group in generators:
        pauli_generators.append([string_to_pauli(n, p) for p in group])
    return pauli_generators

def projector(generators, n_qubits):
    I = qeye(pow(2, n_qubits))
    res = qeye(pow(2, n_qubits))
    dims = [[2]*n_qubits, [2]*n_qubits]
    I.dims = dims
    res.dims = dims
    for gen in generators:
        res *= (I+gen)/2
    return res

def get_proj_eigenstate(projector):
    eigs, vecs = projector.eigenstates(sort='high')
    for n, eig in enumerate(eigs):
        if np.allclose(eig, complex(1)):
            return vecs[n]

def stab_states(n):
    path = os.path.join(os.path.dirname(__file__),
                        str(n)+'_stabs.pkl')
    if os.path.isfile(path):
        print('Loaded')
        return pickle.load(open(path, 'rb'))
    projectors = [projector(g, n) for g in gen_stabiliser_groups(n)]
    if len(projectors) != n_stab(n):
        raise ValueError('Your code is bad and you should feel bad')
    vals = [get_proj_eigenstate(proj)
            for proj in projectors]
    if any([v is None for v in vals]):
        print('Eep')
    pickle.dump(vals, open(path, 'wb'), -1)
    return vals