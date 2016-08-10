"""Function to generate the full set of n-qubit stabiliser groups
using the set of n-qubit paulis with real-valued phase.
Note this set excludes the n-fold identity, and -1*n-fold identity."""
from itertools import combinations
from qutip import commutator, qeye, Qobj, sigmax, sigmay, sigmaz, tensor

import numpy as np

zer_obj = Qobj(np.zeros((4,4)), dims=[[2,2], [2,2]])
I = qeye(2)
X = sigmax()
Y = sigmay()
Z = sigmaz()

__all__ = ['gen_stabiliser_groups', 'BinarySubspace', 'string_to_pauli']

def xor(a,b):
    return (a|b)&~(a&b)

class BinarySubspace(object):
    """Set-like class for bitarray objects to generate a closed subspace."""
    def __init__(self, *data):
        self.order = 0
        self._items = []
        self.generators = []
        for val in data:
            if not isinstance(val, np.ndarray):
                raise ValueError('This class works for numpy arrays only!')
            self.add(val)
    
    def __contains__(self, it):
        for _el in self._items:
            if np.array_equal(_el, it):
                return True
        return False

    def __iter__(self):
        for item in self._items:
            yield item

    def _generate(self, obj):
        for item in self._items:
            new = xor(item, obj)
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
            if np.array_equal(obj, _el):
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
        if not x and not z:
            pauli_chain.append(I)
        elif x and z:
            pauli_chain.append(Y)
        elif x and not z:
            pauli_chain.append(X)
        else:
            pauli_chain.append(Z) 
    if phase:
        return -1 * tensor(pauli_chain)
    return tensor(pauli_chain)

def test_commutivity(n, bits1, bits2):
    if np.array_equal(bits1[:-1], bits2[:-1]): ##Additional check, removes groups
    ## Supposedly 'stabilised' by e.g IX, -IX....
        return False
    p1, p2 = string_to_pauli(n, bits1), string_to_pauli(n, bits2)
    return commutator(p1,p2) == zer_obj

def find_generators(n, bitstrings):
    subspaces = []
    for group in combinations(bitstrings, n):
        if len(group) == 2:
            if not test_commutivity(n, group[0], group[1]):
                continue
        if len(group) > 2:
            if not all([test_commutivity(n, pair[0], pair[1]) 
                        for pair in combinations(group, 2)]): 
                continue
        candidate = BinarySubspace(*group)
        if any([candidate == space for space in subspaces]):
            continue
        subspaces.append(candidate)

    return [tuple(subspace.generators) for subspace in subspaces]

def gen_stabiliser_groups(n):
    bitstrings = []
    for i in range(2, pow(2, 2*n+1)): #We ignore the strings for 0 and 1 as these correspond to I^{n} and -I^{n}.
        bin_string = bin(i)[2:] #strip the 0b from the string
        bin_string = '0'*(2*n+1 - len(bin_string)) + bin_string
        a = np.array([b == '1' for b in bin_string])
        bitstrings.append(a)
    print('Found {} binary strings'.format(len(bitstrings)))
    generators = find_generators(n, bitstrings)
    # print(generators)
    print('Found {} unique generators'.format(len(generators)))
    pauli_generators = []
    for group in generators:
        pauli_generators.append([string_to_pauli(n, p) for p in group])
    return pauli_generators