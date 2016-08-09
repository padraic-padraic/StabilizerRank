"""Function to generate the full set of n-qubit stabiliser groups
using the set of n-qubit paulis with real-valued phase.
Note this set excludes the n-fold identity, and -1*n-fold identity."""
from bitarray import bitarray
from itertools import combinations
from math import pow
from qutip import qeye, sigmax, sigmaz, tensor

I = qeye(2)
X = sigmax()
Z = sigmaz()

class BinarySubspace(object):
    """Set-like class for bitarray objects to generate a closed subspace."""
    self._items = []
    self.generators = []
    def __init__(self, *data):
        self.order = 0
        for val in data:
            if not isinstance(val, bitarray):
                raise ValueError('This class works for bitarrays only!')
            self.add(val)
    
    def __contains__(self, it):
        for _el in self._items:
            if all(it&_el):
                return True
        return False

    def __iter__(self):
        for item in self._items:
            yield item

    def _generate(self, obj):
        for item in self._items:
            new = item^obj
            if new in self:
                continue
            else:
                self.append(new)

    def __eq__(self, other):
        return all([_el in other for _el in self._items])

    def _append(self, obj):
        self.order +=1
        self.generators.append(obj)
        self._items.append(obj)
        self._generate(obj)

    def add(self, obj):
        for _el in self._items:
            if obj == _el:
                return self
        self._append(obj)
        return self

def string_to_pauli(n, bits):
    phase = bits.pop(-1)
    pauli_chain = []
    for bit in zip(bits[:n], bits[n:]):
        if ~x & ~z:
            pauli_chain.append(I)
        elif x&z:
            pauli_chain.append(X*Z)
        elif x&~z:
            pauli_chain.append(X)
        else:
            pauli_chain.append(Z) 
        if phase:
            pauli_chain[-1] *= -1
    return tensor(pauli_chain)

def find_generators(n, bitstrings):
    subspaces = []
    for pair in combinations(bitstrings, 2):
        candidate = BinarySubspace(pair[0], pair[1])
        for subspace in subspaces:
            if candidate == subspace:
                continue
        subspaces.append(candidate)
    return [tuple(subspace.generators) for subspace in subspaces]

def gen_stabiliser_groups(n):
    bitstrings = []
    for i in range(2, pow(2,n)): #We ignore the strings for 0 and 1 as these correspond to I^{n} and -I^{n}.
        bin_string = bin(i)[2:] #strip the 0b from the string
        a = bitarray(n-len(bin_string)) #leftpad with 0s
        a.extend(bin_string)
        bitstrings.append(a)
    generators = find_generators(n, bitstrings)
    pauli_pairs = []
    for pair in generators:
        p1, p2 = pair
        pauli_pairs.append((string_to_pauli(p1), string_to_pauli(p2)))