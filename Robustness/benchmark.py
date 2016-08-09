import numpy as np
import qutip as qt

from cmath import exp as cexp
from itertools import product, combinations
from math import acos, cos, sin, sqrt
from Utils import SL0, n_stab

## Magic state definitions
pi = 2*acos(0.)
beta = acos(1/sqrt(3)) / 2
I, X, Y, Z = qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()
H = (cos(pi/8)*qt.basis(2,0)+sin(pi/8)*qt.basis(2,1)).unit()
F = (cos(beta)*qt.basis(2,0) + cexp(1j*pi/4)*sin(beta)*qt.basis(2,1)).unit()
plus = (qt.basis(2,0)+qt.basis(2,1)).unit()
S = Z.sqrtm()
T = S.sqrtm()
magic_states = [H, S*H, Z*H, S.dag()*H, X*H, S*X*H, Z*X*H, S.dag()*X*H, T*plus,
                S*T*plus, S.dag()*T.dag()*plus,T.dag()*plus]
face_states = [F]# S*F, Z*F, S.dag()*F, X*F, S*X*F, S*Z*F, S.dag()*X*F]

##

## Helpful funcitons
def commutator(A,B):
    return A*B - B*A

def test_equality(g1,g2):
    return all([any([np.allclose(_g1, _g2) for _g2 in g2]) for _g1 in g1])

def projector(group):
    I = group[0]
    return ((I+group[1])/2) * ((I+group[2])/2)
##

print('Generating two qubit operators')
base_elements = [I,X,Y,Z]
two_qubit_ops = []
for pair in product(base_elements, base_elements):
    two_qubit_ops.append(np.matrix(qt.tensor(pair[0], pair[1]).full(), 
                                   dtype=np.complex_))
two_qubit_ops = two_qubit_ops[1:] #Truncate the II operator
for n in range(len(two_qubit_ops)):
    two_qubit_ops.append(-1*two_qubit_ops[n])
print("Found " + str(len(two_qubit_ops)))
null = np.zeros((4,4), dtype=np.complex_)
print('Generating stabiliser groups')
gen_pairs = []
for n, op in enumerate(two_qubit_ops):
    for i in range(len(two_qubit_ops)):
        if i == n:
            continue
        if np.allclose(-1*op,two_qubit_ops[i]):
            continue
        if np.allclose(commutator(op, two_qubit_ops[i]) , null):
            gen_pairs.append((op, two_qubit_ops[i]))
stab_groups = []
for pair in gen_pairs:
    candidate = [np.eye(4), pair[0], pair[1], pair[0]*pair[1]]
    if stab_groups:
        if any([test_equality(candidate, group) for group in stab_groups]):
            continue
        else:
            stab_groups.append(candidate)
    else:
        stab_groups.append(candidate)
print("Found " + str(len(stab_groups)))
print("Generating stabiliser states")
stab_projectors = []
stab_states = []
for group in stab_groups:
    stab_projectors.append(qt.Qobj(projector(group), dims=[[2,2],[2,2]]))
    eigs, vecs = np.linalg.eig(projector(group))
    for i in range(eigs.size):
        if np.allclose(eigs[i], complex(1)):
            stab_states.append(qt.Qobj(vecs[:,i], dims=[[2],[1]]))

basis = [np.matrix('1;0;0;0'), np.matrix('0;1;0;0'), 
         np.matrix('0;0;1;0'), np.matrix('0;0;0;1')]
A = np.matrix(np.zeros((4, 60), dtype=np.complex_))
for i in range(4):
    for j, state in enumerate(stab_states):
        A[i, j] = np.sum(basis[i].H * np.matrix(state.full()))
b = np.matrix(np.zeros((4,1), dtype=np.complex_))

norm_0 = np.zeros(500)
truthy = []
print("Benchmarking 500 random states")
for i in range(500):
    state = qt.rand_ket(4, dims=[[2,2],[1,1]])
    for j in range(4):
        b[j] = np.sum(basis[j].H*np.matrix(state.full()))
    norm_0[i] = np.count_nonzero(SL0(A,b, 1e-12))
print("Average L0 norm is " + str(np.mean(norm_0)))
norm_0 = np.zeros(len(magic_states))
for i, state in enumerate(magic_states):
    for j in range(4):
        b[j] = np.sum(basis[j].H*qt.tensor([state,state]).full())
    norm_0[i] = np.count_nonzero(SL0(A,b,1e-12))
print("Norm for edge type states is "+str(np.mean(norm_0)))
norm_0 = np.zeros(len(face_states))
for i, state in enumerate(face_states):
    for j in range(4):
        b[j] = np.sum(basis[j].H*qt.tensor([state,state]).full())
    norm_0[i] = np.count_nonzero(SL0(A,b,1e-12))
print("Norm for face type states is "+str(np.mean(norm_0)))
