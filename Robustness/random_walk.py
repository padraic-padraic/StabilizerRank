#! /usr/bin/env python
import argparse
import multiprocessing
import numpy as np

from bitarray import bitarray
from cmath import exp as cexp #Differentiate from real-valued exp
from copy import deepcopy
from math import acos, cos, exp, sin, sqrt
from math import pow as fpow #Differentiate from stdlib pow
from random import randrange, random
from qutip import basis, sigmaz, tensor
from Utils import do_anneal, format_for_output, write_string

BETA = acos(1/sqrt(3)) /2
PI = acos(0.)*2
H = (cos(PI/8)*basis(2,0) + sin(PI/8)*basis(2,1)).unit()
F = (cos(BETA)*basis(2,0) + cexp(1j*PI/4)*sin(BETA)*basis(2,1)).unit()
PLUS = (basis(2,0)+basis(2,1)).unit()
S = sigmaz().sqrtm()
Tmat = S.sqrtm()
T = Tmat * PLUS
RT = Tmat.sqrtm() * PLUS
RRT = Tmat.sqrtm().sqrtm() * PLUS

STATES = {'T':T,
          'RT':RT,
          'RRT':RRT,
          'F':F,
          'H':H}

def run_analysis(state, n, chi):
    job = {'ostring': write_string,
                   'fname':"/Users/padraic/Desktop/"+state+".txt",
                   'n_qubits':n,
                   'target_string':state,
                   'func':do_anneal}
    func_inputs = {'target':tensor([STATES[state]]*n),
                       'n_qubits':n}
    job['func_inputs'] = func_inputs
    job['func_inputs']['chi'] = chi
    res = format_for_output(job)
    res.write()

parse = argparse.ArgumentParser('Dispatch the simulated annealer.')
parse.add_argument('state_string', type=str)
parse.add_argument('n_qubits', type=int)
parse.add_argument('chi', type=int)

if __name__ == '__main__':
    args = parse.parse_args()
    run_analysis(args.state_string, args.n_qubits, args.chi)

