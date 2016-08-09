from math import pow

__all__ = ['A', 'G', 'n_stab']

def n_stab(n):
    res = pow(2.,n)
    for i in range(n):
        res *= (pow(2.,n-i)+1)
    return res
def A(n):
    res = 1
    for i in range(n):
        res *= (pow(2.,n)-pow(2.,i))
    return res
def G(n):
    res = pow(2.,n)
    for i in range(n):
        res *= (pow(4.,n)/pow(2.,i) - pow(2.,i))
    return res