#ifndef MEASUREPAULI_H
#define MEASUREPAULI_H
#include <bianryform.h>
#include <complex.h>

typedef struct pauli {
    int n;
    short m;
    short *zeta;
    short *xi;
} pauli;

long double MeasurePauli(stabiliser *phi, pauli *p);

#endif