#ifndef EXPSUM_H
#define EXPSUM_H
#include <complex.h>
#include <binaryform.h>
complex WSum(quadratic_fm *q, int r, short **dimers, int M, short *monomers, short S);

void NonNullS(quadratic_fm *q, short *set_S, short *order_S);

void PartitionBasis(quadratic_fm *q, short **dimers, short *monomers, short S);

complex ExponentialSum(quadratic_fm *q);
#endif