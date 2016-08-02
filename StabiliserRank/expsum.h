#ifndef EXPSUM_H
#define EXPSUM_H
#include <complex.h>
#include <binaryform.h>
complex WSum(quadratic_fm *q, int r, unsigned short **dimers, int M, unsigned short *monomers, short S);

void NonNullS(quadratic_fm *q, unsigned short *set_S, unsigned short *order_S);

void PartitionBasis(quadratic_fm *q, unsigned short **dimers, unsigned short *monomers, short S);

complex ExponentialSum(quadratic_fm *q);
#endif