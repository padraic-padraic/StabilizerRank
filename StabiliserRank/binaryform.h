#ifndef BINARYFORM_H
#define BINARYFORM_H

typedef struct affine_sp {
    int n;
    int k;
    unsigned short *h;
    unsigned short **G;
    unsigned short **GBar;
} affine_sp;

typedef struct quadratic_fm {
    int k;
    unsigned short Q;
    unsigned short *D;
    unsigned short **K;
} quadratic_fm;

int DeleteFromArray(int len; int target, unsigned *arr);

unsigned short Modulo(unsigned short val, unsigned short base);

unsigned short InnerProduct(int len, unsigned short *x, unsigned short *y);

void LeftMultiply(int len, unsigned short **R, unsigned short *x, unsigned short *target, unsigned short base);

void transpose(int len, unsigned short **mat, unsigned short **target);

void MatrixMultiply(int len, unsigned short **left, unsigned short **right, unsigned short ** target);

unsigned short DBasisChange(int len, int index, unsigned short **J, unsigned short **R);

unsigned short QShiftChange(int len, unsigned short **J, unsigned short *y);

void qfm_ShiftChange(quadratic_fm *q, unsigned short *y);


#endif

