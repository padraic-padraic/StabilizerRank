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

typedef {EMPTY, SAME, SUCCESS} result;

int DeleteFromArray(int len; int target, unsigned *arr);

unsigned short Modulo(unsigned short val, unsigned short base);

unsigned short InnerProduct(int len, unsigned short *x, unsigned short *y);

void LeftMultiply(int len, unsigned short **R, unsigned short *x, unsigned short *target, unsigned short base);

void transpose(int len, unsigned short **mat, unsigned short **target);

void MatrixMultiply(int len, unsigned short **left, unsigned short **right, unsigned short ** target);

int RandomInt(int max);

unsigned short DBasisChange(int len, int index, unsigned short **J, unsigned short **R);

unsigned short qfm_ShiftChange(int len, unsigned short **J, unsigned short *y);

void qfm_ShiftChange(quadratic_fm *q, unsigned short *y);

void qfm_DeleteRow(quadratic_fm *q, int target);

void afp_AddVectors(unsigned short *v1, unsigned short *v2);

void shrink(affine_sp *a, quadratic_fm *q, unsigned short *xi, unsigned short alpha);

void lazy_shrink(affine_sp *a, unsigned short *xi, unsigned short alpha);

#endif

