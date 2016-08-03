#ifndef BINARYFORM_H
#define BINARYFORM_H

typedef struct affine_sp {
    int n;
    int k;
    short *h;
    short **G;
    short **GBar;
} affine_sp;

typedef struct quadratic_fm {
    int k;
    short Q;
    short *D;
    short **K;
} quadratic_fm;

typedef struct stabiliser {
    affine_sp *a;
    quadratic_fm *q;
}

typedef enum {EMPTY, SAME, SUCCESS} shrink_result;

int DeleteFromArray(int len; int target, *arr);

quadratic_fm qfm_Copy(quadratic_fm *q);

affine_sp afp_Copy(affine_sp *a);

short Modulo(short val, short base);

short InnerProduct(int len, short *x, short *y);

void LeftMultiply(int len, short **R, short *x, short *target, short base);

short ** transpose(int len, short **mat, short **target);

void MatrixMultiply(int len, short **left, short **right, short ** target);

int RandomInt(int max);

short DBasisChange(int len, int index, short **J, short **R);

void qfm_ShiftChange(quadratic_fm *q, short *y);

void qfm_BasisChange(quadratic_fm *q, short **R);

void qfm_DeleteRow(quadratic_fm *q, int target);

void AddVectors(short *v1, short *v2);

void SwapVectors(affine_sp *a, int target);

shrink_result shrink(stabiliser *phi, short *xi, short alpha);

shrink_result lazy_shrink(stabiliser *phi, short *xi, short alpha);

void extend(affine_sp *a, short *xi);

#endif

