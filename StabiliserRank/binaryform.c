#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <binaryform.h>

int DeleteFromArray(int len; int target, unsigned *arr)
{
    for (int i = target; i<len-1; i++){
        arr[i]=arr[i+1];
    }
    arr[len] = -1;
    return len-1;
}

unsigned short Modulo(unsigned short val, unsigned short base)
{
    return val & base;
}
unsigned short InnerProduct(int len, unsigned short *x, unsigned short *y)
{
    short res = 0;
    for (int i = 0; i < len; i++){
        res += x[i] & y[i];
    }
    return res;
}

void LeftMultiply(int len, unsigned short **R, unsigned short *x, unsigned short *target, unsigned short base)
{
    unsigned short res[len]
    for (int i = 0; i < len; i++){
        res[i] = Modulo(InnerProduct(len, R[i], x), base);
    }
    for (int i = 0; i < len; i++){target[i] = res[i];}
}
void transpose(int len, unsigned short **mat, unsigned short **target)
{
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            target[i][j] = mat[j][i];
        }
    }
}
void MatrixMultiply(int len, unsigned short **left, unsigned short **right, unsigned short ** target)
{
    unsigned short **trans;
    transpose(len, right, trans);
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            target[i][j] = Modulo(InnerProduct(left[i], trans[j]), 8)
        }
    }
}

unsigned short DBasisChange(int len, int index, unsigned short **J, unsigned short **R)
{
    unsigned short res = 0;
    for (int i=0; i < len-1; i++){
        for (int j=i ; j < len; j++){
            res += J[i][j]*R[index][i]*R[index][j];
        }
    }
    return res;
}

unsigned short QShiftChange(int len, unsigned short **J, unsigned short *y)
{
    unsigned short res = 0;
    for (int i = 0; i < len; i++){
        for (int j=i ; j < len; j++){
            res += J[i][j]*y[i]*y[j];
        }
    }
    return res;
}

void qfm_BasisChange(quadratic_fm *q, unsigned short **R){
    LeftMultiply(q->k, R, q->D, q->D);
    for (int i=0; i < len; i++){
        q->D[i] += Modulo(DBasisChange(len, i, q->J, R), 8);
    }
}