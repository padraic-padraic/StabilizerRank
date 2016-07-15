#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
//Big caveat: 
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

void OuterProduct

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

unsigned short DBasisChange(int len, int index, unsigned short **J, unsigned short **R){
    unsigned short res = 0;
    for (int i=0; i < len; i++){
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

void qfm_ShiftChange(quadratic_fm *q, unsigned short *y)
{
    qfm->Q += Modulo(InnerProduct(qfm->D,y), 8);
    qfm->Q += Modulo(QShiftChange(q->k, q->J, y), 8);
}

complex WSum(quadratic_fm *q, int r, unsigned short **dimers, int M, unsigned short *monomers, int S)
{
     //Exp of Q
    complex res= cexp(I*M_PI/4*q->Q); 
    //Product over Monomers
    for (int i=0; i<M; i++){res*= (1+cexp(I*M_PI/4*q->D[monomers[i]]));}
    //Product over Dimers
    for (int i=0; i<r; i++){
        res *= (1 + cexp(I*M_PI/4*q->D[dimers[i][0]])
                  + cexp(I*M_PI/4*q->D[dimers[i][1]])
                  - cexp(I*M_PI/4*(q->D[dimers[i][0]]+q->D[dimers[i][1]])));
    }
    if (2*r + M != q->k){
        complex res2= cexp(I*M_PI/4*(q->Q+q->D[S])); 
        //Product over Monomers
        for (int i=0; i<M; i++){
            res*= (1+cexp(I*M_PI/4*
                            (q->D[monomers[i]] +
                             q->J[monomers[i]][S])
                         )
                  );
        }
        //Product over Dimers
        for (int i=0; i<r; i++){
            res *= (1 + cexp(I*M_PI/4*
                                (q->D[dimers[i]][0] +
                                 q->J[dimers[i][0]][S])
                            )
                      + cexp(I*M_PI/4*
                                (q->D[dimers[i][1]]) + 
                                 q->J[dimers[i][1]][S]
                            )
                      - cexp(I*M_PI/4*
                                (q->J[dimers[i][0]][S] +
                                 q->J[dimers[i][1]][S] +
                                 q->D[dimers[i]][0] +
                                 q->D[dimers[i]][1])
                            )
                    );
        }
        return res+res2;
    } else {
        return res;        
    }
}

complex ExponentialSum(quadratic_fm *q)
{
    return WSum(q, r, M, dimers, monomers, S)
}

int main(){
    //Do nothing
    char a = 1;
    char b = 2;
    printf("%d\n",a^b);//Should equal 3
}