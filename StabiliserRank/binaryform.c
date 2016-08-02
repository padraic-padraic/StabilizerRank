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
quadratic_fm qfm_Copy(quadratic_fm *q){
    quadratic_fm copy;
    copy.k = q->k;
    copy.Q = q->Q;
    copy.D = (unsigned short *)malloc(copy.k*sizeof(unsigned short));
    copy.J = (unsigned short **)malloc(copy.k*sizeof(unsigned short *));
    for (int i = 0; i<copy.k; i++){
        copy.D[i] = q->D[i];
        copy.J[i] = (unsigned short*)malloc(copy.k*sizeof(unsigned short));
        for (int j = 0; j<copy.k; j++){
            copy.J[i][j] = q->J[i][j];
        }
    }
    return copy
}


affine_sp afp_Copy(affine_sp *a){
    affine_sp copy;
    copy.n = a->n;
    copy.k = a->k;
    copy.h = (unsigned short *)malloc(copy.n*sizeof(unsigned short));
    copy.G = (unsigned short **)malloc(copy.n*sizeof(unsigned short *));
    copy.GBar = (unsigned short **)malloc(copy.n*sizeof(unsigned short *));
    for (int i = 0; i < copy.n; i++){
        copy.h[i] = a->h[i];
        copy.G[i] = (unsigned short *)malloc(copy.n*sizeof(unsigned short));
        copy.GBar[i] = (unsigned short *)malloc(copy.n*sizeof(unsigned short));
        for (int j=0; i < copy.n; j++){
            copy.G[i][j] = a->G[i][j];
            copy.GBar[i][j] = a->GBar[i][j]
        }
    }
    return copy;
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

unsigned short qfm_ShiftChange(int len, unsigned short **J, unsigned short *y)
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

void qfm_DeleteIndex(quadratic_fm *q, int target)
{
    for (int i = target; i < q->k -1; i++){
        q->D[i] = q->D[i+1];
    }
    free(q->D[k]);
    for (int i = 0; i < q->k; i++){free(q->J[i][q->k]);}
    free(q->J[q->k]);
    q->k -= 1;
}

void afp_AddVectors(int len, unsigned short *v1, unsigned short *v2)
{
    for (int i = 0; i < len; i++){
        v1[i] = Modulo(v1[i] + v2[i], 2);
    }
}

void afp_SwapVectors(affine_sp *a, int target)
{
    unsigned short scratch_space;
    for (int i = 0; i < a->n; i++){
        scratch_space = a->G[target][i];
        a->G[target][i] = a->G[a->k][i];
        a->G[a->k][i] = scratch_space;
        scratch_space = a->GBar[target][i];
        a->GBar[target][i] = a->GBar[a->k][i]
        a->GBar[a->k][i] = scratch_space;
    }
}

int RandomInt(unsigned short *len)
{    
    int target = (int)((double)rand() / ((double)RAND_MAX + 1) * (*order_S)); //Snippet taken from http://c-faq.com/lib/randrange.html
    return target;
}

result shrink(stabiliser *phi, unsigned short *xi, unsigned short alpha)
{
    affine_sp *a = phi->a;
    quadratic_fm *q = phi->q;
    unsigned short S[a->k];
    unsigned short order_S = 0;
    for (int i = 0; i < a->k; i++){
        if (Modulo(InnerProduct(a->k, xi, a->G[i]), 2) == alpha){
            S[order_S] = i;
            order_S++;
        }
    }
    unsigned short beta = Modulo(InnerProduct(a->k, xi, a->h), 2);
    if (S[0] == 0 && beta == 1){
        return EMPTY;
    } else if (S[0] == 0 && beta ==1){
        return SAME;
    } else {
        int target = RandomInt(&order_S);
        for (int i = 0; i < order_S; i++){
            if (i==target){continue;}
            afp_AddVectors(a->G[i], a->G[target]);
            afp_AddVectors(a->GBar[target], a->GBar[i]);
        }
        unsigned short **R = (unsigned short **)calloc(q->k, sizeof(unsigned short*));
        for (int i = 0; i <q->k;i++){
            R[i]=(unsigned short *)calloc(q->k, sizeof(unsigned short));
            R[i][i] = 1;
        }
        for (int i=0; i < order_S; i++){
            if (i==target){continue;}
            //Set the correct element of R to 1
            R[S[i]][S[target]] = 1;
        }
        qfm_BasisChange(q, R);
        //Reset the transformation matrix
        for (int i =0; i < q->k; i++){
            for (int j=0; j < q->k; j++){R[i][j] = 0;}
        }
    }
    afp_SwapVectors(a, target);
    for (int i = 0; i < q->k; i++){
        //Build basis swap matrix
        if (i==target){
            R[i][q->k] = 1;
        } else if (i == q->k=1){
            R[i][target] = 1;
        } else {
            R[i][i] = 1;
        }
    }
    qfm_BasisChange(q, R);
    if (beta != 0){
        afp_AddVectors(a->h, a->G[a->k]);
        qfm_ShiftChange(q, a->G[a->k]);
    }
    a->k-=1;
    return SUCCESS;
}

result lazy_shrink(stabiliser *phi, unsigned short *xi, unsigned short alpha)
{
    affine_sp *a = phi->a;
    unsigned short S[a->k];
    unsigned short order_S = 0;
    for (int i = 0; i < a->k; i++){
        if (Modulo(InnerProduct(phi->a.k, xi, a->G[i]), 2) == alpha){
            S[order_S] = i;
            order_S++;
        }
    }
    unsigned short beta = Modulo(InnerProduct(phi->a.k, xi, a->h), 2);
    if (S[0] == 0 && beta == 1){
        return EMPTY;
    } else if (S[0] == 0 && beta ==1){
        return SAME;
    } else {
        int target = RandomInt(&order_S);
        for (int i = 0; i < order_S; i++){
            if (i==target){continue;}
            afp_AddVectors(a->G[i], a->G[target]);
            afp_AddVectors(a->GBar[target], a->GBar[i]);
        }
    }
    afp_SwapVectors(a, target);
    if (beta != 0){
        afp_AddVectors(a->h, a->G[a->k]);
    }
    a->k-=1;
    return SUCCESS;
}