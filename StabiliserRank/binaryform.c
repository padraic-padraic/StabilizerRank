#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <binaryform.h>

int DeleteFromArray(int len; int target, *arr)
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
    copy.D = (short *)malloc(copy.k*sizeof(short));
    copy.J = (short **)malloc(copy.k*sizeof(short *));
    for (int i = 0; i<copy.k; i++){
        copy.D[i] = q->D[i];
        copy.J[i] = (short*)malloc(copy.k*sizeof(short));
        for (int j = 0; j<copy.k; j++){
            copy.J[i][j] = q->J[i][j];
        }
    }
    return copy;
}


affine_sp afp_Copy(affine_sp *a){
    affine_sp copy;
    copy.n = a->n;
    copy.k = a->k;
    copy.h = (short *)malloc(copy.n*sizeof(short));
    copy.G = (short **)malloc(copy.n*sizeof(short *));
    copy.GBar = (short **)malloc(copy.n*sizeof(short *));
    for (int i = 0; i < copy.n; i++){
        copy.h[i] = a->h[i];
        copy.G[i] = (short *)malloc(copy.n*sizeof(short));
        copy.GBar[i] = (short *)malloc(copy.n*sizeof(short));
        for (int j=0; i < copy.n; j++){
            copy.G[i][j] = a->G[i][j];
            copy.GBar[i][j] = a->GBar[i][j]
        }
    }
    return copy;
}

short Modulo(short val, short base)
{
    short res = val&base;
    res = res < 0 ? -1*res : res;
    return res;
}
short InnerProduct(int len, short *x, short *y)
{
    short res = 0;
    for (int i = 0; i < len; i++){
        res += x[i] & y[i];
    }
    return res;
}

void LeftMultiply(int len, short **R, short *x, short *target, short base)
{
    short res[len]
    for (int i = 0; i < len; i++){
        res[i] = Modulo(InnerProduct(len, R[i], x), base);
    }
    for (int i = 0; i < len; i++){target[i] = res[i];}
}

short ** transpose(int len, short **mat)
{   short **res = (short **)malloc(len*sizeof(short*));
    for (int i = 0; i <len; i++){res[i] = (short *)malloc(len*sizeof(short));}
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            res[i][j] = mat[j][i];
        }
    }
    return res;
}
void MatrixMultiply(int len, short **left, short **right, short ** target)
{
    short **trans;
    trans = transpose(len, right);
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            target[i][j] = Modulo(InnerProduct(left[i], trans[j]), 8)
        }
    }
} // Double check this? I think it works? Clever!

short DBasisChange(int len, int index, short **J, short **R)
{
    short res = 0;
    for (int i=0; i < len-1; i++){
        for (int j=i ; j < len; j++){
            res += J[i][j]*R[index][i]*R[index][j];
        }
    }
    return res;
}


void qfm_ShiftChange(quadratic_fm *q, short *y)
{
    short res = 0;
    for (int i = 0; i<q->k; i++){res += D[i]*y[i];}
    for (int i = 0; i<q->k; i++){
        for (int j=i; j<q->k; j++){
            res += q->J[i][j]*y[i]*y[j];
        }
    }
    q->Q = Modulo(q->Q+res, 8);
    for (int i=0; i<q->k; i++){
        res = 0;
        for(int j = 0; j< q->k; j++){
            res += q->J[i][j]*y[j];
        }
        q->D[i] = Modulo(q->D[i]+res, 8);
    }
}

void qfm_BasisChange(quadratic_fm *q, short **R)
{
    LeftMultiply(q->k, R, q->D, q->D, 8);
    for (int i=0; i < len; i++){
        q->D[i] += Modulo(DBasisChange(len, i, q->J, R), 8);
    }
    short **res = (short**)calloc(q->k, sizeof(short *))p;
    for (int i = 0; i <q->k; i++){res[i] = (short *)calloc(q->k, sizeof(short));}
    MatrixMultiply(q->k, q->J, transpose(R), res);
    MatrixMultiply(q->k, R, res, q->J);
    for (int i = 0; i<q->k; i++){
        free(res[i]);
    }
    free(res);
}

void qfm_DeleteIndex(quadratic_fm *q, int target) 
{
    for (int i = target; i < q->k-1; i++){
        q->D[i] = q->D[i+1];
        for (int j = 0; j < q->k- 1; j++){
            q->J[i][j] = q->j[i+1][j];
        }
    }
    for (int i = 0; i<q->k; i++){
        for (int j = target; j<q->k-1; j++){q->J[i][j] = q->J[i][j+1];}
    }
    q->D[k] = 0;
    for (int i = 0; i < q->k; i++){q->J[i][q->k] = 0;}
    free(q->J[q->k]);
    q->k -= 1;
}

void AddVectors(int len, short *v1, short *v2)
{
    for (int i = 0; i < len; i++){
        v1[i] = Modulo(v1[i] + v2[i], 2);
    }
}

void SwapVectors(affine_sp *a, int target1, int target2)
{
    short scratch_space;
    for (int i = 0; i < a->n; i++){
        scratch_space = a->G[target1][i];
        a->G[target1][i] = a->G[target2][i];
        a->G[target2][i] = scratch_space;
        scratch_space = a->GBar[target1][i];
        a->GBar[target1][i] = a->GBar[target2][i]
        a->GBar[target2][i] = scratch_space;
    }
}

int RandomIntInRange(short *len)
{    
    int target = (int)((double)rand() / ((double)RAND_MAX + 1) * (*order_S)); //Snippet taken from http://c-faq.com/lib/randrange.html
    return target;
}

shrink_result shrink(stabiliser *phi, short *xi, short alpha)
{
    affine_sp *a = phi->a;
    quadratic_fm *q = phi->q;
    short S[a->k];
    short order_S = 0;
    for (int i = 0; i < a->k; i++){
        if (Modulo(InnerProduct(a->k, xi, a->G[i]), 2) == alpha){
            S[order_S] = i;
            order_S++;
        }
    }
    short beta = Modulo(InnerProduct(a->k, xi, a->h), 2);
    if (S[0] == 0 && beta == 1){
        return EMPTY;
    } else if (S[0] == 0 && beta ==1){
        return SAME;
    } else {
        int target = RandomInt(&order_S);
        for (int i = 0; i < order_S; i++){
            if (i==target){continue;}
            AddVectors(a->G[i], a->G[target]);
            AddVectors(a->GBar[target], a->GBar[i]);
        }
        short **R = (short **)calloc(q->k, sizeof(short*));
        for (int i = 0; i <q->k;i++){
            R[i]=(short *)calloc(q->k, sizeof(short));
            R[i][i] = 1;
        }
        for (int i=0; i < order_S; i++){
            if (i==target){continue;}
            //Set the correct element of R to 1
            R[S[i]][S[target]] = 1;
        }
        qfm_BasisChange(q, R);
    }
    SwapVectors(a, target, a->k);
    //Reset the transformation matrix
        for (int i =0; i < q->k; i++){
            for (int j=0; j < q->k; j++){R[i][j] = 0;}
    }
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
        AddVectors(a->h, a->G[a->k]);
        qfm_ShiftChange(q, a->G[a->k]);
    }
    a->k-=1;
    for (int i = 0; i<q->k; i++){free(R[i]);}
    free(R);
    return SUCCESS;
}

shrink_result lazy_shrink(stabiliser *phi, short *xi, short alpha)
{
    affine_sp *a = phi->a;
    short S[a->k];
    short order_S = 0;
    for (int i = 0; i < a->k; i++){
        if (Modulo(InnerProduct(phi->a.k, xi, a->G[i]), 2) == alpha){
            S[order_S] = i;
            order_S++;
        }
    }
    short beta = Modulo(InnerProduct(phi->a.k, xi, a->h), 2);
    if (S[0] == 0 && beta == 1){
        return EMPTY;
    } else if (S[0] == 0 && beta ==1){
        return SAME;
    } else {
        int target = RandomInt(&order_S);
        for (int i = 0; i < order_S; i++){
            if (i==target){continue;}
            AddVectors(a->G[i], a->G[target]);
            AddVectors(a->GBar[target], a->GBar[i]);
        }
    }
    SwapVectors(a, target, a->k);
    if (beta != 0){
        AddVectors(a->h, a->G[a->k]);
    }
    a->k-=1;
    return SUCCESS;
}

void extend(affine_sp *a, short *xi)
{
    short *S[a->n];
    for (int i = 0; i<n; i++){S[i]=-1;}
    int order_S = 0, order_T = 0;
    for (int i = 0; i < a->n; i++){
        if (InnerProduct(xi, a->GBar[i]) == 1){
            S[order_S++] = i;//Increments order_S and passes the old value to index S
        }
    }
    short *T[order_S];
    for (int i = 0; i<n; i++){T[i]=-1;}
    for (int i = a->k; i<n; i++){
        for (int j = 0; j < order_S; j++){
            if (S[j] == i){
                T[order_T++] = i;
                break;
            }
        }
    }
    if (T[0] == -1){
        return; //Extend has failed
    }
    short i = T[RandInt(order_T)];
    short **R = (short**)malloc(a->n, sizeof(short *));
    for (int i = 0; i<a->n; i++){
        R[i]=(short *)calloc(a->n, sizeof(short));
        R[i][i] = 1;
    }
    for (int j = 0; j<order_S; j++){
        if(S[j] == i){continue;}
        AddVectors(a->GBar[j], a->GBar[i]);
    }
    a->G[i] = xi;
    SwapVectors(a, i, a->k+1);
    a->k+=1;
    return;
}