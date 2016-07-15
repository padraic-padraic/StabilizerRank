#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

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

// void OuterProduct(int len, unsigned short *left, unsigned short *right, unsigned short **target)
// {
//     for (int i = 0; i < len; i++){
//         for (int j = 0; j < len; j++){
//             target[i][j] = 
//         }
//     }
//     return;
// }

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

void qfm_ShiftChange(quadratic_fm *q, unsigned short *y)
{
    qfm->Q += Modulo(InnerProduct(qfm->D,y), 8);
    qfm->Q += Modulo(QShiftChange(q->k, q->J, y), 8);
}

complex WSum(quadratic_fm *q, int r, unsigned short **dimers, int M, unsigned short *monomers, short S)
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
    if (2*r + M < q->k){ //S is defaulted to 0 so this serves as a truth test
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

void NonNullS(quadratic_fm *q, unsigned short *set_S, unsigned short *order_S){
    int target=-1;
    //Prepare the transformation matrix, initially as the Identity
    unsigned short **R = (unsigned short **)calloc(q->k, sizeof(unsigned short*));
    for (int i = 0; i <q->k;i++){
        R[i]=(unsigned short *)calloc(q->k, sizeof(unsigned short));
        R[i][i] = 1;
    }
    target = (int)((double)rand() / ((double)RAND_MAX + 1) * order_S) //Snippet taken from http://c-faq.com/lib/randrange.html
    for (int i=0; i < order_S; i++){
        if (i==target){continue;}
        //Set the correct element of R to 1
        R[set_S[i]][set_S[target]] = 1;
    }
    qfm_BasisChange(q, R);
    S[0] = S[target]
    for (int i =1; i<order_S; i++){S[i]=-1;}
    *order_S = 1;
}

int DeleteFromArray(int len; int target, unsigned *arr){
    for (int i = target; i<len-1; i++){
        arr[i]=arr[i+1];
    }
    arr[len] = -1;
    return len-1;
}

void PartitionBasis(quadratic_fm *q, unsigned short **dimers, unsigned short *monomers, short S){
    short order_E = q->k;
    if (S >-1){order_E--;}
    short *E[order_E];
    short *K[order_E];
    unsigned short **R = (unsigned short *)calloc(q->k*sizeof(unsigned short*));
    for(int i =0; i <q->k; i++){
        R[i]= (unsigned short)calloc(q->k*)*sizeof(unsigned short);
        R[i][i]=1;
    }
    int count = 0, target, target_2, k_count=0, order_M=0, r=0;
    for (short i=0; i<q->k;i++){
        if (i!=S){
            E[count]=i;
            count++;
        }
    }
    while (order_E > 0){
        target = (int)((double)rand() / ((double)RAND_MAX + 1) * order_E);
        for(int i = 0; i<order_E; i++){
            if(i==target){continue;}
            if(q->J[target][E[i]] ==4){
                K[k_count]=E[i];
                k_count++;
            }
        }
        if (k_count==0){
            M[order_M]=(unsigned short)target;
            order_M +=1
            order_E = DeleteFromArray(order_E, target, E);
        }
        else{
            target_2 = K[(int)((double)rand() / ((double)RAND_MAX + 1) * order_K)];
            for (int i=0; i<order_E; i++){
                if (i==target || i ==target_2){continue;}
                //Fill R matrix - TODO
            }
            qfm_BasisChange(q, R);
            dimers[r][0] = E[target];
            dimers[r][1] = E[target_2];
            r++;
            for (int i=0; i<q->k; i++){ //Reset R to identity
                for (int j=0; j<q->k; j++){
                    if (i==j){R[i][j]=1;}
                    else {R[i][j]=0;}
                }
            }
            //Delete a and b from E
            if (target<target_2){
                order_E = DeleteFromArray(order_E, target_2, E);
                order_E = DeleteFromArray(order_E, target, E);
            } else {
                order_E = DeleteFromArray(order_E, target, E);
                order_E = DeleteFromArray(order_E, target_2, E);
            }
        }
        k_count=0;
        for(int i = 0; i <q->k-1; i++){K[i]=-1;}// reset k
    }
    for(int i=0; i<q->k;i++){free(R[i]);}
    free(R);
}

complex ExponentialSum(quadratic_fm *q)
{
    //Init structures to store dimers, monomers
    unsigned short **dimers;
    dimers = (unsigned short **)malloc(q->k*sizeof(unsigned short *));
    for (int i = 0; i < q->k; i++){dimers[i] = (unsigned short*)malloc(2*sizeof(unsigned short));}
    unsigned short *monomers[q->k];
    int r=0, M = 0;
    unsigned short order_S = 0;
    // Calloc ensures memory junk doesn't get thrown in to WSum
    short *set_S = (short *)malloc(q->k*sizeof(int *));
    for (int i=0;i<q->k;i++){S[i]=-1;} //negative values are falsy
    // Populate S
    for (unsigned short i=0; i < q->k; i++){
        if (q->D[i] == 2 || q->D[i] == 6){
            set_S[order_S] = i;
            order_S++;   
        } 
    }
    if (order_S > 1){
        // We need to shrink S down to 1 element and do the associated basis change
        NonNullS(q, set_S, &order_S);
    }
    //Find the dimers and monomers
    PartitionBasis(q, dimers, monomers, S[0]);
    //Do WSum
    complex result = WSum(q, r, M, dimers, monomers, S[0]);
    //Free explicitly allocated arrays
    for (int i = 0; i<len; i++){
        free(dimers[i]);
    }
    free(dimers); 
    free(set_S);
    return result;
}

int main(){
    //Do nothing
    char a = 1;
    char b = 2;
    printf("%d\n",a^b);//Should equal 3
}