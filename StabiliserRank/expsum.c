#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <bianryform.h>
#include <expsum.h>

complex WSum(quadratic_fm *q, int r, short **dimers, int M, short *monomers, short S)
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

void NonNullS(quadratic_fm *q, short *set_S, short *order_S)
{
    int target=-1;
    //Prepare the transformation matrix, initially as the Identity
    short **R = (short **)calloc(q->k, sizeof(short*));
    for (int i = 0; i <q->k;i++){
        R[i]=(short *)calloc(q->k, sizeof(short));
        R[i][i] = 1;
    }
    target = RandomInt(order_S);
    for (int i=0; i < *order_S; i++){
        if (i==target){continue;}
        //Set the correct element of R to 1
        R[set_S[i]][set_S[target]] = 1;
    }
    qfm_BasisChange(q, R);
    S[0] = S[target]
    for (int i =1; i<order_S; i++){S[i]=-1;}
    *order_S = 1;
}

void PartitionBasis(quadratic_fm *q, short **dimers, short *monomers, short S)
{
    short order_E = q->k;
    if (S >-1){order_E--;}
    short *E[order_E];
    short *K[order_E];
    short **R = (short *)calloc(q->k*sizeof(short*));
    for(int i =0; i <q->k; i++){
        R[i]= (short)calloc(q->k*)*sizeof(short);
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
            M[order_M]=(short)target;
            order_M +=1
            order_E = DeleteFromArray(order_E, target, E);
        }
        else{
            target_2 = K[(int)((double)rand() / ((double)RAND_MAX + 1) * order_K)];
            for (int i=0; i<order_E; i++){
                if (i==target || i == target_2){continue;}
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
    short **dimers;
    dimers = (short **)malloc(q->k*sizeof(short *));
    for (int i = 0; i < q->k; i++){dimers[i] = (short*)malloc(2*sizeof(short));}
    short *monomers[q->k];
    int r=0, M = 0;
    short order_S = 0;
    short *set_S = (short *)malloc(q->k*sizeof(int *));
    for (int i=0;i<q->k;i++){S[i]=-1;} //negative values are falsy
    // Populate S
    for (short i=0; i < q->k; i++){
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