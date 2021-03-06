#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <binaryform.h>
#include <expsum.h>
#include <stabilizerproduct.h>

complex stb_InnerProduct(stabiliser *phi1, stabiliser *phi2)
{
    affine_sp a = afp_Copy(phi1->a); //This won't work. I need to define a copy method
    stabiliser stab;
    stab.a = &a;
    stab.q = phi1->q;
    affine_sp *a1 = phi1->a, *a2 = phi2->a;
    quadratic_fm *q1 = phi1->q, *q2 = phi2->q;
    for (int i = a2->k +1; i < a2->n; i++){
        alpha = InnerProduct(a2->h, a2->GBar[i])
        res = shrink(&stab, a2->HBar[i], alpha);
        if (res == EMPTY){
            return (complex) 0;
        }
    }
    short *y = (short *)calloc(a2->k, sizeof(short));
    short *scratch_space = (short *)calloc(a2->n, sizeof(short));
    short res = 0;
    AddVectors(a2->n, scratch_space, a1->h);
    AddVectors(a2->n, scratch_space, a2->h);
    short **R = (short **)calloc(a2->n, sizeof(short*));
    for (int i=0; i<a2->n; i++){R[i] = (short *)calloc(a2->n, sizeof(short));}
    for (int i = 0; i < a2->k; i++){
        y[i] = InnerProduct(scratch_space, a2->Gbar[i]);
        for (int j = 0; j < a.k; j++){
            R[j][i] = InnerProduct(a.G[j], a2.GBar[i]);
        }

    }
    qfm_BasisChange(q2, R);
    qfm_ShiftChange(q2, y);
    for (int i = 0; i<a2->n; i++){a2->h[i] = a1->h[i];}
    //Cleanup
    for(int i = 0; i<q2->k; i++){free(R[i]);}
    free(scratch_space);
    free (R);
    free(y);
    //Find the final quadratic form q
    q.Q = Modulo(q1->Q - q2->Q, 8);
    for (int i = 0; i < q.K; i++){ 
        q.D[i] = Modulo(q1->D[i]-q2->D[i], 8);
        for (int j = 0; j<q.k; j++){
            q.J[i][j] = Modulo(q1->J[i][j] - q2->J[i][j], 8);
        }
    }
    return pow(2, -1*(a1->k - a2->k)/2)*ExponentialSum(q);
}