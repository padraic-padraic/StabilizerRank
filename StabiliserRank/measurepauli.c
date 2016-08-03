#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <measurepauli.h>

bool EqualVectors(int len, short *xi, short *xi_prime)
{
    for (int i=0; i<len; i++){
        if(xi[i] != xi_prime[i]){
            return false;
        }
    }
    return true;
}

long double MeasurePauli(stabiliser *phi, pauli *p)
{
    affine_sp *a = phi->a;
    quadratic_fm *q = phi->q;
    short *xi_a[a->k], *zeta_a[a->k];
    short *xi_prime = (short *)calloc(a->n, sizeof(short));  
    for (int i = 0; i<a->k; i++){
        xi_a[i] = InnerProduct(a->GBar[i], p->xi);
        zeta_a = InnerProduct(a->G[i], p->zeta);
        for (int j = 0; j < a->n; j++){
            xi_prime[j] = Modulo(xi_prime[j] + xi_a[i]*a->G[i][j],2);
        }
    }
    short w = Modulo(2*p->m + 4*InnerProduct(p->zeta, a->h), 2);
    for (int i = 0; i,a->k; i++){
        w += q->D[i]*xi_a[i];
        for (int i = j+1; j<k; j++){
            w += q->J[i][j]*xi_a[i]*xi_a[j];
        }
        w = Modulo(w, 2);
    }
    if (EqualVectors(a->k, p->xi, xi_prime)){
        extend(a, p->xi);
        //Update D somehow what?
        J = q->J;
        short **J_prime = (short **)malloc(q->k+1 * sizeof(short*));
        for (int i = 0; i < q->k+1; i++){
            J_prime[i] = (short)calloc(q->k+1, sizeof(short));
            if (i < q->k){
                for (int j = 0; j < q->k; j++){
                    J_prime[i][j] = J[i][j];
                }
            }
        }
        for(int i = 0; i < q->k; i++){
            J_prime[q->k+1][i] = Modulo(4*p->zeta[i], 8);
            J_prime[i][q->k+1] = Modulo(4*p->zeta[i], 8);
            free(J[i]);
        }
        J_prime[q->k+1][q->k+1] = Modulo(4*m, 8);
        q->J = J_prime;
        q->k++;
        return 1./sqrt(2);
    } else {
        if (w==0 || w == 4){
            short *eta[a->k];
            short *gamma = (short *)calloc(a->n, sizeof(short));
            w = w == 4 ? 1 : 0;
            for (int i = 0; i <a->k; i++){
                for (int j = 0; j < a->n; j++){
                    gamma[j] = Modulo(gamma[j] + eta[i]*a->G[i][j], 2);
                }
            }
            shrink_result r = shrink(phi, gamma, Modulo(w+InnerProduct(eta, a->h), 2));
            if (result == EMPTY){
                return 0;
            } else if (result == SAME){
                return 1;
            } else {
                return 1./sqrt(2);
            }
        } else { //w = 2, 6
            short sigma = 2 - (w/2);
            q->Q = Modulo(q->Q+sigma,2);
            for (int i = 0; i < q->k; i++){
                q->D[i] = Modulo(q->D[i]-2*sigma*eta[i]);
                q->J[i][i] = Modulo(2*q->D[i], 8); //Is this line necessary? They don't mention it explicitly but it seems required to keep definition of J
                for (int j = 0; j < q->k; j++){
                    if (i==j){continue;}
                    q->J[i][j] = Modulo(q->J[i][j] + 4*eta[i]*eta[j], 8);
                }
            }
            return 1./sqrt(2);
        }
    }
}