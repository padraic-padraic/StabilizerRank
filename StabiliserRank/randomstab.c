#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <binaryform.h>

long double Eta(int d, int n)
{
    if (d == 0){return 1;}
    long double res = pow(2, -1*d*(d+1)/2);
    for (int a = 0; a < d; a++){
        res *= (1-pow(2, d-n-a))/ (1-pow(2, -1*a));
    }
    return res;
}

long double Random(void)
{
    return (long double)rand()/(long double) RAND_MAX;
}

int SampleD(int n)
{
    long double sum = 0;
    for (int i = 0; i < n+1; i++){sum+=Eta(i)}
    long double cumulative[n+1];
    for (int i = 0; i < n+1; i++){
        if (i==0){cumulative[i]=Eta(0)/sum;}
        else{cumulative[i] = cumulative[i-1]+(Eta(i)/sum);}
    }
    long dobule p = Random();
    for (int i =0; i < n+1; i++){
        if (p < cumulative[i]){
            return i;
        }
    }
}

unsigned short ** RandomMatrix(int d, int n)
{
    long double r;
    unsigned short **X = (unsigned short **)malloc(d, sizeof(unsigned short *));
    for (int i = 0; i<d; i++){
        X[i] = (unsigned short *)calloc(n, sizeof(unsigned short));
        for (int j = 0; j<n; j++){
            X[i][j] = r<0.5 ? 0 : 1;
        }
    }
    return X;
}

unsigned short ** GenerateX(int d, int n)
{
    while (true){
        unsigned short **X = RandomMatrix(d, n);
        if (Rank(X) == d) {
            return X;
        } else {
            for(int i = 0; i<d;i++){
                free(X[i]);
            }
            free(X);
        }
    }
}

unsigned short * RandomShiftVector(int n)
{
    long double r;
    unsigned short *h = (unsigned short *)calloc(n, sizeof(unsigned short));
    for (int i = 0; i <n; i++){
        r = Random();
        h[i] = r < 0.5 ? 0 : 1;
    }
    return h;
}

unsigned short * RandomD(int n)
{
    unsigned short *D = (unsigned short *)calloc(n, sizeof(unsigned short));
    for (int i = 0; i < n; i++){
        D[i] = Modulo(2* RandInt(4), 8);
    }
    return D;
}

unsigned short ** RandomJ(int n, unsigned short* D)
{
    unsigned short **J = (unsigned short **)malloc(n, sizeof(unsigned short *));
    for (int i =0; i<n; i++){
        J[i] = (unsigned short *)calloc(n, sizeof(unsigned short));
        J[i][i] = Modulo(2*D[i],8);
    }
    return J;
}

stabiliser random_stab(int n)
{
    int d = SampleD(n);
    int k = n-d;
    unsigned short **X = GenerateX(d, n);
    affine_sp a;
    a.n = n;
    a.k = k;
    a.h = (unsigned short*)calloc(a.n, sizeof(unsigned short));
    a.G = (unsigned short **)malloc(a.n*sizeof(unsigned short *));
    a.GBar = (unsigned short **)malloc(a.n*sizeof(unsigned short *));
    for (int i = 0; i<a.n; i++){
        a.G[i] = (unsigned short *)calloc(a.n, sizeof(unsigned short));
        a.GBar[i] = (unsigned short *)calloc(a.n, sizeof(unsigned short));
        a.G[i][i] = 1;
        a.Gbar[i][i] = 1;
    }
    for (int a = 0; a <d; a++){
        lazy_shrink(a, X[a], 0);
    }
    a.h = RandomShiftVector(a.n);
    a.Q = Modulo(RanomInt(8), 8);
    a.D = RandomD(a.n);
    a.J = RandomJ(a.n, a.D);
    //Build struct and return
    stabiliser phi;
    phi.a = &a
    phi.q = &q
    return phi;
}