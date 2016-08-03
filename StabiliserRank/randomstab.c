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

short ** RandomMatrix(int d, int n)
{
    long double r;
    short **X = (short **)malloc(d, sizeof(short *));
    for (int i = 0; i<d; i++){
        X[i] = (short *)calloc(n, sizeof(short));
        for (int j = 0; j<n; j++){
            X[i][j] = r<0.5 ? 0 : 1;
        }
    }
    return X;
}

int Rank(short **X) // Based on methods discussed here https://groups.google.com/forum/?hl=en#!topic/comp.lang.c/6afkseBMcnk as part of the diehard RNG tests
{
    bool found_something;
    int r = 0;
    for (int i = 0; i < d; i++){
        found_something = false;
        for (int j = r; j<n; j++){
            if (X[j][d-i] == 1){
                if (found_something == false){
                    found_something == true;
                    for (int k = 0; k<n; k++){X[j][i] = X[r][i];}
                    r++;
                }
            }
        }
    }
    return r;
}

short ** GenerateX(int d, int n)
{
    while (true){
        short **X = RandomMatrix(d, n);
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

short * RandomShiftVector(int n)
{
    long double r;
    short *h = (short *)calloc(n, sizeof(short));
    for (int i = 0; i <n; i++){
        r = Random();
        h[i] = r < 0.5 ? 0 : 1;
    }
    return h;
}

short * RandomD(int n)
{
    short *D = (short *)calloc(n, sizeof(short));
    for (int i = 0; i < n; i++){
        D[i] = Modulo(2* RandInt(4), 8);
    }
    return D;
}

short ** RandomJ(int n, short* D)
{
    short **J = (short **)malloc(n, sizeof(short *));
    for (int i =0; i<n; i++){
        J[i] = (short *)calloc(n, sizeof(short));
        J[i][i] = Modulo(2*D[i],8);
    }
    //TODO off diagonals
    return J;
}

stabiliser random_stab(int n)
{
    int d = SampleD(n);
    int k = n-d;
    short **X = GenerateX(d, n);
    affine_sp a;
    a.n = n;
    a.k = k;
    a.h = (short*)calloc(a.n, sizeof(short));
    a.G = (short **)malloc(a.n*sizeof(short *));
    a.GBar = (short **)malloc(a.n*sizeof(short *));
    for (int i = 0; i<a.n; i++){
        a.G[i] = (short *)calloc(a.n, sizeof(short));
        a.GBar[i] = (short *)calloc(a.n, sizeof(short));
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