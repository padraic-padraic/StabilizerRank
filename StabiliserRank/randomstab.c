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

unsigned short ** RandomMatrix(unsigned short d, int n)
{
    return;
}

stabiliser random_stab(int n)
{
    int d = SampleD(n);
    int k = n-d;
    unsigned short **X = RandomMatrix(d, n)
    //Build struct and return
    stabiliser phi;
    phi.a = &a
    phi.q = &q
    return phi;
}