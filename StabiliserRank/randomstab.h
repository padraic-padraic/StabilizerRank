#ifndef RANDOMSTAB_H
#define RANDOMSTAB_H
#include <bainaryform.h>

long double Eta(int d, int n);

long double Random(void);

int SampleD(int n);

short ** RandomMatrix(intt d, int n);

short ** GenerateX(int d, int n);

int Rank(short **X);

short * RandomShiftVector(int n);

short * RandomD(int n);

short ** RandomJ(int n, short *D); 

stabiliser random_stab(int n);

#endif