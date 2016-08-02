#ifndef RANDOMSTAB_H
#define RANDOMSTAB_H
#include <bainaryform.h>

long double Eta(int d, int n);

long double Random(void);

int SampleD(int n);

unsigned short ** RandomMatrix(intt d, int n);

unsigned short ** GenerateX(int d, int n);

unsigned short * RandomShiftVector(int n);

unsigned short * RandomD(int n);

unsigned short ** RandomJ(int n, unsigned short *D); 

stabiliser random_stab(int n);

#endif