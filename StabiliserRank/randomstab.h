#ifndef RANDOMSTAB_H
#define RANDOMSTAB_H
#include <bainaryform.h>

long double Eta(int d, int n);

long double Random(void);

int SampleD(int n);

unsigned short ** RandomMatrix(unsigned short d, int n);

stabiliser random_stab(int n);

#endif