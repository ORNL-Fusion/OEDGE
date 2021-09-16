// Random.h : Contains inline randum number generation functions

#ifndef __RANDOM_H__
#define __RANDOM_H__

#define R2IM1		2147483563L
#define R2IM2		2147483399L
#define AM			(1.0/R2IM1)
#define IMM1		(R2IM1 - 1)
#define R2IA1		40014
#define R2IA2		40692
#define R2IQ1		53668
#define R2IQ1INV	(1.0 / R2IQ1)
#define R2IQ2		52774
#define R2IQ2INV	(1.0 / R2IQ2)
#define R2IR1		12211
#define R2IR2		3791
#define NTAB		32
#define NDIV		(1 + IMM1/NTAB)
#define NDIVINV		(1.0 / NDIV)
#define RAN2LOWBOUND 1E-7
#define RAN2SEED	314159265358979323

double ran2(long *idum);

#endif //__RANDOM_H__