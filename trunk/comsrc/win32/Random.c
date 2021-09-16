// Random.c

#include "Random.h"
double ran2(long *idum)
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if(*idum <= 0)
	{
		*idum = (-(*idum) < 1) ? 1 : -(*idum);
		idum2= (*idum);

		for(j = (NTAB + 7); j >= 0; j--)
		{
			k = (*idum) / R2IQ1;
			*idum = R2IA1 * (*idum - (k * R2IQ1)) - k * R2IR1;

			if(*idum < 0) 
				*idum += R2IM1;

			if(j < NTAB)
				iv[j] = *idum;
		}

		iy = iv[0];
	}

	k = (long)((*idum) * R2IQ1INV);
	*idum = R2IA1 * (*idum - k * R2IQ1) - k * R2IR1;
	if(*idum < 0) 
		*idum += R2IM1;

	k = (long)(idum2 * R2IQ2INV);
	idum2 = R2IA2 * (idum2 - k * R2IQ2) - k * R2IR2;
	if(idum2 < 0) 
		idum2 += R2IM2;

	j = (long)(iy * NDIVINV);
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if(iy < 1) 
		iy += IMM1;

	if((temp = AM*iy) > 1.0) 
		return 1.0;
	else
		return temp;
}