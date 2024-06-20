#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//int pi_openacc(int argc, char* argv[])
extern "C" {
        int piopenacc_(void);
        }

int piopenacc_(void)
{
        double niter = 100000;
        double x,y;
        int count=0;
        double z;
        double pi;
        int a = 0;
        time_t t;
        double randomnums[200000];
        double output[100000];
        srand48(((unsigned)time(&t)));
        for (a = 0; a<=2*niter; a++)
                randomnums[a] = (double)drand48();
        int i = 0;
        #pragma acc parallel copy(randomnums[0:200000]) copyout(output[0:100000]) create(x,y,z)
        {
                for (i = 0; i<100000; i++)
                {
                        int temp = i+i;
                        x = randomnums[temp];
                        y = randomnums[temp+1];
                        z = ((x*x)+(y*y));
                        output[i] = z;
                }
                for (i = 0; i<100000; i++)
                {
                        if (output[i]<=1)
                        {
                                ++count;
                        }
                        else
                                continue;
                }
        }
        pi = ((double)count/(double)niter)*4.0;
        printf("Pi openacc: %f\n", pi);

	return 0;
}
