#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
void piserial_(void)
{
    double niter = 100000;
    double x,y;
    int i;
    int count=0;
    double z;
    double pi;
    srand(time(NULL));
    //main loop
    for (i=0; i<niter; ++i)
    {
        //get random points
        x = (double)random()/RAND_MAX;
        y = (double)random()/RAND_MAX;
        z = sqrt((x*x)+(y*y));
        //check to see if point is in unit circle
        if (z<=1)
        {
            ++count;
        }
    }
    //p = 4(m/n)
    pi = ((double)count/(double)niter)*4.0;
    printf("Pi serial: %f\n", pi);
}
