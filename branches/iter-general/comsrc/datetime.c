/* C Subroutines to act as interface for FORTRAN to the 
   C run-time library */

#include <stdio.h>
#include <time.h>
#include <sys/types.h>

int datetest(char *string)
{
time_t tloc;
struct tm *tres;
int  day,month,year;

time(&tloc);

tres = localtime(&tloc);

day = (*tres).tm_mday ;
month = (*tres).tm_mon+1;
year = (*tres).tm_year;

sprintf(string,"%2.2d/%2.2d/%2.2d",month,day,year);

return(1);


}

int clocktest(char *string)
{
time_t tloc;
time(&tloc);
strncpy(string,ctime(&tloc)+11,8);
}

int newsrand(int *seed)
{
char *tmp;
int  tmp2;

static char state[256];

/*tmp = initstate(*seed,state,256);*/

tmp2 = initstate(*seed,state,256);

}


float newrand()
{
static float ran;

ran = (float) (( (double) random() )/(2.147483647e9));

return(ran);

}




