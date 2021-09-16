/* C Subroutines to act as interface for FORTRAN to the 
   C run-time library */


#include <stdio.h>
#include <time.h>
#include <sys/types.h>


void datetest_(char *string)
{
    time_t tloc;
    struct tm *tres;
    int  day,month,year;
    char buf[8];

    time(&tloc);

    tres = localtime(&tloc);

    day = (*tres).tm_mday ;
    month = (*tres).tm_mon+1;
    year = (*tres).tm_year;
    if (year>99) { year = year -100; };

    sprintf(buf,"%2.2d/%2.2d/%2.2d",month,day,year);
    /* IPP/01 - Krieger: sprintf on writing in a string appends an ascii
       zero character to the string. To avoid memory spill over to
       other variables when called from Fortran with a non-zero terminated
       string variable, use strncopy */
    strncpy(string,buf,8);
}


void clocktest_(char *string)
{
    time_t tloc;
    time(&tloc);
    strncpy(string,ctime(&tloc)+11,8);
}


char *initstate();

static long state[64]={
               0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342,
               0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
               0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86,
               0xda672e2a, 0x1588ca88, 0xe369735d, 0x904f35f7,
               0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
               0xde3b81e0, 0xdf0a6fb5, 0xf103bc02, 0x48f340fb,
               0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b,
               0xf5ad9d0e, 0x8999220b, 0x27fb47b9, 0x8a88d77c,
               0x9a219039, 0x32d9d024, 0x9b653182, 0x4da1f242,
               0x7e49554b, 0xb1b1ebb0, 0xdb5c5917, 0x956553fd,
               0x7c2e581a, 0xebad756f, 0xb314e0c7, 0x24433b16,
               0xd367ee2a, 0x1cc8c481, 0xe3637657, 0x912fe537,
               0xd21b8ad1, 0x6ff6f150, 0x611e639b, 0xa2943fac,
               0xde3384e0, 0xaf0a61b5, 0xf133bc00, 0x4af34cfb,
               0x35413793, 0xc122c99a, 0xf5a4bab1, 0x8a48d67b,
               0xf5ed9df1, 0x8949210c, 0x2afb4eb1, 0x8a78d561
               };
/* IPP/01 - Krieger: fixed bus error, which occured because the state
   array was not initialized in newsrand. Now state is initialized
   properly outside newsrand */

void newsrand_(int *seed)
{
    initstate((unsigned)*seed,state,256);
}
