#include <time.h>
#include <stdio.h>
//int enqtim(chtime)
int ENQTIM(chtime)
char chtime[];
{
  int i;
  time_t ltime;
  char* timchr;

  time(&ltime);
  timchr= ctime(&ltime);

  for (i= 0;i<= 7;i++)
    chtime[i]= *(timchr+i+11);
  return 0;
}
