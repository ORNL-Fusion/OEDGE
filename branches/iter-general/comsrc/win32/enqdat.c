#include <time.h>
#include <stdio.h>
//int enqdat(chdate)
int ENQDAT(chdate)
char chdate[];
{  
  time_t ltime;
  char* timchr;

  time(&ltime);
  timchr= ctime(&ltime);

  chdate[0]= timchr[8];
  chdate[1]= timchr[9];
  chdate[2]= '/';
  chdate[3]= timchr[4];
  chdate[4]= timchr[5]-32;
  chdate[5]= timchr[6]-32;
  chdate[6]= '/';
  chdate[7]= timchr[22];
  chdate[8]= timchr[23];
  return 0;
}
