#include <time.h>
#include <stdio.h>
int enqdat(chdate)
char chdate[];
{
  int i;
  long clock;
  char *timchr;

  clock= time(NULL);
  timchr= ctime(&clock);

  chdate[0]= timchr[8];
  chdate[1]= timchr[9];
  chdate[2]= '/';
  chdate[3]= timchr[4];
  chdate[4]= timchr[5]-32;
  chdate[5]= timchr[6]-32;
  chdate[6]= '/';
  chdate[7]= timchr[22];
  chdate[8]= timchr[23];
  return;
}
int enqtim_(chtime)
char chtime[];
{
  int i;
  long clock;
  char *timchr;

  clock= time(NULL);
  timchr= ctime(&clock);

  for (i= 0;i<= 7;i++)
    chtime[i]= *(timchr+i+11);
  return;
}
#include <stdio.h>

int g1pchr_(icommand,icharacter,filename)

char *filename;
int *icommand, *icharacter;

{
  static FILE *fp;

  switch (*icommand)
  {
     case 1:
     fp= fopen(filename,"w+");
     break;

     case 2:
     fclose(fp);
     break;

     case 3:
     putc(*icharacter,fp);
     break;

     case 4:
     *icharacter= getc(fp);
     break;

     default:
     fflush(fp);
  }
  return;
}
int g4getb_(iword,iposn,ibit)

register unsigned *iword,*iposn,*ibit;
{
  *ibit= ((*iword >> (32-*iposn))&1);
  return;
}
int g4getc_(array,iplace,nbitsc,nchrpw,nchar)

char array[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  *nchar= array[*iplace-1]&255;
  return;
}
int g4getk_(iarray,iplace,nbitsc,nchrpw,nchar)

char iarray[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  *nchar= iarray[*iplace-1]&255;
  return;
}
int g4putb_(iword,iposn,ibit)

register unsigned *iword,*iposn,*ibit;
{
  if (*iposn!=1)
    *iword= (((~0<<(33-*iposn))|((1<<(32-*iposn))-1))&*iword)|(*ibit<<(32-*iposn));
  else
    *iword= (((1<<(32-*iposn))-1)&*iword)|(*ibit<<(32-*iposn));
  return;
}
int g4putc_(array,iplace,nbitsc,nchrpw,nchar)

char array[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  array[*iplace-1]= *nchar;
  return;
}
int g4putk_(iarray,iplace,nbitsc,nchrpw,nchar)

char iarray[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  iarray[*iplace-1]= *nchar;
  return;
}













