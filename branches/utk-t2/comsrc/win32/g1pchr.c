#include <stdio.h>

int g1pchr(icommand,icharacter,filename)

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
