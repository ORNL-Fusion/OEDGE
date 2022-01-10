int g4putc(array,iplace,nbitsc,nchrpw,nchar)

char array[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  array[*iplace-1]= *nchar;
  return;
}
