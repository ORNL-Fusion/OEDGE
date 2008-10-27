//int g4getc(array,iplace,nbitsc,nchrpw,nchar)
int G4GETC(array,iplace,nbitsc,nchrpw,nchar)

char array[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  *nchar= array[*iplace-1]&255;
  return 0;
}
