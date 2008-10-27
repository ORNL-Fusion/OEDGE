//int g4putk(iarray,iplace,nbitsc,nchrpw,nchar)
int G4PUTK(iarray,iplace,nbitsc,nchrpw,nchar)

char iarray[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  iarray[*iplace-1]= *nchar;
  return 0;
}
