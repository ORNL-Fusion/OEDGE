//int g4getk(iarray,iplace,nbitsc,nchrpw,nchar)
int G4GETK(iarray,iplace,nbitsc,nchrpw,nchar)

char iarray[];
register int *iplace,*nchar;
int nbitsc,nchrpw;

{
  *nchar= iarray[*iplace-1]&255;
  return 0;
}
