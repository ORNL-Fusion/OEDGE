//int g4getb(iword,iposn,ibit)
int G4GETB(iword,iposn,ibit)

register unsigned *iword,*iposn,*ibit;
{
  *ibit= ((*iword >> (32-*iposn))&1);
  return 0;
}
