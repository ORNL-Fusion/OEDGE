//int g4putb(iword,iposn,ibit)
int G4PUTB(iword,iposn,ibit)

register unsigned *iword,*iposn,*ibit;
{
  if (*iposn!=1)
    *iword= (((~0<<(33-*iposn))|((1<<(32-*iposn))-1))&*iword)|(*ibit<<(32-*iposn));
  else
    *iword= (((1<<(32-*iposn))-1)&*iword)|(*ibit<<(32-*iposn));
  return 0;
}
