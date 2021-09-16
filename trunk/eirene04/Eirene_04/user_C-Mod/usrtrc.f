      SUBROUTINE USRTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      USE CTRIG
      USE COMPRT
      USE CUPD
      IMPLICIT NONE

      INTEGER  IFLAG,ISYM,count
      REAL(DP) XPLO,YPLO,ZPLO

      DATA count /0/
     
      IF (IFLAG.EQ.0) COUNT = COUNT + 1

      WRITE(90,'(I8,3F9.3,3I4,I4,I7,1P,1E10.2,0P,F8.4,I6,7I5,2X,5I5,
     .  1P,E10.2,0P)') 
     .  count,XPLO,YPLO,ZPLO,IFLAG,ISYM,istra,0,
c     .  count,XPLO,YPLO,ZPLO,IFLAG,ISYM,istra,ntrseg,
     .  npanu,e0,weight,nacell,ifpath,iupdte,
     .  ityp,nblock,masurf,msurf,mtsurf,
     .  nrcell,npcell,nblock,nntcll,ntcell,
     .  time

      

      RETURN
      END
