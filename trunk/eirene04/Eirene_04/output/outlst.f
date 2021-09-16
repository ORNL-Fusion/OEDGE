C
C
C
C
C
C        ************
C        *  OUTPUT  *
C        ************
C
C      SUBROUTINE OUTLST
C      SUBROUTINE OUTPLA
C      SUBROUTINE OUTEIR (INDOUT)
C      SUBROUTINE OUTFLX (INDOUT)
C      SUBROUTINE WRREC
C      SUBROUTINE WRSNAP
C      SUBROUTINE GETSCL (ISTRA,FATM,FMOL,FION)
C
C
      SUBROUTINE OUTLST

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLAST
      USE COMXS

      IMPLICIT NONE

      REAL(DP) :: SMMEAN
      INTEGER :: IWR, IRCX, IREL

      call leer(1)
      iwr=0
      do ircx=1,nrcxi
        if (iflrcx(ircx).gt.0) then
          if (iwr.eq.0) then
            WRITE (6,*) 'REJECTION SAMPLING EFFICIENCY IN VELOCX '
            write (6,*) 'IRCX, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLS'
            iwr=1
          endif
          SMMEAN=xcmean(ircx)/(ncmean(ircx)+eps60)
          write (6,*) ircx,SMMEAN,NCMEAN(IRCX)
        endif
      enddo
      call leer(1)
      iwr=0
      do irel=1,nreli
        if (iflrel(irel).gt.0) then
          if (iwr.eq.0) then
            WRITE (6,*) 'REJECTION SAMPLING EFFICIENCY IN VELOEL '
            write (6,*) 'IREL, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLS'
            iwr=1
          endif
          SMMEAN=xemean(irel)/(nemean(irel)+eps60)
          write (6,*) irel,SMMEAN,NEMEAN(IREL)
        endif
      enddo
      call leer(1)
      RETURN
      END
