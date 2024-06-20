      SUBROUTINE EIRENE_OUTLST
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CLAST
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      REAL(DP) :: SMMEAN
      INTEGER :: IWR, IRCX, IREL, IRPI
 
      call EIRENE_leer(1)
      iwr=0
      do ircx=1,nrcxi
        if (iflrcx(ircx).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOCX '
            write (iunout,*)
     .        'IRCX, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLSEIRENE_'
            iwr=1
          endif
          SMMEAN=xcmean(ircx)/(ncmean(ircx)+eps60)
          write (iunout,*) ircx,SMMEAN,NCMEAN(IRCX)
        endif
      enddo
      call EIRENE_leer(1)
      iwr=0
      do irel=1,nreli
        if (iflrel(irel).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOEL '
            write (iunout,*)
     .        'IREL, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLSEIRENE_'
            iwr=1
          endif
          SMMEAN=xemean(irel)/(nemean(irel)+eps60)
          write (iunout,*) irel,SMMEAN,NEMEAN(IREL)
        endif
      enddo
      call EIRENE_leer(1)
      iwr=0
      do irpi=1,nrpii
        if (iflrpi(irpi).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOPI '
            write (iunout,*)
     .        'IRPI, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLSEIRENE_'
            iwr=1
          endif
          SMMEAN=xpmean(irpi)/(npmean(irpi)+eps60)
          write (iunout,*) irpi,SMMEAN,NPMEAN(IRPI)
        endif
      enddo
      call EIRENE_leer(1)
      RETURN
      END
