!pb  18.12.06: COMPUTATION TIME PER PROCESSOR IS SET TO THE MAXIMUM TIME
!pb            THAT IS AVAILABLE
C
      SUBROUTINE PEDIST (XTIM,XX1)
C
C   SUBROUTINE PEDIST CALCULATES THE DISTRIBUTION OF PROCESSORS TO
C   STRATA IN CASE THAT THERE ARE MORE PE'S THAN STRATA
C   DISTRIBUTION OF PE'S IS DONE ACCORDING TO THE DISTRIBUTION OF
C   COMPUTATION TIME.
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPES
      USE COMSOU
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: XTIM(0:NSTRA)
      REAL(DP), INTENT(IN) :: XX1
      REAL(DP) :: TIMPE(0:NSTRA)
      REAL(DP) :: FACP, DELT, SUMTIM, TMEAN
      INTEGER :: IPE, K, I, ISTRA, NPRS_FREE, NPRS_OPT

! calculate mean cpu time per stratum
      sumtim=xtim(nstrai)
      TMEAN=SUMTIM/FLOAT(NPRS)

      WRITE (iunout,*) ' SUMTIM = ',SUMTIM,' MEAN TIME = ',TMEAN

      NPRS_OPT=0
      NPRS_FREE=NPRS
      DO ISTRA=1,NSTRAI
        delt=xtim(istra)-xtim(istra-1)
        IF (delt.GE.1.E-5) THEN
! a stratum that has got computation time gets at least 1 processor
          NPESTR(ISTRA)=1
          NPRS_FREE=NPRS_FREE-1
        ELSE
          NPESTR(ISTRA)=0
        ENDIF
! calculate the optimal number of additional processors according to
! distribution of cpu time done in mcarlo (according to number of particles
! and source strength specified in the input)
        TIMPE(ISTRA)=MAX(delt-TMEAN,1.E-5_DP)/TMEAN
        NPRS_OPT=NPRS_OPT+int(TIMPE(ISTRA))
      ENDDO
      WRITE (iunout,*) ' ISTRA, TIMPE '
      DO ISTRA=1,NSTRAI
        WRITE (iunout,*) ISTRA,TIMPE(ISTRA)
      ENDDO

      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! distribute free processors to strata by their optimal number of processors
      FACP=MIN(1.D0,REAL(NPRS_FREE,KIND(1.D0))/
     .             (REAL(NPRS_OPT,KIND(1.D0))+eps30))
      write (iunout,*) ' facp ',facp
      NPESTR(0)=NPRS
      DO ISTRA=1,NSTRAI
        NPESTR(ISTRA)=NPESTR(ISTRA)+int(TIMPE(ISTRA)*FACP)
        NPRS_FREE=NPRS_FREE-int(TIMPE(ISTRA)*FACP)
      ENDDO
      WRITE (iunout,*) ' NPESTR ',(NPESTR(ISTRA),ISTRA=1,NSTRAI)
      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! if there are still free processors left distribute them to all
! strata with more than tmean cpu time assigned to them using a
! daisy chain mechanism
      ISTRA=0
      DO WHILE (NPRS_FREE.GT.0)
        ISTRA=ISTRA+1
        IF (ISTRA.GT.NSTRAI) ISTRA=1
        IF (TIMPE(ISTRA).GT.1.E-10) THEN
          NPESTR(ISTRA)=NPESTR(ISTRA)+1
          NPRS_FREE=NPRS_FREE-1
        ENDIF
      ENDDO
      WRITE (iunout,*) ' NPESTR '
      WRITE (iunout,'(12I6)') (NPESTR(ISTRA),ISTRA=1,NSTRAI)
      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! assign each processor the number of the stratum it shall work on
      IPE=0
      DO ISTRA=1,NSTRAI
        DO K=1,NPESTR(ISTRA)
          IPE=IPE+1
          NSTRPE(IPE-1)=ISTRA
        ENDDO
      ENDDO
      WRITE (iunout,*) ' IPE, ISTRA '
      WRITE (iunout,'(12I6)') (I,NSTRPE(I),I=0,NPRS-1)

! for each stratum define the number of the first processor
! that does calculations for this stratum
! this is used to determine the groups of processors in the
! accumulation of the results for one stratum
      NPESTA(0)=0
      NPESTA(1)=0
      DO ISTRA=2,NSTRAI
        NPESTA(ISTRA)=NPESTA(ISTRA-1)+NPESTR(ISTRA-1)
      ENDDO
      WRITE (iunout,*) ' NPESTA '
      WRITE (iunout,'(12I6)') (NPESTA(I),I=0,NSTRAI)

      XTIM(1:NSTRAI) = XX1
      CALL MASAGE ('REDEFINED CPU TIME ASSIGNED TO STRATA (SEC) :')
      DO ISTRA=1,NSTRAI
        CALL MASJ1R ('STRATUM, TIME   ',ISTRA,XTIM(ISTRA))
      END DO      

      RETURN
      END
