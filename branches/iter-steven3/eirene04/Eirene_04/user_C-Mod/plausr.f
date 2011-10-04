

      SUBROUTINE PLAUSR
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CSTEP
      USE CGRID
      USE CINIT
      USE COMSOU
      USE CTRIG
      IMPLICIT NONE
      REAL(DP) :: FACTOR, STEP
      INTEGER :: NLINES, ITRI, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, 
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
      INTEGER :: IDEZ
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, DELR, FL
      CHARACTER(100) :: ZEILE

      ZEILE='*   '      
      DO WHILE (ZEILE(1:1) == '*')
         READ (31,'(A100)',END=99) ZEILE
      END DO

      READ (ZEILE,*) NLINES

      CALL ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))
      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0

      DO I=1, NLINES
        READ (31,*) ITRI, ISIDE, FLX, TE, TI, DE
        STRATA:DO ISTRA=1,NSTRAI
          IF (.NOT.NLSRF(ISTRA)) CYCLE
          DO ISRFS=1,NSRFSI(ISTRA)
            ISOR=SORLIM(ISRFS,ISTRA)
            ITEC1=IDEZ(ISOR,1,4)
            ITEC2=IDEZ(ISOR,2,4)
            ITEC3=IDEZ(ISOR,3,4)
            IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) CYCLE
            ISTEP=SORIND(ISRFS,ISTRA)
            IF (ISTEP.EQ.0) THEN
              WRITE (6,*) 'ERROR IN PRIMARY SOURCE DATA '
              WRITE (6,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
              WRITE (6,*) 'NO. ', INSOR(ISRFS,ISTRA),' BUT SORIND.EQ.0.'
              CALL EXIT(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
              CALL MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
              CALL EXIT(1)
            ENDIF
            INDSRF=INSOR(ISRFS,ISTRA)
            IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)
            IF (INMTI(ISIDE,ITRI) == INDSRF) THEN
              IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0._DP
              KSTEP(ISTEP) = KSTEP(ISTEP) + 1
              INOSRC(ISTEP) = ISTRA
              IS1 = ISIDE + 1
              IF (IS1.GT.3) IS1=1
              IF (NSPEZ(ISTRA) <= 0) THEN
                IPLAN(ISTEP)=1
                IPLEN(ISTEP)=NPLSI
              ELSE
                IPLAN(ISTEP)=NSPEZ(ISTRA)
                IPLEN(ISTEP)=NSPEZ(ISTRA)
              END IF
              IRSTEP(ISTEP,KSTEP(ISTEP))=ITRI
              IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
              ITSTEP(ISTEP,KSTEP(ISTEP))=1
              IASTEP(ISTEP,KSTEP(ISTEP))=0
              IBSTEP(ISTEP,KSTEP(ISTEP))=1
              DELR =  SQRT(
     .           (XTRIAN(NECKE(ISIDE,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .           (YTRIAN(NECKE(ISIDE,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
              RRSTEP(ISTEP,KSTEP(ISTEP)+1)=
     .           RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR         
              TESTEP(ISTEP,KSTEP(ISTEP)) = TE    
              DO IPLS=IPLAN(ISTEP), IPLEN(ISTEP)
                TISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = TI
                DISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = DE
                VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITRI)
                VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITRI)
                VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITRI)
                FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
              END DO
              EXIT STRATA
            END IF
          END DO 
        END DO STRATA
      END DO

      DO ISTEP = 1, NSTEP
        IF (KSTEP(ISTEP) > 0) THEN
          NBIN=KSTEP(ISTEP)+1
          FL=STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
          FLUX(INOSRC(ISTEP))=FL
        END IF
      END DO

      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)

 99   RETURN
      END
