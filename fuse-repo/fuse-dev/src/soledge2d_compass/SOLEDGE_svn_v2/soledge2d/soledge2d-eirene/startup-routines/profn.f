C sept.05:  activate this radial profile type also for levgeo=3
C           input radius SEP is relative to flux suface labelling grid RHOSRF
C           rra gt. raa indicates: one more vaccum cell in radial direction.
 
C
      SUBROUTINE EIRENE_PROFN(PRO,PRO0,PROS,P,Q,E,SEP,PROVAC)
C
C  PARABOLIC PROFILE, PLUS EXPONENTIAL DECAY BEYOND RHOSEP
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: PRO0, PROS, P, Q, E, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      REAL(DP) :: EP1R, ELLR, YR, EIRENE_AREAA, ARCA, RHOSEP, RORP, 
     .            FACT, ZR, RDI, AR
      INTEGER  :: J, JM1, EIRENE_LEARCA, NLOCAL
 
      IF (LEVGEO.LE.3) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING CO-ORDINATE AT "SEP"
C
        IF (LEVGEO.EQ.2) THEN
          JM1=EIRENE_LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFN       ')
          AR=EIRENE_AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
        ELSEIF (LEVGEO.EQ.1) THEN
          RHOSEP=SEP
        ELSEIF (LEVGEO.EQ.3) THEN
          RHOSEP=SEP
        ENDIF
      ELSE
        WRITE (iunout,*)
     .  'WARNING: SUBR. PROFN CALLED EIRENE_WITH LEVGEO.GT.3 '
        WRITE (iunout,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      RDI=1./(RHOSEP-RHOSRF(1))
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHOSEP) THEN
          ZR=RHOZNE(J)-RHOSRF(1)
          RORP=ZR*RDI
          IF (RORP.LE.0.D0) THEN
            PRO(J)=PRO0
          ELSE
            FACT=(1.0-RORP**P)**Q
            PRO(J)=(PRO0-PROS)*FACT+PROS
          ENDIF
        ELSE
          PRO(J)=PROS*EXP((RHOSEP-RHOZNE(J))/E)
        ENDIF
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
