C EIRENE07 COMPILATION
C ===== SOURCE: calc_spectrum.f
      SUBROUTINE CALC_SPECTRUM (WT,IND,ISC)

      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMPRT
      USE CUPD
      USE CGRID
      USE CGEOM
      USE CCONA
      USE COMUSR
c slmod begin
      USE CTRIG
c slmod end

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IND, ISC
      REAL(DP), INTENT(IN) :: WT
      INTEGER :: ISPC, I, IS, IC, IRDO, IRD, IAT, IML, IIO, IPL
      REAL(DP) :: ADD, WV, DIST, WTR, SPCVX, SPCVY, SPCVZ, CDYN, EB
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:), CNDYNM(:), CNDYNI(:), 
     .                               CNDYNP(:)
      TYPE(EIRENE_SPECTRUM), POINTER :: P

      IF ((ISC == 0) .AND. (IND .NE. 1)) RETURN

      SELECT CASE (ITYP)
      CASE (0)
        IS = IPHOT
        CDYN = 1._DP
      CASE (1)
        IS = IATM
        IF (.NOT.ALLOCATED(CNDYNA)) THEN
          ALLOCATE (CNDYNA(NATM))
          DO IAT=1,NATMI
            CNDYNA(IAT)=AMUA*RMASSA(IAT)
          END DO
        END IF
        CDYN = CNDYNA(IATM)
      CASE (2)
        IS = IMOL
        IF (.NOT.ALLOCATED(CNDYNM)) THEN
          ALLOCATE (CNDYNM(NMOL))
          DO IML=1,NMOLI
            CNDYNM(IML)=AMUA*RMASSM(IML)
          END DO
        END IF
        CDYN = CNDYNM(IMOL)
      CASE (3)
        IS = IION
        IF (.NOT.ALLOCATED(CNDYNI)) THEN
          ALLOCATE (CNDYNI(NION))
          DO IIO=1,NIONI
            CNDYNI(IIO)=AMUA*RMASSI(IIO)
          END DO
        END IF
        CDYN = CNDYNI(IION)
      CASE (4)
        IS = IPLS
        IF (.NOT.ALLOCATED(CNDYNP)) THEN
          ALLOCATE (CNDYNP(NPLS))
          DO IPL=1,NPLSI
            CNDYNP(IPL)=AMUA*RMASSP(IPL)
          END DO
        END IF
        CDYN = CNDYNP(IPLS)
      END SELECT

      IF (ISC == 0) THEN    ! SURFACE

        DO ISPC=1,NADSPC
          P => ESTIML(ISPC)%PSPC
          IF ((P%ISRFCLL == ISC) .AND.
     .        (P%ISPCSRF == MSURF) .AND.
     .        (P%IPRTYP == ITYP) .AND.
     .        ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN

            SELECT CASE(ESTIML(ISPC)%PSPC%ISPCTYP)
            CASE (1)
              ADD = WT
            CASE (2)
c slmod begin
             stop 'should not be here'
c slmod end
              ADD = WT*E0
            CASE DEFAULT
c slmod begin
             stop 'should not be here'
c sldmod end
              ADD = 0._DP
            END SELECT

            EB = E0

            IF (EB < ESTIML(ISPC)%PSPC%SPCMIN) THEN
              I = 0
            ELSEIF (EB >= ESTIML(ISPC)%PSPC%SPCMAX) THEN
              I = ESTIML(ISPC)%PSPC%NSPC + 1
            ELSE
              I = (EB - ESTIML(ISPC)%PSPC%SPCMIN) *
     .             ESTIML(ISPC)%PSPC%SPCDELI + 1
            END IF
            ESTIML(ISPC)%PSPC%SPC(I) = ESTIML(ISPC)%PSPC%SPC(I) + ADD
            ESTIML(ISPC)%PSPC%ESP_MIN= MIN(ESTIML(ISPC)%PSPC%ESP_MIN,EB)
            ESTIML(ISPC)%PSPC%ESP_MAX= MAX(ESTIML(ISPC)%PSPC%ESP_MAX,EB)
            ESTIML(ISPC)%PSPC%IMETSP = 1
c slmod begin
c            write(0,*) 'E0=',E0
            SCORE_ITYP=ITYP
            SCORE_IS=IS
            SCORE_ISPC=ISPC
            SCORE_I=I 
            SCORE_ADD=ADD
c slmod end

c slmod begin
c            write(0,'(A,3I6,2X,3I6,2X,4I6,2X,1P,E10.2,0P,2F10.3)') 
c     .              'debug: fuck!',
c     .        ESTIML(ISPC)%PSPC%ISPCTYP,i,ispc,
c     .        isc,msurf,ityp,
c     .        IPOLG,
c     .        IPOLGN,MRSURF,necke(ipolgn,mrsurf),   ! *** THESE GIVE THE TRIANGLE ***
c     .        SNGL(ESTIML(ISPC)%PSPC%SPC(I)),
c     .        SNGL(e0),SNGL(wt)
c slmod end


          END IF
        END DO

      ELSE     ! CELL
         
        WV=WEIGHT/VEL
        DO IC=1,NCOU
          DIST=CLPD(IC)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(IC)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)


          DO ISPC=1,NADSPC
            P => ESTIML(ISPC)%PSPC
            IF ((P%ISRFCLL > 0) .AND.
     .          (((P%ISRFCLL == 1).AND.(P%ISPCSRF == IRD)) .OR.      ! scoring cell
     .           ((P%ISRFCLL == 2).AND.(P%ISPCSRF == IRDO))) .AND.   ! geometry cell
     .          (P%IPRTYP == ITYP) .AND.
     .          ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN

              SELECT CASE(ESTIML(ISPC)%PSPC%ISPCTYP)
              CASE (1)
                ADD = WTR
              CASE (2)
                ADD = WTR*E0
c slmod begin
             stop 'should not be here'
c sldmod end
              CASE (3)
                ADD = WTR*VEL*CDYN
c slmod begin
             stop 'should not be here'
c sldmod end
              CASE DEFAULT
c slmod begin
             stop 'should not be here'
c sldmod end
                ADD = 0._DP
              END SELECT
c slmod begin
c... Want velocity spectrum, not energy:
              EB = VEL
c
c              EB = E0
c slmod end
              IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
                SPCVX = ESTIML(ISPC)%PSPC%SPCVX
                SPCVY = ESTIML(ISPC)%PSPC%SPCVY
                SPCVZ = ESTIML(ISPC)%PSPC%SPCVZ
                EB = EB * (SPCVX*VELX+SPCVY*VELY+SPCVZ*VELZ)
              END IF

              IF (EB < ESTIML(ISPC)%PSPC%SPCMIN) THEN
                I = 0
              ELSEIF (EB >= ESTIML(ISPC)%PSPC%SPCMAX) THEN
                I = ESTIML(ISPC)%PSPC%NSPC + 1
              ELSE
                I = (EB - ESTIML(ISPC)%PSPC%SPCMIN) *
     .               ESTIML(ISPC)%PSPC%SPCDELI + 1
              END IF
              ESTIML(ISPC)%PSPC%SPC(I) = ESTIML(ISPC)%PSPC%SPC(I) + ADD
              ESTIML(ISPC)%PSPC%ESP_MIN = 
     .               MIN(ESTIML(ISPC)%PSPC%ESP_MIN,EB)
              ESTIML(ISPC)%PSPC%ESP_MAX = 
     .               MAX(ESTIML(ISPC)%PSPC%ESP_MAX,EB)
              ESTIML(ISPC)%PSPC%IMETSP = 1

c           write(0,'(A,3I6,1P,E10.2,0P,F10.3,P,7E10.2,0P)') 
c     .                'debug: fuck!',
c     .                ESTIML(ISPC)%PSPC%ISPCTYP,i,ispc,
c     .                SNGL(ESTIML(ISPC)%PSPC%SPC(I)),
c     .                weight,
c     .                vel,wv,DIST,wtr,add,e0,eb

            END IF
          END DO
        
        END DO

      END IF

      RETURN
      END SUBROUTINE CALC_SPECTRUM
C ===== SOURCE: getscl.f
C
C
      SUBROUTINE GETSCL (ISTRA,FA,FM,FI,FPH)
C
C  FIND SCALING FACTORS TO ENFORCE PARTICLE BALANCE
C  SIMPLE VERSION: NOT SPLIT BY SPECIES, ONLY BY TYPE
C  MODIFIED JAN/95: INCLUDE SURFACE TALLIES IN MATRIX, NOT IN
C  INHOMOGENITY
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COUTAU
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), INTENT(OUT) :: FA, FM, FI, FPH
      REAL(DP) :: FC(3), P(3,4), B(3)
      REAL(DP) :: DTB1, DTB2, DTA, DETER, DTB3, FNEN
      REAL(DP) :: P11, P12, P13, P21, P22, P23, P31, P32, P33,
     .          B1, B2, B3
      INTEGER :: ICOL, IROW, I, J, I1, I2, J1, J2
      LOGICAL :: LCOLM(3), LROW(3)
C
C
      FC(1)=1.
      FC(2)=1.
      FC(3)=1.
C
C P(..,1)*FC(1)
C P(..,2)*FC(2)
C P(..,3)*FC(3)
C
      P(1,1)=PAATI(0,ISTRA)+POTATI(0,ISTRA)+PRFAAI(0,ISTRA)+
     .       PGENAI(0,ISTRA)
      P(1,2)=PMATI(0,ISTRA)+PRFMAI(0,ISTRA)
      P(1,3)=PIATI(0,ISTRA)+PRFIAI(0,ISTRA)
      P(2,1)=PAMLI(0,ISTRA)+PRFAMI(0,ISTRA)
      P(2,2)=PMMLI(0,ISTRA)+POTMLI(0,ISTRA)+PRFMMI(0,ISTRA)+
     .       PGENMI(0,ISTRA)
      P(2,3)=PIMLI(0,ISTRA)+PRFIMI(0,ISTRA)
      P(3,1)=PAIOI(0,ISTRA)+PRFAII(0,ISTRA)
      P(3,2)=PMIOI(0,ISTRA)+PRFMII(0,ISTRA)
      P(3,3)=PIIOI(0,ISTRA)+POTIOI(0,ISTRA)+PRFIII(0,ISTRA)+
     .       PGENII(0,ISTRA)
C
      B(1)=-(PPATI(0,ISTRA)+WTOTA(0,ISTRA))
      B(2)=-(PPMLI(0,ISTRA)+WTOTM(0,ISTRA))
      B(3)=-(PPIOI(0,ISTRA)+WTOTI(0,ISTRA))
C
      ICOL=0
      IROW=0
      DO 1 I=1,3
        LROW(I)=P(I,1)**2+P(I,2)**2+P(I,3)**2.GT.EPS30
        LCOLM(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2.GT.EPS30
        IF (LROW(I)) IROW=IROW+1
        IF (LCOLM(I)) ICOL=ICOL+1
1     CONTINUE
C
      IF (IROW.EQ.0) THEN
C  NO ROW IS NON ZERO, I.E. NO PARTICLES FOLLOWED
        GOTO 1000
C
      ELSEIF (IROW.EQ.1) THEN
C  ONLY ONE ROW (NO. I) IS NON ZERO, I.E., ONLY ATOMS, ONLY MOLECULES
C                                    OR  ONLY TEST IONS ARE FOLLOWED
         DO 10 I=1,3
           IF (LROW(I)) THEN
             IF (LCOLM(1)) THEN
               FC(1)=(B(1)-P(I,2)-P(I,3))/P(I,1)
             ELSEIF (LCOLM(2)) THEN
               FC(2)=(B(2)-P(I,3))/P(I,2)
             ELSE
C  NOTHING TO BE DONE, ALL FC'S ARE 1.
             ENDIF
           ENDIF
10       CONTINUE
C
C
      ELSEIF (IROW.EQ.2) THEN
C  TWO ROWS ARE NON ZERO
        J1=0
        J2=0
C  DETERMINE THE INDICES FOR THE NON ZERO ROWS
        DO 20 J=1,3
          IF (LROW(J)) THEN
            IF (J1.EQ.0) THEN
              J1=J
            ELSE
              J2=J
            ENDIF
          ENDIF
20      CONTINUE
C
        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NON ZERO
          DO 30 I=1,3
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
30        CONTINUE
C
        ELSE
C  MORE THAN ONE COLUMN IS NON ZERO
          I1=0
          I2=0
          IF (ICOL.EQ.2) THEN
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
            DO 40 I=1,3
              IF (LCOLM(I)) THEN
                IF (I1.EQ.0) THEN
                  I1=I
                ELSE
                  I2=I
                ENDIF
              ENDIF
40          CONTINUE
C
          ELSE
            I1=1
            I2=2
            B(J1)=B(J1)-P(J1,3)
            B(J2)=B(J2)-P(J2,3)
          ENDIF
C
          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSEIF (IROW.EQ.3) THEN
C
        IF (ICOL.EQ.1) THEN
          DO 50 I=1,3
            IF (LCOLM(I)) FC(I)=B(1)/P(1,I)
50        CONTINUE
C
        ELSEIF (ICOL.EQ.2) THEN
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
          DO 60 I=1,3
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
              ENDIF
            ENDIF
60        CONTINUE
C
          FNEN=P(1,I1)*P(2,I2)-P(2,I1)*P(1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(1)*P(2,I2)-B(2)*P(1,I2))/FNEN
            IF (ABS(P(1,I2)).GT.EPS12) THEN
              FC(I2)=(B(1)-P(1,I1)*FC(I1))/P(1,I2)
            ELSEIF (ABS(P(2,I2)).GT.EPS12) THEN
              FC(I2)=(B(2)-P(2,I1)*FC(I1))/P(2,I2)
            ENDIF
          ENDIF
C
        ELSE
C  THE WHOLE MATRIX IS TO BE USED
          P11=P(1,1)
          P21=P(2,1)
          P31=P(3,1)
          P12=P(1,2)
          P22=P(2,2)
          P32=P(3,2)
          P13=P(1,3)
          P23=P(2,3)
          P33=P(3,3)
          B1=B(1)
          B2=B(2)
          B3=B(3)
          dta=deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          fc(1)=dtb1/(dta+1.d-30)
          fc(2)=dtb2/(dta+1.d-30)
          fc(3)=dtb3/(dta+1.d-30)
        ENDIF
      ENDIF
C
1000  CONTINUE

!  FOR THE TIME BEING
      FPH=1._DP

      CALL LEER(1)
      WRITE (iunout,*) 
     .  'EIRENE RECOMMENDED RESCALING OF VOLUME AVERAGED '
      WRITE (iunout,*) 'TALLIES DUE TO STATISTICAL ERRORS IN BALANCE '
      CALL MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .             FC(1),FC(2),FC(3),fph)
      CALL LEER(2)
C
      FA=FC(1)
      FM=FC(2)
      FI=FC(3)


C
      RETURN
      END
C ===== SOURCE: getscl4.f
C
C
      SUBROUTINE GETSCL4 (ISTRA,FA,FM,FI,FPH)
C
C  FIND SCALING FACTORS TO ENFORCE PARTICLE BALANCE
C  SIMPLE VERSION: NOT SPLIT BY SPECIES, ONLY BY TYPE
C  MODIFIED JAN/95: INCLUDE SURFACE TALLIES IN MATRIX, NOT IN
C  INHOMOGENITY
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COUTAU
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), INTENT(OUT) :: FA, FM, FI, FPH
      REAL(DP) :: FC(4), P(4,5), B(4), PP(4,4), FFC(3), BB(3)
      REAL(DP) :: DTB1, DTB2, DTA, DETER, DTB3, FNEN, DTB4, DETER4X4
      REAL(DP) :: P11, P12, P13, P21, P22, P23, P31, P32, P33,
     .            B1, B2, B3
      INTEGER :: ICOL, IROW, I, J, I1, I2, J1, J2, IC, JOUT, IOUT
      LOGICAL :: LCOLM(4), LROW(4)
C
C
      FC(1)=1.
      FC(2)=1.
      FC(3)=1.
      FC(4)=1.
C
C P(..,1)*FC(1)
C P(..,2)*FC(2)
C P(..,3)*FC(3)
C P(..,4)*FC(4)
C
      P(1,1)=PAATI(0,ISTRA)+POTATI(0,ISTRA)+PRFAAI(0,ISTRA)+
     .       PGENAI(0,ISTRA)
      P(1,2)=PMATI(0,ISTRA)+PRFMAI(0,ISTRA)
      P(1,3)=PIATI(0,ISTRA)+PRFIAI(0,ISTRA)
      P(1,4)=0._DP
      P(2,1)=PAMLI(0,ISTRA)+PRFAMI(0,ISTRA)
      P(2,2)=PMMLI(0,ISTRA)+POTMLI(0,ISTRA)+PRFMMI(0,ISTRA)+
     .       PGENMI(0,ISTRA)
      P(2,3)=PIMLI(0,ISTRA)+PRFIMI(0,ISTRA)
      P(2,4)=0._DP
      P(3,1)=PAIOI(0,ISTRA)+PRFAII(0,ISTRA)
      P(3,2)=PMIOI(0,ISTRA)+PRFMII(0,ISTRA)
      P(3,3)=PIIOI(0,ISTRA)+POTIOI(0,ISTRA)+PRFIII(0,ISTRA)+
     .       PGENII(0,ISTRA)
      P(3,4)=0._DP
      P(4,1)=0._DP
      P(4,2)=0._DP
      P(4,3)=0._DP
      P(4,4)=1._DP
C
      B(1)=-(PPATI(0,ISTRA)+WTOTA(0,ISTRA))
      B(2)=-(PPMLI(0,ISTRA)+WTOTM(0,ISTRA))
      B(3)=-(PPIOI(0,ISTRA)+WTOTI(0,ISTRA))
      B(4)=1._DP
C
      ICOL=0
      IROW=0
      DO 1 I=1,4
        LROW(I)=P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,4)**2.GT.EPS30
        LCOLM(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2+P(4,I)**2.GT.EPS30
        IF (LROW(I)) IROW=IROW+1
        IF (LCOLM(I)) ICOL=ICOL+1
1     CONTINUE
C
      IF (IROW.EQ.0) THEN
C  NO ROW IS NON ZERO, I.E. NO PARTICLES FOLLOWED
        GOTO 1000
C
      ELSEIF (IROW.EQ.1) THEN
C  ONLY ONE ROW (NO. I) IS NON ZERO, I.E., ONLY ATOMS, ONLY MOLECULES
C                                    OR  ONLY TEST IONS ARE FOLLOWED
         DO 10 I=1,4
           IF (LROW(I)) THEN
             IF (LCOLM(1)) THEN
               FC(1)=(B(I)-P(I,2)-P(I,3)-P(I,4))/P(I,1)
             ELSEIF (LCOLM(2)) THEN
               FC(2)=(B(I)-P(I,3)-P(I,4))/P(I,2)
             ELSEIF (LCOLM(3)) THEN
               FC(3)=(B(I)-P(I,4))/P(I,3)
             ELSE
               FC(4)=B(I)/P(I,4)
             ENDIF
           ENDIF
10       CONTINUE
C
C
      ELSEIF (IROW.EQ.2) THEN
C  TWO ROWS ARE NON ZERO
        J1=0
        J2=0
C  DETERMINE THE INDICES FOR THE NON ZERO ROWS
        DO 20 J=1,4
          IF (LROW(J)) THEN
            IF (J1.EQ.0) THEN
              J1=J
            ELSE
              J2=J
            ENDIF
          ENDIF
20      CONTINUE
C
        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NON ZERO
          DO 30 I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
30        CONTINUE
C
        ELSE
C  MORE THAN ONE COLUMN IS NON ZERO
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE FIRST TWO NON ZERO COLUMNS
          DO 40 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
40        CONTINUE
C
! nur spalten rechts von spalte i2 koennen noch werte enthalten
          DO I=I2+1,4
            B(J1)=B(J1)-P(J1,I)
            B(J2)=B(J2)-P(J2,I)
          END DO
C
          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSEIF (IROW.EQ.3) THEN
C  THREE ROWS ARE NON ZERO
C
C  DETERMINE THE ROW WHICH IS COMPLETELY ZERO
        JOUT = 0
        DO J=1,4
          IF (.NOT.LROW(J)) JOUT=J
        END DO

        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NON ZERO
          J1=1
          IF (JOUT.EQ.1) J1=2
          DO  I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
          END DO

        ELSEIF (ICOL.EQ.2) THEN
C  TWO COLUMNS ARE NON ZERO
          J1=1
          IF (JOUT.EQ.J1) J1=J1+1
          J2=J1+1
          IF (JOUT.EQ.J2) J2=J2+1

C  DETERMINE THE INDICES FOR THE FIRST TWO NON ZERO COLUMNS
          I1=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
          END DO

          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF

        ELSE

C  AT LEAST THREE COLUMNS ARE NON ZERO
          IF (ICOL.EQ.4) THEN
            DO J=1,4
              B(J)=B(J)-P(J,4)
            END DO
            IOUT=4
          ELSE
            IOUT=0
            DO I=1,4
              IF (.NOT.LCOLM(I)) IOUT=I
            END DO
          END IF

          J1=0
          DO J=1,4
            IF (J.EQ.JOUT) CYCLE
            J1=J1+1
            I1=0
            BB(J1)=B(J)
            DO I=1,4
              IF (I.EQ.IOUT) CYCLE
              I1=I1+1
              PP(J1,I1)=P(J,I)
            END DO
          END DO

          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=BB(1)
          B2=BB(2)
          B3=BB(3)
          dta=deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          Ffc(1)=dtb1/(dta+1.d-30)
          Ffc(2)=dtb2/(dta+1.d-30)
          Ffc(3)=dtb3/(dta+1.d-30)

          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
              IF (IC == 3) EXIT
            END IF
          END DO
        END IF
C
C
      ELSEIF (IROW.EQ.4) THEN
C  ALL ROWS ARE NON ZERO
C
        IF (ICOL.EQ.1) THEN
          DO 50 I=1,4
            IF (LCOLM(I)) FC(I)=B(I)/P(1,I)
50        CONTINUE
C
        ELSEIF (ICOL.EQ.2) THEN
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
          DO 60 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
              ENDIF
            ENDIF
60        CONTINUE
C
          FNEN=P(1,I1)*P(2,I2)-P(2,I1)*P(1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(1)*P(2,I2)-B(2)*P(1,I2))/FNEN
            IF (ABS(P(1,I2)).GT.EPS12) THEN
              FC(I2)=(B(1)-P(1,I1)*FC(I1))/P(1,I2)
            ELSEIF (ABS(P(2,I2)).GT.EPS12) THEN
              FC(I2)=(B(2)-P(2,I1)*FC(I1))/P(2,I2)
            ENDIF
          ENDIF
C
        ELSEIF (ICOL.EQ.3) THEN

          IC=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              PP(1:3,IC) = P(1:3,I)
            END IF
          END DO

          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=B(1)
          B2=B(2)
          B3=B(3)
          dta=deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          Ffc(1)=dtb1/(dta+1.d-30)
          Ffc(2)=dtb2/(dta+1.d-30)
          Ffc(3)=dtb3/(dta+1.d-30)

          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
            END IF
          END DO

        ELSE
C  THE WHOLE MATRIX IS TO BE USED
          pp(1:4,1:4) = p(1:4,1:4)
          dta=deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,1) = b(1:4)
          dtb1=deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,2) = b(1:4)
          dtb2=deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,3) = b(1:4)
          dtb3=deter4x4(pp)

          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,4) = b(1:4)
          dtb4=deter4x4(pp)

          fc(1)=dtb1/(dta+1.d-30)
          fc(2)=dtb2/(dta+1.d-30)
          fc(3)=dtb3/(dta+1.d-30)
          fc(4)=dtb4/(dta+1.d-30)
        ENDIF
      ENDIF
C
1000  CONTINUE

!  FOR THE TIME BEING

      CALL LEER(1)
      WRITE (iunout,*) 
     .  'EIRENE RECOMMENDED RESCALING OF VOLUME AVERAGED '
      WRITE (iunout,*) 'TALLIES DUE TO STATISTICAL ERRORS IN BALANCE '
      CALL MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .             FC(1),FC(2),FC(3),FC(4))
      CALL LEER(2)
C
      FA=FC(1)
      FM=FC(2)
      FI=FC(3)
      FPH=FC(4)


C
      RETURN
      END
C ===== SOURCE: statis.f
C  may 06:  bug fix: SDC initialized to zero
C
      SUBROUTINE STATIS

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CGRID
      USE CSDVI
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX
      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT
      REAL(DP), ALLOCATABLE, SAVE :: VECTOR(:), VECTRC(:,:),
     .          SD(:),   SDC(:,:)
      INTEGER, ALLOCATABLE, SAVE :: IADD(:),    IGFF(:),
     .                              IADDW(:),   IGFFW(:),
     .                              IADDC(:,:), IGFFC(:,:),
     .                              IND(:,:),   IIND(:),    INDSS(:,:)
      REAL(DP) :: D1, DS1, D2S, DS2, DSA, DD22, DA1, DD11, D2, DD12,
     .          ZFLUXQ, DS, ZNM, SD2S, SD2, SG2, SG, DA, D, DD, DA2,
     .          D2S11, D2S22, D2S12, SG12, SG1, DSA1, DSA2,
     .          SAV, SD1S, SD1, XNM
      INTEGER :: ISCO2, NR1, NP2, NT3, INP, IGF, IC,
     .           I, IRU, IIN, J, IR, IGS,
     .           ICO,  IS, ITL2, ISCO1, ITL1, IGE,
     .           NSYM, IGI, ITL, ISCO, NSYH, J1, J2, IP, IG, IT
      INTEGER, SAVE :: NSB, NRW
C
C
      ENTRY STATS0

      IMETCL = 0
      NCLMT = 0
      NCLMTS = 0
      LMETSP = .FALSE.

      IF (NSIGI.EQ.0) RETURN
C
      IF (.NOT.ALLOCATED(IADD)) THEN
        AllOCATE (IADD(NSD))
        AllOCATE (IGFF(NSD))
        AllOCATE (IADDW(NSDW))
        AllOCATE (IGFFW(NSDW))
        AllOCATE (IADDC(2,NCV))
        AllOCATE (IGFFC(2,NCV))
        AllOCATE (IND(NRTAL,8))
        AllOCATE (IIND(NRTAL))
        AllOCATE (INDSS(NRTAL,8))
        AllOCATE (VECTOR(MAX(NRTAL,NLMPGS)))
        AllOCATE (VECTRC(2,NRTAL))
        AllOCATE (SD(0:(MAX(NRTAL,NLMPGS))))
        AllOCATE (SDC(2,0:NRTAL))
        SD=0._DP
        SD2=0._DP
        SDC=0._DP
      END IF

C  FILL IIND, INDSS ARRAYS FOR THOSE "AVERAGE" CELLS, TO WHICH "REAL"
C  CELL IR ALSO CONTRIBUTES
C  IIND: HOW MANY CELLS
C  INDSS: WHICH CELLS
      CALL INDTAL(IND,NRTAL,NR1TAL,NP2TAL,NT3TAL,NBMLT)
      DO IR=1,NSBOX_TAL
        IIND(IR)=0
        IIN=0
        DO J=1,8
          IF (IND(IR,J).NE.0) THEN
            IIND(IR)=IIND(IR)+1
            IIN=IIN+1
            INDSS(IR,IIN)=J
          ENDIF
        ENDDO
      ENDDO
C
      DO 101 J=1,NSIGVI
        IGFF(J)=NFIRST(IIH(J))
        IADD(J)=NADDV(IIH(J))
101   CONTINUE
C
      DO 102 J=1,NSIGCI
        IGFFC(1,J)=NFIRST(IIHC(1,J))
        IGFFC(2,J)=NFIRST(IIHC(2,J))
        IADDC(1,J)=NADDV(IIHC(1,J))
        IADDC(2,J)=NADDV(IIHC(2,J))
102   CONTINUE
C
      DO 108 J=1,NSIGSI
        IGFFW(J)=NFRSTW(IIHW(J))
        IADDW(J)=NADDW(IIHW(J))
108   CONTINUE
C
      RETURN
C
      ENTRY STATS1(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
C
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
      NRW=NSIN

      IF ((NSIGVI > 0) .OR. (NSIGCI > 0)) THEN
C  HISTORY HAS TOUCHED NCLMT CELLS.
C  IT CONTRIBUTES TO NCLMTS CELLS (AVERAGES)
C  THESE CELL NUMBERS ARE STORED HERE ON ICLMT ARRAY
        NCLMTS = NCLMT
        DO I=1,NCLMT
          IR = ICLMT(I)
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            IF (IMETCL(IRU) == 0) THEN
              NCLMTS = NCLMTS+1
              IMETCL(IRU) = NCLMTS
              ICLMT(NCLMTS) = IRU
            END IF
          END DO
        END DO
      END IF
C
C
      IF (NSIGVI.EQ.0) GOTO 1020
C
C  ARE THERE CONTRIBUTIONS TO THE REQUESTED TALLY FROM THIS HISTORY?
      DO 1012 IC=1,NSIGVI
        INP=IADD(IC)
        IGF=IGFF(IC)
        IGS=IGH(IC)
        ITL=IIH(IC)
        IF (NSPAN(ITL) == 0) THEN
          ISCO = 1
        ELSE
          ISCO = 0
          IF (IGS == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL):NSPEN(ITL))) ) ISCO = 1
          ELSE
            IF (LMETSP(NSPAN(ITL)+IGS-1)) ISCO = 1
          END IF
        END IF
        IF (ISCO == 0) GOTO 1012
C
        IF (.NOT.LP.AND..NOT.LT) GOTO 1005
C  USE SYMMETRY IN POLOIDAL/Y AND/OR TOROIDAL/Z COORDINATE
        IF (IGS.LE.0) THEN
          IGI=1
          IGE=IGF
        ELSE
          IGI=IGS
          IGE=IGS
        ENDIF
        IF (LP) THEN
          NSYM=NP2
          NSYH=(NSYM-1)/2
          DO 1003 IG=IGI,IGE
          DO 1003 IR=1,NR1
          DO 1003 IT=1,NT3
          DO 1003 IP=1,NSYH
                J1=IR+((IT-1)*NP2+IP-1)*NR1
                J2=IR+((IT-1)*NP2+NSYM-IP-1)*NR1
                SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                ESTIMV(INP+IG,J1)=SAV
                ESTIMV(INP+IG,J2)=SAV
1003      CONTINUE
        ENDIF
        IF (LT) THEN
          NSYM=NT3
          NSYH=(NSYM-1)/2
          DO 1004 IG=IGI,IGE
          DO 1004 IR=1,NR1
          DO 1004 IP=1,NP2
          DO 1004 IT=1,NSYH
                J1=IR+((IT-1)*NP2+IP-1)*NR1
                J2=IR+((NSYM-IT-1)*NP2+IP-1)*NR1
                SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                ESTIMV(INP+IG,J1)=SAV
                ESTIMV(INP+IG,J2)=SAV
1004      CONTINUE
        ENDIF
1005    CONTINUE
C
C  FILL ARRAY VECTOR, EITHER FOR SUM OVER SPECIES OR FOR INDIVIDUAL SPECIES
C  VECTOR IS FILLED ONLY FOR THOSE CELLS,
C  WHICH HAVE BEEN TOUCHED BY THIS HISTROY
C  VECTOR CONTAINS THE SUM FROM ALL HISTORIES IN THESE CELLS UP TO THE
C  PRESENT HISTORY
        IF (IGS.NE.0) THEN
           DO ICO = 1,NCLMT
             IR = ICLMT(ICO)
             VECTOR(ICO)=ESTIMV(INP+IGS,IR)
           END DO
        ELSE
          DO 1014 ICO=1,NCLMT
            VECTOR(ICO)=0.
1014      CONTINUE
          DO 1015 IS=1,IGF
          DO 1015 ICO=1,NCLMT
            IR = ICLMT(ICO)
            VECTOR(ICO)=VECTOR(ICO)+ESTIMV(INP+IS,IR)
1015      CONTINUE
        ENDIF

        SD1S = 0.D0
C  FILL ARRAY SD WITH THE INDIVIDUAL CONTRIBUTION FROM THIS HISTROY,
C  IN EACH CELL THAT HAS BEEN TOUCHED BY THIS HISTORY
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTOR(ICO)-SDVIA(IC,IR)
          SD1S=SD1S+SD1
          SDVIA(IC,IR)=VECTOR(ICO)
          SD(IR) = SD1
C  FILL SD ALSO FOR OTHER "CELL", TO WHICH CELL IR CONTRIBUTES
C  I.E., AVERAGES OVER COORDINATES OR OVER THE ENTIRE COMPUTATIONAL DOMAIN
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SD(IR)
          SIGMA(IC,IR)=SIGMA(IC,IR)+SD1*SD1
          SD(IR)=0._DP
        END DO
        SGMS(IC)=SGMS(IC)+SD1S*SD1S
1012  CONTINUE
C
C
1020  CONTINUE
      IF (NSIGSI.EQ.0) GOTO 1030
      DO 1022 IC=1,NSIGSI
        INP=IADDW(IC)
        IGF=IGFFW(IC)
        IGS=IGHW(IC)
        IF (IGS.NE.0) THEN
          DO 1023 IR=1,NRW
            VECTOR(IR)=ESTIMS(INP+IGS,IR)
1023      CONTINUE
        ELSE
          DO 1024 IR=1,NRW
            VECTOR(IR)=0.
1024      CONTINUE
          DO 1025 IS=1,IGF
          DO 1025 IR=1,NRW
            VECTOR(IR)=VECTOR(IR)+ESTIMS(INP+IS,IR)
1025      CONTINUE
        ENDIF
C
        SD1S=0.
        DO 1021 IR=1,NRW
          SD1=VECTOR(IR)-SDVIAW(IC,IR)
          SD1S=SD1S+SD1
          SIGMAW(IC,IR)=SIGMAW(IC,IR)+SD1*SD1
          SDVIAW(IC,IR)=VECTOR(IR)
1021    CONTINUE
        SGMWS(IC)=SGMWS(IC)+SD1S*SD1S
1022  CONTINUE
C
1030  CONTINUE
C
      IF (NSIGCI.EQ.0) GOTO 1050
C
      DO 1032 IC=1,NSIGCI
        ITL1=IIHC(1,IC)
        ITL2=IIHC(2,IC)
        IF (NSPAN(ITL1) == 0) THEN
          ISCO1 = 1
        ELSE
          ISCO1 = 0
          IF (IGHC(1,IC) == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL1):NSPEN(ITL1))) ) ISCO1 = 1
          ELSE
            IF (LMETSP(NSPAN(ITL1)+IGHC(1,IC)-1)) ISCO1 = 1
          END IF
        END IF
        IF (NSPAN(ITL2) == 0) THEN
          ISCO2 = 1
        ELSE
          ISCO2 = 0
          IF (IGHC(2,IC) == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL2):NSPEN(ITL2))) ) ISCO2 = 1
          ELSE
            IF (LMETSP(NSPAN(ITL2)+IGHC(2,IC)-1)) ISCO2 = 1
          END IF
        END IF
        IF (ISCO1+ISCO2 == 0) GOTO 1032
C
        DO 1037 I=1,2
          INP=IADDC(I,IC)
          IGF=IGFFC(I,IC)
          IGS=IGHC(I,IC)
C
          IF (.NOT.LP.AND..NOT.LT) GOTO 1035
          IF (IGS.LE.0) THEN
            IGI=1
            IGE=IGF
          ELSE
            IGI=IGS
            IGE=IGS
          ENDIF
          IF (LP) THEN
            NSYM=NP2
            NSYH=(NSYM-1)/2
            DO 1033 IG=IGI,IGE
            DO 1033 IR=1,NR1
            DO 1033 IT=1,NT3
            DO 1033 IP=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((IT-1)*NP2+NSYM-IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
1033        CONTINUE
          ENDIF
          IF (LT) THEN
            NSYM=NT3
            NSYH=(NSYM-1)/2
            DO 1034 IG=IGI,IGE
            DO 1034 IR=1,NR1
            DO 1034 IP=1,NP2
            DO 1034 IT=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((NSYM-IT-1)*NP2+IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
1034        CONTINUE
          ENDIF
1035      CONTINUE
C
          IF (IGS.NE.0) THEN
            DO ICO = 1,NCLMT
              IR = ICLMT(ICO)
              VECTRC(I,ICO)=ESTIMV(INP+IGS,IR)
            END DO
          ELSE
            DO 1044 ICO=1,NCLMT
              VECTRC(I,ICO)=0.
1044        CONTINUE
            DO 1045 IS=1,IGF
            DO 1045 ICO=1,NCLMT
              IR = ICLMT(ICO)
              VECTRC(I,ICO)=VECTRC(I,ICO)+ESTIMV(INP+IS,IR)
1045        CONTINUE
          ENDIF
1037    CONTINUE
C
C
        SD1S = 0.D0
        SD2S = 0.D0
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTRC(1,ICO)-SDVIAC(1,IC,IR)
          SD2 = VECTRC(2,ICO)-SDVIAC(2,IC,IR)
          SD1S=SD1S+SD1
          SD2S=SD2S+SD2
          SDVIAC(1,IC,IR)=VECTRC(1,ICO)
          SDVIAC(2,IC,IR)=VECTRC(2,ICO)
          SDC(1,IR) = SD1
          SDC(2,IR) = SD2
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SDC(1,IRU)=SDC(1,IRU)+SD1
            SDC(2,IRU)=SDC(2,IRU)+SD2
          END DO
        END DO
C
        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SDC(1,IR)
          SD2=SDC(2,IR)
          SIGMAC(0,IC,IR)=SIGMAC(0,IC,IR)+SD1*SD2
          SIGMAC(1,IC,IR)=SIGMAC(1,IC,IR)+SD1*SD1
          SIGMAC(2,IC,IR)=SIGMAC(2,IC,IR)+SD2*SD2
          SDC(1:2,IR)=0._DP
        END DO
        SGMCS(0,IC)=SGMCS(0,IC)+SD1S*SD2S
        SGMCS(1,IC)=SGMCS(1,IC)+SD1S*SD1S
        SGMCS(2,IC)=SGMCS(2,IC)+SD2S*SD2S
1032  CONTINUE
C
C
1050  CONTINUE
      RETURN
C
      ENTRY STATS2(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.D0) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NSIGVI.EQ.0) GOTO 2200
C
      DO 2112 IC=1,NSIGVI
        INP=IADD(IC)
        IGF=IGFF(IC)
        IGS=IGH(IC)
        IF (IGS.NE.0) THEN
          DO 2113 IR=1,NSB
            VECTOR(IR)=ESTIMV(INP+IGS,IR)
2113      CONTINUE
        ELSE
          DO 2114 IR=1,NSB
            VECTOR(IR)=0.
2114      CONTINUE
          DO 2115 IS=1,IGF
          DO 2115 IR=1,NSB
            VECTOR(IR)=VECTOR(IR)+ESTIMV(INP+IS,IR)
2115      CONTINUE
        ENDIF
C
        DS=0.
        DO 2011 IR=1,NSB
          SD1=VECTOR(IR)
          DS=DS+SD1
          DO 2016 IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
2016      CONTINUE
2011    CONTINUE
C
        DO 2111 IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA(IC,IR)-DD/XN)
C RELATIV STANDARD DEVIATION FOR CURRENT STRATUM
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA(IC,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV(IC,IR)=STV(IC,IR)+SG2*ZFLUXQ/XNM/XN
          EE(IC,IR)=EE(IC,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
2111    CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS(IC)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS(IC)=SG*FSIG
C
        STVS(IC)=STVS(IC)+SG2*ZFLUXQ/XNM/XN
        EES(IC)=EES(IC)+DS*ZFLUX/XN
2112  CONTINUE
C
2200  CONTINUE
      IF (NSIGSI.EQ.0) GOTO 2300
      DO 2212 IC=1,NSIGSI
        INP=IADDW(IC)
        IGF=IGFFW(IC)
        IGS=IGHW(IC)
        DS=0.
        IF (IGS.NE.0) THEN
          DO 2213 IR=1,NRW
            VECTOR(IR)=ESTIMS(INP+IGS,IR)
2213      CONTINUE
        ELSE
          DO 2214 IR=1,NRW
            VECTOR(IR)=0.
2214      CONTINUE
          DO 2215 IS=1,IGF
          DO 2215 IR=1,NRW
            VECTOR(IR)=VECTOR(IR)+ESTIMS(INP+IS,IR)
2215      CONTINUE
        ENDIF
        DO 2211 IR=1,NRW
          D=VECTOR(IR)
          DS=DS+D
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMAW(IC,IR)-DD/XN)
C RELATIV STANDARD DEVIATION FOR CURRENT STRATUM
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMAW(IC,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STVW(IC,IR)=STVW(IC,IR)+SG2*ZFLUXQ/XNM/XN
          FF(IC,IR)=FF(IC,IR)+D*ZFLUX/XN
2211    CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMWS(IC)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMWS(IC)=SG*FSIG
C
        STVWS(IC)=STVWS(IC)+SG2*ZFLUXQ/XNM/XN
        FFS(IC)=FFS(IC)+DS*ZFLUX/XN
2212  CONTINUE
C
2300  CONTINUE
C
      IF (NSIGCI.EQ.0) GOTO 2400
C
      DO 2312 IC=1,NSIGCI
        DO 2317 I=1,2
          INP=IADDC(I,IC)
          IGF=IGFFC(I,IC)
          IGS=IGHC(I,IC)
          IF (IGS.NE.0) THEN
            DO 2313 IR=1,NSB
              VECTRC(I,IR)=ESTIMV(INP+IGS,IR)
2313        CONTINUE
          ELSE
            DO 2314 IR=1,NSB
              VECTRC(I,IR)=0.
2314        CONTINUE
            DO 2315 IS=1,IGF
            DO 2315 IR=1,NSB
              VECTRC(I,IR)=VECTRC(I,IR)+ESTIMV(INP+IS,IR)
2315        CONTINUE
          ENDIF
2317    CONTINUE
C
        DS1=0.
        DS2=0.
        DO 2311 IR=1,NSB
          SD1=VECTRC(1,IR)
          SD2=VECTRC(2,IR)
          DS1=DS1+SD1
          DS2=DS2+SD2
          DO 2316 IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SDC(1,IRU)=SDC(1,IRU)+SD1
            SDC(2,IRU)=SDC(2,IRU)+SD2
2316      CONTINUE
2311    CONTINUE
        DO 2411 IR=1,NSB
          D1=SDC(1,IR)
          D2=SDC(2,IR)
          DD12=D1*D2
          DD11=D1*D1
          DD22=D2*D2
          DA1=ABS(D1)
          DA2=ABS(D2)
          SG12=         SIGMAC(0,IC,IR)-DD12/XN
          SG1 =MAX(0._DP,SIGMAC(1,IC,IR)-DD11/XN)
          SG2 =MAX(0._DP,SIGMAC(2,IC,IR)-DD22/XN)
C ABSOLUTE STANDARD DEVIATION AND COVARIANCES
          SIGMAC(0,IC,IR)=SG12/XNM/XN
          SIGMAC(1,IC,IR)=SQRT(SG1/XNM/XN)
          SIGMAC(2,IC,IR)=SQRT(SG2/XNM/XN)
          SDC(1:2,IR)=0._DP
2411    CONTINUE
        D2S12=DS1*DS2
        D2S11=DS1*DS1
        D2S22=DS2*DS2
        DSA1=ABS(DS1)
        DSA2=ABS(DS2)
        SG12=         SGMCS(0,IC)-D2S12/XN
        SG1 =MAX(0._DP,SGMCS(1,IC)-D2S11/XN)
        SG2 =MAX(0._DP,SGMCS(2,IC)-D2S22/XN)
        SGMCS(0,IC)=SG12/XNM/XN
        SGMCS(1,IC)=SQRT(SG1/XNM/XN)
        SGMCS(2,IC)=SQRT(SG2/XNM/XN)
2312  CONTINUE
C
2400  RETURN

      ENTRY STATS3

      IF (ALLOCATED(IADD)) THEN
         DEAllOCATE (IADD)
         DEAllOCATE (IGFF)
         DEAllOCATE (IADDW)
         DEAllOCATE (IGFFW)
         DEAllOCATE (IADDC)
         DEAllOCATE (IGFFC)
         DEAllOCATE (IND)
         DEAllOCATE (IIND)
         DEAllOCATE (INDSS)
         DEAllOCATE (VECTOR)
         DEAllOCATE (VECTRC)
         DEAllOCATE (SD)
         DEAllOCATE (SDC)
      END IF

      RETURN
      END
C ===== SOURCE: statis_spc.f
C
C
      SUBROUTINE STATIS_SPC

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CSDVI
      USE CSDVI_COP
      USE COUTAU
      USE COMSOU

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX
      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT

      REAL(DP), ALLOCATABLE :: SD(:)
      REAL(DP) :: XNM, DS, ZFLUXQ, D2S, SG,
     .            DSA, DD, D, SG2, DA, SD1, SD1S
      INTEGER :: NSB, NP2, NR1, NT3, I, ISPC
C
      SAVE
C
      ENTRY STATS0_SPC

      IF (NADSPC > 0) THEN
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%PSPC%IMETSP = 0
        END DO
      END IF
C
      RETURN

C
      ENTRY STATS1_SPC(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
C
C
      IF (NADSPC.EQ.0) RETURN
C
C
C  STATISTICS FOR SPECTRA
      DO ISPC=1,NADSPC
C
        IF (ESTIML(ISPC)%PSPC%IMETSP > 0) THEN
          SD1S=0.
          ALLOCATE (SD(0:ESTIML(ISPC)%PSPC%NSPC+1))
          SD = 0.D0
          DO I = 0,ESTIML(ISPC)%PSPC%NSPC+1
            SD1=ESTIML(ISPC)%PSPC%SPC(I)-ESTIML(ISPC)%PSPC%SDV(I)
            SD1S=SD1S+SD1
            ESTIML(ISPC)%PSPC%SDV(I)=ESTIML(ISPC)%PSPC%SPC(I)
            SD(I) = SD1
          END DO

          DO I = 0,ESTIML(ISPC)%PSPC%NSPC+1
            SD1=SD(I)
            ESTIML(ISPC)%PSPC%SGM(I)=ESTIML(ISPC)%PSPC%SGM(I)+SD1*SD1
          END DO
          ESTIML(ISPC)%PSPC%SGMS=ESTIML(ISPC)%PSPC%SGMS+SD1S*SD1S
          DEALLOCATE (SD)
        END IF
      END DO
C
C
      RETURN
C
      ENTRY STATS2_SPC(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NADSPC.EQ.0) RETURN
C
C  STATISTICS FOR MOMENTUM SOURCES
      DO ISPC=1,NADSPC
C
        ALLOCATE (SD(0:ESTIML(ISPC)%PSPC%NSPC+1))

        SD=ESTIML(ISPC)%PSPC%SPC
        DS=SUM(SD)

        DO I=0,ESTIML(ISPC)%PSPC%NSPC+1
          D=SD(I)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,ESTIML(ISPC)%PSPC%SGM(I)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          ESTIML(ISPC)%PSPC%SGM(I)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
! STV
          IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
            SMESTL(ISPC)%PSPC%SGM(I)=SMESTL(ISPC)%PSPC%SGM(I)+
     .                               SG2*ZFLUXQ/XNM/XN
! EE
            SMESTL(ISPC)%PSPC%SDV(I)=SMESTL(ISPC)%PSPC%SDV(I)+D*ZFLUX/XN
          END IF
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,ESTIML(ISPC)%PSPC%SGMS-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        ESTIML(ISPC)%PSPC%SGMS=SG*FSIG
C
        IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
          SMESTL(ISPC)%PSPC%STVS=SMESTL(ISPC)%PSPC%STVS+
     .                           SG2*ZFLUXQ/XNM/XN
          SMESTL(ISPC)%PSPC%EES =SMESTL(ISPC)%PSPC%EES+DS*ZFLUX/XN
        END IF

        DEALLOCATE (SD)
      END DO
C
2200  CONTINUE
      RETURN
      END
C ===== SOURCE: store.f
C
      SUBROUTINE STORE(IFLAG)

C  PURPOSE: STORE TRAJECTORIES
C  TO BE WRITTEN

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      INTEGER :: IFIRST, IFLAG
      DATA IFIRST/0/

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        OPEN (UNIT=16,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      ENDIF
      GOTO (10,20,30,40,50,60,70,80,90),IFLAG
      WRITE (iunout,*) 'IFLAG OUT OF RANGE IN SUBR. STORE '
      WRITE (iunout,*) 'EXIT CALLED '
      CALL EXIT_OWN(1)
C  LOCATE
10    CONTINUE
20    CONTINUE
C  IONIZATION
30    CONTINUE
40    CONTINUE
50    CONTINUE
C  SURFACE
60    CONTINUE
70    CONTINUE
80    CONTINUE
C  NEW CELL
90    CONTINUE
      RETURN
      END
C ===== SOURCE: update.f
C  27.6.05 updphot: iadd removed
C  21.01.06 photon background for test atoms: removed
C  18.04.06 test ions and atoms: syncronized
C           bug fix: V0_para  --> parmom_0 for elastic momentum source
C                                 contribution from atoms.
C
      SUBROUTINE UPDATE
C
C ESTIMATORS ARE UPDATED FOR EACH TRACK TAKING T/VEL SEC.
C T (CM) IS STORED ON CLPD ARRAY FOR ONE OR MORE CELLS, THAT HAVE
C BEEN CROSSED WITHOUT COLLISION.
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CUPD
      USE CGRID
      USE CSPEZ
      USE CGEOM
      USE COMPRT
      USE CSDVI
      USE COMXS
      USE CCONA
      USE PHOTON

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN OUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                            XSTORV2(NSTORV,N2ND+N3RD)
      INTEGER, INTENT(IN OUT) :: IFLAG
      REAL(DP) :: WTRSIG, DIST, WTR, WTRE0, WV, VELQ, CNDYNPH, WTRV,
     .            V0_PARB, PARMOM_0, P
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:), CNDYNM(:), CNDYNI(:)
      INTEGER :: IRD,  I, IRDO, INUM,
     .           IPL, IAT, IA,
     .           IM,  IIO, IP, IML, II, NPBGK,
     .           IBGK, I1, I2, IPLV, IPLSV
C SECONDARY SPECIES IDENTIFIERS
      INTEGER ::  IAT1,IAT2,IML1,IML2,IIO1,IIO2,IPH1,IPH2,IPL1,IPL2
C EL PROCESSES
      INTEGER ::      IAEL,IREL
      INTEGER ::      IMEL
      INTEGER ::      IIEL
C CX PROCESSES
      INTEGER ::      IACX,IRCX
      INTEGER ::      IMCX
      INTEGER ::      IICX
C OT PROCESSES
      INTEGER ::      IAOT,IROT,UPDF
C PI PROCESSES (to be done)
      INTEGER ::      IAPI,IRPI
C EI PROCESSES
      INTEGER ::      IAEI,IREI
      INTEGER ::      IMEI
      INTEGER ::      IIEI
C
C  ESTIMATORS FOR ATOMS
C
      ENTRY UPDATM (XSTOR2,XSTORV2,IFLAG)
C
      WV=WEIGHT/VEL
      NPBGK=NPBGKA(IATM)
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
      IF ((NPBGK.GT.0).AND.LBGKV)
     .   CALL UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)

      IF (IUPDTE == 2) RETURN

      IF (.NOT.ALLOCATED(CNDYNA)) THEN
        ALLOCATE (CNDYNA(NATM))
        DO IAT=1,NATMI
          CNDYNA(IAT)=AMUA*RMASSA(IAT)
        END DO
      END IF
          
      VELQ=VEL*VEL
C
      DO 51 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTR*VEL*CNDYNA(IATM)
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENA) EDENA(IATM,IRD)=EDENA(IATM,IRD)+WTRE0
        IF (LPDENA) PDENA(IATM,IRD)=PDENA(IATM,IRD)+WTR
        IF (LEDENA.OR.LPDENA) LMETSP(NSPH+IATM)=.TRUE.

        IF (LVXDENA) VXDENA(IATM,IRD)=VXDENA(IATM,IRD)+WTRV*VELX
        IF (LVYDENA) VYDENA(IATM,IRD)=VYDENA(IATM,IRD)+WTRV*VELY
        IF (LVZDENA) VZDENA(IATM,IRD)=VZDENA(IATM,IRD)+WTRV*VELZ
        IF (LVXDENA.OR.LVYDENA.OR.LVZDENA) LMETSP(NSPH+IATM)=.TRUE. 
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 51
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)-WTRSIG
        IF (LEAAT) EAAT(IRD)     =EAAT(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGACX(IATM,0,0).EQ.0) GOTO 43
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 44  IACX=1,NACXI(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IPLS=LGACX(IATM,IACX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LPAAT.AND.IESTCX(IRCX,1).NE.0) THEN
            PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            IF (LPAPL) THEN
              PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT1,IRD)= PAAT(IAT1,IRD)+WTRSIG
                LMETSP(NSPH+IAT1)=.TRUE.
              END IF
C           ELSEIF (N1STX(IRCX,1).EQ.2) THEN
C             IML1=N1STX(IRCX,2)
C             LOGMOL(IML1,ISTRA)=.TRUE.
C             IF (LPAML) THEN
C               PAML(IML1,IRD)= PAML(IML1,IRD)+WTRSIG
C               LMETSP(NSPA+IML1)=.TRUE.
C             END IF
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO1,IRD)= PAIO(IIO1,IRD)+WTRSIG
                LMETSP(NSPAM+IIO1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL1,IRD)= PAPL(IPL1,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL1)=.TRUE.
              END IF
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT2,IRD)= PAAT(IAT2,IRD)+WTRSIG
                LMETSP(NSPH+IAT2)=.TRUE.
              END IF
C           ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
C             IML2=N2NDX(IRCX,2)
C             LOGMOL(IML2,ISTRA)=.TRUE.
C             IF (LPAML) THEN
C               PAML(IML2,IRD)= PAML(IML2,IRD)+WTRSIG
C               LMETSP(NSPA+IML2)=.TRUE.
C             END IF
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO2,IRD)= PAIO(IIO2,IRD)+WTRSIG
                LMETSP(NSPAM+IIO2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL2,IRD)= PAPL(IPL2,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL2)=.TRUE.
              END IF
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LEA) THEN
            IF (LEAAT.AND.IESTCX(IRCX,3).NE.0) THEN
              EAAT(IRD)   = EAAT(IRD) + WTRSIG*E0
            ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEAPL) EAPL(IRD)   = EAPL(IRD) - WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
              IF (N1STX(IRCX,1).EQ.1) THEN
                IAT1=N1STX(IRCX,2)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LEAAT) EAAT(IRD) = EAAT(IRD) + WTRSIG*ESIGCX(IRCX,1)
C             ELSEIF (N1STX(IRCX,1).EQ.2) THEN
C               IML1=N1STX(IRCX,2)
C               LOGMOL(IML1,ISTRA)=.TRUE.
C               IF (LEAML) EAML(IRD) = EAML(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.3) THEN
                IIO1=N1STX(IRCX,2)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LEAIO) EAIO(IRD) = EAIO(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LEAPL) EAPL(IRD) = EAPL(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
              IF (N2NDX(IRCX,1).EQ.1) THEN
                IAT2=N2NDX(IRCX,2)
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LEAAT) EAAT(IRD) = EAAT(IRD) + WTRSIG*E0
C             ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
C               IML2=N2NDX(IRCX,2)
C               LOGMOL(IML2,ISTRA)=.TRUE.
C               IF (LEAML) EAML(IRD) = EAML(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
                IIO2=N2NDX(IRCX,2)
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LEAIO) EAIO(IRD) = EAIO(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LEAPL) EAPL(IRD) = EAPL(IRD) + WTRSIG*E0
              ENDIF
            ENDIF
          ENDIF
C
44      CONTINUE
43      CONTINUE
C
C  ELASTIC NEUTRAL BULK-ION COLLISION CONTRIBUTION
C
        IF (LGAEL(IATM,0,0).EQ.0) GOTO 60
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 61  IAEL=1,NAELI(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IPLS=LGAEL(IATM,IAEL,1)
C  DO NOT UPDATE BGK TALLIES HERE
          IBGK=NPBGKP(IPLS,1)
          IF (IBGK.NE.0) GOTO 61
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVEL(IREL)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTEL(IREL,1).NE.0) THEN
            IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C  UPDATE TRACKLENGTH ESTIMATOR
C           IF (LPAPL) THEN
C             PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
C             PAPL(IPLS,IRD)=PAPL(IPLS,IRD)+WTRSIG
C             LMETSP(NSPAMI+IPLS)=.TRUE.
C           END IF
            IF (LPAAT) THEN
              PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
              LMETSP(NSPH+IATM)=.TRUE.
            END IF
          ENDIF
C
          IF (LEA) THEN
            IF (IESTEL(IREL,3).NE.0) THEN
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
            ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEAPL) EAPL(IRD)=EAPL(IRD)-WTRSIG*ESIGEL(IREL,1)
C
C  FIRST SECONDARY: = INCIDENT ION. REMAINS SAME PARTICLE BY DEFAULT
              IF (LEAPL) EAPL(IRD)=EAPL(IRD)+WTRSIG*E0
C  SECOND SECONDARY: = INCIDENT ATOM. REMAINS SAME PARTICLE BY DEFAULT
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEL(IREL,1)
            ENDIF
          ENDIF
C
61      CONTINUE
60      CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGAEI(IATM,0).EQ.0) GOTO 57
C
        DO 55 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          IF (SIGVEI(IREI).LE.0.D0) GOTO 55
C
          WTRSIG=WTR*SIGVEI(IREI)
C
C  COLLISION ESTIMATOR:
C  NOT AVAILABLE, HENCE: NO NEED FOR COMPENSATION HERE
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          IF (LPAEL) PAEL(IRD)=PAEL(IRD)+WTRSIG*PELDS(IREI)
C
C  POST COLLISION CONTRIBUTIONS
          DO IA=1,IPATDS(IREI,0)
            IAT=IPATDS(IREI,IA)
            LOGATM(IAT,ISTRA)=.TRUE.
            IF (LPAAT) THEN
              PAAT(IAT,IRD)=PAAT(IAT,IRD)+PATDS(IREI,IAT)*WTRSIG
              LMETSP(NSPH+IAT)=.TRUE.
            END IF
          END DO
          DO IM=1,IPMLDS(IREI,0)
            IML=IPMLDS(IREI,IM)
            LOGMOL(IML,ISTRA)=.TRUE.
            IF (LPAML) THEN
              PAML(IML,IRD)=PAML(IML,IRD)+PMLDS(IREI,IML)*WTRSIG
              LMETSP(NSPA+IML)=.TRUE.
            END IF
          END DO
          DO II=1,IPIODS(IREI,0)
            IIO=IPIODS(IREI,II)
            LOGION(IIO,ISTRA)=.TRUE.
            IF (LPAIO) THEN
              PAIO(IIO,IRD)=PAIO(IIO,IRD)+PIODS(IREI,IIO)*WTRSIG
              LMETSP(NSPAM+IIO)=.TRUE.
            END IF
          END DO
          DO IP=1,IPPLDS(IREI,0)
            IPL=IPPLDS(IREI,IP)
            LOGPLS(IPL,ISTRA)=.TRUE.
            IF (LPAPL) THEN
              PAPL(IPL,IRD)=PAPL(IPL,IRD)+PPLDS(IREI,IPL)*WTRSIG
              LMETSP(NSPAMI+IPL)=.TRUE.
            END IF
          END DO
C
          IF ((IESTEI(IREI,3).EQ.0) .AND. LEAEL) 
     .      EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGEI(IREI,5)

          IF (LEA) THEN
            IF (IESTEI(IREI,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
C
            ELSE
C 
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEI(IREI,1)
              IF (LEAML) EAML(IRD)=EAML(IRD)+WTRSIG*ESIGEI(IREI,2)
              IF (LEAIO) EAIO(IRD)=EAIO(IRD)+WTRSIG*ESIGEI(IREI,3)
              IF (LEAPL) EAPL(IRD)=EAPL(IRD)+WTRSIG*ESIGEI(IREI,4)
C 
            ENDIF
          ENDIF
55      CONTINUE
57      CONTINUE
C
C  PLASMA ION IMPACT CONTRIBUTION
C
        IF (LGAPI(IATM,0,0).EQ.0) GOTO 59
        DO 58  IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          IPLS=LGAPI(IATM,IAPI,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.

          WTRSIG=WTR*SIGVPI(IRPI)
C
C  PRE COLLISION BULK ION CONTRIBUTION, ASSUME: INCIDENT ION IS LOST
C
          IF (LPAPL) THEN
            PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
          IF (LEAPL) EAPL(IRD)     =EAPL(IRD)     +WTRSIG*ESIGPI(IRPI,1)
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION
C
          IF (LEAEL) EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGPI(IRPI,2)
          IF (LPAEL) PAEL(IRD)=PAEL(IRD)+WTRSIG*PELPI(IRPI)
C
C SO NICHT  EAIO(IRD)     = EAIO(IRD)     +WTRSIG*E0
          DO II=1,IPIOPI(IRPI,0)
            IIO=IPIOPI(IRPI,II)
            LOGION(IIO,ISTRA)=.TRUE.
            IF (LPAIO) THEN
              PAIO(IIO,IRD)= PAIO(IIO,IRD)+WTRSIG*PIOPI(IRPI,IIO)
              LMETSP(NSPAM+IIO)=.TRUE.
            END IF
          ENDDO
C SO NICHT  EAPL(IRD)     = EAPL(IRD)     +WTRSIG*E0
          DO IP=1,IPPLPI(IRPI,0)
            IPL=IPPLPI(IRPI,IP)
            LOGPLS(IPL,ISTRA)=.TRUE.
            IF (LPAPL) THEN
              PAPL(IPL,IRD)= PAPL(IPL,IRD)+WTRSIG*PPLPI(IRPI,IPL)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END IF
          ENDDO
58      CONTINUE
59      CONTINUE
C
C.........................................
C
C   MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C.........................................
C
C
C  CONTRIBUTIONS FROM ATOMS
C
C
        IF (LMAPL) THEN
          V0_PARB=VEL*(VELX*BXIN(IRDO)+VELY*BYIN(IRDO)+VELZ*BZIN(IRDO))
          PARMOM_0=V0_PARB*CNDYNA(IATM)
C     
          IF (LGACX(IATM,0,0).EQ.0) GOTO 159
          DO 156 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRDO,IPLS)) GOTO 156
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 156
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHA
C     
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            MAPL(IPLS,IRD)=MAPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*PARMOM(IPL,IRDO)
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLV,IRDO))
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
156       CONTINUE
159       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 161 IAEI=1,NAEII(IATM)
            IREI=LGAEI(IATM,IAEI)
            IF (PPLDS(IREI,0).GT.0) THEN
              DO 162 IPL=1,NPLSI
                P=PPLDS(IREI,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IREI)*P
C  NEW BULK ION IPL
                  IPLV=MPLSV(IPL)
                  MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                          SIGN(1._DP,BVIN(IPLV,IRDO))
                  LMETSP(NSPAMI+IPL)=.TRUE.
                ENDIF
162           CONTINUE
            ENDIF
161       CONTINUE
C
C  ION IMPACT IONIZATION CONTRIBUTION: NOT INCLUDED
C
C
C  ELASTIC CONTRIBUTION FROM ATOMS
C
          IF (LGAEL(IATM,0,0).EQ.0) GOTO 180
C  DEFAULT TRACKLENGTH ESTIMATOR (BGK APPROXIMATION)
          DO 181 IAEL=1,NAELI(IATM)
            IREL=LGAEL(IATM,IAEL,0)
            IPLS=LGAEL(IATM,IAEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 181
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 181
C
            WTRSIG=WTR*SIGVEL(IREL)
C     
            MAPL(IPLS,IRD)=MAPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IPL2=IPLS
            MAPL(IPL2,IRD)=MAPL(IPL2,IRD)+WTRSIG*PARMOM_0*
     .                     SIGN(1._DP,BVIN(IPLSV,IRDO))
            LMETSP(NSPAMI+IPL2)=.TRUE.
181       CONTINUE
180       CONTINUE
C
        END IF
51    CONTINUE
      RETURN
C
C
C  ESTIMATORS FOR MOLECULES
C
      ENTRY UPDMOL (XSTOR2,XSTORV2,IFLAG)
C
      WV=WEIGHT/VEL
      NPBGK=NPBGKM(IMOL)
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
      IF ((NPBGK.GT.0).AND.LBGKV)
     .   CALL UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)

      IF (IUPDTE == 2) RETURN

      IF (.NOT.ALLOCATED(CNDYNM)) THEN
        ALLOCATE (CNDYNM(NMOL))
        DO IML=1,NMOLI
          CNDYNM(IML)=AMUA*RMASSM(IML)
        END DO
      END IF
C
      VELQ=VEL*VEL
C
      DO 71 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTR*VEL*CNDYNM(IMOL)
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENM) EDENM(IMOL,IRD)=EDENM(IMOL,IRD)+WTRE0
        IF (LPDENM) PDENM(IMOL,IRD)=PDENM(IMOL,IRD)+WTR
        IF (LEDENM.OR.LPDENM) LMETSP(NSPA+IMOL)=.TRUE.

        IF (LVXDENM) VXDENM(IMOL,IRD)=VXDENM(IMOL,IRD)+WTRV*VELX
        IF (LVYDENM) VYDENM(IMOL,IRD)=VYDENM(IMOL,IRD)+WTRV*VELY
        IF (LVZDENM) VZDENM(IMOL,IRD)=VZDENM(IMOL,IRD)+WTRV*VELZ
        IF (LVXDENM.OR.LVYDENM.OR.LVZDENM) LMETSP(NSPA+IMOL)=.TRUE. 
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 71
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        IF (LPMML) PMML(IMOL,IRD)=PMML(IMOL,IRD)-WTRSIG
        IF (LEMML) EMML(IRD)     =EMML(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGMCX(IMOL,0,0).EQ.0) GOTO 79
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 76  IMCX=1,NMCXI(IMOL)
          IRCX=LGMCX(IMOL,IMCX,0)
          IPLS=LGMCX(IMOL,IMCX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LPMML.AND.(IESTCX(IRCX,1).NE.0)) THEN
            PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            IF (LPMPL) THEN
              PMPL(IPLS,IRD)=PMPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              IF (LPMAT) THEN
                PMAT(IAT1,IRD)= PMAT(IAT1,IRD)+WTRSIG
                LMETSP(NSPH+IAT1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              IF (LPMML) THEN
                PMML(IML1,IRD)= PMML(IML1,IRD)+WTRSIG
                LMETSP(NSPA+IML1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              IF (LPMIO) THEN
                PMIO(IIO1,IRD)= PMIO(IIO1,IRD)+WTRSIG
                LMETSP(NSPAM+IIO1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              IF (LPMPL) THEN
                PMPL(IPL1,IRD)= PMPL(IPL1,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL1)=.TRUE.
              END IF
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              IF (LPMAT) THEN
                PMAT(IAT2,IRD)= PMAT(IAT2,IRD)+WTRSIG
                LMETSP(NSPH+IAT2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              IF (LPMML) THEN
                PMML(IML2,IRD)= PMML(IML2,IRD)+WTRSIG
                LMETSP(NSPA+IML2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              IF (LPMIO) THEN
                PMIO(IIO2,IRD)= PMIO(IIO2,IRD)+WTRSIG
                LMETSP(NSPAM+IIO2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              IF (LPMPL) THEN
                PMPL(IPL2,IRD)= PMPL(IPL2,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL2)=.TRUE.
              END IF
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LEM) THEN
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEMML) EMML(IRD) = EMML(IRD) + WTRSIG*E0
            ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEMPL) EMPL(IRD) = EMPL(IRD) - WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
              IF (N1STX(IRCX,1).EQ.1) THEN
                IAT1=N1STX(IRCX,2)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LEMAT) EMAT(IRD) = EMAT(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.2) THEN
                IML1=N1STX(IRCX,2)
                LOGMOL(IML1,ISTRA)=.TRUE.
                IF (LEMML) EMML(IRD) = EMML(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.3) THEN
                IIO1=N1STX(IRCX,2)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LEMIO) EMIO(IRD) = EMIO(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LEMPL) EMPL(IRD) = EMPL(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ENDIF
C  SECOND SECONDARY: PREVIOUS MOLECULE IMOL
              IF (N2NDX(IRCX,1).EQ.1) THEN
                IAT2=N2NDX(IRCX,2)
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LEMAT) EMAT(IRD) = EMAT(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
                IML2=N2NDX(IRCX,2)
                LOGMOL(IML2,ISTRA)=.TRUE.
                IF (LEMML) EMML(IRD) = EMML(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
                IIO2=N2NDX(IRCX,2)
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LEMIO) EMIO(IRD) = EMIO(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LEMPL) EMPL(IRD) = EMPL(IRD) + WTRSIG*E0
              ENDIF
            ENDIF
          ENDIF
C
76      CONTINUE
79      CONTINUE
C
C  ELASTIC NEUTRAL BULK-ION COLLISION CONTRIBUTION
C
        IF (LGMEL(IMOL,0,0).EQ.0) GOTO 80
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 81  IMEL=1,NMELI(IMOL)
          IREL=LGMEL(IMOL,IMEL,0)
          IPLS=LGMEL(IMOL,IMEL,1)
C  DO NOT UPDATE BGK TALLIES HERE
          IBGK=NPBGKP(IPLS,1)
          IF (IBGK.NE.0) GOTO 81
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVEL(IREL)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LPMML.AND.(IESTEL(IREL,1)).NE.0) THEN
            PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
          ELSE
C  UPDATE TRACKLENGTH ESTIMATOR
C           IF (LPMPL) THEN
C           PMPL(IPLS,IRD)=PMPL(IPLS,IRD)-WTRSIG
C           PMPL(IPLS,IRD)=PMPL(IPLS,IRD)+WTRSIG
C           LMETSP(NSPAMI+IPLS)=.TRUE.
C           END IF
            IF (LPMML) THEN
              PMML(IMOL,IRD)=PMML(IMOL,IRD)+WTRSIG
              LMETSP(NSPA+IMOL)=.TRUE.
            END IF
          ENDIF
C
          IF (LEM) THEN
            IF (IESTEL(IREL,3).NE.0) THEN
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
              IF (LEMML) EMML(IRD)=EMML(IRD)+WTRSIG*E0
            ELSE
              IF (LEMPL) THEN
                EMPL(IRD)=EMPL(IRD)-WTRSIG*ESIGEL(IREL,1)
                EMPL(IRD)=EMPL(IRD)+WTRSIG*E0
              END IF
              IF (LEMML) EMML(IRD)=EMML(IRD)+WTRSIG*ESIGEL(IREL,1)
            ENDIF
          ENDIF
C
81      CONTINUE
80      CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGMEI(IMOL,0).EQ.0) GOTO 100
C
        DO 90 IMEI=1,NMDSI(IMOL)
          IREI=LGMEI(IMOL,IMEI)
          IF (SIGVEI(IREI).LE.0.D0) GOTO 90
C
          WTRSIG=WTR*SIGVEI(IREI)
C
C  COLLISION ESTIMATOR:
C  NOT AVAILABLE, HENCE: NO NEED FOR COMPENSATION HERE
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          IF (LPMEL) PMEL(IRD)=PMEL(IRD)+WTRSIG*PELDS(IREI)
C
C  POST COLLISION CONTRIBUTIONS
          DO IA=1,IPATDS(IREI,0)
            IAT=IPATDS(IREI,IA)
            LOGATM(IAT,ISTRA)=.TRUE.
            IF (LPMAT) THEN
              PMAT(IAT,IRD)=PMAT(IAT,IRD)+PATDS(IREI,IAT)*WTRSIG
              LMETSP(NSPH+IAT)=.TRUE.
            END IF
          END DO
          DO IM=1,IPMLDS(IREI,0)
            IML=IPMLDS(IREI,IM)
            LOGMOL(IML,ISTRA)=.TRUE.
            IF (LPMML) THEN
              PMML(IML,IRD)=PMML(IML,IRD)+PMLDS(IREI,IML)*WTRSIG
              LMETSP(NSPA+IML)=.TRUE.
            END IF
          END DO
          DO II=1,IPIODS(IREI,0)
            IIO=IPIODS(IREI,II)
            LOGION(IIO,ISTRA)=.TRUE.
            IF (LPMIO) THEN
              PMIO(IIO,IRD)=PMIO(IIO,IRD)+PIODS(IREI,IIO)*WTRSIG
              LMETSP(NSPAM+IIO)=.TRUE.
            END IF
          END DO
          DO IP=1,IPPLDS(IREI,0)
            IPL=IPPLDS(IREI,IP)
            LOGPLS(IPL,ISTRA)=.TRUE.
            IF (LPMPL) THEN
              PMPL(IPL,IRD)=PMPL(IPL,IRD)+PPLDS(IREI,IPL)*WTRSIG
              LMETSP(NSPAMI+IPL)=.TRUE.
            END IF
          END DO
C
          IF ((IESTEI(IREI,3).EQ.0) .AND. LEMEL)
     .      EMEL(IRD)=EMEL(IRD)+WTRSIG*ESIGEI(IREI,5)

          IF (LEM) THEN
            IF (IESTEI(IREI,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
              IF (LEMML) EMML(IRD)=EMML(IRD)+WTRSIG*E0
C
            ELSE
C  
              IF (LEMAT) EMAT(IRD)=EMAT(IRD)+WTRSIG*ESIGEI(IREI,1)
              IF (LEMML) EMML(IRD)=EMML(IRD)+WTRSIG*ESIGEI(IREI,2)
              IF (LEMIO) EMIO(IRD)=EMIO(IRD)+WTRSIG*ESIGEI(IREI,3)
              IF (LEMPL) EMPL(IRD)=EMPL(IRD)+WTRSIG*ESIGEI(IREI,4)
C
            ENDIF
          ENDIF
90      CONTINUE
100     CONTINUE
C
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C  CONTRIBUTIONS FROM MOLECULES
C
        IF (LMMPL) THEN
          V0_PARB=VEL*(VELX*BXIN(IRDO)+VELY*BYIN(IRDO)+VELZ*BZIN(IRDO))
          PARMOM_0=V0_PARB*CNDYNM(IMOL)
C
          IF (LGMCX(IMOL,0,0).EQ.0) GOTO 590
          DO 560 IMCX=1,NMCXI(IMOL)
            IRCX=LGMCX(IMOL,IMCX,0)
            IPLS=LGMCX(IMOL,IMCX,1)
            IF (LGVAC(IRDO,IPLS)) GOTO 560
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 560
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHM
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            MMPL(IPLS,IRD)=MMPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              MMPL(IPL,IRD)=MMPL(IPL,IRD)+WTRSIG*PARMOM(IPL,IRDO)
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              MMPL(IPL,IRD)=MMPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLV,IRDO))
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
560       CONTINUE
590       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 610 IMEI=1,NMDSI(IMOL)
            IREI=LGMEI(IMOL,IMEI)
            IF (PPLDS(IREI,0).GT.0) THEN
              DO 620 IPL=1,NPLSI
                P=PPLDS(IREI,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IREI)*P
C  NEW BULK ION IPL
                  IPLV=MPLSV(IPL)
                  MMPL(IPL,IRD)=MMPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                          SIGN(1._DP,BVIN(IPLV,IRDO))
                  LMETSP(NSPAMI+IPL)=.TRUE.
                ENDIF
620           CONTINUE
            ENDIF
610       CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM MOLECULES
C
C
          IF (LGMEL(IMOL,0,0).EQ.0) GOTO 800
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 810 IMEL=1,NMELI(IMOL)
            IREL=LGMEL(IMOL,IMEL,0)
            IPLS=LGMEL(IMOL,IMEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C     
            IF (IBGK.NE.0) GOTO 810
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 810
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            MMPL(IPLS,IRD)=MMPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IPL2=IPLS
            MMPL(IPL2,IRD)=MMPL(IPL2,IRD)+WTRSIG*PARMOM_0*
     .                     SIGN(1._DP,BVIN(IPLSV,IRDO))
            LMETSP(NSPAMI+IPL2)=.TRUE.
810       CONTINUE
800       CONTINUE

        END IF

71    CONTINUE
      RETURN
C
C
C  ESTIMATORS FOR TEST IONS
C
      ENTRY UPDION (XSTOR2,XSTORV2,IFLAG)
C
      WV=WEIGHT/VEL
      NPBGK=NPBGKI(IION)
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
      IF ((NPBGK.GT.0).AND.LBGKV)
     .   CALL UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)

      IF (IUPDTE == 2) RETURN

      IF (.NOT.ALLOCATED(CNDYNI)) THEN
        ALLOCATE (CNDYNI(NION))
        DO IIO=1,NIONI
          CNDYNI(IIO)=AMUA*RMASSI(IIO)
        END DO
      END IF
C
      VELQ=VEL*VEL
C
      DO 111 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTR*VEL*CNDYNI(IION)
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENI) EDENI(IION,IRD)=EDENI(IION,IRD)+WTRE0
        IF (LPDENI) PDENI(IION,IRD)=PDENI(IION,IRD)+WTR
        IF (LEDENI.OR.LPDENI) LMETSP(NSPAM+IION)=.TRUE.

        IF (LVXDENI) VXDENI(IION,IRD)=VXDENI(IION,IRD)+WTRV*VELX
        IF (LVYDENI) VYDENI(IION,IRD)=VYDENI(IION,IRD)+WTRV*VELY
        IF (LVZDENI) VZDENI(IION,IRD)=VZDENI(IION,IRD)+WTRV*VELZ
        IF (LVXDENI.OR.LVYDENI.OR.LVZDENI) LMETSP(NSPAM+IION)=.TRUE. 
C
C    ESTIMATORS FOR SOURCES AND SINKS
C    NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C    POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 111
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        IF (LPIIO) PIIO(IION,IRD)=PIIO(IION,IRD)-WTRSIG
        IF (LEIIO) EIIO(IRD)     =EIIO(IRD)     -WTRSIG*E0
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGICX(IION,0,0).EQ.0) GOTO 119
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 116  IICX=1,NICXI(IION)
          IRCX=LGICX(IION,IICX,0)
          IPLS=LGICX(IION,IICX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LPIIO.AND.IESTCX(IRCX,1).NE.0) THEN
            PIIO(IION,IRD)=PIIO(IION,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            IF (LPIPL) THEN
              PIPL(IPLS,IRD)=PIPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              IF (LPIAT) THEN
                PIAT(IAT1,IRD)= PIAT(IAT1,IRD)+WTRSIG
                LMETSP(NSPH+IAT1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              IF (LPIML) THEN
                PIML(IML1,IRD)= PIML(IML1,IRD)+WTRSIG
                LMETSP(NSPA+IML1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              IF (LPIIO) THEN
                PIIO(IIO1,IRD)= PIIO(IIO1,IRD)+WTRSIG
                LMETSP(NSPAM+IIO1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              IF (LPIPL) THEN
                PIPL(IPL1,IRD)= PIPL(IPL1,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL1)=.TRUE.
              END IF
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              IF (LPIAT) THEN
                PIAT(IAT2,IRD)= PIAT(IAT2,IRD)+WTRSIG
                LMETSP(NSPH+IAT2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              IF (LPIML) THEN
                PIML(IML2,IRD)= PIML(IML2,IRD)+WTRSIG
                LMETSP(NSPA+IML2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              IF (LPIIO) THEN
                PIIO(IIO2,IRD)= PIIO(IIO2,IRD)+WTRSIG
                LMETSP(NSPAM+IIO2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              IF (LPIPL) THEN
                PIPL(IPL2,IRD)= PIPL(IPL2,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL2)=.TRUE.
              END IF
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LEIO) THEN
            IF (LEIIO.AND.(IESTCX(IRCX,3).NE.0)) THEN
              EIIO(IRD) = EIIO(IRD) + WTRSIG*E0
            ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEIPL) EIPL(IRD)   = EIPL(IRD) - WTRSIG*ESIGCX(IRCX,1)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
              IF (N1STX(IRCX,1).EQ.1) THEN
                IAT1=N1STX(IRCX,2)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LEIAT) EIAT(IRD) = EIAT(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.2) THEN
                IML1=N1STX(IRCX,2)
                LOGMOL(IML1,ISTRA)=.TRUE.
                IF (LEIML) EIML(IRD) = EIML(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.3) THEN
                IIO1=N1STX(IRCX,2)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LEIIO) EIIO(IRD) = EIIO(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LEIPL) EIPL(IRD) = EIPL(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
              IF (N2NDX(IRCX,1).EQ.1) THEN
                IAT2=N2NDX(IRCX,2)
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LEIAT) EIAT(IRD) = EIAT(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
                IML2=N2NDX(IRCX,2)
                LOGMOL(IML2,ISTRA)=.TRUE.
                IF (LEIML) EIML(IRD) = EIML(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
                IIO2=N2NDX(IRCX,2)
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LEIIO) EIIO(IRD) = EIIO(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LEIPL) EIPL(IRD) = EIPL(IRD) + WTRSIG*E0
              ENDIF
            ENDIF
          ENDIF
C
116     CONTINUE
119     CONTINUE
C
C  ELECTRON IMPACT COLLISION CONTRIBUTION
C
        IF (LGIEI(IION,0).EQ.0) GOTO 130
C
        DO 120 IIEI=1,NIDSI(IION)
          IREI=LGIEI(IION,IIEI)
          IF (SIGVEI(IREI).LE.0.D0) GOTO 120
C
          WTRSIG=WTR*SIGVEI(IREI)
C
C  COLLISION ESTIMATOR:
C  NOT AVAILABLE, HENCE: NO NEED FOR COMPENSATION HERE
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
          IF (LPIEL) PIEL(IRD)=PIEL(IRD)+WTRSIG*PELDS(IREI)
C
C  POST COLLISION CONTRIBUTIONS
          DO IA=1,IPATDS(IREI,0)
            IAT=IPATDS(IREI,IA)
            LOGATM(IAT,ISTRA)=.TRUE.
            IF (LPIAT) THEN
              PIAT(IAT,IRD)=PIAT(IAT,IRD)+PATDS(IREI,IAT)*WTRSIG
              LMETSP(NSPH+IAT)=.TRUE.
            END IF
          END DO
          DO IM=1,IPMLDS(IREI,0)
            IML=IPMLDS(IREI,IM)
            LOGMOL(IML,ISTRA)=.TRUE.
            IF (LPIML) THEN
              PIML(IML,IRD)=PIML(IML,IRD)+PMLDS(IREI,IML)*WTRSIG
              LMETSP(NSPA+IML)=.TRUE.
            END IF
          END DO
          DO II=1,IPIODS(IREI,0)
            IIO=IPIODS(IREI,II)
            LOGION(IIO,ISTRA)=.TRUE.
            IF (LPIIO) THEN
              PIIO(IIO,IRD)=PIIO(IIO,IRD)+PIODS(IREI,IIO)*WTRSIG
              LMETSP(NSPAM+IIO)=.TRUE.
            END IF
          END DO
          DO IP=1,IPPLDS(IREI,0)
            IPL=IPPLDS(IREI,IP)
            LOGPLS(IPL,ISTRA)=.TRUE.
            IF (LPIPL) THEN
              PIPL(IPL,IRD)=PIPL(IPL,IRD)+PPLDS(IREI,IPL)*WTRSIG
              LMETSP(NSPAMI+IPL)=.TRUE.
            END IF
          END DO
C
          IF ((IESTEI(IREI,3).EQ.0) .AND. LEIEL) 
     .      EIEL(IRD)=EIEL(IRD)+WTRSIG*ESIGEI(IREI,5)

          IF (LEIO) THEN
            IF (IESTEI(IREI,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
              IF (LEIIO) EIIO(IRD)=EIIO(IRD)+WTRSIG*E0
C
            ELSE
C
              IF (LEIAT) EIAT(IRD)=EIAT(IRD)+WTRSIG*ESIGEI(IREI,1)
              IF (LEIML) EIML(IRD)=EIML(IRD)+WTRSIG*ESIGEI(IREI,2)
              IF (LEIIO) EIIO(IRD)=EIIO(IRD)+WTRSIG*ESIGEI(IREI,3)
              IF (LEIPL) EIPL(IRD)=EIPL(IRD)+WTRSIG*ESIGEI(IREI,4)
C     
            ENDIF
          ENDIF
120     CONTINUE
130     CONTINUE
C
C  PLASMA ION IMPACT CONTRIBUTION
C
C  NOT AVAILABLE FOR TEST IONS
C
C.........................................
C
C   MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C.........................................
C
C
C  CONTRIBUTIONS FROM TEST IONS
C
C
        IF (LMIPL) THEN
          V0_PARB=VEL*(VELX*BXIN(IRDO)+VELY*BYIN(IRDO)+VELZ*BZIN(IRDO))
          PARMOM_0=V0_PARB*CNDYNI(IION)
C
          IF (LGICX(IION,0,0).EQ.0) GOTO 5900
          DO 5600 IICX=1,NICXI(IION)
            IRCX=LGICX(IION,IICX,0)
            IPLS=LGICX(IION,IICX,1)
            IF (LGVAC(IRDO,IPLS)) GOTO 5600
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 5600
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHI
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST
            MIPL(IPLS,IRD)=MIPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              MIPL(IPL,IRD)=MIPL(IPL,IRD)+WTRSIG*PARMOM(IPL,IRDO)
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              IPLV=MPLSV(IPL)
              MIPL(IPL,IRD)=MIPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,BVIN(IPLV,IRDO))
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
5600      CONTINUE
5900      CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 6100 IIEI=1,NIDSI(IION)
            IREI=LGIEI(IION,IIEI)
            IF (PPLDS(IREI,0).GT.0) THEN
              DO 6200 IPL=1,NPLSI
                P=PPLDS(IREI,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IREI)*P
C  NEW BULK ION IPL
                  IPLV=MPLSV(IPL)
                  MIPL(IPL,IRD)=MIPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                          SIGN(1._DP,BVIN(IPLV,IRDO))
                  LMETSP(NSPAMI+IPL)=.TRUE.
                ENDIF
6200          CONTINUE
            ENDIF
6100      CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM TEST IONS
C
          IF (LGIEL(IION,0,0).EQ.0) GOTO 8000
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 8100 IIEL=1,NIELI(IION)
            IREL=LGIEL(IION,IIEL,0)
            IPLS=LGIEL(IION,IIEL,1)
            IPLSV=MPLSV(IPLS)
            IBGK=NPBGKP(IPLS,1)
C     
            IF (IBGK.NE.0) GOTO 8100
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 8100
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            MIPL(IPLS,IRD)=MIPL(IPLS,IRD)-WTRSIG*PARMOM(IPLS,IRDO)
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IPL2=IPLS
            MIPL(IPL2,IRD)=MIPL(IPL2,IRD)+WTRSIG*PARMOM_0*
     .                     SIGN(1._DP,BVIN(IPLSV,IRDO))
            LMETSP(NSPAMI+IPL2)=.TRUE.
8100      CONTINUE
8000      CONTINUE

        END IF
111   CONTINUE
      RETURN
C
C
C  ESTIMATORS FOR PHOTONS
C
      ENTRY UPDPHOT (XSTOR2,XSTORV2,IFLAG)
C
      WV=WEIGHT/VEL
c     NPBGK=NPBGKPH(IPHOT)
C
      IF (NADVI.GT.0) CALL UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
cdr generalise flag NPBGK to mean: model collision term for bi-linear collision
cdr only in this case: set backgound radiation intensity profiles.
cdr   IF (NPBGK.GT.0) CALL ....
C
      IF (IUPDTE == 2) RETURN

      CNDYNPH = EV_TO_ERG/CLIGHT
      VELQ=VEL*VEL
C
      DO 131 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTRE0*VEL*CNDYNPH
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENPH) EDENPH(IPHOT,IRD)=EDENPH(IPHOT,IRD)+WTRE0
        IF (LPDENPH) PDENPH(IPHOT,IRD)=PDENPH(IPHOT,IRD)+WTR
        IF (LEDENPH.OR.LPDENPH) LMETSP(IPHOT)=.TRUE.

        IF (LVXDENPH) VXDENPH(IPHOT,IRD)=VXDENPH(IPHOT,IRD)+WTRV*VELX
        IF (LVYDENPH) VYDENPH(IPHOT,IRD)=VYDENPH(IPHOT,IRD)+WTRV*VELY
        IF (LVZDENPH) VZDENPH(IPHOT,IRD)=VZDENPH(IPHOT,IRD)+WTRV*VELZ
        IF (LVXDENPH.OR.LVYDENPH.OR.LVZDENPH) LMETSP(IPHOT)=.TRUE. 
C
C  ESTIMATORS FOR SOURCES AND SINKS
C  NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C  POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 131
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES ARE LOST
C
        IF ((LAST_EVENT%IFLAG == 1) .AND.
     .      (LAST_EVENT%NCELL == IRD)) THEN

! collision estimator for first cell ("brick") along the track
! in case of a collision sample 1 (the whole weight)
! in case of no collision sample 0
          IF ((IFLAG == 4).OR.(IFLAG == 5)) THEN
            IF (LPPHPHT) PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)-WEIGHT
            IF (LEPHPHT) EPHPHT(IRD)      =EPHPHT(IRD)      -WEIGHT*E0
          ELSE
!  nothing to be done, sample a 0
          END IF

        ELSE
          WTRSIG=WTR*(SIGTOT-SIGBGK)
          IF (LPPHPHT) PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)-WTRSIG
          IF (LEPHPHT) EPHPHT(IRD)      =EPHPHT(IRD)      -WTRSIG*E0
        END IF
C
C  OTHER (OT) CONTRIBUTION for photons
cdr  in analogy to CX processes. better: analogy to PI ? wg. stim emission ?
C
        IF (PHV_LGPHOT(IPHOT,0,0).EQ.0) GOTO 133
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO IAOT=1,PHV_NPHOTI(IPHOT)
          IROT=PHV_LGPHOT(IPHOT,IAOT,0)
cdr  check irot=0 in initialisation phase only once
cdr       if(irot == 0) cycle
cdr
          IPLS=PHV_LGPHOT(IPHOT,IAOT,1)
          UPDF=PHV_LGPHOT(IPHOT,IAOT,4)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVOT(IROT)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (PHV_IESTOTph(iphot,IROT,1).NE.0) THEN
            IF (LPPHPHT) PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)+WTRSIG
cdr         if(updf==1) PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)+WTRSIG !prob. wrong
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
cdr  if(ipls > 0) then
cdr  do this check in initialisation, only once
            IF (LPPHPL) THEN
              PPHPL(IPLS,IRD)=PPHPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C
cdr-test    goto 133
!pb         IF (PHV_N1STOTph(iphot,IROT,3).NE.0) THEN
!pb PHV_N1STOTph(iphot,IROT,3) does not include bulk
C  FIRST SECONDARY:
              IF (PHV_N1STOTph(iphot,IROT,1).EQ.1) THEN
                IAT1=PHV_N1STOTph(iphot,IROT,2)
                INUM=PHV_N1STOTph(iphot,IROT,3)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LPPHAT) THEN
                  PPHAT(IAT1,IRD)= PPHAT(IAT1,IRD)+WTRSIG*INUM
                  LMETSP(NSPH+IAT1)=.TRUE.
                END IF
              ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.2) THEN
                IML1=PHV_N1STOTph(iphot,IROT,2)
                LOGMOL(IML1,ISTRA)=.TRUE.
                IF (LPPHML) THEN
                  PPHML(IML1,IRD)= PPHML(IML1,IRD)+WTRSIG
                  LMETSP(NSPA+IML1)=.TRUE.
                END IF
              ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.3) THEN
                IIO1=PHV_N1STOTph(iphot,IROT,2)
                INUM=PHV_N1STOTph(iphot,IROT,3)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LPPHIO) THEN
                  PPHIO(IIO1,IRD)= PPHIO(IIO1,IRD)+WTRSIG*INUM
                  LMETSP(NSPAM+IIO1)=.TRUE.
                END IF
              ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.4) THEN
                IPL1=PHV_N1STOTph(iphot,IROT,2)
C               INUM=PHV_N1STOTph(iphot,IROT,3)
                INUM=1
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LPPHPL) THEN
                  PPHPL(IPL1,IRD)= PPHPL(IPL1,IRD)+WTRSIG*INUM
csw added updf check (stim.em)
cdr  stim emission: am besten: 2 secondaries in group 1. hier jedoch:
cdr  dazu PI-process vervollstandigen.
cdr  test photon + ipls --> ipl2 (bulk particle)
cdr  if this bulk is a photon (same as iphot), dann noch eins raus auf tally
                  LMETSP(NSPAMI+IPL1)=.TRUE.
                END IF
csw added branch
              ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.0) then
                iph1=phv_n1stotph(iphot,irot,2)
                inum=phv_n1stotph(iphot,irot,3)
cdr  test iph1 > 0 only once, in initialisation. here: removed
                if ((inum > 0) .and. (iph1 > 0)) then
                  logphot(iph1,istra)=.true.
                  IF (LPPHPHT) THEN
                    PPHPHT(iph1,ird)=PPHPHT(iph1,ird)+wtrsig*inum
                    LMETSP(IPH1)=.TRUE.
                  END IF
                end if
              ENDIF
!pb         END IF

!pb         IF (PHV_N2NDOTph(iphot,IROT,3).NE.0) THEN
!pb PHV_N2NDOTph(iphot,IROT,3) does not include bulk
C  SECOND SECONDARY:
              IF (PHV_N2NDOTph(iphot,IROT,1).EQ.1) THEN
                IAT2=PHV_N2NDOTph(iphot,IROT,2)
                INUM=PHV_N2NDOTph(iphot,IROT,3)
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LPPHAT) THEN
                  PPHAT(IAT2,IRD)= PPHAT(IAT2,IRD)+WTRSIG*INUM
                  LMETSP(NSPH+IAT2)=.TRUE.
                END IF
              ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.2) THEN
                IML2=PHV_N2NDOTph(iphot,IROT,2)
                LOGMOL(IML2,ISTRA)=.TRUE.
                IF (LPPHML) THEN
                  PPHML(IML2,IRD)= PPHML(IML2,IRD)+WTRSIG
                  LMETSP(NSPA+IML2)=.TRUE.
                END IF
              ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.3) THEN
                IIO2=PHV_N2NDOTph(iphot,IROT,2)
                INUM=PHV_N2NDOTph(iphot,IROT,3)
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LPPHIO) THEN
                  PPHIO(IIO2,IRD)= PPHIO(IIO2,IRD)+WTRSIG*INUM
                  LMETSP(NSPAM+IIO2)=.TRUE.
                END IF
              ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.4) THEN
                IPL2=PHV_N2NDOTph(iphot,IROT,2)
C               INUM=PHV_N2NDOTph(iphot,IROT,3)
                INUM=1
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LPPHPL) THEN
                  PPHPL(IPL2,IRD)= PPHPL(IPL2,IRD)+WTRSIG*INUM
csw added updf check (stim.em)
cdr  stim emission: am besten: 2 secondaries in group 2.
cdr  dazu PI-process vervollstandigen.
cdr  test photon + ipls --> ipl2 (bulk particle)
cdr  if this bulk is a photon (same as iphot), dann noch eins rauf auf tally
                  LMETSP(NSPAMI+IPL2)=.TRUE.
                END IF
csw added branch
              ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.0) then
                iph2=phv_n2ndotph(iphot,irot,2)
                inum=phv_n2ndotph(iphot,irot,3)
cdr test iph2 > 0 removed, to be done only once in initialisation
                if ((inum > 0) .and. (iph2 > 0)) then
                  logphot(iph2,istra)=.true.
                  IF (LPPHPHT) THEN
                    PPHPHT(iph2,ird)=PPHPHT(iph2,ird)+wtrsig*inum
                    LMETSP(iph2)=.true.
                  END IF
                END IF
              ENDIF
!pb         ENDIF
          ENDIF
C
C  PARTICLE ESTIMATORS DONE. NEXT: ENERGY ESTIMATORS
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LEPHPHT.AND.(PHV_IESTOTph(iphot,IROT,3).NE.0)) THEN
            EPHPHT(IRD)=EPHPHT(IRD)          +WTRSIG*E0
cdr         if(updf==1) EPHPHT(IRD)=EPHPHT(IRD)+WTRSIG*E0 ! verm. falsch
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
cdr if(ipls > 0) then : this check only once in initialisation. removed
cdr         IF (LEPHEL) EPHPL(IRD)    =EPHPL(IRD)  -WTRSIG*E0  vermutl. falsch
cdr  gibt es schon ESIGOT ? ist dann IROT das richtige argument ?
cdr         IF (LEPHPL) EPHPL(IRD)    =EPHPL(IRD)  -WTRSIG*ESIGOT(IROT)
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)

!dr       IF (PHV_N1STOTph(iphot,IROT,3).NE.0) THEN
!dr PHV_N1STOTph(iphot,IROT,3) does not include bulk
C  FIRST SECONDARY:
cdr: erstmal raus. vermutlich muss hier die energie des
cdr  neuen schwerteilchens stehen, nicht die E0 des IPHOT.
cdr         IF (PHV_N1STOTph(iphot,IROT,1).EQ.1) THEN
cdr           IAT1=PHV_N1STOTph(iphot,IROT,2)
cdr           LOGATM(IAT1,ISTRA)=.TRUE.
cdr           IF (LEPHAT) EPHAT(IRD)     = EPHAT(IRD)  +WTRSIG*E0
C           ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.2) THEN
C             IML1=PHV_N1STOTph(iphot,IROT,2)
C             LOGMOL(IML1,ISTRA)=.TRUE.
C             IF (LEPHML) EPHML(IRD)     = EPHML(IRD)     +WTRSIG*E0
cdr         ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.3) THEN
cdr           IIO1=PHV_N1STOTph(iphot,IROT,2)
cdr           LOGION(IIO1,ISTRA)=.TRUE.
cdr           IF (LEPHIO) EPHIO(IRD)     = EPHIO(IRD)     +WTRSIG*E0
cdr         ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.4) THEN
cdr           IPL1=PHV_N1STOTph(iphot,IROT,2)
cdr           LOGPLS(IPL1,ISTRA)=.TRUE.
cdr           IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)   +WTRSIG*E0
csw added updf check (stim.em)
cdr           if(LEPHPL.AND.(updf==1.and.phv_is_plsphot(ipl1)>0)) then
cdr              EPHPL(IRD)     = EPHPL(IRD)+WTRSIG*E0
cdr           endif
csw added branch
cdr: photon scattering, e0_in = e0_out, also: delta scattering,
cdr  mit eventuell teilchenmultiplikation
cdr         ELSEIF (PHV_N1STOTph(iphot,IROT,1).EQ.0) THEN
              IF (PHV_N1STOTph(iphot,IROT,1).EQ.0) THEN
                IPH1=PHV_N1STOTph(iphot,IROT,2)
                INUM=PHV_N1STOTph(iphot,IROT,2)
cdr  if(iph1 > 0) then abfrage hier raus, nur in initialisation
                LOGPHOT(IPH1,ISTRA)=.TRUE.
                IF (LEPHPHT) EPHPHT(IRD)=EPHPHT(IRD) +WTRSIG*E0*INUM
              ENDIF
!dr         ENDIF

!dr     IF (PHV_N2NDOTph(iphot,IROT,3).NE.0) THEN
!dr PHV_N2NDOTph(iphot,IROT,3) does not include bulk
C  SECOND SECONDARY:
cdr : gleiches problem wie oben E0_out ? bei test photon impact mit E0 ?
cdr         IF (PHV_N2NDOTph(iphot,IROT,1).EQ.1) THEN
cdr           IAT2=PHV_N2NDOTph(iphot,IROT,2)
cdr           LOGATM(IAT2,ISTRA)=.TRUE.
cdr           IF (LEPHAT) EPHAT(IRD)     = EPHAT(IRD)     +WTRSIG*E0
C           ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.2) THEN
C             IML2=PHV_N2NDOTph(iphot,IROT,2)
C             LOGMOL(IML2,ISTRA)=.TRUE.
C             IF (LEPHML) EPHML(IRD)     = EPHML(IRD)     +WTRSIG*E0
cdr         ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.3) THEN
cdr           IIO2=PHV_N2NDOTph(iphot,IROT,2)
cdr           LOGION(IIO2,ISTRA)=.TRUE.
cdr           IF (LEPHIO) EPHIO(IRD)     = EPHIO(IRD)     +WTRSIG*E0
cdr         ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.4) THEN
cdr           IPL2=PHV_N2NDOTph(iphot,IROT,2)
cdr           LOGPLS(IPL2,ISTRA)=.TRUE.
cdr           IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)     +WTRSIG*E0
csw added updf check (stim.em)
cdr           if(updf==1.and.phv_is_plsphot(ipl2)>0) then
cdr              IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)     +WTRSIG*E0
cdr           endif
csw added branch
cdr         ELSEIF (PHV_N2NDOTph(iphot,IROT,1).EQ.0) THEN
              IF (PHV_N2NDOTph(iphot,IROT,1).EQ.0) THEN
                IPH2=PHV_N2NDOTph(iphot,IROT,2)
                INUM=PHV_N2NDOTph(iphot,IROT,3)
cdr if(iph2 > 0) then  ! dieser test nur in initialisation phase
                if ((inum > 0) .and. (iph2 > 0)) then
                  LOGPHOT(IPH2,ISTRA)=.TRUE.
                  IF (LEPHPHT) EPHPHT(IRD) = EPHPHT(IRD)+WTRSIG*E0*INUM
                end if
              ENDIF
!dr         ENDIF
          ENDIF
C
       ENDDO
133    CONTINUE
131    CONTINUE
      RETURN
      END
C ===== SOURCE: update_surface.f

      SUBROUTINE UPDATE_SURFACE (ITOLD)

      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMPRT
      USE CSPEZ

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITOLD

      IF (MSURF .LE. 0) RETURN

      IF (ITYP.EQ.0) THEN
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        IF (ITOLD.EQ.0) THEN
          IF (LPRFPHPHT)
     .      PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WEIGHT
          IF (LERFPHPHT)
     .      ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFPHPHT)
     .        PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WEIGHT
            IF (LERFPHPHT)
     .        ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)
     .                               +E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WEIGHT
          IF (LERFAAT)
     .      ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAAT)
     .        PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WEIGHT
            IF (LERFAAT)
     .        ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMAT) PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WEIGHT
          IF (LERFMAT)
     .      ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMAT)
     .        PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WEIGHT
            IF (LERFMAT)
     .        ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIAT) PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WEIGHT
          IF (LERFIAT)
     .      ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIAT)
     .        PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WEIGHT
            IF (LERFIAT)
     .        ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPAT)
     .      PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WEIGHT
          IF (LERFPAT)
     .      ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFPAT)
     .        PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WEIGHT
            IF (LERFPAT)
     .        ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAML) PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WEIGHT
          IF (LERFAML)
     .      ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAML)
     .        PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WEIGHT
            IF (LERFAML)
     .        ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WEIGHT
          IF (LERFMML)
     .      ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMML)
     .        PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WEIGHT
            IF (LERFMML)
     .        ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIML) PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WEIGHT
          IF (LERFIML)
     .      ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIML)
     .        PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WEIGHT
            IF (LERFIML)
     .        ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPML)
     .      PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WEIGHT
          IF (LERFPML)
     .      ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFPML)
     .        PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WEIGHT
            IF (LERFPML)
     .        ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAIO) PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WEIGHT
          IF (LERFAIO)
     .      ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAIO)
     .        PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WEIGHT
            IF (LERFAIO)
     .        ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMIO) PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WEIGHT
          IF (LERFMIO)
     .      ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMIO)
     .        PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WEIGHT
            IF (LERFMIO)
     .        ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WEIGHT
          IF (LERFIIO)
     .      ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIIO)
     .        PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WEIGHT
            IF (LERFIIO)
     .        ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPIO)
     .      PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WEIGHT
          IF (LERFPIO)
     .      ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFPIO)
     .        PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WEIGHT
            IF (LERFPIO)
     .        ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ENDIF
      
      RETURN
      END SUBROUTINE UPDATE_SURFACE
