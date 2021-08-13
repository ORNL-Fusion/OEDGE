CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  TEXTOR  *  (HELIUM INCLUDED)
C           ************
C
      SUBROUTINE PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CGRID
      USE CGEOM
      USE CINIT
      USE CCONA
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: NTR, I, J, LL
      REAL(DP), ALLOCATABLE, SAVE :: PLAS(:,:)
      REAL(DP) :: bnorm
      INTEGER, SAVE :: INDAR(9)=(/ (0, i=1,9) /) 
      
      IF (.NOT.ALLOCATED(PLAS)) THEN

        LL=LEN_TRIM(CASENAME)

        FILENAME=CASENAME(1:LL) // '.plasma'
        OPEN (UNIT=31,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')

        ZEILE='*   '      
        DO WHILE (ZEILE(1:1) == '*')
           READ (31,'(A100)') ZEILE
        END DO

        READ (ZEILE,*) NTR

        IF (NTR /= NR1ST-1) THEN
          WRITE (6,*) ' WRONG NUMBER OF TRIANGLES IN PLASMA FILE'
          WRITE (6,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
          CALL EXIT(1)
        END IF

!        ALLOCATE (PLAS(2+4*NPLS,NRAD))
        ALLOCATE (PLAS(9,NRAD))
        PLAS=0._DP
        DO I=1, NTR
          READ (31,*) J, PLAS(:,I)
          bnorm=sqrt(plas(7,i)**2 + plas(8,i)**2 + plas(9,I)**2)
          if (bnorm < eps12) then
             plas(7,i) = 0._dp
             plas(8,i) = 0._dp
             plas(9,i) = 1._dp
          else
             plas(7:9,i) = plas(7:9,i)/bnorm
          end if
        END DO
        plas(9,ntr+1:nrad) = 1._dp
      END IF
      
      IF (INDX <= 1) THEN
        PRO(1:N) = PLAS(INDX+1,1:N)
        INDAR(INDX+1) = INDAR(INDX+1)+1
      ELSEIF (INDX == 1+1*NPLS) THEN
        IF (INDAR(3) == 0) THEN
           PRO(1:N) = PLAS(3,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(3) = INDAR(3) + 1
      ELSEIF (INDX == 1+2*NPLS) THEN
        IF (INDAR(4) == 0) THEN
           PRO(1:N) = PLAS(4,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(4) = INDAR(4) + 1
      ELSEIF (INDX == 1+3*NPLS) THEN
        IF (INDAR(5) == 0) THEN
           PRO(1:N) = PLAS(5,1:N)     
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(5) = INDAR(5) + 1
      ELSEIF (INDX == 1+4*NPLS) THEN
        IF (INDAR(6) == 0) THEN
           PRO(1:N) = PLAS(6,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(6) = INDAR(6) + 1
      ELSEIF (INDX == 1+5*NPLS) THEN
        IF (INDAR(7) == 0) THEN
           PRO(1:N) = PLAS(7,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(7) = INDAR(7) + 1
      ELSEIF (INDX == 2+5*NPLS) THEN
        IF (INDAR(8) == 0) THEN
           PRO(1:N) = PLAS(8,1:N)
        ELSE
           PRO(1:N) = 0._DP
        END IF
        INDAR(8) = INDAR(8) + 1
      ELSEIF (INDX == 3+5*NPLS) THEN
        IF (INDAR(9) == 0) THEN
           PRO(1:N) = PLAS(9,1:N)
        ELSE
           PRO(1:N) = 1._DP
        END IF
        INDAR(9) = INDAR(9) + 1
      ELSE
        WRITE (6,*) ' PROUSR: NO DATA PROVIDED FOR INDEX ',INDX
        PRO(1:N) = 0._DP
      END IF

      RETURN
      END
