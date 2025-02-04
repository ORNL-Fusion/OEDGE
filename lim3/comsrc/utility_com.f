C     -*-Fortran-*-
C
C
C
C  *********************************************************************
C  *  PRQ:  PRINTS A REAL*8 NUMBER                                    *
C  *********************************************************************
C
      SUBROUTINE PRQ (NAME, R)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      real*8      R
      IF (ABS(R).LT.0.1.OR.ABS(R).GE.1000.0) THEN
        WRitE (datunit,'(1X,A,1P,G15.8)') NAME,R
      ELSEIF (ABS(R).LT.1.0) THEN
        WRITE (datunit,'(1X,A,F15.8)') NAME,R
      ELSE
        WRITE (datunit,'(1X,A,F15.8)') NAME,R
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRQ2: PRINTS TWO REAL*8 NUMBERS                                 *
C  *********************************************************************
C
      SUBROUTINE PRQ2 (NAME, R1, R2)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      REAL*8      R1,R2
      IF (ABS(R1).LT.0.01.OR.ABS(R1).GE.1000.0.OR.
     >    ABS(R2).LT.0.01.OR.ABS(R2).GE.1000.0) THEN
        WRITE (datunit,'(1X,A,1P,G11.3,1x,g11.3)') NAME,R1,R2
      ELSEIF (ABS(R1).LT.1.0.AND.ABS(R2).LT.1.0) THEN
        WRITE (datunit,'(1X,A,F8.4,3X,F8.4)') NAME,R1,R2
      ELSE
        WRITE (datunit,'(1X,A,F8.3,3X,F8.3)') NAME,R1,R2
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRR0:  PRINTS A REAL NUMBER                                      *
C  *********************************************************************
C
      SUBROUTINE PRR0(NAME, R)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      REAL      R
      IF ((ABS(R).NE.0.0.AND.ABS(R).LT.0.001).OR.ABS(R).GE.1.0E+06) THEN
        WRITE (datunit,'(1X,A,1P,G11.3)') NAME,R
      ELSEIF (ABS(R).NE.0.0E+00) THEN
        WRITE (datunit,'(1X,A,F10.3)') NAME,R
      ELSE
        WRITE (datunit,'(1X,A,A)') NAME,'     0    '
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRR:  PRINTS A REAL NUMBER                                       *
C  *********************************************************************
C
      SUBROUTINE PRR (NAME, R)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      REAL      R
      IF (ABS(R).LT.0.1.OR.ABS(R).GE.1000.0) THEN
        WRITE (datunit,'(1X,A,1P,G11.3)') NAME,R
      ELSEIF (ABS(R).LT.1.0) THEN
        WRITE (datunit,'(1X,A,F8.4)') NAME,R
      ELSE
        WRITE (datunit,'(1X,A,F8.3)') NAME,R
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRR2: PRINTS TWO REAL NUMBERS                                    *
C  *********************************************************************
C
      SUBROUTINE PRR2 (NAME, R1, R2)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      REAL      R1,R2
      IF (ABS(R1).LT.0.1.OR.ABS(R1).GE.1000.0.OR.
     >    ABS(R2).LT.0.1.OR.ABS(R2).GE.1000.0) THEN
        WRITE (datunit,'(1X,A,1P,G11.3,1x,g11.3)') NAME,R1,R2
      ELSEIF (ABS(R1).LT.10.0.AND.ABS(R2).LT.10.0) THEN
        WRITE (datunit,'(1X,A,F8.4,3X,F8.4)') NAME,R1,R2
      ELSE
        WRITE (datunit,'(1X,A,F8.3,3X,F8.3)') NAME,R1,R2
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRR3: PRINTS THREE REAL NUMBERS                                  *
C  *********************************************************************
C
      SUBROUTINE PRR3 (NAME, R1, R2 , R3 )
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      REAL      R1,R2,r3
      IF (ABS(R1).LT.0.1.OR.ABS(R1).GE.1000.0.OR.
     >    ABS(R2).LT.0.1.OR.ABS(R2).GE.1000.0.or.
     >    ABS(R3).LT.0.1.OR.ABS(R3).GE.1000.0) THEN
        WRITE (datunit,'(1X,A,1P,3(G10.3,1x))') NAME,R1,R2,R3
      ELSEIF (ABS(R1).LT.1.0.AND.ABS(R2).LT.1.0.AND.ABS(R3).LT.1.0) THEN
        WRITE (datunit,'(1X,A,3(F10.4,1x))') NAME,R1,R2,R3
      ELSE
        WRITE (datunit,'(1X,A,3(F10.4,1x))') NAME,R1,R2,R3
      ENDIF
      RETURN
      END
C
C  *********************************************************************
C  *  PRRMATDIV: PRINTS A 2-DIMIENSIONAL REAL ARRAY                    *
C  *********************************************************************
C
      SUBROUTINE PRRMATDIV(A,IDIMA,IDIM1,IDIM2,IWT,TIT)
      implicit none
c
c      IMPLICIT REAL (A-H,O-Z)
C PRINT OF REAL MATRIX A
C INPUT
C -------
C BY ARGUMENT-LIST:
C    A        - MATRIX OF DIMENSION A(IDIMA,IDIM2)
C               FIRST DIMENSION IS OCCUPIED ONLY WITH IDIM1 ELEMENTS
C    IDIMA    - LEADING DIMENSION OF A
C    IDIM1    - NUMBER OF ROWS OF A
C    IDIM2    - NUMBER OF COLUMNS OF A
C    IWT      - OUTPUT-CHANNEL FOR PRINTOUT
C    TIT      - CHARACTER STRING FOR TITLE
C=======================================================================
      integer idima,idim1,idim2,iwt
      real A(IDIMA,IDIM2)
      CHARACTER*(*) TIT
c
c     Local variables
c
      integer inum
      DATA INUM/9/
C INUM = NUMBER OF COLUMNS ON PAGE
c
      integer i1,i2,in,nr,j,n,n1,i,ii

C
      WRITE(IWT,50)TIT
C
      IF(IDIM2.GT.1)GOTO 20
         I1=1
         I2=0
         NR=IDIM1
         IN=(IDIM1-1)/INUM+1
         DO 10 II=1,IN
            I2=I2+INUM
            IF(NR.LT.INUM)I2=I2-INUM+NR
            NR=NR-INUM
            WRITE(IWT,80)(I,I=I1,I2)
            WRITE(IWT,80)
            WRITE(IWT,70)(A(I,1),I=I1,I2)
   10       I1=I1+INUM
         GOTO 90
C
   20 N1=1
      N=0
      NR=IDIM2
      IN=(IDIM2-1)/INUM +1
      DO 40 II=1,IN
         N=N+INUM
         IF (NR .LT. INUM ) N=N-INUM+NR
         NR=NR-INUM
         WRITE(IWT,80)(J,J=N1,N)
         WRITE(IWT,80)
         DO 30 I=1,IDIM1
   30       WRITE(IWT,60)I,(A(I,J),J=N1,N)
         N1=N1+INUM
   40    CONTINUE
C
   50 FORMAT(////,1X,A/1X,132('-'))
   60 FORMAT(1X,I4,   2X,1P,9E13.5)
   70 FORMAT(7X,1P,9E13.5)
   80 FORMAT(/,11X, 9(I4,9X))
   90 RETURN
      END
C
C  *********************************************************************
C  *  PRI:  PRINTS AN INTEGER                                          *
C  *********************************************************************
C
      SUBROUTINE PRI (NAME, I)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      INTEGER   I
      if (i.gt.9.99e6) then 
         WRITE (datunit,'(1X,A,I12)') NAME,I
      else
         WRITE (datunit,'(1X,A,I7)') NAME,I
      endif
      RETURN
      END
C
C  *********************************************************************
C  *  PRI2: PRINTS TWO INTEGERS                                        *
C  *********************************************************************
C
      SUBROUTINE PRI2 (NAME, I1, I2)
      use mod_params
      implicit none
c      include 'params' 
      CHARACTER NAME*(*)
      INTEGER   I1,I2
      if (i1.gt.9.999e6.or.i2.gt.9.99e6) then 
         WRITE (datunit,'(1X,A,I12,4X,I12)') NAME,I1,I2
      else
         WRITE (datunit,'(1X,A,I7,4X,I7)') NAME,I1,I2
      endif
      RETURN
      END

C
C  *********************************************************************
C  *  PRC:  PRINTS A CHARACTER STRING                                  *
C  *********************************************************************
C
      SUBROUTINE PRC(STRING)
      use mod_params
      implicit none
c      include 'params' 
c
      integer stringlen,lenstr
      external lenstr
      CHARACTER STRING*(*)
c
      stringlen = lenstr(string)
c
      WRITE (datunit,'(1X,A)') STRING(1:stringlen)
      RETURN
      END
C
C  *********************************************************************
C  *  PRB:  PRINTS A BLANK LINE                                        *
C  *********************************************************************
C
      SUBROUTINE PRB
      use mod_params
      implicit none
c      include 'params' 
      WRITE (datunit,'(1X)')
      RETURN
      END
C
C  *********************************************************************
C  *  PRP:  PRINTS A PAGE THROW                                        *
C  *********************************************************************
C
      SUBROUTINE PRP
      use mod_params
      implicit none
c      include 'params' 
      WRITE (datunit,'(''1'')')
      RETURN
      END
C
C  *********************************************************************
C  *  RINOUT: READS IN / WRITES OUT AN UNFORMATTED ARRAY OF REALS.     *
C  *  THE ARRAYS ARE READ/WRITTEN ON CHANNEL 8, TO A DATASET WITH      *
C  *  ATTRIBUTES BLKSIZE=6160, RECFM=VBS, LREC=6160, TRKS=(20,20)      *
C  *  OPT(1:1) SHOULD BE 'R' OR 'W', AND OPT(3:8) IS THE ARRAY NAME    *
C  *  (USED IN WRITE STATEMENT AT END OF ROUTINE).                     *
C  *                                                                   *
C  *          CHRIS FARRELL    MARCH 1989                              *
C  *********************************************************************
C
      SUBROUTINE RINOUT (OPT,RARRAY,N)
      implicit none
      INTEGER I,J,N,IBLOCK,ierr,len,lenstr
      external lenstr
      CHARACTER OPT*(*)
      REAL RARRAY(N)
      DATA IBLOCK /1500/
C
      IF     (OPT(1:1).EQ.'R') THEN
        DO 100 I = 1, N, IBLOCK
          READ  (8,ERR=300,iostat=ierr)
     >          (RARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (8,ERR=400,iostat=ierr)
     >          (RARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF
 
      len = lenstr(opt)

      WRITE (6,9001) OPT(3:len),REAL(4*N)
c      WRITE (0,9001) OPT(3:len),REAL(4*N)

      return

  300 Write (0,*) 'ERROR READING: ',OPT,' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',OPT,' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',OPT,' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',OPT,' : ERROR=',ierr
      return

 9001 FORMAT(1X,'RINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
      END
C
C  *********************************************************************
C  *  DINOUTU: READS IN / WRITES OUT AN UNFORMATTED ARRAY OF REALS.
C  *  THE ARRAYS ARE READ/WRITTEN ON CHANNEL IONUM, TO A DATASET WITH  *
C  *  ATTRIBUTES BLKSIZE=6160, RECFM=VBS, LREC=6160, TRKS=(20,20)      *
C  *  OPT(1:1) SHOULD BE 'R' OR 'W', AND OPT(3:8) IS THE ARRAY NAME    *
C  *  (USED IN WRITE STATEMENT AT END OF ROUTINE).                     *
C  *                                                                   *
C  *          CHRIS FARRELL    MARCH 1989                              *
C  *********************************************************************
C
      SUBROUTINE DINOUTU (OPT,DARRAY,N,IONUM)
      implicit none
      INTEGER I,J,N,IBLOCK,IONUM,ierr
      CHARACTER OPT*(*)
      REAL*8 DARRAY(N)
      DATA IBLOCK /750/
C
      IF     (OPT(1:1).EQ.'R') THEN
        DO 100 I = 1, N, IBLOCK
          READ  (IONUM,ERR=300,iostat=ierr)
     >          (DARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (IONUM,ERR=400,iostat=ierr)
     >          (DARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF

      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return
C      IF (4*N.GT.10000) WRITE (6,9001) OPT(3:8),REAL(4*N)
C 9001 FORMAT(1X,'RINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
C
      RETURN
      END
C
C  *********************************************************************
C  *  DINOUT: WRITE ONLY ROUTINE, CONVERTS D.P. TO REAL WHEN WRITING.  *
C  *********************************************************************
C
      SUBROUTINE DINOUT (OPT,DARRAY,N)
      implicit none
      INTEGER I,J,N,IBLOCK,ierr
      CHARACTER OPT*(*)
      DOUBLE PRECISION DARRAY(N)
      DATA IBLOCK /1500/
C
      IF     (OPT(1:1).EQ.'R') THEN
        WRITE (6,*) ' DINOUT: ERROR!  USE ONLY FOR WRITING!'
        STOP
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (8,ERR=400,iostat=ierr)
     >          (SNGL(DARRAY(J)),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF
      IF (4*N.GT.10000) WRITE (6,9001) OPT(3:len_trim(opt)),REAL(4*N)

      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return


 9001 FORMAT(1X,'DINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
      END
C
C  *********************************************************************
C  *  R8INOUT: READ/WRITE ROUTINE FOR R*8 - NO CONVERSION              *
C  *********************************************************************
C
      SUBROUTINE R8INOUT (OPT,DARRAY,N)
      implicit none
      INTEGER I,J,N,IBLOCK,ierr
      CHARACTER OPT*(*)
      real*8 DARRAY(N)
      DATA IBLOCK /1500/
C
      IF     (OPT(1:1).EQ.'R') THEN
        DO 100 I = 1, N, IBLOCK
          READ  (8,ERR=300,iostat=ierr)
     >          (DARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (8,ERR=400,iostat=ierr)
     >          (DARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF
      IF (8*N.GT.10000) WRITE (6,9001) OPT(3:len_trim(opt)),REAL(8*N)

      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return


 9001 FORMAT(1X,'R8INOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
      END
C
C  *********************************************************************
C  *  IINOUT: WRITE / READ INTEGER ARRAY,  SIMILAR TO RINOUT.          *
C  *********************************************************************
C
      SUBROUTINE IINOUT (OPT,IARRAY,N)
      implicit none
      INTEGER I,J,N,IBLOCK,IARRAY(N),ierr
      CHARACTER OPT*(*)
      DATA IBLOCK /1500/
C
      IF     (OPT(1:1).EQ.'R') THEN
        DO 100 I = 1, N, IBLOCK
          READ  (8,ERR=300,iostat=ierr)
     >          (IARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        DO 200 I = 1, N, IBLOCK
          WRITE (8,ERR=400,iostat=ierr)
     >          (IARRAY(J),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE
      ENDIF
      WRITE (6,9001) OPT(3:len_trim(opt)),REAL(4*N)
      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return

 9001 FORMAT(1X,'IINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
      END
C
C
C  *********************************************************************
C  *  IINOUT2: WRITE / READ INTEGER ARRAY,  SIMILAR TO IINOUT.         *
C  *           Except done by elements                                 *
C  *********************************************************************
C
      SUBROUTINE IINOUT2 (OPT,IARRAY,M,N,L,U)
      implicit none
      CHARACTER OPT*(*)
      INTEGER M,N,l,u,ierr
      integer IARRAY(L,U)
c
      integer tot,s1,s2,iblock,i,j,k,cnt
      parameter (IBLOCK=1500)
      integer tmparray(iblock)
C
      tot = m * n
      IF     (OPT(1:1).EQ.'R') THEN
        s1 = 1
        s2 = 1
        DO 100 I = 1, tot, IBLOCK
          READ (8,ERR=300,iostat=ierr)
     >         (tmpARRAY(J),J=1,MIN(tot-(i-1)*iblock,IBLOCK-1))
          cnt = 0
          do 110 k = s1, m
             do 110 j = s2,n
                cnt = cnt + 1
                iarray(k,j) = tmparray(cnt)
                if (cnt.eq.iblock-1) goto 115
 110      continue
c
 115      continue
          s1 = k
          s2 = j+1
          if (s2.eq.n+1) then
             s2 = 1
             s1 = k+1
          endif
c
  100   CONTINUE
      ELSEIF (OPT(1:1).EQ.'W') THEN
        s1 = 1
        s2 = 1
        DO 200 I = 1, tot, IBLOCK
          cnt = 0
          do 210 k = s1, m
             do 210 j = s2,n
                cnt = cnt + 1
                tmparray(cnt) = iarray(k,j)
                if (cnt.eq.iblock-1) goto 215
 210      continue
c
 215      continue
          s1 = k
          s2 = j+1
          if (s2.eq.n+1) then
             s2 = 1
             s1 = k+1
          endif
          WRITE (8,ERR=400,iostat=ierr)
     >          (tmpARRAY(J),J=1,MIN(tot-(i-1)*iblock,IBLOCK-1))
c
 200    CONTINUE
      ENDIF
      WRITE (6,9001) OPT(3:len_trim(opt)),REAL(4*N)
      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (6,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return
 9001 FORMAT(1X,'IINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
      END



C
C
C
      SUBROUTINE FITTER (N1,X1,F1,N2,X2,F2,METHOD)
      implicit none
      INTEGER N1,N2
      REAL    X1(N1),F1(N1),X2(N2),F2(N2)
      CHARACTER*(*) METHOD
C
C  *********************************************************************
C  *                                                                   *
C  *  FITTER:   THIS ROUTINE INTERPOLATES BETWEEN A SET OF GIVEN DATA  *
C  *  POINTS AND EXTRAPOLATES OFF THE ENDS OF THE RANGE.  SEVERAL      *
C  *  METHODS ARE USED ACCORDING TO THE VALUE OF "METHOD":-            *
C  *                                                                   *
C  *  "SPLINE"  THIS OPTION USES THE HARWELL CUBIC SPLINE GENERATORS   *
C  *  TO MODEL A SYSTEM OF A FEW DATA POINTS (N1,X1,F1), FROM WHICH WE *
C  *  INTERPOLATE / EXTRAPOLATE TO THE REQUIRED NEW SYSTEM OF MORE     *
C  *  DATAPOINTS (N2,X2,F2).  ALL ARGUMENTS ARE INPUT EXCEPT F2, WHICH *
C  *  IS CALCULATED IN THIS ROUTINE.                                   *
C  *                                                                   *
C  *  "LINEAR"  LINEAR INTERPOLATION                                   *
C  *                                                                   *
C  *  CHRIS FARRELL  (HUNTERSKIL)  SEPT 88                             *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER I1,I2
      REAL FDASH1(1001),WORK(3003),TG01B
C-----------------------------------------------------------------------
      IF     (METHOD.EQ.'SPLINE') THEN
        IF (N1.GT.1001) GOTO 1001
        IF (N1.GT.1) THEN
          CALL TB04A (N1,X1,F1,FDASH1,WORK)
          IF (WORK(1).GT.1.E-6) GOTO 1001
        ENDIF
        DO 100 I2 = 1, N2
          IF     (X2(I2).LE.X1(1)) THEN
            F2(I2) = F1(1)
          ELSEIF (X2(I2).GE.X1(N1)) THEN
            F2(I2) = F1(N1)
          ELSE
            F2(I2) = TG01B (-1,N1,X1,F1,FDASH1,X2(I2))
          ENDIF
  100   CONTINUE
C-----------------------------------------------------------------------
      ELSEIF (METHOD.EQ.'LINEAR') THEN
        DO 210 I2 = 1, N2
          IF     (X2(I2).LE.X1(1)) THEN
            F2(I2) = F1(1)
C           WRITE (6,'(1X,F7.4,'' ='',F7.4)') F2(I2),F1(1)
          ELSEIF (X2(I2).GE.X1(N1)) THEN
            F2(I2) = F1(N1)
C           WRITE (6,'(1X,F7.4,'' ='',F7.4)') F2(I2),F1(N1)
          ELSE
            I1 = 2
  200       IF (X2(I2).GT.X1(I1)) THEN
              I1 = I1 + 1
              GOTO 200
            ENDIF
            F2(I2) = F1(I1) - (X1(I1)-X2(I2)) * (F1(I1)-F1(I1-1)) /
     >                                          (X1(I1)-X1(I1-1))
C           WRITE (6,9001) F2(I2),F1(I1),X1(I1),X2(I2),F1(I1),F1(I1-1),
C    >                                                 X1(I1),X1(I1-1)
          ENDIF
  210   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
      RETURN
C
 1001 CALL PRC ('FITTER:  CUBIC SPLINE ERROR.  CHECK N1 <= 1001')
      STOP
C
 9001 FORMAT(1X,F7.4,' =',F7.4,' - (',F7.4,' -',F7.4,') * (',F7.4,
     >  ' -',F7.4,') / (',F7.4,' -',F7.4,')')
      END






C
C  *********************************************************************
C  *  IZERO:  ZEROES AN INTEGER ARRAY ...  L.D.HORTON    FEB 1994      *
C  *********************************************************************
C
      SUBROUTINE IZERO (IARRAY, N)
      implicit none
      INTEGER I,N
      INTEGER IARRAY(N)
      DO 100 I = 1, N
        IARRAY(I) = 0
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  RZERO:  ZEROES A REAL ARRAY  ...     C.M.FARRELL   NOV 1987      *
C  *********************************************************************
C
      SUBROUTINE RZERO (RARRAY, N)
      implicit none
      INTEGER I,N
      REAL RARRAY(N)
      DO 100 I = 1, N
        RARRAY(I) = 0.0
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  DZERO:  ZEROES A D.P. ARRAY  ...     C.M.FARRELL   FEB 1988      *
C  *********************************************************************
C
      SUBROUTINE DZERO (DARRAY, N)
      implicit none
      INTEGER I,N
      DOUBLE PRECISION DARRAY(N)
      DO 100 I = 1, N
        DARRAY(I) = 0.0D0
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  QZERO:  ZEROES AN EXTENDED PRECISION ARRAY ... D. ELDER MAR 1995 *
C  *********************************************************************
C
      SUBROUTINE QZERO (QARRAY, N)
      implicit none
      INTEGER I,N
      REAL*8 QARRAY(N)
      DO 100 I = 1, N
        QARRAY(I) = 0.0
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  RINIT:  INITIALISES A REAL ARRAY  ...   G.J.RADFORD   JUN 1993   *
C  *********************************************************************
C
      SUBROUTINE RINIT (RARRAY, N, A)
      implicit none
      INTEGER I,N
      REAL RARRAY(N), A
      DO 100 I = 1, N
        RARRAY(I) = A
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  IINIT:  INITIALISES AN INTEGER ARRAY  ...                        *
C  *********************************************************************
C
      SUBROUTINE IINIT (IARRAY, N, A)
      implicit none
      INTEGER I,N
      integer IARRAY(N), A
      DO 100 I = 1, N
        IARRAY(I) = A
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  DINIT:  INITIALISES A D.P. ARRAY  ... G.J.RADFORD   JUN 1993     *
C  *********************************************************************
C
      SUBROUTINE DINIT (DARRAY, N, A)
      implicit none
      INTEGER I,N
      DOUBLE PRECISION DARRAY(N), A
      DO 100 I = 1, N
        DARRAY(I) = A
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *  QINIT:  INITIALISES AN EXTENDED PRECISION ARRAY                  *
C  *********************************************************************
C
      SUBROUTINE QINIT (QARRAY, N, A)
      implicit none
      INTEGER I,N
      REAL*8 QARRAY(N), A
      DO 100 I = 1, N
        QARRAY(I) = A
  100 CONTINUE
      RETURN
      END
C
C  *********************************************************************
C  *                                                                   *
C  *  IPOS   : FINDS NEAREST HIGHER VALUE IN RS ARRAY TO GIVEN R.      *
C  *           RS ARRAY MUST BE IN ASCENDING ORDER                     *
C  *                                                                   *
C  *  CHRIS FAR RELL    FEBRUARY 1989
C  *                                                                   *
C  *********************************************************************
C
      INTEGER FUNCTION IPOS (R, RS, NRS)
      implicit none
      INTEGER NRS,ILOW,IMID
      REAL    R,RS(NRS)
C
c     NRS = 0 is an error condition - however, it appears that LIM
c     sometimes does this when calculating time points in cases where
c     the case being run is not time dependent so NTS=0. In any case, 
c     IPOS should return some value in error cases - so IPOS will be 
c     set to 1 initially. A fix has been added to LIM setting NTS to 1. 
c
      if (nrs.eq.0) then 
         ipos = 1 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' IPOS ERROR:'//
     >            ' NUMBER OF ELEMENTS IS ZERO ',nrs,r
         WRITE (0,'(a,i6,3(1x,g12.5))') ' IPOS ERROR:'//
     >            ' NUMBER OF ELEMENTS IS ZERO',nrs,r
         return
      elseif (RS(1).GT.RS(NRS)) then 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' IPOS ERROR: DESCENDING ORDER',
     >                  nrs,r,rs(1),rs(nrs)
      endif
C
      ILOW = 0
      IPOS = NRS
      IF (NRS.EQ.1) RETURN
100   CONTINUE
      IMID = (IPOS + ILOW) / 2
      IF (R.GT.RS(IMID)) THEN
        ILOW = IMID
      ELSE
        IPOS = IMID
      ENDIF
      IF (IPOS-ILOW.GT.1) GOTO 100
C
      RETURN
      END

C
C  *********************************************************************
C  *                                                                   *
C  *  IPOSQ  : FINDS NEAREST HIGHER VALUE IN RS ARRAY TO GIVEN R.      *
C  *           RS ARRAY MUST BE IN ASCENDING ORDER
C  *           R, RS ARE REAL*8
C  *                                                                   *
C  *********************************************************************
C
      INTEGER FUNCTION IPOSQ (R, RS, NRS)
      implicit none
      INTEGER NRS,ILOW,IMID
      REAL*8    R,RS(NRS)
C
c     NRS = 0 is an error condition - however, it appears that LIM
c     sometimes does this when calculating time points in cases where
c     the case being run is not time dependent so NTS=0. In any case, 
c     IPOSQ should return some value in error cases - so IPOSQ will be 
c     set to 1 initially. A fix has been added to LIM setting NTS to 1. 
c

      if (nrs.eq.0) then 
         iposq = 1 
         WRITE (6,'(a,i6,3(1x,g12.5))') 'IPOSQ ERROR:'//
     >            ' NUMBER OF ELEMENTS IS ZERO ',nrs,r
         WRITE (0,'(a,i6,3(1x,g12.5))') 'IPOSQ ERROR:'//
     >            ' NUMBER OF ELEMENTS IS ZERO',nrs,r
         return
      elseif (RS(1).GT.RS(NRS)) then 
         WRITE (6,'(a,i6,3(1x,g12.5))') 'IPOSQ ERROR: DESCENDING ORDER',
     >                  nrs,r,rs(1),rs(nrs)
      endif
C
      ILOW = 0
      IPOSQ = NRS
      IF (NRS.EQ.1) RETURN
100   CONTINUE
      IMID = (IPOSQ + ILOW) / 2
      IF (R.GT.RS(IMID)) THEN
        ILOW = IMID
      ELSE
        IPOSQ = IMID
      ENDIF
      IF (IPOSQ-ILOW.GT.1) GOTO 100
C
      RETURN
      END
      
C
C  *********************************************************************
C  *                                                                   *
C  *  JPOS   : FINDS NEAREST HIGHER VALUE IN RS ARRAY TO GIVEN R.      *
C  *           RS ARRAY MUST BE IN DESCENDING ORDER                    *
C  *                                                                   *
C  *  CHRIS FARRELL    FEBRUARY 1989                                   *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER FUNCTION JPOS (R, RS, NRS)
      implicit none
      INTEGER NRS,ILOW,IMID
      REAL    R,RS(NRS)
C
      if (nrs.eq.0) then 
         jpos = 1 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' JPOS ERROR:'//
     >            ' NUMBER OF ELEMENTS IS ZERO',
     >                  nrs,r,rs(1),rs(nrs)
         return
      elseif (RS(1).GT.RS(NRS)) then 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' JPOS ERROR: ASCENDING ORDER',
     >                  nrs,r,rs(1),rs(nrs)
      endif
C
      ILOW = 1
      JPOS = NRS + 1
100   CONTINUE
      IMID = (JPOS + ILOW) / 2
      IF (R.LE.RS(IMID)) THEN
        ILOW = IMID
      ELSE
        JPOS = IMID
      ENDIF
      IF (JPOS-ILOW.GT.1) GOTO 100
      JPOS = JPOS - 1
C
      RETURN
      END
C
C
C
      INTEGER FUNCTION LENSTR (ASTR)
      implicit none
      CHARACTER*(*) ASTR
C
C  *********************************************************************
C  *                                                                   *
C  *  LENSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS.                                     *
C  *                                                                   *
C  *********************************************************************
C
      integer i
c
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            LENSTR = I
            RETURN
         ENDIF
   10 CONTINUE
      LENSTR = 1
      RETURN
      END
C
C
C
      INTEGER FUNCTION EXTSTR (ASTR,start)
      implicit none
      CHARACTER*(*) ASTR
      integer start,i
c
      extstr = 0
C
C  *********************************************************************
C  *                                                                   *
C  *  EXTSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS AND THE EFFECTIVE STARTING POINT BY  *
c  *          IGNORING LEADING BLANKS.                                 *
C  *                                                                   *
C  *********************************************************************
C
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            extSTR = I
            goto 20
         ENDIF
   10 CONTINUE

 20   if (extstr.gt.1) then
         do i = 1, extstr
            if (astr(I:I).ne.' ') then
               start = i
               return
            endif
         end do
      else
         extSTR = 1
         start  = 1
      endif
c
      RETURN
      END
c
c ======================================================================
c
c
c
c
      SUBROUTINE ER(routine,message,*)
      use error_handling
      use mod_params
      !use mod_slcom
      IMPLICIT none

      CHARACTER routine*(*),message*(*)

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      call errmsg(routine,message)
      call errmsg(routine,message,slout)

c      WRITE(0    ,'(4A)') ' ERROR ',routine,': ',message
c      WRITE(SLOUT,'(4A)') ' ERROR ',routine,': ',message

      RETURN 1
      END


      subroutine intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect,flag)
      implicit none   
      real*8 ra,za,rb,zb,r1,z1,r2,z2,rint,zint
      integer sect,flag
c
c     Calculates intersection of two lines and sets flag if the
c     intersection is due to lines being parallel. 
c
c     FLAG:
c
c     0 - normal intersection
c 
c     1 - lines colinear horizontal
c
c     2 - lines colinear vertical
c
c     3 - lines colinear parallel
c
c
c     SECT:
c
c     Instead of just true and false for the intersection - this 
c     code has been generalized to 4 return values.
c
c     0 - intersection does not occur between end-points or no intersection
c
c     1 - intersection between the lines occurs between both sets of 
c         specified end points 
c
c     2 - intersection between the lines occurs between the first
c         set of end-points but not the second
c
c     3 - intersection between the lines occurs between the second
c         set of end-points but not the first
c
      real*8 ma,m1,ba,b1
      real*8 rdsta,rdstb
      logical verta,vert1 
      real*8 eps
      parameter (eps=1.0d-8)
c
      integer warnings
      data warnings/0/
c
      logical :: debug = .false.

c
      verta = .false.
      vert1 = .false.

      rint = 0.0
      zint = 0.0

      flag = 0
      sect = 0
c
c     Define slopes of lines
c 
      if (ra.eq.rb) then 
         verta = .true.
         ma = 0.0  
         ba = ra
      else
         ma = (zb-za)/(rb-ra)
         ba = za - ma * ra
      endif
c
      if (r1.eq.r2) then 
         vert1 = .true.
         m1 = 0.0  
         b1 = r1
      else
         m1 = (z2-z1)/(r2-r1)
         b1 = z1 - m1 * r1
      endif
c
c     Debug:
c
      if (debug) then 
         write(6,'(a,8g20.12)') 'IS2DP:',ra,za,rb,zb,r1,z1,r2,z2
         write(6,'(a,8g20.12)') '      ',ma,ba,m1,b1
         write(6,'(a,6l6)')     '      ',verta,vert1,ba.eq.b1,ma.eq.m1,
     >                     abs(ba-b1).lt.eps,abs(ma-m1).lt.eps
      endif
c
c
c     Find intersection 
c
      if (verta.and.vert1) then 
c
c        Line segments may overlap
c        Error - set sect to false and issue error message
c
c        Do nothing for parallel case
c
         if (abs(ba-b1).lt.eps) then 
c
            warnings = warnings + 1
c
            write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                     ' COLINEAR-VERTICAL',warnings

            if (debug) then 
               write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                     ' COLINEAR-VERTICAL',warnings
            endif 

            write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                  r2,z2,ma,m1,ba,b1
            
c           set flag to 1 for colinear vertical
            flag = 1
c
c           Check for intersection region and return the mid-point of 
c           the overlap as the intersection point. 
c
            call find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,
     >                                      ma,ba,rint,zint,0)

c
         endif
c
      elseif (verta) then 
c
c        Line A is vertical - line 1 is not -> ra = rb = rint
c
         rint = ra
         zint = m1 * rint + b1         
c
      elseif (vert1) then  
c
c        Line 1 is vertical - line A is not -> r1 = r2 = rint
c
         rint = r1
         zint = ma * rint + ba        
c
      else
c
c        Neiither line vertical - check for parallel or numerically parallel lines
c  
         if (abs(ma-m1).lt.eps) then 
c         if (ma.eq.m1) then 
c
c           Check for the same line 
c
            if (abs(ba-b1).lt.eps) then 
c            if (ba.eq.b1) then 
c
               warnings = warnings+1
c
               if (ma.eq.0.0) then 
                  write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-HORIZONTAL',warnings
                 if (debug) then 
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-HORIZONTAL',warnings
                 endif

                 write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                                r2,z2,ma,m1,ba,b1
c                set flag to 2 for colinear horizontal
                 flag = 2
               else
                  write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-PARALLEL',warnings
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-PARALLEL',warnings
                 write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                                r2,z2,ma,m1,ba,b1
c                set flag to 3 for colinear parallel
                 flag = 3
               endif
c
c              Check for intersection region and return the mid-point of 
c              the overlap as the intersection point. 
c
               call find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,
     >                                      ma,ba,rint,zint,1)
c
            endif  
c
c        Calculate intersection 
c
         else

            rint = (b1-ba)/(ma-m1)
            zint = ma*rint + ba

         endif 
c
      endif
c
c     Now determine if rint,zint lies inside both of the rectangles
c     defined by the two pairs of end points. 
c 
      if (
     >  (dabs(dabs(rint-ra)+dabs(rint-rb)-dabs(ra-rb)).le.eps).and.
     >  (dabs(dabs(zint-za)+dabs(zint-zb)-dabs(za-zb)).le.eps).and. 
     >  (dabs(dabs(rint-r1)+dabs(rint-r2)-dabs(r1-r2)).le.eps).and.
     >  (dabs(dabs(zint-z1)+dabs(zint-z2)-dabs(z1-z2)).le.eps)) then 
           sect = 1
      elseif (
     >  (dabs(dabs(rint-ra)+dabs(rint-rb)-dabs(ra-rb)).le.eps).and.
     >  (dabs(dabs(zint-za)+dabs(zint-zb)-dabs(za-zb)).le.eps)) then 
           sect = 2
      elseif (
     >  (dabs(dabs(rint-r1)+dabs(rint-r2)-dabs(r1-r2)).le.eps).and.
     >  (dabs(dabs(zint-z1)+dabs(zint-z2)-dabs(z1-z2)).le.eps)) then 
           sect = 3 
      endif 
c
c     Check for the case of a numerical interection where the point 
c     is exceptionally close to the wall. 
c
      if (sect.eq.3) then 
         rdsta = sqrt((ra-rint)**2+(za-zint)**2)
         rdstb = sqrt((rb-rint)**2+(zb-zint)**2)
c 
c        If the point is close enough to wall then count it as the intersection
c
         if (rdsta.lt.eps.or.rdstb.lt.eps) then 
            sect=1 
         endif
c
c        Print notification if this code triggers
c
         if (debug) then
            write(6,'(a,2g20.12,2l6)') '      ',rdsta,rdstb,
     >                                 rdsta.lt.eps,rdstb.lt.eps
         endif

c
      endif 


c
c      write(6,'(a,2g18.10,i8)') '      ',rint,zint,sect
c
c
c      write (6,'(a,6l4,1p,10(g14.7))') 'DEBUG I2A:',verta,vert1,
c     >  ((abs(rint-ra)+abs(rint-rb)-abs(ra-rb)).lt.eps),
c     >  ((abs(zint-za)+abs(zint-zb)-abs(za-zb)).lt.eps),
c     >  ((abs(rint-r1)+abs(rint-r2)-abs(r1-r2)).lt.eps),
c     >  ((abs(zint-z1)+abs(zint-z2)-abs(z1-z2)).lt.eps),
c     >  (abs(rint-ra)+abs(rint-rb)-abs(ra-rb)),
c     >  (abs(zint-za)+abs(zint-zb)-abs(za-zb)),
c     >  (abs(rint-r1)+abs(rint-r2)-abs(r1-r2)),
c     >  (abs(zint-z1)+abs(zint-z2)-abs(z1-z2))
c
c      write (6,'(a,1p,10(g14.7))') 'DEBUG I2B:',ra,rint,rb,za,zint,zb
c      write (6,'(a,1p,10(g14.7))') 'DEBUG I2C:',r1,rint,r2,z1,zint,z2
c      write (6,'(a,1p,10(g14.7))') 'DEBUG I2C:',ma,ba,m1,b1
c

c
      return 
      end
c
c
c
      subroutine find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,
     >                                      ma,ba,rint,zint,
     >                                      vert)
      implicit none
      real*8 ra,za,rb,zb,r1,z1,r2,z2,ma,ba,rint,zint
      integer vert,sect

c
c     The lines are known to be co-linear.
c     vert = 0 - vertical lines
c     vert = 1 - not vertical lines
c
c     The code returns the center of the overlap region of the 2 line
c     segments if it exists
c
      real*8 zstart,zend,rstart,rend
c
c     vertical lines
c
      if (vert.eq.0) then 
c
c        Organize by Z-coordinate
c
         zstart=max(min(za,zb),min(z1,z2)) 
         zend  =min(max(za,zb),max(z1,z2))
c
c        No overlap
c
         if (zend.lt.zstart) then 
c            sect = 0
            rint = 0.0
            zint = 0.0
c
c        Get center of overlap region
c
         else
c            sect = 1
            zint = (zstart+zend)/2.0
            rint = ra
         endif
c   
c     Base analysis on R - and use line equation for Z. 
c
      else
c
         rstart=max(min(ra,rb),min(r1,r2)) 
         rend  =min(max(ra,rb),max(r1,r2))
c
c
c        No overlap
c
         if (rend.lt.rstart) then 
c            sect = 0
            rint = 0.0
            zint = 0.0
c
c        Get center of overlap region
c
         else
c
c            sect = 1
c
            rint = (rstart+rend)/2.0
            zint =  rint * ma + ba
c
         endif
c
      endif 
c
c      write(6,'(a,10g18.10,i6)') 'WARNING: INTSECT2:'//
c     >                    ' INTERSECTION REGION IS CO-LINEAR:',
c     >                      ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect
c
      return
      end
c
c
c
      REAL FUNCTION ATAN2C (ARGZ,ARGR)
      use mod_params
      implicit none
      REAL ARGZ,ARGR
C     INCLUDE "PARAMS"
c      include 'params'
C
C     THIS ACTS AS AN ERROR-CHECKING FRONT-END TO THE ATAN2
C     IMPLICIT FUNCTION. IT RETURNS APPROPRIATE ANGLE VALUES FOR
C     EITHER OF THE ARGUMENTS EQUAL TO ZERO AND RETURNS A
C     ZERO VALUE IF BOTH ARGUMENTS ARE EQUAL TO ZERO. SINCE
C     THE TANGENT IS UNDEFINED IN THIS CASE.
C
C     D. ELDER  SEPTEMBER 1992
C
      IF (ARGZ.EQ.0.0) THEN
         IF (ARGR.GT.0.0) THEN
            ATAN2C = 0.0
         ELSEIF (ARGR.LT.0.0) THEN
            ATAN2C = PI
         ELSE
            ATAN2C = 0.0
         ENDIF
      ELSEIF (ARGR.EQ.0.0) THEN
         IF (ARGZ.GT.0.0) THEN
            ATAN2C = PI /2.0
         ELSEIF (ARGZ.LT.0.0) THEN
            ATAN2C = - PI /2.0
         ELSE
            ATAN2C = 0.0
         ENDIF
      ELSE
         ATAN2C = ATAN2(ARGZ,ARGR)
      ENDIF
      RETURN
      END
c
c
c
      REAL*8 FUNCTION DATAN2C (ARGZ,ARGR)
      use mod_params
      implicit none
      REAL*8 ARGZ,ARGR
C     INCLUDE "PARAMS"
c      include 'params'
C
C     THIS ACTS AS AN ERROR-CHECKING FRONT-END TO THE ATAN2
C     IMPLICIT FUNCTION. IT RETURNS APPROPRIATE ANGLE VALUES FOR
C     EITHER OF THE ARGUMENTS EQUAL TO ZERO AND RETURNS A
C     ZERO VALUE IF BOTH ARGUMENTS ARE EQUAL TO ZERO. SINCE
C     THE TANGENT IS UNDEFINED IN THIS CASE.
C
C     D. ELDER  SEPTEMBER 1992
C
      IF (ARGZ.EQ.0.0) THEN
         IF (ARGR.GT.0.0) THEN
            DATAN2C = 0.0
         ELSEIF (ARGR.LT.0.0) THEN
            DATAN2C = PI
         ELSE
            DATAN2C = 0.0
         ENDIF
      ELSEIF (ARGR.EQ.0.0) THEN
         IF (ARGZ.GT.0.0) THEN
            DATAN2C = PI /2.0
         ELSEIF (ARGZ.LT.0.0) THEN
            DATAN2C = - PI /2.0
         ELSE
            DATAN2C = 0.0
         ENDIF
      ELSE
         DATAN2C = ATAN2(ARGZ,ARGR)
      ENDIF
      RETURN
      END
c
c
c
      subroutine find_free_unit_number(unit)
      implicit none
      integer unit
c
c     FIND_FREE_UNIT_NUMBER:
c
c     This routine scans through unit numbers looking for one that
c     is not currently in use. This number is returned. This code
c     is based on the assumption that any unit numbers returned will
c     be used before this routine is called again asking for another 
c     number - otherwise it will likely return the previous value.
c
      integer test_unit
      logical unit_open

      test_unit = 10
      unit_open = .true.

      ! Check for unit number assignment.  
      Do While (Unit_open)
         test_unit=test_unit + 1
         Inquire (Unit = test_unit, Opened = Unit_open)
      End Do

      unit = test_unit

      return
      end
