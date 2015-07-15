C     -*-Fortran-*-
C
C
      SUBROUTINE RDRAR(RS,NRS, MAXNRS, RMIN, RMAX, ASCEND, NAME, IERR)
      use error_handling
      implicit none
      INTEGER   NRS, MAXNRS
      REAL      RS(MAXNRS), RMIN, RMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
      LOGICAL   ASCEND
C
C  *********************************************************************
C  *                                                                   *
C  *  RDRAR:   ROUTINE READS A REAL VARIABLE LENGTH 1D ARRAY           *
C  *    AND CHECKS THAT IT IS SORTED IN ASCENDING ORDER,               *
C  *    THAT ALL THE VALUES ARE WITHIN A RANGE AND THAT NO             *
C  *    TWO VALUES ARE EQUAL.                                          *
C  *    ANY OF THESE RULES BEING BROKEN RESULTS IN THE                 *
C  *    OFFENDING VALUE BEING REMOVED FROM THE ARRAY, AN               *
C  *    ERROR MESSAGE BEING OUTPUT AND AN ERROR FLAG BEING             *
C  *    SET.                                                           *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  RS     : ARRAY VALUES RETURNED IN THIS                           *
C  *  NRS    : NUMBER OF VALUES READ RETURNED IN THIS                  *
C  *  MAXNRS : MAXIMUM NUMBER OF VALUES TO BE READ                     *
C  *  RMIN   : MINIMUM ALLOWED VALUE (EXCLUSIVE)                       *
C  *  RMAX   : MAXIMUM ALLOWED VALUE (EXCLUSIVE)                       *
C  *  ASCEND : INDICATES WHETHER ASCENDING ORDER CHECK IS TO BE MADE   *
C  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL   JAN 1988                                       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N
      REAL      RLAST,R
C
C---- READ IN AND TEST SIZE OF ARRAY
C
      NRS = 0
      CALL RDC (COMENT, NAME, IERR)
      CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
      IF (N.EQ.0) RETURN
C
C---- READ IN AND TEST EACH ARRAY VALUE
C
      RLAST = RMIN
      IR = 1
   50 CONTINUE
      MESAGE = 'END OF FILE ON UNIT 5'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRAR'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      WRITE (MESAGE,'(A,I5,A)') 'EXPECTING',N,' REALS, ONE PER LINE'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) R
C
      IF (ASCEND) THEN
        WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        IF (R.LT.RLAST) GOTO 9999
      ENDIF
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      IF (R.GT.RMAX) GOTO 9999
      RS(IR) = R
      IR = IR + 1
      IF (IR.LE.N) GOTO 50
C
C---- SET UP NUMBER OF VALID VALUES READ
C
      NRS = IR - 1
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDRAR-READ ERROR',name//mesage)
         call errmsg('RDRAR-LAST LINE ',buffer)

      else

         call dbgmsg('RDRAR-READ ERROR',name//mesage)
         call dbgmsg('RDRAR-LAST LINE ',buffer)

      endif


c      WRITE (7,'(1X,2A,3(/1X,A))')
c     > 'RDRAR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
C
C
C
c slmod begin 
c
c This routine is based on RDRARN but has been modified to read an array
c of integers, rather than reals. -March, 2001
c
      SUBROUTINE RDIARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,
     >                                                   NFS,NAME,IERR)
      use error_handling
      implicit none
      INTEGER   NRS, MAXNRS, NFS
      INTEGER   RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
c      REAL      RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
      LOGICAL   ASCEND
C
C  *********************************************************************
C  *                                                                   *
C  *  RDIARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
C  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
C  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
C  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
C  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
C  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
C  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
C  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
C  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
C  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
C  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
C  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
C  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N,I
      INTEGER   RLAST,R,F(10)
c      REAL      RLAST,R,F(10)
c
C
C---- READ IN AND TEST SIZE OF ARRAY
C
      NRS = 0
      CALL RDC (COMENT, NAME, IERR)
      CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
      IF (N.EQ.0) RETURN
C
C---- READ IN AND TEST EACH SET OF ARRAY VALUES
C
      RLAST = RMIN
      IR = 1
   50 CONTINUE
      MESAGE = 'END OF FILE ON UNIT 5'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDIARN'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
c      IF (BUFFER(2:2).EQ.'*') THEN
c        CALL ReadUnstructuredInput(BUFFER)
c        GOTO 100
c      ENDIF
c slmod end
C
      WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,
     >  ' REALS,',NFS+1,' PER LINE'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
      IF (ASCEND) THEN
        WRITE (MESAGE,'(I11,A,I11)') R,' LESS THAN PREV/MIN',RLAST
        IF (R.LT.RLAST) GOTO 9999
      ENDIF
      WRITE (MESAGE,'(I11,A,I11)') R,' MORE THAN MAXIMUM',RMAX
      IF (R.GT.RMAX) GOTO 9999
      RS(IR,1) = R
C
      DO 120 I = 1, NFS
        WRITE (MESAGE,'(I11,A,I11)') F(I),' LESS THAN MINIMUM',FMIN
        IF (F(I).LT.FMIN) GOTO 9999
        WRITE (MESAGE,'(I11,A,I11)') F(I),' MORE THAN MAXIMUM',FMAX
        IF (F(I).GT.FMAX) GOTO 9999
        RS(IR,1+I) = F(I)
  120 CONTINUE
C
C---- SET UP NUMBER OF VALID VALUES READ
C
      IR = IR + 1
      IF (IR.LE.N) GOTO 50
      NRS = IR - 1
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDIARN-READ ERROR',name//mesage)
         call errmsg('RDIARN-LAST LINE ',buffer)

      else
         call dbgmsg('RDIARN-READ ERROR',name//mesage)
         call dbgmsg('RDIARN-LAST LINE ',buffer)
      endif


c      WRITE (7,'(1X,2A,3(/1X,A))')
c     > 'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
c      RETURN
      END
c slmod end
      SUBROUTINE RDRARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,
     >                                                   NFS,NAME,IERR)
      implicit none
      INTEGER   NRS, MAXNRS, NFS
      REAL      RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
      LOGICAL   ASCEND
C
C  *********************************************************************
C  *                                                                   *
C  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
C  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
C  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
C  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
C  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
C  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
C  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
C  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
C  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
C  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
C  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
C  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
C  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N,I
c slmod begin
      REAL      RLAST,R,F(20)
c
c      REAL      RLAST,R,F(10)
c slmod end
C
C---- READ IN AND TEST SIZE OF ARRAY
C
      NRS = 0
c slmod begin
c...  Required for Intel compiler, but likely I can find a flag that initializes all variables:
      IERR = 0
c slmod end
      CALL RDC (COMENT, NAME, IERR)
      CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)

      IF (N.EQ.0) RETURN
C
C---- READ IN AND TEST EACH SET OF ARRAY VALUES
C
      RLAST = RMIN
      IR = 1
   50 CONTINUE
      MESAGE = 'END OF FILE ON UNIT 5'
c slmod begin
c... This 72 character limit has always been annoying, and I can't see
c    any reason not to increase it since BUFFER*512 is declared
c    in READER:
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
c
c  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9999,END=9999) BUFFER
c slmod end
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,
     >  ' REALS,',NFS+1,' PER LINE'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
      IF (ASCEND) THEN
        WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        IF (R.LT.RLAST) GOTO 9999
      ENDIF
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      IF (R.GT.RMAX) GOTO 9999
      RS(IR,1) = R
C
      DO 120 I = 1, NFS
        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
        IF (F(I).LT.FMIN) GOTO 9999
        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
        IF (F(I).GT.FMAX) GOTO 9999
        RS(IR,1+I) = F(I)
  120 CONTINUE
C
C---- SET UP NUMBER OF VALID VALUES READ
C
      IR = IR + 1
      IF (IR.LE.N) GOTO 50
      NRS = IR - 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     > 'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',
     > BUFFER(1:72)
      RETURN
      END
C
C
C
      SUBROUTINE RDR(R, TSTMIN, RMIN, TSTMAX, RMAX, NAME, IERR)
      use error_handling
      implicit none
      REAL      R, RMIN, RMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDR:  ROUTINE READS IN A REAL.
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  R      : VALUE                                                   *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  RMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  RMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R = 0.0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & REAL NUMBER'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,R
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS LESS THAN MINIMUM ', RMIN
      IF (TSTMIN.AND.(R.LT.RMIN)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS MORE THAN MAXIMUM ', RMAX
      IF (TSTMAX.AND.(R.GT.RMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDR-READ ERROR',name//mesage)
         call errmsg('RDR-LAST LINE ',buffer)

      else
         call dbgmsg('RDR-READ ERROR',name//mesage)
         call dbgmsg('RDR-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
c
c
c
      SUBROUTINE RDR2(R1, R2,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
      use error_handling
      IMPLICIT  NONE
      REAL      R1, R2, RMIN, RMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDR2:  ROUTINE READS IN 2 REALS.
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  R1,R2  : VALUES                                                  *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  RMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  RMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR2'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & 2 REAL NUMBERS'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2
C
c     Test minimums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999
c
c     Test Maximums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDR2-READ ERROR',name//mesage)
         call errmsg('RDR2-LAST LINE ',buffer)

      else
         call dbgmsg('RDR2-READ ERROR',name//mesage)
         call dbgmsg('RDR2-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDR2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDR3(R1, R2, R3,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
      use error_handling
      IMPLICIT  NONE
      REAL      R1, R2, R3, RMIN, RMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDR3:  ROUTINE READS IN 3 REALS.
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  R1,R2,R3 : VALUES                                                *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  RMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  RMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
      R3 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR3'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & 3 REAL NUMBERS'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2,R3
C
c     Test minimums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R3,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R3.LT.RMIN)) GOTO 9999
C
c     Test Maximums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R3,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R3.GT.RMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDR3-READ ERROR',name//mesage)
         call errmsg('RDR3-LAST LINE ',buffer)

      else
         call dbgmsg('RDR3-READ ERROR',name//mesage)
         call dbgmsg('RDR3-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDR3: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDQAR(RS,NRS, MAXNRS, RMIN, RMAX, ASCEND, NAME, IERR)
      use error_handling
      implicit none
      INTEGER   NRS, MAXNRS
      real*8   RS(MAXNRS), RMIN, RMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
      LOGICAL   ASCEND
C
C  *********************************************************************
C  *                                                                   *
C  *  RDQAR:   ROUTINE READS A REAL VARIABLE LENGTH 1D ARRAY           *
C  *    AND CHECKS THAT IT IS SORTED IN ASCENDING ORDER,               *
C  *    THAT ALL THE VALUES ARE WITHIN A RANGE AND THAT NO             *
C  *    TWO VALUES ARE EQUAL.                                          *
C  *    ANY OF THESE RULES BEING BROKEN RESULTS IN THE                 *
C  *    OFFENDING VALUE BEING REMOVED FROM THE ARRAY, AN               *
C  *    ERROR MESSAGE BEING OUTPUT AND AN ERROR FLAG BEING             *
C  *    SET.                                                           *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  RS     : ARRAY VALUES RETURNED IN THIS                           *
C  *  NRS    : NUMBER OF VALUES READ RETURNED IN THIS                  *
C  *  MAXNRS : MAXIMUM NUMBER OF VALUES TO BE READ                     *
C  *  RMIN   : MINIMUM ALLOWED VALUE (EXCLUSIVE)                       *
C  *  RMAX   : MAXIMUM ALLOWED VALUE (EXCLUSIVE)                       *
C  *  ASCEND : INDICATES WHETHER ASCENDING ORDER CHECK IS TO BE MADE   *
C  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL   JAN 1988                                       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      INCLUDE   'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N
      real*8   RLAST,R
C
C---- READ IN AND TEST SIZE OF ARRAY
C
      NRS = 0
      CALL RDC (COMENT, NAME, IERR)
      CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
      IF (N.EQ.0) RETURN
C
C---- READ IN AND TEST EACH ARRAY VALUE
C
      RLAST = RMIN
      IR = 1
   50 CONTINUE
      MESAGE = 'END OF FILE ON UNIT 5'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      WRITE (MESAGE,'(A,I5,A)') 'EXPECTING',N,' REALS, ONE PER LINE'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) R
C
      IF (ASCEND) THEN
        WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        IF (R.LT.RLAST) GOTO 9999
      ENDIF
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      IF (R.GT.RMAX) GOTO 9999
      RS(IR) = R
      IR = IR + 1
      IF (IR.LE.N) GOTO 50
C
C---- SET UP NUMBER OF VALID VALUES READ
C
      NRS = IR - 1
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDQAR-READ ERROR',name//mesage)
         call errmsg('RDQAR-LAST LINE ',buffer)

      else
         call dbgmsg('RDQAR-READ ERROR',name//mesage)
         call dbgmsg('RDQAR-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     > 'RDQAR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
c
      RETURN
      END
C
C
C
      SUBROUTINE RDQARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,
     >                                                   NFS,NAME,IERR)
      use error_handling
      implicit none
      INTEGER   NRS, MAXNRS, NFS
      real*8  RS(MAXNRS,1+NFS)
      real    RMIN, RMAX,FMIN,FMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
      LOGICAL   ASCEND
C
C  *********************************************************************
C  *                                                                   *
C  *  RDQARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)   *
C  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
C  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
C  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
C  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
C  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
C  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
C  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
C  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
C  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
C  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
C  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
C  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
C  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
C  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      INCLUDE   'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N,I
      real*8   RLAST,R,F(10)
C
C---- READ IN AND TEST SIZE OF ARRAY
C
      NRS = 0
      CALL RDC (COMENT, NAME, IERR)
      CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
      IF (N.EQ.0) RETURN
C
C---- READ IN AND TEST EACH SET OF ARRAY VALUES
C
      RLAST = RMIN
      IR = 1
   50 CONTINUE
      MESAGE = 'END OF FILE ON UNIT 5'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQARN'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,
     >  ' REALS,',NFS+1,' PER LINE'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
      IF (ASCEND) THEN
        WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        IF (R.LT.RLAST) GOTO 9999
      ENDIF
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      IF (R.GT.RMAX) GOTO 9999
      RS(IR,1) = R
C
      DO 120 I = 1, NFS
        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
        IF (F(I).LT.FMIN) GOTO 9999
        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
        IF (F(I).GT.FMAX) GOTO 9999
        RS(IR,1+I) = F(I)
  120 CONTINUE
C
C---- SET UP NUMBER OF VALID VALUES READ
C
      IR = IR + 1
      IF (IR.LE.N) GOTO 50
      NRS = IR - 1
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDQARN-READ ERROR',name//mesage)
         call errmsg('RDQARN-LAST LINE ',buffer)

      else
         call dbgmsg('RDQARN-READ ERROR',name//mesage)
         call dbgmsg('RDQARN-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     > 'RDQARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDQ(R, TSTMIN, RMIN, TSTMAX, RMAX, NAME, IERR)
      use error_handling
      implicit none
      real*8  R, RMIN, RMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDQ:  ROUTINE READS IN A REAL*8                                 *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  R      : VALUE                                                   *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  RMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  RMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      INCLUDE   'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R = 0.0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & REAL NUMBER'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,R
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS LESS THAN MINIMUM ', RMIN
      IF (TSTMIN.AND.(R.LT.RMIN)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS MORE THAN MAXIMUM ', RMAX
      IF (TSTMAX.AND.(R.GT.RMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDQ-READ ERROR',name//mesage)
         call errmsg('RDQ-LAST LINE ',buffer)

      else
         call dbgmsg('RDQ-READ ERROR',name//mesage)
         call dbgmsg('RDQ-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDQ: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
c
c
c
      SUBROUTINE RDQ2(R1, R2,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
      use error_handling
      IMPLICIT  NONE
      REAL*8      R1, R2, RMIN, RMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDQ2:  ROUTINE READS IN 2 REAL*8.                               *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  R1,R2  : VALUES                                                  *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  RMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  RMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ2'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & 2 REAL NUMBERS'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2
C
c     Test minimums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
      IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999
c
c     Test Maximums
c
      WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999
C
      WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
      IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDQ2-READ ERROR',name//mesage)
         call errmsg('RDQ2-LAST LINE ',buffer)

      else
         call dbgmsg('RDQ2-READ ERROR',name//mesage)
         call dbgmsg('RDQ2-LAST LINE ',buffer)
      endif

c
c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDQ2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
c
      RETURN
      END

C
C
C
      SUBROUTINE RDI(I, TSTMIN, IMIN, TSTMAX, IMAX, NAME, IERR)
      use error_handling
      implicit none
      INTEGER   I, IMIN, IMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDI:  THIS ROUTINE READS IN AN INTEGER.                          *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  I      : VALUE                                                   *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  IMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  IMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL   JAN 1988                                       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
c slmod begin
      CHARACTER COMENT*128,MESAGE*72
c
c      CHARACTER COMENT*72,MESAGE*72
c slmod end
C
      I = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & IN TEGER'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,I
C
      WRITE (MESAGE,'(I11,A,I11)') I,' IS LESS THAN MINIMUM ', IMIN
      IF (TSTMIN.AND.(I.LT.IMIN)) GOTO 9999
C
      WRITE (MESAGE,'(I11,A,I11)') I,' IS MORE THAN MAXIMUM ', IMAX
      IF (TSTMAX.AND.(I.GT.IMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDI-READ ERROR',name//mesage)
         call errmsg('RDI-LAST LINE ',buffer)

      else
         call dbgmsg('RDI-READ ERROR',name//mesage)
         call dbgmsg('RDI-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDI: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
c
c
c
      SUBROUTINE RDI2(I1, I2, TSTMIN, IMIN, TSTMAX, IMAX, NAME, IERR)
      use error_handling
      implicit none
      INTEGER   I1, I2, IMIN, IMAX
      LOGICAL   TSTMIN, TSTMAX
      CHARACTER NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDI:  THIS ROUTINE READS IN TWO INTEGERS.                        *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  I      : VALUE                                                   *
C  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
C  *  IMIN   : MINIMUM VALUE                                           *
C  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
C  *  IMAX   : MAXIMUM VALUE                                           *
C  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *    CHRIS FARRELL   JAN 1988                                       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      I1 = 0
      I2 = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI2'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & 2 INTEGERS'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,I1,I2
C
c     Test Min
c
      WRITE (MESAGE,'(I11,A,I11)') I1,' IS LESS THAN MINIMUM ', IMIN
      IF (TSTMIN.AND.(I1.LT.IMIN)) GOTO 9999
C
      WRITE (MESAGE,'(I11,A,I11)') I2,' IS LESS THAN MINIMUM ', IMIN
      IF (TSTMIN.AND.(I2.LT.IMIN)) GOTO 9999
c
c     Test Max
C
      WRITE (MESAGE,'(I11,A,I11)') I1,' IS MORE THAN MAXIMUM ', IMAX
      IF (TSTMAX.AND.(I1.GT.IMAX)) GOTO 9999
C
      WRITE (MESAGE,'(I11,A,I11)') I2,' IS MORE THAN MAXIMUM ', IMAX
      IF (TSTMAX.AND.(I2.GT.IMAX)) GOTO 9999
C
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDI2-READ ERROR',name//mesage)
         call errmsg('RDI2-LAST LINE ',buffer)

      else
         call dbgmsg('RDI2-READ ERROR',name//mesage)
         call dbgmsg('RDI2-LAST LINE ',buffer)
      endif
c
c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDI2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
c
      RETURN
      END
C
C
C
      SUBROUTINE RDC(STRING, NAME, IERR)
      use error_handling
      implicit none
      CHARACTER STRING*(*), NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDC:  THIS ROUTINE READS IN A CHARACTER STRING.                  *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  STRING : CHARACTER STRING                                        *
C  *  NAME   : NAME OF CHARACTER STRING                                *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      STRING = ' '
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

c     Feb/2008 - jde - changed all buffer reads to A256 from A72
c                    - this one left at 512 for reading title - other
c                      entries could be increased to 512 if needed
c                      buffer is 512 - * specifier can not be used
c                      since the input text contains quoted character
c                      strings
  100 IF (IBUF.EQ.0) READ (5,'(A512)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDC'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'EXPECTING COMMENT & CHARACTER STRING'
      IBUF = 0
      READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDC-READ ERROR',name//mesage)
         call errmsg('RDC-LAST LINE ',buffer)

      else
         call dbgmsg('RDC-READ ERROR',name//mesage)
         call dbgmsg('RDC-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
C
C
C
      SUBROUTINE RDBUFFER(STRING, NAME, IERR)
      use error_handling
      implicit none
      CHARACTER STRING*(*), NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDBUFFER:  THIS ROUTINE READS IN AND RETURNS THE ENTIRE INPUT    *
C  *             BUFFER SO IT CAN BE PROCESSED BY THE CALLING ROUTINE  *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  STRING : CHARACTER STRING                                        *
C  *  NAME   : NAME OF CHARACTER STRING                                *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      STRING = ' '
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

c     Feb/2008 - jde - changed all buffer reads to A256 from A72
c                    - this one left at 512 for reading title - other
c                      entries could be increased to 512 if needed
c                      buffer is 512 - * specifier can not be used
c                      since the input text contains quoted character
c                      strings
  100 IF (IBUF.EQ.0) READ (5,'(A512)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDBUFFER'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      MESAGE = 'RETURNS BUFFER'
      IBUF = 0
      string = trim(buffer)
c      READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDBUFFER-READ ERROR',name//mesage)
         call errmsg('RDBUFFER-LAST LINE ',buffer)

      else
         call dbgmsg('RDBUFFER-READ ERROR',name//mesage)
         call dbgmsg('RDBUFFER-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
C
C
C
      SUBROUTINE RDBUFFERX(STRING, NAME, IERR)
      use error_handling
      implicit none
      CHARACTER STRING*(*), NAME*(*)
      INTEGER   IERR
C
C  *********************************************************************
C  *                                                                   *
C  *  RDBUFFERX:  THIS ROUTINE READS IN AND RETURNS THE ENTIRE INPUT   *
C  *             BUFFER SO IT CAN BE PROCESSED BY THE CALLING ROUTINE  *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *    THIS ROUTINE REQUIRES THE FIRST CHARACTER IN THE BUFFER TO BE  *
C  *    "X" OR IT ISSUES AN ERROE MESSAGE AND SETS THE ERROR FLAG      *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  STRING : CHARACTER STRING                                        *
C  *  NAME   : NAME OF CHARACTER STRING                                *
C  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      STRING = ' '
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

c     Feb/2008 - jde - changed all buffer reads to A256 from A72
c                    - this one left at 512 for reading title - other
c                      entries could be increased to 512 if needed
c                      buffer is 512 - * specifier can not be used
c                      since the input text contains quoted character
c                      strings
  100 IF (IBUF.EQ.0) READ (5,'(A512)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDBUFFERX'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      IF (BUFFER(2:2).ne.'X'.and.buffer(2:2).ne.'x') then 
         MESAGE = 'FIRST CHARACTER IN BUFFER IS NOT "X"'
         string = trim(buffer)
         goto 9999
      endif

      MESAGE = 'RETURNS BUFFER'
      IBUF = 0
      string = trim(buffer)
c      READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
      RETURN
C
 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDBUFFERX-READ ERROR',name//mesage)
         call errmsg('RDBUFFERX-LAST LINE ',buffer)

      else
         call dbgmsg('RDBUFFERX-READ ERROR',name//mesage)
         call dbgmsg('RDBUFFERX-LAST LINE ',buffer)
      endif

c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      RETURN
      END
c
C  *********************************************************************
C  *  RDCAR:  READS A SERIES OF CHARACTER STRINGS                      *
C  *********************************************************************
c
      SUBROUTINE RDCAR(STRINGS,nstrings,maxstrings,NAME, IERR)
      use error_handling
      IMPLICIT  none
      INTEGER   IERR,nstrings,maxstrings
      CHARACTER*(*) STRINGS(maxstrings), NAME
C
C  *********************************************************************
C  *                                                                   *
C  *  RDCAR:  THIS ROUTINE READS IN A SERIES OF CHARACTER STRINGS.     *
C  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
C  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
C  *                                                                   *
C  *      PARAMETERS -                                                 *
C  *  STRINGS   : CHARACTER STRING ARRAY                               *
c  *  NSTRINGS  : NUMBER OF STRINGS READ                               *
c  *  MAXSTRINGS: MAXIMUM NUMBER OF STRINGS TO READ                    *
C  *  NAME      : NAME OF CHARACTER STRING                             *
C  *  IERR      : SET TO 1 IF AN ERROR FOUND                           *
C  *                                                                   *
C  *     CHRIS FARRELL    JAN 1988                                     *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c
      nstrings = 0
c
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (5,'(A256)',ERR=9999,END=9999) BUFFER
      write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDCAR'
c
      IF (BUFFER(1:1).EQ.'$'.or.buffer(1:1).eq.'c'.or.
     >    buffer(1:1).eq.'C') GOTO 100
c slmod begin
      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
        CALL ReadUnstructuredInput(BUFFER)
        GOTO 100
      ENDIF
c slmod end
C
      nstrings = nstrings + 1
c
      if (nstrings.gt.maxstrings) goto 9999
c
      strings(nstrings) = ' '
c
      MESAGE = 'EXPECTING CHARACTER STRING'
      IBUF = 0
c
      strings(nstrings) = buffer
c
c      READ (BUFFER,*,ERR=9999,END=9999) STRINGS(nstrings)
c
      goto 100
c
c
 9999 continue

      if (nstrings.eq.0) then

         if (ierr.eq.0) then 

            IERR = 1

            call errmsg('RDCAR-READ ERROR',name//mesage)
            call errmsg('RDCAR-LAST LINE ',buffer)
 
         else
            call dbgmsg('RDCAR-READ ERROR',name//mesage)
            call dbgmsg('RDCAR-LAST LINE ',buffer)
         endif


c        WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDCAR: ERROR READING - NO PIN DATA',nstrings,NAME,MESAGE,
c     >  'LAST LINE READ :-',BUFFER

      else

        WRITE (9,'(1X,2A,3(/1X,A))')
     >  'RDCAR: FINISHED READING ',nstrings,NAME,MESAGE,
     >  'LAST LINE READ :-',BUFFER

      endif
c
      RETURN
      END
C
C  *********************************************************************
C  *  PRQ:  PRINTS A REAL*8 NUMBER                                    *
C  *********************************************************************
C
      SUBROUTINE PRQ (NAME, R)
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
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
      implicit none
      include 'params' 
      WRITE (datunit,'(1X)')
      RETURN
      END
C
C  *********************************************************************
C  *  PRP:  PRINTS A PAGE THROW                                        *
C  *********************************************************************
C
      SUBROUTINE PRP
      implicit none
      include 'params' 
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
     >            ' NUMBER OF ELEMENTS IS ZERO',
c slmod begin
     >                  nrs,r,rs(1)
c     >                  nrs,r,rs(1),rs(nrs)
c slmod end
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
      IMPLICIT none

      CHARACTER routine*(*),message*(*)

      INCLUDE 'params'
      INCLUDE 'slcom'

      call errmsg(routine,message)
      call errmsg(routine,message,slout)

c      WRITE(0    ,'(4A)') ' ERROR ',routine,': ',message
c      WRITE(SLOUT,'(4A)') ' ERROR ',routine,': ',message

      RETURN 1
      END


      subroutine intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect)
      implicit none   
      real*8 ra,za,rb,zb,r1,z1,r2,z2,rint,zint
      integer sect
c
c     Calculates intersection of two lines and sets flag if the
c     intersection is between the end points of both line segments.
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
            write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                     ' COLINEAR-VERTICAL',warnings
            write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                  r2,z2,ma,m1,ba,b1
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
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-HORIZONTAL',warnings
                 write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                                r2,z2,ma,m1,ba,b1
               else
                  write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-PARALLEL',warnings
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//
     >                         ' COLINEAR-PARALLEL',warnings
                 write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,
     >                                                r2,z2,ma,m1,ba,b1
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
      implicit none
      REAL ARGZ,ARGR
C     INCLUDE "PARAMS"
      include 'params'
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
      implicit none
      REAL*8 ARGZ,ARGR
C     INCLUDE "PARAMS"
      include 'params'
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
      subroutine set_bcomponents(br,bz,bt)
      implicit none
      include 'params'
      include 'cgeom'
c
      real br(maxnks,maxnrs)
      real bz(maxnks,maxnrs)
      real bt(maxnks,maxnrs)
c
c     This routine calculates the direction of the magnetic field vector
c     for each cell on the grid - assuming a cylindrical geometry. 
c
c     This code does not calculate the absolute magnitude of the B field only the 
c     normalized direction vector.
c
c     First need to calculate the cell axis vector for each cell on the 
c     grid. From these - caclulate the BT vector that scales with the 
c     other two and the bratio for the cell. Finally, renormalize the 
c     magnitude of the components so that the bfield vector has a magnitude
c     of 1. 
c
c     Start the bfield calculations at the ik=1 target for simplicity
c       
      integer ik,ir
      real btot
      real btmp
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            br(ik,ir) = krb(ik,ir)-krb(ik-1,ir)
            bz(ik,ir) = kzb(ik,ir)-kzb(ik-1,ir)
c
c           Temporarily set bt to contain Btot - which is 
c           Bpol * kbfs
c
            btot = kbfs(ik,ir) * sqrt(br(ik,ir)**2 + bz(ik,ir)**2)
c
c           Now use the bratio to calculate btotal and thus btoroidal
c
c           bratio = Bpol/B
c           kbfs   = 1/bratio = B/Bpol
c     
c           B = kbfs * Bpol
c
c           Bt = sqrt(B**2 - Br**2 - Bz**2) 
c     
            btmp = btot**2 - br(ik,ir)**2 -bz(ik,ir)**2
            if (btmp.ge.0.0) then 
               bt(ik,ir) = sqrt(btmp)
            else
               write(6,'(a,2i6,10(1x,g12.5))') 
     >              'BFIELD CALCULATION ERROR:',
     >              ik,ir,btot,br(ik,ir),bz(ik,ir),kbfs(ik,ir),btmp
               bt(ik,ir) = 0.0
            endif
c
c           Normalize the b-vector to magnitude 1
c
            if (btot.gt.0.0) then 
               br(ik,ir) = br(ik,ir)/btot
               bz(ik,ir) = bz(ik,ir)/btot
               bt(ik,ir) = bt(ik,ir)/btot
            endif
c
            write(6,'(a,2i6,10(1x,g12.5))') 'BFIELD:',ik,ir,
     >          kbfs(ik,ir),br(ik,ir),bz(ik,ir),bt(ik,ir),btot,
     >        sqrt(br(ik,ir)**2+bz(ik,ir)**2+bt(ik,ir)**2)

         end do
      end do
c
         

      return
      end
