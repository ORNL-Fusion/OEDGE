C     -*-Fortran-*-
C
C
      SUBROUTINE RDRAR(RS,NRS, MAXNRS, RMIN, RMAX, ASCEND, NAME, IERR)
      use error_handling
      use mod_io_units
      use mod_reader
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
c     include 'reader'
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
c                    - added buff_format to common to make buffer size changes easy
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRAR'
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
         call errmsg('RDRAR-LAST LINE ',trim(buffer))

      else

         call dbgmsg('RDRAR-READ ERROR',name//mesage)
         call dbgmsg('RDRAR-LAST LINE ',trim(buffer))

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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
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
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDIARN'
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
         call errmsg('RDIARN-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDIARN-READ ERROR',name//mesage)
         call dbgmsg('RDIARN-LAST LINE ',trim(buffer))
      endif


c      WRITE (7,'(1X,2A,3(/1X,A))')
c     > 'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
c      RETURN
      END
c slmod end
      SUBROUTINE RDRARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,
     >                                                   NFS,NAME,IERR)
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
      INTEGER   IR,N,I
c slmod begin
      REAL      RLAST,R,F(20)
c
c      REAL      RLAST,R,F(10)
c slmod end
      integer :: tmp_ierr
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
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
c
c  100 IF (IBUF.EQ.0) READ (stdin,'(A72)',ERR=9999,END=9999) BUFFER
c slmod end
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'
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
c
c     jde - This hack is needed to maintain input file compatibility
c           when I added functionality to the bgplasopt function. 
c
      if (trim(name).eq. 'SET OF BG PLASMA OPTIONS BY RING') then
         READ (BUFFER,*,iostat=tmp_ierr,END=9999) R,(F(I),I=1,NFS)
         if (tmp_ierr.ne.0) then
c
c           Read the old size of the array  
c            
            WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(9)*N,
     >            ' REALS,',9,' PER LINE'
            READ (BUFFER,*,err=9999,END=9999) R,(F(I),I=1,8)
c
c           Put in placeholders if the rest was read correctly
c            
            do i = 9,12
               f(i) = -1.0
            end do
c
         endif   
c
      else         
         READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
      endif 
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
     > trim(buffer)
      WRITE (0,'(1X,2A,3(/1X,A))')
     > 'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',
     > trim(buffer)
      RETURN
      END
C
C
C
      SUBROUTINE RDR(R, TSTMIN, RMIN, TSTMAX, RMAX, NAME, IERR)
      use error_handling
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R = 0.0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR'
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
         call errmsg('RDR-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDR-READ ERROR',name//mesage)
         call dbgmsg('RDR-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR2'
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
         call errmsg('RDR2-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDR2-READ ERROR',name//mesage)
         call dbgmsg('RDR2-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
      R3 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR3'
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
         call errmsg('RDR3-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDR3-READ ERROR',name//mesage)
         call dbgmsg('RDR3-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     INCLUDE   'reader'
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
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
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
         call errmsg('RDQAR-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDQAR-READ ERROR',name//mesage)
         call dbgmsg('RDQAR-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     INCLUDE   'reader'
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
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQARN'
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
         call errmsg('RDQARN-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDQARN-READ ERROR',name//mesage)
         call dbgmsg('RDQARN-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     INCLUDE   'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R = 0.0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
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
         call errmsg('RDQ-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDQ-READ ERROR',name//mesage)
         call dbgmsg('RDQ-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      R1 = 0.0
      R2 = 0.0
c
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ2'
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
         call errmsg('RDQ2-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDQ2-READ ERROR',name//mesage)
         call dbgmsg('RDQ2-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
c slmod begin
      CHARACTER COMENT*128,MESAGE*72
c
c      CHARACTER COMENT*72,MESAGE*72
c slmod end
C     
      
c      write(0,*) 'RDI:',stdin,i,imin,imax,tstmin,tstmax,ibuf
c      write(0,*)  'RDI2:',len(buffer),':',trim(buffer),':',name

      I = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
 100  IF (IBUF.EQ.0) READ (stdin,buff_format,
     >    iostat=ierr,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI'
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

      write(0,*) 'ERROR EXIT:',ierr

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDI-READ ERROR',name//mesage)
         call errmsg('RDI-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDI-READ ERROR',name//mesage)
         call dbgmsg('RDI-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER COMENT*72,MESAGE*72
C
      I1 = 0
      I2 = 0
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI2'
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
         call errmsg('RDI2-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDI2-READ ERROR',name//mesage)
         call dbgmsg('RDI2-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
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
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDC'
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
         call errmsg('RDC-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDC-READ ERROR',name//mesage)
         call dbgmsg('RDC-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER MESAGE*72
!      CHARACTER COMENT*72,MESAGE*72
C
      STRING = ' '
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

c     Feb/2008 - jde - changed all buffer reads to A256 from A72
c                    - this one left at 512 for reading title - other
c                      entries could be increased to 512 if needed
c                      buffer is 512 - * specifier can not be used
c                      since the input text contains quoted character
c                      strings
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)')
     >              name,':',trim(buffer),'RDBUFFER'
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
         call errmsg('RDBUFFER-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDBUFFER-READ ERROR',name//mesage)
         call dbgmsg('RDBUFFER-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER MESAGE*72
!      CHARACTER COMENT*72,MESAGE*72
C
      STRING = ' '
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

c     Feb/2008 - jde - changed all buffer reads to A256 from A72
c                    - this one left at 512 for reading title - other
c                      entries could be increased to 512 if needed
c                      buffer is 512 - * specifier can not be used
c                      since the input text contains quoted character
c                      strings
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)')
     >             name,':',trim(buffer),'RDBUFFERX'
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
         call errmsg('RDBUFFERX-LAST LINE ',trim(buffer))

      else
         call dbgmsg('RDBUFFERX-READ ERROR',name//mesage)
         call dbgmsg('RDBUFFERX-LAST LINE ',trim(buffer))
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
      use mod_io_units
      use mod_reader
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
c     include 'reader'
      CHARACTER MESAGE*72
!      CHARACTER COMENT*72,MESAGE*72
C
      MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
c
      nstrings = 0
c
c     Feb/2008 - jde - changed all buffer reads to A256 from A72
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDCAR'
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
            call errmsg('RDCAR-LAST LINE ',trim(buffer))
 
         else
            call dbgmsg('RDCAR-READ ERROR',name//mesage)
            call dbgmsg('RDCAR-LAST LINE ',trim(buffer))
         endif


c        WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDCAR: ERROR READING - NO PIN DATA',nstrings,NAME,MESAGE,
c     >  'LAST LINE READ :-',BUFFER

      else

        WRITE (echout,'(1X,2A,3(/1X,A))')
     >  'RDCAR: FINISHED READING ',nstrings,NAME,MESAGE,
     >  'LAST LINE READ :-',trim(BUFFER)

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
      use mod_io_units
      implicit none
c     include 'params' 
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
      !use mod_params
      use mod_io_units
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
!     use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
      WRITE (datunit,'(1X)')
      RETURN
      END
C
C  *********************************************************************
C  *  PRP:  PRINTS A PAGE THROW                                        *
C  *********************************************************************
C
      SUBROUTINE PRP
      use mod_io_units
      !use mod_params
      implicit none
c     include 'params' 
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
      use mod_io_units
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

      WRITE (stddbg,9001) OPT(3:len),REAL(4*N)
c      WRITE (0,9001) OPT(3:len),REAL(4*N)

      return

  300 Write (0,*) 'ERROR READING: ',OPT,' : ERROR=',ierr
      write (stddbg,*) 'ERROR READING: ',OPT,' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',OPT,' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',OPT,' : ERROR=',ierr
      return

 9001 FORMAT(1X,'RINOUT: SIZE OF ',A6,' =',-6P,F9.4,' MB')
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
      use mod_io_units
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

      WRITE (stddbg,9001) OPT(3:len_trim(opt)),REAL(4*N)
      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return
C      IF (4*N.GT.10000) WRITE (stddbg,9001) OPT(3:8),REAL(4*N)
 9001 FORMAT(1X,'DINOUTU: SIZE OF ',A6,' =',-6P,F6.2,' MB')
C
      RETURN
      END
C
C  *********************************************************************
C  *  DINOUT: WRITE ONLY ROUTINE, CONVERTS D.P. TO REAL WHEN WRITING.  *
C  *********************************************************************
C
      SUBROUTINE DINOUT (OPT,DARRAY,N)
      use mod_io_units
      implicit none
      INTEGER I,J,N,IBLOCK,ierr
      CHARACTER OPT*(*)
      DOUBLE PRECISION DARRAY(N)
      DATA IBLOCK /1500/
      real,allocatable :: tmparray(:)
C     
      IF     (OPT(1:1).EQ.'R') THEN
        WRITE (stddbg,*) ' DINOUT: ERROR!  USE ONLY FOR WRITING!'
        STOP
      ELSEIF (OPT(1:1).EQ.'W') THEN

c     jdemod - conversion to real causes a floating point exception
c              when the double precision value is too large to fit in a real           
c     
c     The maximum real value is about 1e38 - to avoid this issue all values
c     > 1.0e38 are converted to 1e38.
c
c
      allocate(tmparray(n))
      tmparray = 0.0
      do i = 1,n
         if (abs(darray(i)).gt.1d38) then
            tmparray(i) = sngl(sign(1d38,darray(i)))
         elseif (abs(darray(i)).lt.1d-37) then
            tmparray(i) = 0.0
         else
            tmparray(i)=sngl(darray(i))
         endif
      end do

c      do i=1,n
c         write(stddbg,'(a,i8,10(1x,g12.5))') 'tmparray:',i,
c     >         tmparray(i),darray(i)
c      end do
      
         
      DO 200 I = 1, N, IBLOCK
c
c           WRITE (8,ERR=400,iostat=ierr)
c     >          (SNGL(DARRAY(J)),J=I,MIN(N,I+IBLOCK-1))

           WRITE (8,ERR=400,iostat=ierr)
     >          (tmpARRAY(J),J=I,MIN(N,I+IBLOCK-1))
c           WRITE (8,ERR=400,iostat=ierr)
c     >          (SNGL(tmpARRAY(J)),J=I,MIN(N,I+IBLOCK-1))
  200   CONTINUE

        if (allocated(tmparray)) deallocate(tmparray)
      ENDIF

      WRITE (stddbg,9001) OPT(3:len_trim(opt)),REAL(4*N)
c      IF (4*N.GT.10000) WRITE (stddbg,9001) OPT(3:len_trim(opt)),REAL(4*N)

      return

!  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
!      write (stddbg,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
!      ! Deallocate temp storage if it has been allocated
!      if (allocated(tmparray)) deallocate(tmparray)
!      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      ! Deallocate temp storage if it has been allocated
      if (allocated(tmparray)) deallocate(tmparray)
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
      use mod_io_units
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
      IF (8*N.GT.10000) WRITE (stddbg,9001)
     >                   OPT(3:len_trim(opt)),REAL(8*N)

      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
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
      use mod_io_units
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
      WRITE (stddbg,9001) OPT(3:len_trim(opt)),REAL(4*N)
      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
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
      use mod_io_units
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
      WRITE (stddbg,9001) OPT(3:len_trim(opt)),REAL(4*N)
      return

  300 Write (0,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR READING: ',trim(OPT),' : ERROR=',ierr
      return

  400 Write (0,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      write (stddbg,*) 'ERROR WRITING: ',trim(OPT),' : ERROR=',ierr
      return
 9001 FORMAT(1X,'IINOUT: SIZE OF ',A6,' =',-6P,F6.2,' MB')
      RETURN
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




c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadIR(line,ival,rval,imin,imax,tag)
      use mod_io_units
      !use mod_params
      !use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*)
      INTEGER fp,ival,imin,imax
      REAL    rval

      INTEGER i
      REAL    r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i,r

      IF (i.LT.imin.OR.i.GT.imax)
     .  CALL ER('ReadI','Out of bounds: '//line,*99)

      ival = i
      rval = r

      WRITE(STDDBG,'(A)') line
      WRITE(STDDBG,'(5X,2A,I4,1P,E10.2)') tag,' = ',ival,rval

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
 99   WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(STDERR,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadI(line,ival,imin,imax,tag)
      use mod_io_units

      !use mod_params
      !use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*)
      INTEGER ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i

      IF (i.LT.imin.OR.i.GT.imax) then 

        write (0,*)  'READI:ERROR:',i,imin,imax 
        CALL ER('ReadI','Out of bounds: '//line,*99)

      endif

      ival = i

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
99    WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(STDERR,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadC(line,cval,tag)
      use mod_io_units
      !use mod_params
      !use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*),cval*(*)

!      INTEGER fp,ival,imin,imax

!      INTEGER i
      CHARACTER comment*72

      !integer erout1

      !erout1 = 0

      READ (line,*,ERR=98,END=98) comment,cval

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,A)') tag,' = ',cval

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
      WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(STDERR,'(5X,2A)')    'CVAL = ''',cval,''''
      STOP 'READC'
      END
c
c ======================================================================
c
c
c
      SUBROUTINE Read2I(line,ival1,ival2,imin,imax,tag)
      use mod_io_units

      !use mod_params
      !use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*)
      INTEGER ival1,ival2,imin,imax
!      INTEGER fp,ival1,ival2,imin,imax

      INTEGER i1,i2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i1,i2

      IF (i1.LT.imin.OR.i1.GT.imax.OR.
     .    i2.LT.imin.OR.i2.GT.imax)
     .  CALL ER('Read2I','Out of bounds: '//line,*99)

      ival1 = i1
      ival2 = i2

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival1
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival2

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
99    WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE ReadR(line,rval,rmin,rmax,tag)
      use mod_io_units

      !use mod_params
      !use mod_slcom
      IMPLICIT none

      CHARACTER line*(*),tag*(*)
      REAL rval,rmin,rmax

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      REAL r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r

      IF (r.LT.rmin.OR.r.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval = r

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(2A,G10.3)') tag,' = ',rval

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
99    WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(STDERR,'(5X,A,3G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r,rval,rmin,rmax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE Read2R(line,rval1,rval2,rmin,rmax,tag)
      use mod_io_units

      !use mod_params
      !use mod_slcom
      IMPLICIT none

      CHARACTER line*(*),tag*(*)
      REAL rval1,rval2,rmin,rmax

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      REAL r1,r2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r1,r2

      IF (r1.LT.rmin.OR.r1.GT.rmax.OR.
     .    r2.LT.rmin.OR.r2.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval1 = r1
      rval2 = r2

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(2A,2G10.3)') tag,' = ',rval1,rval2

      RETURN
98    WRITE(STDERR,*) 'Problem reading unstructured input'
99    WRITE(STDERR,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(STDERR,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(STDERR,'(5X,A,6G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r1,r2,rval1,rval2,rmin,rmax
      STOP
      END

c
c ======================================================================
c
c subroutine: ISet
c
c
      SUBROUTINE ISet(array,ndim,val)
      implicit none   

      INTEGER ndim
      INTEGER array(ndim),val

      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

      RETURN
      END

c
c ======================================================================
c
c subroutine: RSet
c
c
      SUBROUTINE RSet(array,ndim,val)
      implicit none 
      INTEGER ndim
      REAL array(ndim),val

      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

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
      use mod_io_units
      !use mod_params
      !use mod_slcom
      IMPLICIT none

      CHARACTER routine*(*),message*(*)

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      call errmsg(routine,message)
      call errmsg(routine,message,stddbg)

c      WRITE(0    ,'(4A)') ' ERROR ',routine,': ',message
c      WRITE(SLOUT,'(4A)') ' ERROR ',routine,': ',message

      RETURN 1
      END

C
C
C
      SUBROUTINE ADASRD(YEAR,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)     
c      SUBROUTINE ADASRD(YEAR,YEARDF,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)
C
C  READ THE REQUESTED RATE COEFFICIENT FROM THE ADAS MASTER ELEMENT
C  FILES:
C        ICLASS = 1: RECOMBINATION RATE COEFFICIENT
C                 2: IONISATION RATE COEFFICIENT
C                 3: CHARGE EXCHANGE RECOMBINATION COEFFICIENT
C                 4: POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLU
C                 5: POWER COEFFICIENT FOR LINE RADIATION
C                 6: POWER COEFFICIENT FOR CHARGE EXCHANGE
C                 added by Krieger, IPP 5/95
C                 7: PHOTON EMISSIVITY FOR DIAGNOSTIC LINES
C  THIS ROUTINE USES THE STANDARD ADAS EXTRACTION ROUTINE D2DATA AND
C  REALLY ONLY PROVIDES A 'CLEAN' INTERFACE, TAKING CARE OF CHANGES
C  IN UNITS AND IN PRECISION OF VARIABLES.  IF THE REQUESTED DATA
C  DOESN'T EXIST (IFAIL=1 RETURNED FROM D2DATA) THE PROGRAM IS STOPPED.
C
      !use mod_params
      use mod_cadas2
      IMPLICIT NONE
C     INCLUDE   "PARAMS"
c     include    'params'
C     INCLUDE   "CADAS2"
c     include    'cadas2'
C
      CHARACTER*2 YEAR
      character*80 class
      INTEGER IZ0, IZ1, ICLASS, NPTS,len,lenstr
      external lenstr 
      REAL TE(NPTS), NE(NPTS), COEF(NPTS)
c
c      logical lintrp(maxpts)
C
      INTEGER J
!      INTEGER I, J
C
C
c     additional diagnostic output; Krieger IPP/97
c
      class = ' '
      if (iclass.eq.1) then
        class = 'RECOMBINATION RATE COEFFICIENT'
      else if (iclass.eq.2) then
        class = 'IONISATION RATE COEFFICIENT'
      else if (iclass.eq.3) then
        class = 'CHARGE EXCHANGE RECOMBINATION COEFFICIENT'
      else if (iclass.eq.4) then
        class = 'POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLUNG'
      else if (iclass.eq.5) then
        class = 'POWER COEFFICIENT FOR LINE RADIATION'
      else if (iclass.eq.6) then
        class = 'POWER COEFFICIENT FOR CHARGE EXCHANGE'
      else if (iclass.eq.7) then
        class = 'PHOTON EMISSIVITY FOR DIAGNOSTIC LINES'
      endif
c
      IEVCUT = 0
C
      DO J = 1, NPTS
        DTEV(J) = DBLE(ALOG10(TE(J)))
        DDENS(J) = DBLE(ALOG10(NE(J)*1.0E-6))
      ENDDO

      CALL D2DATA(YEAR, YEAR, TITLF, IFAIL,
     >            IZ0, IZ1, ICLASS, NPTS, IEVCUT,
     >            MAXADS, ITMAXD, IDMAXD, IZMAXD,
     >            DTEV, DDENS,
     >            DTEVD, DDENSD, DRCOFD, ZDATA,
     >            DRCOFI
c     >            , LINTRP
     >            )
c
      IF (IFAIL.EQ.1) THEN
        len = lenstr(class) 
        WRITE(6,1000) IZ0, IZ1,YEAR
        write(6,1001) iclass, class(1:len)
        WRITE(7,1000) IZ0, IZ1,YEAR
        write(7,1001) iclass, class(1:len)
        WRITE(0,1000) IZ0, IZ1,YEAR
        write(0,1001) iclass, class(1:len)
        STOP
      ENDIF
C
C  EXTRAPOLATED VALUES ARE RETURNED AS ZERO!
C
      DO J = 1, NPTS
        IF (DRCOFI(J).NE.0.0) THEN
          COEF(J) = 10.**SNGL(DRCOFI(J)) * 1.0E-6
        ELSE
          COEF(J) = 0.0
c
          len = lenstr(class)
          write (6,'(a,4i4,3(1x,g12.5))') 
     >           'Extrapolated coefficient:'//class(1:len),
     >            j,iz0,iz1,npts,te(j),ne(j)
c
        ENDIF

      ENDDO
C
 1000 FORMAT(' ERROR READING REQUESTED ATOMIC DATA!',/,
     >       ' MASTER ELEMENT FILE FOR NUCLEAR CHARGE ',I2,
     >       ' AND ION CHARGE ',I2,
     >       ' WAS NOT FOUND IN YEAR ',A2)
 1001 format(' CLASS ',i2,1x,a)
C
      RETURN
      END
