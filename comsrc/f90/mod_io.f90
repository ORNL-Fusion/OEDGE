module mod_io

  use mod_params
  use mod_reader
  use mod_io_units
  use allocate_arrays
  
  
  implicit none


  ! Static or pre-allocated inputs
  interface divrd

     module procedure rdr,rdr2,rdr3,rdrar,rdrarn,&
          rdi,rdi2,rdiarn, &
          rdq,rdq2,rdqar,rdqarn, &
          rdc,rdcar, &
          readir, readi, readc, read2i, readr, read2r
  end interface divrd
     
  ! allocatable inputs
  interface divrda
     module procedure rdrarn_alloc,rdrarn_alloc_check,rdiarn_alloc,rdrar_alloc
  end interface divrda
     

contains

  



  ! =======================================================
  ! 
  ! Original READ routines
  !
  ! =======================================================


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
    !
    !  *********************************************************************
    !  *                                                                   *
    !  *  RDRAR:   ROUTINE READS A REAL VARIABLE LENGTH 1D ARRAY           *
    !  *    AND CHECKS THAT IT IS SORTED IN ASCENDING ORDER,               *
    !  *    THAT ALL THE VALUES ARE WITHIN A RANGE AND THAT NO             *
    !  *    TWO VALUES ARE EQUAL.                                          *
    !  *    ANY OF THESE RULES BEING BROKEN RESULTS IN THE                 *
    !  *    OFFENDING VALUE BEING REMOVED FROM THE ARRAY, AN               *
    !  *    ERROR MESSAGE BEING OUTPUT AND AN ERROR FLAG BEING             *
    !  *    SET.                                                           *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  RS     : ARRAY VALUES RETURNED IN THIS                           *
    !  *  NRS    : NUMBER OF VALUES READ RETURNED IN THIS                  *
    !  *  MAXNRS : MAXIMUM NUMBER OF VALUES TO BE READ                     *
    !  *  RMIN   : MINIMUM ALLOWED VALUE (EXCLUSIVE)                       *
    !  *  RMAX   : MAXIMUM ALLOWED VALUE (EXCLUSIVE)                       *
    !  *  ASCEND : INDICATES WHETHER ASCENDING ORDER CHECK IS TO BE MADE   *
    !  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL   JAN 1988                                       *
    !  *                                                                   *
    !  *********************************************************************


    CHARACTER COMENT*72,MESAGE*72
    INTEGER   IR,N
    REAL      RLAST,R
    !
    !---- READ IN AND TEST SIZE OF ARRAY
    !
    NRS = 0
    CALL RDC (COMENT, NAME, IERR)
    CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
    IF (N.EQ.0) RETURN
    !
    !---- READ IN AND TEST EACH ARRAY VALUE
    !
    RLAST = RMIN
    IR = 1
50  CONTINUE
    MESAGE = 'END OF FILE ON UNIT 5'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
    !                    - added buff_format to common to make buffer size changes easy
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRAR'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end
    !
    WRITE (MESAGE,'(A,I5,A)') 'EXPECTING',N,' REALS, ONE PER LINE'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) R
    !
    IF (ASCEND) THEN
       WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
       IF (R.LT.RLAST) GOTO 9999
    ENDIF
    !
    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
    IF (R.GT.RMAX) GOTO 9999
    RS(IR) = R
    IR = IR + 1
    IF (IR.LE.N) GOTO 50
    !
    !---- SET UP NUMBER OF VALID VALUES READ
    !
    NRS = IR - 1
    RETURN
    !
9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDRAR-READ ERROR',name//mesage)
       call errmsg('RDRAR-LAST LINE ',trim(buffer))

    else

       call dbgmsg('RDRAR-READ ERROR',name//mesage)
       call dbgmsg('RDRAR-LAST LINE ',trim(buffer))

    endif


    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     > 'RDRAR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDRAR


  !
  !
  !
  ! slmod begin 
  !
  ! This routine is based on RDRARN but has been modified to read an array
  ! of integers, rather than reals. -March, 2001
  !
  SUBROUTINE RDIARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,NFS,NAME,IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    INTEGER   NRS, MAXNRS, NFS
    INTEGER   RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
    !      REAL      RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
    CHARACTER NAME*(*)
    INTEGER   IERR
    LOGICAL   ASCEND

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDIARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
    !  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
    !  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
    !  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
    !  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
    !  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
    !  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
    !  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
    !  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
    !  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
    !  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
    !  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
    !  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72
    INTEGER   IR,N,I
    INTEGER   RLAST,R,F(10)
    !
    !
    !---- READ IN AND TEST SIZE OF ARRAY
    !
    NRS = 0
    CALL RDC (COMENT, NAME, IERR)
    CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
    IF (N.EQ.0) RETURN
    !
    !---- READ IN AND TEST EACH SET OF ARRAY VALUES
    !
    RLAST = RMIN
    IR = 1
50  CONTINUE
    MESAGE = 'END OF FILE ON UNIT 5'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDIARN'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    !      IF (BUFFER(2:2).EQ.'*') THEN
    !        CALL ReadUnstructuredInput(BUFFER)
    !        GOTO 100
    !      ENDIF
    ! slmod end
    !
    WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' INTS,',NFS+1,' PER LINE'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
    IF (ASCEND) THEN
       WRITE (MESAGE,'(I11,A,I11)') R,' LESS THAN PREV/MIN',RLAST
       IF (R.LT.RLAST) GOTO 9999
    ENDIF
    WRITE (MESAGE,'(I11,A,I11)') R,' MORE THAN MAXIMUM',RMAX
    IF (R.GT.RMAX) GOTO 9999
    RS(IR,1) = R

    DO I = 1, NFS
       WRITE (MESAGE,'(I11,A,I11)') F(I),' LESS THAN MINIMUM',FMIN
       IF (F(I).LT.FMIN) GOTO 9999
       WRITE (MESAGE,'(I11,A,I11)') F(I),' MORE THAN MAXIMUM',FMAX
       IF (F(I).GT.FMAX) GOTO 9999
       RS(IR,1+I) = F(I)
    end do
    !
    !---- SET UP NUMBER OF VALID VALUES READ
    !
    IR = IR + 1
    IF (IR.LE.N) GOTO 50
    NRS = IR - 1
    RETURN
9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDIARN-READ ERROR',name//mesage)
       call errmsg('RDIARN-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDIARN-READ ERROR',name//mesage)
       call dbgmsg('RDIARN-LAST LINE ',trim(buffer))
    endif


    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     > 'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
    !      RETURN
  END SUBROUTINE RDIARN


  SUBROUTINE RDRARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,NFS,NAME,IERR)
    use mod_io_units
    use mod_reader
    implicit none
    INTEGER   NRS, MAXNRS, NFS
    REAL      RS(MAXNRS,1+NFS), RMIN, RMAX,FMIN,FMAX
    CHARACTER NAME*(*)
    INTEGER   IERR
    LOGICAL   ASCEND

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
    !  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
    !  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
    !  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
    !  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
    !  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
    !  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
    !  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
    !  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
    !  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
    !  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
    !  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
    !  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72
    INTEGER   IR,N,I
    REAL      RLAST,R,F(20)
    integer :: tmp_ierr
    !     
    !---- READ IN AND TEST SIZE OF ARRAY
    !
    NRS = 0
    IERR = 0

    CALL RDC (COMENT, NAME, IERR)
    CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)

    IF (N.EQ.0) RETURN
    !
    !---- READ IN AND TEST EACH SET OF ARRAY VALUES
    !
    RLAST = RMIN
    IR = 1
50  CONTINUE
    MESAGE = 'END OF FILE ON UNIT 5'
    ! slmod begin
    !... This 72 character limit has always been annoying, and I can't see
    !    any reason not to increase it since BUFFER*512 is declared
    !    in READER:
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    !
    !  100 IF (IBUF.EQ.0) READ (stdin,'(A72)',ERR=9999,END=9999) BUFFER
    ! slmod end
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end
    !
    WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'
    IBUF = 0
    !
    !     jde - This hack is needed to maintain input file compatibility
    !           when I added functionality to the bgplasopt function. 
    !
    if (trim(name).eq. 'SET OF BG PLASMA OPTIONS BY RING') then
       READ (BUFFER,*,iostat=tmp_ierr,END=9999) R,(F(I),I=1,NFS)
       if (tmp_ierr.ne.0) then
          !
          !           Read the old size of the array  
          !            
          WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(9)*N,' REALS,',9,' PER LINE'
          READ (BUFFER,*,err=9999,END=9999) R,(F(I),I=1,8)
          !
          !           Put in placeholders if the rest was read correctly
          !            
          do i = 9,12
             f(i) = -1.0
          end do

       endif

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

    DO I = 1, NFS
       WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
       IF (F(I).LT.FMIN) GOTO 9999
       WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
       IF (F(I).GT.FMAX) GOTO 9999
       RS(IR,1+I) = F(I)
    end do
    !
    !---- SET UP NUMBER OF VALID VALUES READ
    !
    IR = IR + 1
    IF (IR.LE.N) GOTO 50
    NRS = IR - 1
    RETURN

9999 IERR = 1
    WRITE (7,'(1X,2A,3(/1X,A))')'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(buffer)
    WRITE (stderr,'(1X,2A,3(/1X,A))')'RDRARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(buffer)
    RETURN
  END SUBROUTINE RDRARN



  SUBROUTINE RDR(R, TSTMIN, RMIN, TSTMAX, RMAX, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    REAL      R, RMIN, RMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDR:  ROUTINE READS IN A REAL.
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  R      : VALUE                                                   *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  RMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  RMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    R = 0.0
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end
    !
    MESAGE = 'EXPECTING COMMENT & REAL NUMBER'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,R

    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS LESS THAN MINIMUM ', RMIN
    IF (TSTMIN.AND.(R.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS MORE THAN MAXIMUM ', RMAX
    IF (TSTMAX.AND.(R.GT.RMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDR-READ ERROR',name//mesage)
       call errmsg('RDR-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDR-READ ERROR',name//mesage)
       call dbgmsg('RDR-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDR




  SUBROUTINE RDR2(R1, R2,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    IMPLICIT  NONE
    REAL      R1, R2, RMIN, RMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR
    !  *********************************************************************
    !  *                                                                   *
    !  *  RDR2:  ROUTINE READS IN 2 REALS.
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  R1,R2  : VALUES                                                  *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  RMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  RMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    R1 = 0.0
    R2 = 0.0

    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR2'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & 2 REAL NUMBERS'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2

    !     Test minimums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999

    !     Test Maximums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDR2-READ ERROR',name//mesage)
       call errmsg('RDR2-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDR2-READ ERROR',name//mesage)
       call dbgmsg('RDR2-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDR2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
    RETURN
  END SUBROUTINE RDR2



  SUBROUTINE RDR3(R1, R2, R3,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    IMPLICIT  NONE
    REAL      R1, R2, R3, RMIN, RMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDR3:  ROUTINE READS IN 3 REALS.
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  R1,R2,R3 : VALUES                                                *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  RMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  RMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    R1 = 0.0
    R2 = 0.0
    R3 = 0.0

    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDR3'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & 3 REAL NUMBERS'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2,R3

    !     Test minimums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R3,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R3.LT.RMIN)) GOTO 9999

    !     Test Maximums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R3,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R3.GT.RMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDR3-READ ERROR',name//mesage)
       call errmsg('RDR3-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDR3-READ ERROR',name//mesage)
       call dbgmsg('RDR3-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDR3: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
    RETURN
  END SUBROUTINE RDR3



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

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDQAR:   ROUTINE READS A REAL VARIABLE LENGTH 1D ARRAY           *
    !  *    AND CHECKS THAT IT IS SORTED IN ASCENDING ORDER,               *
    !  *    THAT ALL THE VALUES ARE WITHIN A RANGE AND THAT NO             *
    !  *    TWO VALUES ARE EQUAL.                                          *
    !  *    ANY OF THESE RULES BEING BROKEN RESULTS IN THE                 *
    !  *    OFFENDING VALUE BEING REMOVED FROM THE ARRAY, AN               *
    !  *    ERROR MESSAGE BEING OUTPUT AND AN ERROR FLAG BEING             *
    !  *    SET.                                                           *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  RS     : ARRAY VALUES RETURNED IN THIS                           *
    !  *  NRS    : NUMBER OF VALUES READ RETURNED IN THIS                  *
    !  *  MAXNRS : MAXIMUM NUMBER OF VALUES TO BE READ                     *
    !  *  RMIN   : MINIMUM ALLOWED VALUE (EXCLUSIVE)                       *
    !  *  RMAX   : MAXIMUM ALLOWED VALUE (EXCLUSIVE)                       *
    !  *  ASCEND : INDICATES WHETHER ASCENDING ORDER CHECK IS TO BE MADE   *
    !  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL   JAN 1988                                       *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72
    INTEGER   IR,N
    real*8   RLAST,R

    !---- READ IN AND TEST SIZE OF ARRAY

    NRS = 0
    CALL RDC (COMENT, NAME, IERR)
    CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
    IF (N.EQ.0) RETURN

    !---- READ IN AND TEST EACH ARRAY VALUE

    RLAST = RMIN
    IR = 1
50  CONTINUE
    MESAGE = 'END OF FILE ON UNIT 5'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    WRITE (MESAGE,'(A,I5,A)') 'EXPECTING',N,' REALS, ONE PER LINE'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) R

    IF (ASCEND) THEN
       WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
       IF (R.LT.RLAST) GOTO 9999
    ENDIF

    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
    IF (R.GT.RMAX) GOTO 9999
    RS(IR) = R
    IR = IR + 1
    IF (IR.LE.N) GOTO 50

    !---- SET UP NUMBER OF VALID VALUES READ

    NRS = IR - 1
    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDQAR-READ ERROR',name//mesage)
       call errmsg('RDQAR-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDQAR-READ ERROR',name//mesage)
       call dbgmsg('RDQAR-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     > 'RDQAR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDQAR



  SUBROUTINE RDQARN (RS,NRS,MAXNRS,RMIN,RMAX,ASCEND,FMIN,FMAX,NFS,NAME,IERR)
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

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDQARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)   *
    !  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
    !  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
    !  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
    !  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
    !  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
    !  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
    !  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).  MAX OF 10 FUNCTIONS ALLOWED*
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
    !  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
    !  *  MAXNRS : MAXIMUM NUMBER OF SETS OF VALUES TO BE READ             *
    !  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
    !  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
    !  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
    !  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
    !  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72
    INTEGER   IR,N,I
    real*8   RLAST,R,F(10)

    !---- READ IN AND TEST SIZE OF ARRAY

    NRS = 0
    CALL RDC (COMENT, NAME, IERR)
    CALL RDI (N, .TRUE., 0 ,.TRUE., MAXNRS, NAME, IERR)
    IF (N.EQ.0) RETURN

    !---- READ IN AND TEST EACH SET OF ARRAY VALUES

    RLAST = RMIN
    IR = 1
50  CONTINUE
    MESAGE = 'END OF FILE ON UNIT 5'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQARN'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
    IF (ASCEND) THEN
       WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
       IF (R.LT.RLAST) GOTO 9999
    ENDIF
    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
    IF (R.GT.RMAX) GOTO 9999
    RS(IR,1) = R

    DO I = 1, NFS
       WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
       IF (F(I).LT.FMIN) GOTO 9999
       WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
       IF (F(I).GT.FMAX) GOTO 9999
       RS(IR,1+I) = F(I)
    end do

    !---- SET UP NUMBER OF VALID VALUES READ

    IR = IR + 1
    IF (IR.LE.N) GOTO 50
    NRS = IR - 1
    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDQARN-READ ERROR',name//mesage)
       call errmsg('RDQARN-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDQARN-READ ERROR',name//mesage)
       call dbgmsg('RDQARN-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     > 'RDQARN: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
    RETURN
  END SUBROUTINE RDQARN



  SUBROUTINE RDQ(R, TSTMIN, RMIN, TSTMAX, RMAX, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    real*8  R, RMIN, RMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDQ:  ROUTINE READS IN A REAL*8                                 *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  R      : VALUE                                                   *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  RMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  RMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    R = 0.0
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & REAL NUMBER'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,R

    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS LESS THAN MINIMUM ', RMIN
    IF (TSTMIN.AND.(R.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R,' IS MORE THAN MAXIMUM ', RMAX
    IF (TSTMAX.AND.(R.GT.RMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDQ-READ ERROR',name//mesage)
       call errmsg('RDQ-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDQ-READ ERROR',name//mesage)
       call dbgmsg('RDQ-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDQ: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDQ



  SUBROUTINE RDQ2(R1, R2,TSTMIN,RMIN,TSTMAX,RMAX,NAME,IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    IMPLICIT  NONE
    REAL*8      R1, R2, RMIN, RMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDQ2:  ROUTINE READS IN 2 REAL*8.                               *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  R1,R2  : VALUES                                                  *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  RMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  RMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    R1 = 0.0
    R2 = 0.0

    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDQ2'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & 2 REAL NUMBERS'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,R1,R2

    !     Test minimums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R1.LT.RMIN)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS LESS THAN MINIMUM ',RMIN
    IF (TSTMIN.AND.(R2.LT.RMIN)) GOTO 9999

    !     Test Maximums

    WRITE (MESAGE,'(G11.4,A,G11.4)') R1,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R1.GT.RMAX)) GOTO 9999

    WRITE (MESAGE,'(G11.4,A,G11.4)') R2,' IS MORE THAN MAXIMUM ',RMAX
    IF (TSTMAX.AND.(R2.GT.RMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDQ2-READ ERROR',name//mesage)
       call errmsg('RDQ2-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDQ2-READ ERROR',name//mesage)
       call dbgmsg('RDQ2-LAST LINE ',trim(buffer))
    endif


    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDQ2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDQ2




  SUBROUTINE RDI(I, TSTMIN, IMIN, TSTMAX, IMAX, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    INTEGER   I, IMIN, IMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDI:  THIS ROUTINE READS IN AN INTEGER.                          *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  I      : VALUE                                                   *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  IMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  IMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL   JAN 1988                                       *
    !  *                                                                   *
    !  *********************************************************************

    ! slmod begin
    CHARACTER COMENT*128,MESAGE*72

    !      CHARACTER COMENT*72,MESAGE*72
    ! slmod end

    I = 0
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & IN TEGER'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,I

    WRITE (MESAGE,'(I11,A,I11)') I,' IS LESS THAN MINIMUM ', IMIN
    IF (TSTMIN.AND.(I.LT.IMIN)) GOTO 9999

    WRITE (MESAGE,'(I11,A,I11)') I,' IS MORE THAN MAXIMUM ', IMAX
    IF (TSTMAX.AND.(I.GT.IMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDI-READ ERROR',name//mesage)
       call errmsg('RDI-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDI-READ ERROR',name//mesage)
       call dbgmsg('RDI-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDI: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDI



  SUBROUTINE RDI2(I1, I2, TSTMIN, IMIN, TSTMAX, IMAX, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    INTEGER   I1, I2, IMIN, IMAX
    LOGICAL   TSTMIN, TSTMAX
    CHARACTER NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDI:  THIS ROUTINE READS IN TWO INTEGERS.                        *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  I      : VALUE                                                   *
    !  *  TSTMIN : IF TRUE VALUE MUST NOT BE LESS THAN RMIN                *
    !  *  IMIN   : MINIMUM VALUE                                           *
    !  *  TSTMAX : IF TRUE VALUE MUST NOT BE GREATER THAN RMAX             *
    !  *  IMAX   : MAXIMUM VALUE                                           *
    !  *  NAME   : NAME OF VARIABLE (FOR ERROR MESSAGES)                   *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *    CHRIS FARRELL   JAN 1988                                       *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    I1 = 0
    I2 = 0
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'
    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDI2'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & 2 INTEGERS'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,I1,I2

    !     Test Min

    WRITE (MESAGE,'(I11,A,I11)') I1,' IS LESS THAN MINIMUM ', IMIN
    IF (TSTMIN.AND.(I1.LT.IMIN)) GOTO 9999

    WRITE (MESAGE,'(I11,A,I11)') I2,' IS LESS THAN MINIMUM ', IMIN
    IF (TSTMIN.AND.(I2.LT.IMIN)) GOTO 9999

    !     Test Max

    WRITE (MESAGE,'(I11,A,I11)') I1,' IS MORE THAN MAXIMUM ', IMAX
    IF (TSTMAX.AND.(I1.GT.IMAX)) GOTO 9999

    WRITE (MESAGE,'(I11,A,I11)') I2,' IS MORE THAN MAXIMUM ', IMAX
    IF (TSTMAX.AND.(I2.GT.IMAX)) GOTO 9999

    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDI2-READ ERROR',name//mesage)
       call errmsg('RDI2-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDI2-READ ERROR',name//mesage)
       call dbgmsg('RDI2-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDI2: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDI2



  SUBROUTINE RDC(STRING, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    CHARACTER STRING*(*), NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDC:  THIS ROUTINE READS IN A CHARACTER STRING.                  *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  STRING : CHARACTER STRING                                        *
    !  *  NAME   : NAME OF CHARACTER STRING                                *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    STRING = ' '
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
    !                    - this one left at 512 for reading title - other
    !                      entries could be increased to 512 if needed
    !                      buffer is 512 - * specifier can not be used
    !                      since the input text contains quoted character
    !                      strings
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDC'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'EXPECTING COMMENT & CHARACTER STRING'
    IBUF = 0
    READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDC-READ ERROR',name//mesage)
       call errmsg('RDC-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDC-READ ERROR',name//mesage)
       call dbgmsg('RDC-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDC



  SUBROUTINE RDBUFFER(STRING, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    CHARACTER STRING*(*), NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDBUFFER:  THIS ROUTINE READS IN AND RETURNS THE ENTIRE INPUT    *
    !  *             BUFFER SO IT CAN BE PROCESSED BY THE CALLING ROUTINE  *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  STRING : CHARACTER STRING                                        *
    !  *  NAME   : NAME OF CHARACTER STRING                                *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    STRING = ' '
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
    !                    - this one left at 512 for reading title - other
    !                      entries could be increased to 512 if needed
    !                      buffer is 512 - * specifier can not be used
    !                      since the input text contains quoted character
    !                      strings
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDBUFFER'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    MESAGE = 'RETURNS BUFFER'
    IBUF = 0
    string = trim(buffer)
    !      READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDBUFFER-READ ERROR',name//mesage)
       call errmsg('RDBUFFER-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDBUFFER-READ ERROR',name//mesage)
       call dbgmsg('RDBUFFER-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDBUFFER



  SUBROUTINE RDBUFFERX(STRING, NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    implicit none
    CHARACTER STRING*(*), NAME*(*)
    INTEGER   IERR

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDBUFFERX:  THIS ROUTINE READS IN AND RETURNS THE ENTIRE INPUT   *
    !  *             BUFFER SO IT CAN BE PROCESSED BY THE CALLING ROUTINE  *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *    THIS ROUTINE REQUIRES THE FIRST CHARACTER IN THE BUFFER TO BE  *
    !  *    "X" OR IT ISSUES AN ERROE MESSAGE AND SETS THE ERROR FLAG      *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  STRING : CHARACTER STRING                                        *
    !  *  NAME   : NAME OF CHARACTER STRING                                *
    !  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    STRING = ' '
    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
    !                    - this one left at 512 for reading title - other
    !                      entries could be increased to 512 if needed
    !                      buffer is 512 - * specifier can not be used
    !                      since the input text contains quoted character
    !                      strings
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDBUFFERX'
    IF (BUFFER(1:1).EQ.'$') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    IF (BUFFER(2:2).ne.'X'.and.buffer(2:2).ne.'x') then 
       MESAGE = 'FIRST CHARACTER IN BUFFER IS NOT "X"'
       string = trim(buffer)
       goto 9999
    endif

    MESAGE = 'RETURNS BUFFER'
    IBUF = 0
    string = trim(buffer)
    !      READ (BUFFER,*,ERR=9999,END=9999) COMENT,STRING
    RETURN

9999 continue

    if (ierr.eq.0) then 

       IERR = 1

       call errmsg('RDBUFFERX-READ ERROR',name//mesage)
       call errmsg('RDBUFFERX-LAST LINE ',trim(buffer))

    else
       call dbgmsg('RDBUFFERX-READ ERROR',name//mesage)
       call dbgmsg('RDBUFFERX-LAST LINE ',trim(buffer))
    endif

    !      WRITE (7,'(1X,2A,3(/1X,A))')
    !     >  'RDC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

    RETURN
  END SUBROUTINE RDBUFFERX





  !  *********************************************************************
  !  *  RDCAR:  READS A SERIES OF CHARACTER STRINGS                      *
  !  *********************************************************************

  SUBROUTINE RDCAR(STRINGS,nstrings,maxstrings,NAME, IERR)
    use error_handling
    use mod_io_units
    use mod_reader
    IMPLICIT  none
    INTEGER   IERR,nstrings,maxstrings
    CHARACTER*(*) STRINGS(maxstrings), NAME

    !  *********************************************************************
    !  *                                                                   *
    !  *  RDCAR:  THIS ROUTINE READS IN A SERIES OF CHARACTER STRINGS.     *
    !  *    IF IT IS INVALID OR AN ERROR OCCURS THEN AN ERROR FLAG         *
    !  *    IS SET AND AN ERROR MESSAGE IS OUTPUT.                         *
    !  *                                                                   *
    !  *      PARAMETERS -                                                 *
    !  *  STRINGS   : CHARACTER STRING ARRAY                               *
    !  *  NSTRINGS  : NUMBER OF STRINGS READ                               *
    !  *  MAXSTRINGS: MAXIMUM NUMBER OF STRINGS TO READ                    *
    !  *  NAME      : NAME OF CHARACTER STRING                             *
    !  *  IERR      : SET TO 1 IF AN ERROR FOUND                           *
    !  *                                                                   *
    !  *     CHRIS FARRELL    JAN 1988                                     *
    !  *                                                                   *
    !  *********************************************************************

    CHARACTER COMENT*72,MESAGE*72

    MESAGE = 'PROBLEM WITH UNIT 5.  END OF FILE?'

    nstrings = 0

    !     Feb/2008 - jde - changed all buffer reads to A256 from A72
100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9999,END=9999) BUFFER
    write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDCAR'

    IF (BUFFER(1:1).EQ.'$'.or.buffer(1:1).eq.'c'.or.buffer(1:1).eq.'C') GOTO 100
    ! slmod begin
    IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
       CALL ReadUnstructuredInput(BUFFER)
       GOTO 100
    ENDIF
    ! slmod end

    nstrings = nstrings + 1

    if (nstrings.gt.maxstrings) goto 9999

    strings(nstrings) = ' '

    MESAGE = 'EXPECTING CHARACTER STRING'
    IBUF = 0

    strings(nstrings) = buffer

    !      READ (BUFFER,*,ERR=9999,END=9999) STRINGS(nstrings)

    goto 100


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


       !        WRITE (7,'(1X,2A,3(/1X,A))')
       !     >  'RDCAR: ERROR READING - NO PIN DATA',nstrings,NAME,MESAGE,
       !     >  'LAST LINE READ :-',BUFFER

    else

       WRITE (echout,'(1X,2A,3(/1X,A))')'RDCAR: FINISHED READING ',nstrings,NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)

    endif

    RETURN
  END SUBROUTINE RDCAR















  ! ======================================================================
  !
  ! Original Unstructured Input read routines
  !
  ! ======================================================================



  SUBROUTINE ReadIR(line,ival,rval,imin,imax,tag)
    use mod_params
    !use mod_slcom
    IMPLICIT none

    CHARACTER line*(*),tag*(*)
    INTEGER fp,ival,imin,imax
    REAL    rval

    INTEGER i
    REAL    r
    CHARACTER comment*72

    READ (line,*,ERR=98,END=98) comment,i,r

    !only check imax IF imax>imin - always check imin
    if (i.lt.imin.or.(i.gt.imax.and.imin.lt.imax)) call er('ReadI','Out of bounds: '//line,*99)

    ival = i
    rval = r

    WRITE(SLOUT,'(A)') line
    WRITE(SLOUT,'(5X,2A,I4,1P,E10.2)') tag,' = ',ival,rval

    RETURN
98  WRITE(EROUT,*) 'Problem reading unstructured input'
99  WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
    WRITE(EROUT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
    STOP
  END SUBROUTINE ReadIR


  ! ======================================================================



  SUBROUTINE ReadI(line,ival,imin,imax,tag)

    use mod_params
    !use mod_slcom
    IMPLICIT none

    !     INCLUDE 'params'
    !     INCLUDE 'slcom'

    CHARACTER line*(*),tag*(*)
    INTEGER fp,ival,imin,imax

    INTEGER i
    CHARACTER comment*72

    READ (line,*,ERR=98,END=98) comment,i

    IF (i.LT.imin.OR.i.GT.imax) then 

       write (0,*)  'READI:ERROR:',i,imin,imax 
       CALL ER('ReadI','Out of bounds: '//line,*99)

    endif

    ival = i

    WRITE(SLOUT,'(A)')        line
    WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival

    RETURN
98  WRITE(EROUT,*) 'Problem reading unstructured input'
99  WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
    WRITE(EROUT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
    STOP
  END SUBROUTINE ReadI


  ! ======================================================================



  SUBROUTINE ReadC(line,cval,tag)

    use mod_params
    !use mod_slcom
    IMPLICIT none

    !     INCLUDE 'params'
    !     INCLUDE 'slcom'

    CHARACTER line*(*),tag*(*),cval*(*)
    INTEGER fp,ival,imin,imax

    INTEGER i
    CHARACTER comment*72

    integer erout1

    erout1 = 0

    READ (line,*,ERR=98,END=98) comment,cval

    WRITE(SLOUT,'(A)')        line
    WRITE(SLOUT,'(5X,2A,A)') tag,' = ',cval

    RETURN
98  WRITE(EROUT1,*) 'Problem reading unstructured input'
99  WRITE(EROUT1,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT1,'(5X,2A)')    'TAG  = ''',tag,''''
    WRITE(EROUT1,'(5X,2A)')    'CVAL = ''',cval,''''
    STOP 'READC'
  END SUBROUTINE ReadC

  ! ======================================================================



  SUBROUTINE Read2I(line,ival1,ival2,imin,imax,tag)

    use mod_params
    !use mod_slcom
    IMPLICIT none

    !     INCLUDE 'params'
    !     INCLUDE 'slcom'

    CHARACTER line*(*),tag*(*)
    INTEGER fp,ival1,ival2,imin,imax

    INTEGER i1,i2
    CHARACTER comment*72

    READ (line,*,ERR=98,END=98) comment,i1,i2

    IF (i1.LT.imin.OR.i1.GT.imax.OR.i2.LT.imin.OR.i2.GT.imax) CALL ER('Read2I','Out of bounds: '//line,*99)

    ival1 = i1
    ival2 = i2

    WRITE(SLOUT,'(A)')        line
    WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival1
    WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival2

    RETURN
98  WRITE(EROUT,*) 'Problem reading unstructured input'
99  WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
    STOP
  END SUBROUTINE Read2I

  ! ======================================================================



  SUBROUTINE ReadR(line,rval,rmin,rmax,tag)

    use mod_params
    !use mod_slcom
    IMPLICIT none

    CHARACTER line*(*),tag*(*)
    REAL rval,rmin,rmax

    !     INCLUDE 'params'
    !     INCLUDE 'slcom'

    REAL r
    CHARACTER comment*72

    READ (line,*,ERR=98,END=98) comment,r

    IF (r.LT.rmin.OR.r.GT.rmax) CALL ER('ReadR','Out of bounds: '//line,*99)

    rval = r

    WRITE(SLOUT,'(A)')        line
    WRITE(SLOUT,'(2A,G10.3)') tag,' = ',rval

    RETURN
98  WRITE(EROUT,*) 'Problem reading unstructured input'
99  WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
    WRITE(EROUT,'(5X,A,3G10.3)') 'R,RVAL,RMIN,RMAX = ',r,rval,rmin,rmax
    STOP
  END SUBROUTINE ReadR


  ! ======================================================================



  SUBROUTINE Read2R(line,rval1,rval2,rmin,rmax,tag)

    use mod_params
    !use mod_slcom
    IMPLICIT none

    CHARACTER line*(*),tag*(*)
    REAL rval1,rval2,rmin,rmax

    !     INCLUDE 'params'
    !     INCLUDE 'slcom'

    REAL r1,r2
    CHARACTER comment*72

    READ (line,*,ERR=98,END=98) comment,r1,r2

    IF (r1.LT.rmin.OR.r1.GT.rmax.OR.r2.LT.rmin.OR.r2.GT.rmax) CALL ER('ReadR','Out of bounds: '//line,*99)

    rval1 = r1
    rval2 = r2

    WRITE(SLOUT,'(A)')        line
    WRITE(SLOUT,'(2A,2G10.3)') tag,' = ',rval1,rval2

    RETURN
98  WRITE(EROUT,*) 'Problem reading unstructured input'
99  WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
    WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
    WRITE(EROUT,'(5X,A,6G10.3)') 'R,RVAL,RMIN,RMAX = ',r1,r2,rval1,rval2,rmin,rmax
    STOP
  END SUBROUTINE Read2R



  
  ! ======================================================
  !
  ! Allocatable inputs
  !
  ! ======================================================


  
        subroutine rdrarn_alloc(rs,nrs,rmin,rmax,ascend,fmin,fmax,nfs,name,ierr)
          use mod_reader
          use mod_io_units
          use allocate_arrays
          implicit none
          integer :: nrs,nfs,ierr
          logical :: ascend
          real :: rmin,rmax,fmin,fmax
          character*(*) :: name
          real,allocatable,intent(out) :: rs(:,:)

!
!  *********************************************************************
!  *                                                                   *
!  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
!  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
!  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
!  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
!  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
!  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
!  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
!  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).                             *
!  *                                                                   *
!  *      PARAMETERS -                                                 *
!  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
!  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
!  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
!  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
!  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
!  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
!  *                                                                   *
!  *********************************************************************
!
      !include 'reader'

      character COMENT*72,MESAGE*72
      character*10 :: rform,sform
      integer ::  IR,N,I,testmax

      integer :: bufflen,nform

      real    ::  RLAST,R
      real,allocatable :: f(:)
      
      !write(0,*) 'Called RDRARN_ALLOC:',nrs,nfs,allocated(rs)

!
!---- READ IN AND TEST SIZE OF ARRAY
!
      NRS = 0
      IERR = 0
      ! turn off maximum array size checking since it is allocated
      testmax = 100000

      call RDC (COMENT, NAME, IERR)
      ! jdemod - allow a value of 0 to be specified in which case the array
      !          won't be allocated
      call RDI (N, .true., 0 ,.false., testmax, NAME, IERR)

      !write(0,*) 'RDRARN_ALLOC:',n

      if (N.eq.0) return

      !
      ! Data needs to be read - so allocate temp array - allocate_array handles deallocation if needed
      !
      call allocate_array(f,nfs+1,'rdrarn_alloc - local function',ierr)

      !
      ! Allocate storage array
      !

      call allocate_array(rs,n,nfs+1,'rdrarn_alloc - allocate input variable',ierr)

      !
      ! Assign arrays to read data of appropriate buffer size
      !

      bufflen = len(buffer)
      nform=int(log10(real(bufflen))+1)
      write(sform,'(A,I1,A)') '(A,I',nform,'A)'
      write(rform,sform) '(A',bufflen,')'
      

!
!---- READ IN AND TEST EACH SET OF ARRAY VALUES
!
      RLAST = RMIN
      IR = 1

   50 continue

      MESAGE = 'END OF FILE ON UNIT 5'

  100 if (IBUF.eq.0) read (stdin,rform,ERR=9999,end=9999) BUFFER


      ! write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'


      if (BUFFER(1:1).eq.'$') goto 100
      if (BUFFER(2:2).eq.'*'.or.BUFFER(2:2).eq.'{') then
        call ReadUnstructuredInput(BUFFER)
        goto 100
      endif

      write (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'

      IBUF = 0
      read (BUFFER,*,ERR=9999,end=9999) R,(F(I),I=1,NFS)

      if (ASCEND) then
        write (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        if (R.lt.RLAST) goto 9999
      endif

      write (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      if (R.gt.RMAX) goto 9999

      RS(IR,1) = R

      do I = 1, NFS
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
        if (F(I).lt.FMIN) goto 9999
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
        if (F(I).gt.FMAX) goto 9999
        RS(IR,1+I) = F(I)
      end do

!
!---- SET UP NUMBER OF VALID VALUES READ
!
      IR = IR + 1
      if (IR.le.N) goto 50
      NRS = IR - 1

      ! deallocate local function 
      if (allocated(f)) deallocate(f)
      
      return

 9999 IERR = 1
      write (7,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      write (0,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      ! clear storage used if error condition encountered
      if (allocated(f)) deallocate(f)
      if (allocated(rs)) deallocate(rs)

      return

    end subroutine rdrarn_alloc


        subroutine rdrarn_alloc_check(rs,nrs,maxnrs,rmin,rmax,ascend,fmin,fmax,nfs,name,ierr)
          use mod_reader
          use mod_io_units
          use allocate_arrays
          implicit none
          integer :: nrs,nfs,ierr,maxnrs
          logical :: ascend
          real :: rmin,rmax,fmin,fmax
          character*(*) :: name
          real,allocatable,intent(out) :: rs(:,:)

!
!  *********************************************************************
!  *                                                                   *
!  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
!  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
!  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
!  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
!  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
!  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
!  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
!  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).                             *
!  *                                                                   *
!  *      PARAMETERS -                                                 *
!  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
!  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
!  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
!  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
!  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
!  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
!  *                                                                   *
!  *********************************************************************
!
      !include 'reader'

      character COMENT*72,MESAGE*72
      character*10 :: rform,sform
      integer ::  IR,N,I

      integer :: bufflen,nform

      real    ::  RLAST,R
      real,allocatable :: f(:)
      
      !write(0,*) 'Called RDRARN_ALLOC:',nrs,nfs,allocated(rs)

!
!---- READ IN AND TEST SIZE OF ARRAY
!
      NRS = 0
      IERR = 0

      call RDC (COMENT, NAME, IERR)
      ! jdemod - allow a value of 0 to be specified in which case the array
      !          won't be allocated
      call RDI (N, .true., 0 ,.true., maxnrs, NAME, IERR)
      if (N.eq.0) return

      !
      ! Data needs to be read - so allocate temp array - allocate_array handles deallocation if needed
      !
      call allocate_array(f,nfs+1,'rdrarn_alloc - local function',ierr)
      !
      ! Allocate output array
      !
      call allocate_array(rs,n,nfs+1,'rdrarn_alloc - allocate input variable',ierr)

      !
      ! Assign arrays to read data of appropriate buffer size
      ! This calculates a format string that creates a format '(a<len>)' where len is the length of the buffer
      !
      bufflen = len(buffer)
      nform=int(log10(real(bufflen))+1)
      write(sform,'(A,I1,A)') '(A,I',nform,'A)'
      write(rform,sform) '(A',bufflen,')'
      

!
!---- READ IN AND TEST EACH SET OF ARRAY VALUES
!
      RLAST = RMIN
      IR = 1

   50 continue

      MESAGE = 'END OF FILE ON UNIT 5'

  100 if (IBUF.eq.0) read (stdin,rform,ERR=9999,end=9999) BUFFER


      ! write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'

      if (BUFFER(1:1).eq.'$') goto 100
      if (BUFFER(2:2).eq.'*'.or.BUFFER(2:2).eq.'{') then
        call ReadUnstructuredInput(BUFFER)
        goto 100
      endif

      write (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'

      IBUF = 0
      read (BUFFER,*,ERR=9999,end=9999) R,(F(I),I=1,NFS)

      if (ASCEND) then
        write (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        if (R.lt.RLAST) goto 9999
      endif

      write (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      if (R.gt.RMAX) goto 9999

      RS(IR,1) = R

      do I = 1, NFS
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
        if (F(I).lt.FMIN) goto 9999
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
        if (F(I).gt.FMAX) goto 9999
        RS(IR,1+I) = F(I)
      end do

!
!---- SET UP NUMBER OF VALID VALUES READ
!
      IR = IR + 1
      if (IR.le.N) goto 50
      NRS = IR - 1

      ! deallocate local function 
      if (allocated(f)) deallocate(f)
      
      return

 9999 IERR = 1
      write (7,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      write (0,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      ! clear storage used if error condition encountered
      if (allocated(f)) deallocate(f)
      if (allocated(rs)) deallocate(rs)

      return

    end subroutine rdrarn_alloc_check


    
    subroutine rdiarn_alloc(rs,nrs,rmin,rmax,ascend,fmin,fmax,nfs,name,ierr)
          use mod_reader
          use mod_io_units
          use allocate_arrays
          implicit none
          integer :: nrs,nfs,ierr
          logical :: ascend
          real :: rmin,rmax,fmin,fmax
          character*(*) :: name
          integer,allocatable :: rs(:,:)

!
!  *********************************************************************
!  *                                                                   *
!  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
!  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
!  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
!  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
!  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
!  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
!  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
!  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).                             *
!  *                                                                   *
!  *      PARAMETERS -                                                 *
!  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
!  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
!  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
!  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
!  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
!  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
!  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
!  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
!  *                                                                   *
!  *********************************************************************
!
      !include 'reader'

      character COMENT*72,MESAGE*72
      integer ::  IR,N,I,testmax

      character*10 :: rform,sform
      integer :: bufflen,nform

      real    ::  RLAST,R
      real,allocatable :: f(:)
      
      !write(0,*) 'Called RDRARN_ALLOC:',nrs,nfs,allocated(rs)

!
!---- READ IN AND TEST SIZE OF ARRAY
!
      NRS = 0
      IERR = 0
      ! turn off maximum array size checking since it is allocated
      testmax = 100000

      ! calls to rdc not needed since header line has the unstructured input tag
      call RDC (COMENT, NAME, IERR)
      ! jdemod - allow a value of 0 to be specified in which case the array
      !          won't be allocated
      call RDI (N, .true., 0 ,.false., testmax, NAME, IERR)

      !write(0,*) 'RDRARN_ALLOC:',n

      if (N.eq.0) return

      !
      ! Data needs to be read - so allocate temp array - allocate_array handles deallocation if needed
      !
      call allocate_array(f,nfs+1,'rdrarn_alloc - local function',ierr)

      !
      ! Allocate storage array
      !

      call allocate_array(rs,n,nfs+1,'rdrarn_alloc - allocate input variable',ierr)

      !
      ! Assign arrays to read data of appropriate buffer size
      !

      bufflen = len(buffer)
      nform=int(log10(real(bufflen))+1)
      write(sform,'(A,I1,A)') '(A,I',nform,'A)'
      write(rform,sform) '(A',bufflen,')'
      

!
!---- READ IN AND TEST EACH SET OF ARRAY VALUES
!
      RLAST = RMIN
      IR = 1

   50 continue

      MESAGE = 'END OF FILE ON UNIT 5'

  100 if (IBUF.eq.0) read (stdin,rform,ERR=9999,end=9999) BUFFER


      ! write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'


      if (BUFFER(1:1).eq.'$') goto 100
      if (BUFFER(2:2).eq.'*'.or.BUFFER(2:2).eq.'{') then
        call ReadUnstructuredInput(BUFFER)
        goto 100
      endif

      write (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'

      IBUF = 0
      read (BUFFER,*,ERR=9999,end=9999) R,(F(I),I=1,NFS)

      if (ASCEND) then
        write (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        if (R.lt.RLAST) goto 9999
      endif

      write (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      if (R.gt.RMAX) goto 9999

      RS(IR,1) = R

      do I = 1, NFS
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
        if (F(I).lt.FMIN) goto 9999
        write (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
        if (F(I).gt.FMAX) goto 9999
        RS(IR,1+I) = F(I)
      end do

!
!---- SET UP NUMBER OF VALID VALUES READ
!
      IR = IR + 1
      if (IR.le.N) goto 50
      NRS = IR - 1

      ! deallocate local function 
      if (allocated(f)) deallocate(f)
      
      return

 9999 IERR = 1
      write (0,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      write (7,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      ! clear storage used if error condition encountered
      if (allocated(f)) deallocate(f)
      if (allocated(rs)) deallocate(rs)

      return

    end subroutine rdiarn_alloc



    
    subroutine RDRAR_alloc(RS,NRS, MAXNRS, RMIN, RMAX, ASCEND, NAME, IERR)
      use error_handling
      use mod_reader
      use mod_io_units
      use allocate_arrays
      implicit none
      integer   NRS, MAXNRS
      real      RMIN, RMAX
      character NAME*(*)
      integer   IERR
      logical   ASCEND
      real, allocatable :: rs(:)
      
!  *********************************************************************
!  *                                                                   *
!  *  RDRAR:   ROUTINE READS A REAL VARIABLE LENGTH 1D ARRAY           *
!  *    AND CHECKS THAT IT IS SORTED IN ASCENDING ORDER,               *
!  *    THAT ALL THE VALUES ARE WITHIN A RANGE AND THAT NO             *
!  *    TWO VALUES ARE EQUAL.                                          *
!  *    ANY OF THESE RULES BEING BROKEN RESULTS IN THE                 *
!  *    OFFENDING VALUE BEING REMOVED FROM THE ARRAY, AN               *
!  *    ERROR MESSAGE BEING OUTPUT AND AN ERROR FLAG BEING             *
!  *    SET.                                                           *
!  *                                                                   *
!  *      PARAMETERS -                                                 *
!  *  RS     : ARRAY VALUES RETURNED IN THIS                           *
!  *  NRS    : NUMBER OF VALUES READ RETURNED IN THIS                  *
!  *  MAXNRS : MAXIMUM NUMBER OF VALUES TO BE READ                     *
!  *  RMIN   : MINIMUM ALLOWED VALUE (EXCLUSIVE)                       *
!  *  RMAX   : MAXIMUM ALLOWED VALUE (EXCLUSIVE)                       *
!  *  ASCEND : INDICATES WHETHER ASCENDING ORDER CHECK IS TO BE MADE   *
!  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
!  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
!  *                                                                   *
!  *    CHRIS FARRELL   JAN 1988                                       *
!  *                                                                   *
!  *********************************************************************
!
!     INCLUDE   "READER"
!      include 'reader'
      character COMENT*72,MESAGE*72
      integer   IR,N
      real      RLAST,R
!
!---- READ IN AND TEST SIZE OF ARRAY
!
      NRS = 0
      call RDC (COMENT, NAME, IERR)
      ! turn off maximum array size checking since it is allocated
      call RDI (N, .true., 0 ,.true., MAXNRS, NAME, IERR)
      if (N.eq.0) return

      !
      ! allocate rs
      !
      call allocate_array(rs,n,'rdrar_alloc input array',ierr)

!
!---- READ IN AND TEST EACH ARRAY VALUE
!
      RLAST = RMIN
      IR = 1
   50 continue
      MESAGE = 'END OF FILE ON UNIT 5'
!     Feb/2008 - jde - changed all buffer reads to A256 from A72
!                    - added buff_format to common to make buffer size changes easy
  100 if (IBUF.eq.0) read(STDIN,buff_format,ERR=9999,end=9999) BUFFER
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRAR_ALLOC'
      if (BUFFER(1:1).eq.'$') goto 100
! slmod begin
      if (BUFFER(2:2).eq.'*'.or.BUFFER(2:2).eq.'{') then
        call ReadUnstructuredInput(BUFFER)
        goto 100
      endif
! slmod end
!
      write (MESAGE,'(A,I5,A)') 'EXPECTING',N,' REALS, ONE PER LINE'
      IBUF = 0
      read (BUFFER,*,ERR=9999,end=9999) R
!
      if (ASCEND) then
        write (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
        if (R.lt.RLAST) goto 9999
      endif
!
      write (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
      if (R.gt.RMAX) goto 9999
      RS(IR) = R
      IR = IR + 1
      if (IR.le.N) goto 50
!
!---- SET UP NUMBER OF VALID VALUES READ
!
      NRS = IR - 1
      return

 9999 continue

      if (ierr.eq.0) then 

         IERR = 1

         call errmsg('RDRAR-READ ERROR',name//mesage)
         call errmsg('RDRAR-LAST LINE ',trim(buffer))

         !if (allocated(rs)) deallocate(rs)
      else

         call dbgmsg('RDRAR-READ ERROR',name//mesage)
         call dbgmsg('RDRAR-LAST LINE ',trim(buffer))

      endif


!      WRITE (7,'(1X,2A,3(/1X,A))')
!     > 'RDRAR: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER

      return
    end subroutine RDRAR_ALLOC



    ! ===================================
    !
    ! Other input routines
    !
    ! ===================================

      SUBROUTINE RDG1 (GRAPH,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,ISELD,IERR)
!      use mod_io_units
!      use mod_io
      implicit none
      INTEGER   ISELE,ISELR,ISELX,ISELD,IERR,ADASYR
      CHARACTER GRAPH*(*), ADASID*(*),ADASEX*(*)

!  *********************************************************************
!  *                                                                   *
!  *  RDG1 : READ IN SELECTOR SWITCHES FOR ADAS PLRP CALCULATIONS      *
!  *                                                                   *
!  *********************************************************************

!     INCLUDE   "READER"
!     include 'reader'
      CHARACTER MESAGE*72

      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG1'
      IF (BUFFER(1:1).EQ.'$') GOTO 100

!     Feature Only useful in OUT

!     jdemod - Added so that global plot modifiers could be read from
!              anywhere. 
!     - not needed when plotting line profile data in DIVIMP
!     
!      IF (BUFFER(2:2).EQ.'#') THEN
!        CALL Read_AdditionalPlotData(BUFFER)
!        GOTO 100
!      ENDIF

!      write(0,'(a,8i5)')
!     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx

      MESAGE = 'EXPECTING 2 CHAR, 1 INT, 1 CHAR  AND 4 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,ISELD

!      write(0,'(a,8i5)')
!     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx

!      write(0,'(3a)')
!     >  'RDG1:',buffer,':'
!      write(0,'(3a)')
!     >  'RDG1:',graph,':'
!      write(0,'(3a)')
!     >  'RDG1:',adasid,':'
!      write(0,'(3a)')
!     >  'RDG1:',adasex,':'


      RETURN

 9998 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))') 'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      WRITE (7,'(1X,A,4(/1X,A))') 'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      RETURN

 9999 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))') 'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      WRITE (7,'(1X,A,4(/1X,A))') 'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      RETURN
    END SUBROUTINE RDG1




      SUBROUTINE RD_lp_los (GRAPH,lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_instrument_width,lp_bin_width,ierr)
!      use mod_io_units
!      use mod_reader
        implicit none
      INTEGER   IERR
      real lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_instrument_width,lp_bin_width
      CHARACTER GRAPH*(*)

!  *********************************************************************
!  *                                                                   *
!  *  RD_LP_LOS : LOS DEFINITION FOR LINE PROFILE CALCULATION          *
!  *                                                                   *
!  *********************************************************************

      CHARACTER MESAGE*72

      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (stdin,buff_format,ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RD_LP'
      IF (BUFFER(1:1).EQ.'$') GOTO 100

!     jdemod - Added so that global plot modifiers could be read from
!              anywhere. 
!     - not needed when calculating line profile data in DIVIMP
!     
!      IF (BUFFER(2:2).EQ.'#') THEN
!        CALL Read_AdditionalPlotData(BUFFER)
!        GOTO 100
!      ENDIF

      MESAGE = 'EXPECTING 1 CHAR, 6 REALS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_instrument_width,lp_bin_width

      RETURN

 9998 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))') 'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      WRITE (7,'(1X,A,4(/1X,A))') 'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      RETURN

 9999 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))') 'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      WRITE (7,'(1X,A,4(/1X,A))') 'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',trim(BUFFER)
      RETURN
    END SUBROUTINE RD_lp_los



    
end module mod_io
