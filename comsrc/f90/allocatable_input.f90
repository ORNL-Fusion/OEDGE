module allocatable_input
      implicit none

      public

      interface allocate_array_input
         module procedure rdrarn_alloc,rdrarn_alloc_check,rdrar_alloc,rdiarn_alloc,rdiar_alloc
      end interface allocate_array_input
      
      ! This module provides the functionality required to read in and allocate an input quantity at the same time

      contains


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


    
    subroutine RDiAR_alloc(RS,NRS, MAXNRS, RMIN, RMAX, ASCEND, NAME, IERR)
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
      integer, allocatable :: rs(:)
      
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
      call RDI (N, .true., 0 ,.false., 0, NAME, IERR)
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
      write(echout,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDIAR_ALLOC'
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
    end subroutine RDiAR_ALLOC



  end module allocatable_input
