module allocatable_input_data
  implicit none

  public

  !
  ! The purpose of this module is to accumulate declarations of 
  ! allocatable data as well as routines to read and manipulate it
  !



  ! Unstructured Input data
  !
  ! TAG Q44 - core plasma profiles
  !
  integer :: ncoreprofile
  real,allocatable :: coreprofile(:,:)



contains



  !        subroutine rdrarn_alloc(rs,nrs,rmin,rmax,ascend,fmin,fmax,nfs,name,ierr)
  !          use mod_reader
  !          implicit none
  !          integer :: nrs,nfs,ierr
  !          logical :: ascend
  !          real :: rmin,rmax,fmin,fmax
  !          character*(*) :: name
  !          real,allocatable :: rs(:,:)
  !
  !!
  !!  *********************************************************************
  !!  *                                                                   *
  !!  *  RDRARN:  ROUTINE READS IN A SET OF VALUES (X, F1(X),F2(X),...)
  !!  *    WHERE X REPRESENTS AN X POSITION AND FI(X) A FUNCTION VALUE    *
  !!  *    AT THAT X POSITION.  THE X VALUES MUST BE IN ASCENDING ORDER,  *
  !!  *    WITH NO TWO VALUES BEING EQUAL, AND MUST LIE WITHIN THE GIVEN  *
  !!  *    RANGE RMIN TO RMAX.  THE FUNCTION VALUES MUST LIE WITHIN THE   *
  !!  *    RANGE FMIN TO FMAX.  THE QUANTITY NFS GIVES THE NUMBER OF      *
  !!  *    FUNCTIONS GIVEN ALONGSIDE THE X POSITIONS  (EG NFS=1, JUST     *
  !!  *    ONE FUNCTION, NFS=2, 2 FUNCTIONS).                             *
  !!  *                                                                   *
  !!  *      PARAMETERS -                                                 *
  !!  *  RS     : ARRAY: VALUES RETURNED IN THIS                          *
  !!  *  NRS    : NUMBER OF SETS OF VALUES READ RETURNED IN THIS          *
  !!  *  RMIN   : MINIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
  !!  *  RMAX   : MAXIMUM ALLOWED X VALUE (EXCLUSIVE)                     *
  !!  *  ASCEND : INDICATES WHETHER X VALS CHECKED FOR ASCENDING ORDER    *
  !!  *  NFS    : NUMBER OF FUNCTION VALUES TO BE READ AFTER EACH X VALUE *
  !!  *  FMIN   : MINIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
  !!  *  FMAX   : MAXIMUM ALLOWED F VALUE (EXCLUSIVE)                     *
  !!  *  NAME   : NAME OF ARRAY (FOR ERROR MESSAGES)                      *
  !!  *  IERR   : SET TO 1 IF AN ERROR FOUND                              *
  !!  *                                                                   *
  !!  *********************************************************************
  !!
  !      !include 'reader'
  !
  !      CHARACTER COMENT*72,MESAGE*72
  !      character*10 :: rform,sform
  !      INTEGER ::  IR,N,I,testmax
  !
  !      integer :: bufflen,nform
  !
  !      REAL    ::  RLAST,R
  !      real,allocatable :: f(:)
  !      
  !      !write(0,*) 'Called RDRARN_ALLOC:',nrs,nfs,allocated(rs)
  !
  !!
  !!---- READ IN AND TEST SIZE OF ARRAY
  !!
  !      NRS = 0
  !      IERR = 0
  !      ! turn off maximum array size checking since it is allocated
  !      testmax = 100000
  !
  !      CALL RDC (COMENT, NAME, IERR)
  !      ! jdemod - allow a value of 0 to be specified in which case the array
  !      !          won't be allocated
  !      CALL RDI (N, .TRUE., 0 ,.FALSE., testmax, NAME, IERR)
  !
  !      !write(0,*) 'RDRARN_ALLOC:',n
  !
  !      IF (N.EQ.0) RETURN
  !
  !      !
  !      ! Data needs to be read - so allocate temp array
  !      !
  !      allocate(f(nfs+1))
  !
  !      !
  !      ! Allocate storage array
  !      !
  !
  !      allocate(rs(n,nfs+1))
  !
  !      !
  !      ! Assign arrays to read data of appropriate buffer size
  !      !
  !
  !      bufflen = len(buffer)
  !      nform=int(log10(real(bufflen))+1)
  !      write(sform,'(A,I1,A)') '(A,I',nform,'A)'
  !      write(rform,sform) '(A',bufflen,')'
  !      
  !
  !!
  !!---- READ IN AND TEST EACH SET OF ARRAY VALUES
  !!
  !      RLAST = RMIN
  !      IR = 1
  !
  !   50 CONTINUE
  !
  !      MESAGE = 'END OF FILE ON UNIT 5'
  !
  !  100 IF (IBUF.EQ.0) READ (5,rform,ERR=9999,END=9999) BUFFER
  !
  !
  !      ! write(9,'(1x,a20,a1,a,1x,a6)') name,':',trim(buffer),'RDRARN'
  !
  !
  !      IF (BUFFER(1:1).EQ.'$') GOTO 100
  !      IF (BUFFER(2:2).EQ.'*'.OR.BUFFER(2:2).EQ.'{') THEN
  !        CALL ReadUnstructuredInput(BUFFER)
  !        GOTO 100
  !      ENDIF
  !
  !      WRITE (MESAGE,'(A,I5,A,I2,A)') 'EXPECTING',(NFS+1)*N,' REALS,',NFS+1,' PER LINE'
  !
  !      IBUF = 0
  !      READ (BUFFER,*,ERR=9999,END=9999) R,(F(I),I=1,NFS)
  !
  !      IF (ASCEND) THEN
  !        WRITE (MESAGE,'(G11.4,A,G11.4)') R,' LESS THAN PREV/MIN',RLAST
  !        IF (R.LT.RLAST) GOTO 9999
  !      ENDIF
  !
  !      WRITE (MESAGE,'(G11.4,A,G11.4)') R,' MORE THAN MAXIMUM',RMAX
  !      IF (R.GT.RMAX) GOTO 9999
  !
  !      RS(IR,1) = R
  !
  !      DO I = 1, NFS
  !        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' LESS THAN MINIMUM',FMIN
  !        IF (F(I).LT.FMIN) GOTO 9999
  !        WRITE (MESAGE,'(G11.4,A,G11.4)') F(I),' MORE THAN MAXIMUM',FMAX
  !        IF (F(I).GT.FMAX) GOTO 9999
  !        RS(IR,1+I) = F(I)
  !      end do
  !
  !!
  !!---- SET UP NUMBER OF VALID VALUES READ
  !!
  !      IR = IR + 1
  !      IF (IR.LE.N) GOTO 50
  !      NRS = IR - 1
  !      RETURN
  !
  ! 9999 IERR = 1
  !      WRITE (7,'(1X,2A,3(/1X,A))') 'RDRARN_ALLOC: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',trim(BUFFER)
  !      if (allocated(rs)) deallocate(rs)
  !
  !      RETURN
  !
  !    END subroutine rdrarn_alloc


end module allocatable_input_data
