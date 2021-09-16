!     -*-Fortran-*-
!
! add a special section for simple stores, like single entries 
!    {HEADER INFO FOR ALL ENTRIES...} (can be made robust?)
!    {DATA ITEM} {_TAG_} units data_type value 
! text header explaining file including comment characters
! tags repeated every 50 lines or so
! clean up ugly repetative buffer loads, use buffer array
! automatic detection of how large the array is that is being passed
! more intelligent initial allocation and expansion of allocation
! remove screen dumping, or hide behind 'debug' logical
! find out why the density profile looks funny with u-lin-0002a 
! change routine names to inXX
!
!
!     CALL OpenInterface...
!
!     DO iz = 1, cion
!       DO ir = 2, nrs
!         IF (idring(ir).EQ.BOUNDARY) CYCLE
!         CALL PutData(ddlims(1,ir,iz),nks(ir),'ADAMS_DATA1','s-1 m-3')
!         CALL PutData(vel   (1,ir,iz),nks(ir),'ADAMS_DATA2','s-1 m-3')
!         CALL PutData(temp  (1,ir,iz),nks(ir),'ADAMS_DATA3','s-1 m-3')
!       ENDDO
!     ENDDP
!

!     ==================================================================
!     ==================================================================
!
      MODULE MOD_INTERFACE_TAGS
      IMPLICIT none
      PUBLIC

!...  Data types: 
      INTEGER, PARAMETER, PUBLIC :: DTY_B = 1,     ! byte
     .                              DTY_I = 2,     ! integer (4 byte)
     .                              DTY_R = 3,     ! single precision real (4 byte)
     .                              DTY_D = 4      ! double precision real (8 byte)    


!...  Define some preset data and unit strings:   *** REMOVE THESE? ***

      INTEGER, PARAMETER, PUBLIC :: E_T = 001,    ! Electron temperature
     .                              E_N = 002     ! Electron density

      INTEGER, PARAMETER, PUBLIC :: F_T_1 = 003,  !
     .                              F_N_1 = 004   !

      INTEGER, PARAMETER, PUBLIC :: N_T_1 = 101,  ! 
     .                              N_N_1 = 102   !  

      INTEGER, PARAMETER, PUBLIC :: I_T_1 = 201,  ! CII temperature
     .                              I_N_1 = 202   !     density
  

      END MODULE MOD_INTERFACE_TAGS
!
!     ==================================================================
!     ==================================================================
!
      MODULE MOD_INTERFACE
      USE mod_interface_tags
      IMPLICIT none
      PRIVATE


      INTERFACE inPutData

        MODULE PROCEDURE PutDataI ,PutDataR ,PutDataD,
     .                   PutDatumI,PutDatumR,PutDatumD

      END INTERFACE

      INTERFACE inGetData

        MODULE PROCEDURE GetDataI ,GetDataR ,GetDataD,
     .                   GetDatumI,GetDatumR,GetDatumD

      END INTERFACE


!...  Definitions:
!     ==================================================================

!
!     ------------------------------------------------------------------
      TYPE :: type_interface
        REAL           :: version
        CHARACTER(512) :: file_name
        INTEGER        :: file_pointer
        LOGICAL        :: file_open
        LOGICAL        :: file_stream
      ENDTYPE type_interface
!
!     ------------------------------------------------------------------
      TYPE :: type_handle
        REAL           :: version
        CHARACTER(512) :: id
        CHARACTER(512) :: date_time
        CHARACTER(128) :: device
        CHARACTER(128) :: code
        CHARACTER(256) :: case
        INTEGER        :: shot
        REAL           :: time
      ENDTYPE type_handle
!
!     ------------------------------------------------------------------
      TYPE :: type_data
!       Descriptors:
        TYPE(type_handle) :: handle
!       Data:
        REAL           :: version
        CHARACTER(256) :: tag
        CHARACTER(256) :: units
        INTEGER        :: type   ! Variable type of data
        INTEGER        :: n      ! Number of array elements stored
        INTEGER        :: nmax   ! Maximum number of array elements that can be stored... (not sure this is the way to go)
        CHARACTER(1), ALLOCATABLE :: bdata(:)
        INTEGER*4,    ALLOCATABLE :: idata(:)
        REAL*4,       ALLOCATABLE :: rdata(:)
        REAL*8,       ALLOCATABLE :: ddata(:)
      ENDTYPE type_data


      INTEGER :: MAXNDAT = 300,   
     .           MAXNTAG = 100,
     .           MAXNMAX = 100000, ! 100 KB allocated for each data item
     .           MAXJUMP = 100E+6  ! 100 MB is the max jump in memory allocation per item


!..  .Declarations:
!     ==================================================================

      PUBLIC :: inOpenInterface, inCloseInterface, inPutData , inGetData, inGetDataSize !, GetHandle

      TYPE(type_interface) :: interface

      INTEGER :: ndat
      TYPE(type_data), ALLOCATABLE :: dat(:)


!...  Data types: 
      INTEGER, PARAMETER, PUBLIC :: ITF_READ  = 1,  ! open interface for reading data
     .                              ITF_WRITE = 2   ! open interface for writing data

!...  Routines:
!     ==================================================================
!
      CONTAINS
!
!     ---------------------------------------------------------------------
!     Public: 
!     ---------------------------------------------------------------------
!
!...eventually want to specify whether or not the data is to be dumped to a file, and what the unit # is...also some default max. src array size?
      SUBROUTINE inOpenInterface(file_name,io_select)
      IMPLICIT none

      INTEGER  , INTENT(IN) :: io_select
      CHARACTER, INTENT(IN) :: file_name*(*)

!     Blank file name string (just to be sure):
      WRITE(interface%file_name,'(512X)')  

      interface%file_name    = TRIM(file_name)
      interface%file_open    = .FALSE.
      interface%file_pointer = 99

      ndat = 0
      IF (.NOT.ALLOCATED(dat)) ALLOCATE(dat(MAXNDAT))

!...  Just to be sure:
      dat(1:MAXNDAT)%n    = 0
      dat(1:MAXNDAT)%nmax = 0

!...  Return handle:


      RETURN
 99   STOP
      END SUBROUTINE inOpenInterface
!
!     ---------------------------------------------------------------------
!
      SUBROUTINE inCloseInterface
      IMPLICIT none

      IF (.TRUE.) CALL WriteInterfaceData

!     Clean up:
      IF (ALLOCATED(dat)) DEALLOCATE(dat)

!     Close file stream:
      IF (interface%file_open) THEN
        CLOSE(interface%file_pointer)
        interface%file_open = .FALSE.
      ENDIF

      RETURN
 99   STOP
      END SUBROUTINE inCloseInterface
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDatumI(src,dtag,units) 
      IMPLICIT none

!...  Input:
      INTEGER,  INTENT(IN) :: src
      CHARACTER,INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      INTEGER src_single(1)

      src_single(1) = src
      CALL PutDataI(src_single,dtag,units)

      RETURN
 99   STOP
      END SUBROUTINE PutDatumI
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDataI(src,dtag,units) 
      IMPLICIT none

!...  Input:
      INTEGER,   INTENT(IN) :: src(:)
      CHARACTER, INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      INTEGER             :: idat,nsrc,n,nmax
      REAL*4, ALLOCATABLE :: store(:)

      nsrc = SIZE(src)

      CALL FindDataIndex(nsrc,DTY_I,dtag,units,idat)

      n    = dat(idat)%n
      nmax = dat(idat)%nmax

!...  Check if more memory needs to be allocated:

c      write(0,*) ' max  ' ,nsrc,n,nmax,nmax-n

      IF (nsrc.GT.nmax-n) THEN
        IF (.NOT.ALLOCATED(dat(idat)%idata)) THEN
!         Fresh start:
          nmax = MAX(MAXNMAX,nsrc)
          ALLOCATE(dat(idat)%idata(nmax))
        ELSE
!         Store current data but deallocate, then reallocate:
          ALLOCATE(store(n))
          store(1:n) = dat(idat)%idata(1:n)    ! Can get awkward for large meshes... maybe 10% increase each time...?
          DEALLOCATE(dat(idat)%idata)
          nmax = MIN(MAXJUMP,2*nmax) ! n + nsrc  ! Need a feature that dumps everything to disk if the demands get too large...
          ALLOCATE(dat(idat)%idata(nmax))        ! Also need to be able to store data more efficiently, eliminating
          dat(idat)%idata(1:n) = store(1:n)      ! blanks in the data file...
          DEALLOCATE(store)
        ENDIF
      ENDIF

!...  Store data:
c      write(0,*) ' mod_i' ,n,nsrc,idat,' '//TRIM(dtag)
c      write(0,*) ' mod_i' ,src(1:nsrc)
c      write(0,*) '      ' ,SIZE(dat,1)
c      write(0,*) '      ' ,SIZE(dat(1)%idata,1)
c      write(0,*) '      ' ,dat(1)%idata(1)

      dat(idat)%nmax = nmax
      dat(idat)%idata(n+1:n+nsrc) = src(1:nsrc)
      dat(idat)%n = dat(idat)%n + nsrc

      RETURN
 99   STOP
      END SUBROUTINE PutDataI
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDatumR(src,dtag,units) 
      IMPLICIT none

!...  Input:
      REAL,     INTENT(IN) :: src
      CHARACTER,INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      REAL src_single(1)

      src_single(1) = src
      CALL PutDataR(src_single,dtag,units)

      RETURN
 99   STOP
      END SUBROUTINE PutDatumR
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDatumD(src,dtag,units) 
      IMPLICIT none

!...  Input:
      REAL*8,   INTENT(IN) :: src
      CHARACTER,INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      REAL*8 src_single(1)

      src_single(1) = src
      CALL PutDataD(src_single,dtag,units)

      RETURN
 99   STOP
      END SUBROUTINE PutDatumD
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDataR(src,dtag,units)
      IMPLICIT none

!...  Input:
      REAL*4   , INTENT(IN) :: src(:)
      CHARACTER, INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      INTEGER             :: idat,nsrc,n,nmax
      REAL*4, ALLOCATABLE :: store(:)

      nsrc = SIZE(src)
      
      CALL FindDataIndex(nsrc,DTY_R,dtag,units,idat)

      n    = dat(idat)%n
      nmax = dat(idat)%nmax
!...  Check if more memory needs to be allocated:
      IF (nsrc.GT.nmax-n) THEN
        IF (.NOT.ALLOCATED(dat(idat)%rdata)) THEN
!         Fresh start:
          nmax = MAX(MAXNMAX,nsrc)
          ALLOCATE(dat(idat)%rdata(nmax))
        ELSE
!         Store current data but deallocate, then reallocate:
          ALLOCATE(store(n))
          store(1:n) = dat(idat)%rdata(1:n) 
          DEALLOCATE(dat(idat)%rdata)
          nmax = MIN(MAXJUMP,2*nmax) ! n + nsrc
          ALLOCATE(dat(idat)%rdata(nmax))
          dat(idat)%rdata(1:n) = store(1:n)
          DEALLOCATE(store)
        ENDIF
      ENDIF
!...  Store data:
      dat(idat)%nmax = nmax
      dat(idat)%rdata(n+1:n+nsrc) = src(1:nsrc)
      dat(idat)%n = dat(idat)%n + nsrc
      RETURN
 99   STOP
      END SUBROUTINE PutDataR
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE PutDataD(src,dtag,units)
      IMPLICIT none

!...  Input:
      REAL*8   , INTENT(IN) :: src(:)
      CHARACTER, INTENT(IN) :: dtag*(*),units*(*)
!...  Locals:
      INTEGER             :: idat,nsrc,n,nmax
      REAL*8, ALLOCATABLE :: store(:)

      nsrc = SIZE(src)
      
      CALL FindDataIndex(nsrc,DTY_D,dtag,units,idat)

      n    = dat(idat)%n
      nmax = dat(idat)%nmax

!...  Check if more memory needs to be allocated:
      IF (nsrc.GT.nmax-n) THEN
        IF (.NOT.ALLOCATED(dat(idat)%ddata)) THEN
!         Fresh start:
          nmax = MAX(MAXNMAX,nsrc)
          ALLOCATE(dat(idat)%ddata(nmax))
        ELSE
!         Store current data but deallocate, then reallocate:
          ALLOCATE(store(n))
          store(1:n) = dat(idat)%ddata(1:n) 
          DEALLOCATE(dat(idat)%ddata)
          nmax = MIN(MAXJUMP,2*nmax) ! n + nsrc  -- 
          ALLOCATE(dat(idat)%ddata(nmax))
          dat(idat)%ddata(1:n) = store(1:n)
          DEALLOCATE(store)
        ENDIF
      ENDIF

!...  Store data:
      dat(idat)%nmax = nmax
      dat(idat)%ddata(n+1:n+nsrc) = src(1:nsrc)
      dat(idat)%n = dat(idat)%n + nsrc

      RETURN
 99   STOP
      END SUBROUTINE PutDataD
!
!     -----------------------------------------------------------------------
!     Private:
!     -----------------------------------------------------------------------
!
      SUBROUTINE FindDataIndex(nsrc,dtype,dtag,units,idat)
      IMPLICIT none

!...  Input:
      INTEGER   :: nsrc,dtype,idat
      CHARACTER :: dtag*(*),units*(*)
!...  Locals:
      INTEGER   :: i1

!...  Check if data item with specificed DTAG already exists:
      idat = 0
      DO i1 = 1, ndat
        IF (TRIM(dtag).EQ.TRIM(dat(i1)%tag)) THEN
          IF (idat.EQ.0) THEN
            idat = i1
          ELSE
            CALL InterfaceError('PutData','Duplicate CTAG found',*99)
          ENDIF
        ENDIF
      ENDDO

      IF (idat.EQ.0) THEN
!...    Setup data item and allocate some memory:

        IF (ndat+1.GT.MAXNDAT)           !... in future, just expand the array...
     .    CALL InterfaceError('PutData','Storage limit exceeded',*99)

        ndat = ndat + 1
        idat = ndat
        dat(idat)%version = 0.01
        dat(idat)%tag     = TRIM(dtag)
        dat(idat)%units   = TRIM(units)
        dat(idat)%type    = dtype
        dat(idat)%n       = 0

c        WRITE(0,*) 'TAG:',TRIM(dtag)//'<'
      ELSE
!...    Check that everything is consistent with previsous definition: 
      ENDIF


      RETURN
 99   STOP
      END SUBROUTINE FindDataIndex
!
!     ---------------------------------------------------------------------
!
      SUBROUTINE WriteInterfaceData
      IMPLICIT none

      INTEGER       fp,i,j,k,n,cw(MAXNDAT),nmax,pw,count
      CHARACTER*512 fname,buffer(10)

!     Number of columns in output file:
      INTEGER, PARAMETER :: ncol = 4

      fp    = interface%file_pointer
      fname = interface%file_name

      IF (.NOT.interface%file_open) 
     .  OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE')

      interface%file_open = .TRUE.

      SELECTCASE (1)
        CASE (1)
!         ASCII, columns:
          pw = 22

          WRITE(fp,'(F4.2)') 1.00  ! Version

          WRITE(fp,'(A,I6)') '{FILE BLOCK FORMAT}',1 
          WRITE(fp,'(A,I6)') '{FILE INDENT}      ',pw
          WRITE(fp,'(A)')    '*'
          WRITE(fp,'(A,I6)') '{FILE INDEX}       ',ndat
          WRITE(fp,10) '*'

          DO i = 1, ndat, ncol
            WRITE(fp,10) '*'
            WRITE(fp,'(A,I6)') '{FILE COLUMNS}',MIN(i+ncol-1,ndat)-i+1
            nmax = -1
            cw = -1 
c...        Determine the character width of each data column:
            DO j = i, MIN(i+ncol-1,ndat)
              SELECTCASE (dat(j)%type)
c                CASE (DTY_B)
                CASE (DTY_I)
                  cw(j) = MAX(11,LEN_TRIM(dat(j)%tag)+4)  ! ...proper maximum...?
                CASE (DTY_R)
                  cw(j) = MAX(14,LEN_TRIM(dat(j)%tag)+4)
                CASE (DTY_D)
                  cw(j) = MAX(22,LEN_TRIM(dat(j)%tag)+4)
                CASE DEFAULT
                  CALL InterfaceError
     .              ('WriteInterfaceData','Unknown data type',*99)
              ENDSELECT
              cw(j) = MAX(cw(j),LEN_TRIM(dat(j)%units))
              nmax = MAX(nmax,dat(j)%n)
            ENDDO

!...        Collect data header information:
            WRITE(buffer(1:10),'(512X)')              
            WRITE(buffer(1),10) '{FILE COLUMN WIDTH}'
            WRITE(buffer(2),10) '{DATA TAG}'
            WRITE(buffer(3),10) '{DATA UNITS}'         
            WRITE(buffer(4),10) '{DATA TYPE}'         
            WRITE(buffer(5),10) '{DATA N}'         
10          FORMAT(A)
c
c           jdemod - I am putting in Karl's change for now so he can compile on the SUN
c                    machines - however, this code writes an interface file so it may be
c                    very sensitive to file format. Need to get Steve to check since this
c                    is his code. 
c
c           IPP/08 Krieger - SUN f95 compiler rejects format(I) specifiers
c           without explicit width (which is an unofficial f77 extension)
c           Put in 12 for long integers and hope that's ok
11          FORMAT(I11)
c11          FORMAT(I)
c
            n = pw    
            DO j = i, MIN(i+ncol-1,ndat)
c             WRITE(0,*) n,n+cw(j),cw(j)+4
             WRITE(buffer(1)(n+1:n+cw(j)),11) cw(j)+4
             WRITE(buffer(2)(n+1:n+cw(j)),10) '{'//TRIM(dat(j)%tag)//'}'
             WRITE(buffer(3)(n+1:n+cw(j)),10) TRIM(dat(j)%units)
             WRITE(buffer(4)(n+1:n+cw(j)),11) dat(j)%type
             WRITE(buffer(5)(n+1:n+cw(j)),11) dat(j)%n
             n = n + cw(j) + 4
            ENDDO
c...        Write header information:
            WRITE(fp,10) TRIM(buffer(1))
            WRITE(fp,10) '*'
            DO j = 2, 5
              WRITE(fp,10) TRIM(buffer(j))
            ENDDO
c...        List data:
            WRITE(fp,10) '{DATA VALUES}'         
            count = 0
            WRITE(buffer(2)(1:pw),'(A,256X)') '*'
            DO k = 1, nmax
              count = count + 1
              WRITE(buffer(1),'(512X)')              
              WRITE(buffer(1),'(I9)') k
              n = pw    
              DO j = i, MIN(i+ncol-1,ndat)
                IF (dat(j)%n.GE.k) THEN
                  SELECTCASE (dat(j)%type)
c                    CASE (DTY_B)
                    CASE (DTY_I)
                      WRITE(buffer(1)(n+1:n+cw(j)),12) dat(j)%idata(k)
12                    FORMAT(I11)
                    CASE (DTY_R)
                      WRITE(buffer(1)(n+1:n+cw(j)),13) dat(j)%rdata(k)
13                    FORMAT(1P,E14.6,0P)
                    CASE (DTY_D)
                      WRITE(buffer(1)(n+1:n+cw(j)),14) dat(j)%ddata(k)
14                    FORMAT(1P,E22.14,0P)
                    CASE DEFAULT
                      CALL InterfaceError
     .                  ('WriteInterfaceData','Unknown data type',*99)
                  ENDSELECT
                ELSE
                  WRITE(buffer(1)(n+1:n+cw(j)),'(128X)')
                ENDIF
                n = n + cw(j) + 4 
              ENDDO
              IF (MOD(count,50).EQ.0) WRITE(fp,10) TRIM(buffer(2))  ! Write the list of tags every 50 lines, 
              WRITE(fp,10) TRIM(buffer(1))                          ! for easy reference...
            ENDDO

          ENDDO

          WRITE(fp,10) '*'
          WRITE(fp,10) '{FILE END}'
          WRITE(fp,10) '*'

        CASE DEFAULT
          CALL InterfaceError('WriteInterfaceData','Unknown file '//
     .                        'format option',*99)
      ENDSELECT


      RETURN
 99   STOP
      END SUBROUTINE WriteInterfaceData
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE inGetDataSize
      IMPLICIT none
      RETURN
 99   STOP
      END SUBROUTINE inGetDataSize
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDataI(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      INTEGER  , INTENT(OUT) :: src(:)
      INTEGER :: nsrc
      nsrc = SIZE(src)
      src(1:nsrc) = -1
      RETURN
 99   STOP
      END SUBROUTINE GetDataI
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDatumI(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      INTEGER  , INTENT(OUT) :: src
      src = -1
      RETURN
 99   STOP
      END SUBROUTINE GetDatumI
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDataR(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      REAL     , INTENT(OUT) :: src(:)
      INTEGER :: nsrc
      nsrc = SIZE(src)
      src(1:nsrc) = -1.0
      RETURN
 99   STOP
      END SUBROUTINE GetDataR
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDatumR(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      REAL     , INTENT(OUT) :: src
      src = -1.0
      RETURN
 99   STOP
      END SUBROUTINE GetDatumR
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDataD(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      REAL*8   , INTENT(OUT) :: src(:)
      INTEGER :: nsrc
      nsrc = SIZE(src)
      src(1:nsrc) = -1.0D0
      RETURN
 99   STOP
      END SUBROUTINE GetDataD
!
!   -----------------------------------------------------------------------
!
      SUBROUTINE GetDatumD(src,dtag)
      IMPLICIT none
      CHARACTER, INTENT(IN)  :: dtag*(*)
      REAL*8   , INTENT(OUT) :: src
      src = -1.0D0
      RETURN
 99   STOP
      END SUBROUTINE GetDatumD
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
      SUBROUTINE InterfaceError(char1,char2,*)
      IMPLICIT none

      CHARACTER :: char1*(*), char2*(*)

      WRITE(0,'(A)') 'ERROR: '// char1(1:LEN_TRIM(char1))//', '//   
     .                           char2(1:LEN_TRIM(char2))

      RETURN 1
 99   STOP
      END SUBROUTINE InterfaceError

!     ==================================================================
      END MODULE MOD_INTERFACE


