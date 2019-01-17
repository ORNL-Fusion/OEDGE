module mod_cadas2

  use mod_params


!C                                                                       
!      COMMON /CADAS2/ DTEV, DDENS, DTEVD, DDENSD, DRCOFD, ZDATA,        
!     >                DRCOFI, TITLF,                                    
!     >                IFAIL, IEVCUT, ITMAXD, IDMAXD, IZMAXD             
!C                                                                       
!      REAL*8 DTEV(MAXNXS), DDENS(MAXNXS), DTEVD(MAXADS),                
!     >       DDENSD(MAXADS), DRCOFD(MAXADS,MAXADS,MAXADS),              
!     >       ZDATA(MAXADS), DRCOFI(MAXNXS)                              
!      CHARACTER TITLF*80                                                
!      INTEGER IFAIL, IEVCUT, ITMAXD, IDMAXD, IZMAXD                     

  implicit none
  private

      REAL*8,allocatable,public:: DTEV(:), DDENS(:), DTEVD(:),&                
            DDENSD(:), DRCOFD(:,:,:),&              
            ZDATA(:), DRCOFI(:)                              
      CHARACTER,public:: TITLF*80                                                
      INTEGER,public:: IFAIL, IEVCUT, ITMAXD, IDMAXD, IZMAXD                     


  public :: allocate_mod_cadas2, deallocate_mod_cadas2


contains

  subroutine allocate_mod_cadas2
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)
    call allocate_array(DDENS ,maxnxs,'DDENS',ierr)
    call allocate_array(DTEVD ,maxads,'DTEVD',ierr)
    call allocate_array(DDENSD,maxads,'DDENSD',ierr)
    call allocate_array(DRCOFD,maxads,maxads,maxads,'DRCOFD',ierr)
    call allocate_array(DRCOFI,maxnxs,'DRCOFI',ierr)
    call allocate_array(ZDATA ,maxads,'ZDATA',ierr)

  end subroutine allocate_mod_cadas2


  subroutine deallocate_mod_cadas2
    use mod_params
    use allocate_arrays
    implicit none

    !write(0,*) '2A'
    deallocate(DTEV  )
    !write(0,*) '2B'
    deallocate(DDENS )
    !write(0,*) '2C'
    deallocate(DTEVD )
    !write(0,*) '2D'
    deallocate(DDENSD)
    !write(0,*) '2E'
    deallocate(DRCOFD)
    !write(0,*) '2F'
    deallocate(DRCOFI)
    !write(0,*) '2G'
    deallocate(ZDATA )
    !write(0,*) '2H'

  end subroutine deallocate_mod_cadas2



end module mod_cadas2
