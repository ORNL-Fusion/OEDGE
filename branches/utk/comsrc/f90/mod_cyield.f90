module mod_cyield


!  c     -*-Fortran-*-
!C                                                                       
!      COMMON /CYIELD/ CETH,CETF,CQ,CIDATA,NTARS,flux_frac                               
!      REAL            CETH(7,12),CETF(7,12),CQ(7,12),                   
!     >                flux_frac
!      integer         ntars
!      LOGICAL         CIDATA(7,12)                                      


  implicit none
  private


      REAL,public::   CETH(7,12),CETF(7,12),CQ(7,12),flux_frac
      integer,public::         ntars
      LOGICAL,public::         CIDATA(7,12)                                      

  
  public :: allocate_mod_cyield, deallocate_mod_cyield


contains

  subroutine allocate_mod_cyield
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_cyield


  subroutine deallocate_mod_cyield
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_cyield



end module mod_cyield
