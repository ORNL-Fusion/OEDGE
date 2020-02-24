module mod_zommv

  use mod_params

!C                                                                               
!      COMMON /ZOMMV/  ZICABS,ZIFABS,ZILABS,ZISABS,ZRTABS,ZRVABS,ZTBS,           
!     >                ZRAVAV,ZICAY0,ZICAAW,ZICA2L                               
!      REAL            ZICABS(MAXIZS), ZIFABS(MAXIZS), ZILABS(MAXIZS)            
!      REAL            ZISABS(MAXIZS), ZICAY0, ZICAAW, ZICA2L                    
!      REAL            ZRTABS(MAXIZS), ZRVABS(MAXIZS), ZRAVAV(MAXIZS)            
!      REAL            ZTBS(MAXIZS)                                              


  implicit none
  private

      REAL,public:: ZICABS(MAXIZS), ZIFABS(MAXIZS), ZILABS(MAXIZS)            
      REAL,public:: ZISABS(MAXIZS), ZICAY0, ZICAAW, ZICA2L                    
      REAL,public:: ZRTABS(MAXIZS), ZRVABS(MAXIZS), ZRAVAV(MAXIZS)            
      REAL,public:: ZTBS(MAXIZS)                                              

 
  public :: allocate_mod_zommv, deallocate_mod_zommv


contains

  subroutine allocate_mod_zommv
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_zommv


  subroutine deallocate_mod_zommv
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_zommv



end module mod_zommv
