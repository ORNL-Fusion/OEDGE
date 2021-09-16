module mod_printr

  use mod_params

!  C                                                                               
!      COMMON /PRINTR/ IQY0,  IQYL8, IQYL4, IQYL2, IQYL,  IQY3L2,IQY7L4,         
!     >                IQY158,IQY2L, IQYS2, IQY3S2,IQY5S2,IQXA,  IQXA2,          
!     >                IQXA4, IQXIN, IQXOUT,IQXAW4,IQXAW2,IQXAW, IQXFAC,         
!     >                IY0,   IYL8,  IYL4,  IYL2,  IYL,   IY3L2, IY7L4,          
!     >                IY158, IY2L,  IYS2,  IY3S2, IY5S2, IXA,   IXA2,           
!     >                IXA4,  IXIN,  IXOUT, IXAW4, IXAW2, IXAW,  IXFAC,          
!     >                IQY5S, IQY15S,IY0LT                                             
!      INTEGER         IQY0,  IQYL8, IQYL4, IQYL2, IQYL,  IQY3L2,IQY7L4,         
!     >                IQY158,IQY2L, IQYS2, IQY3S2,IQY5S2,IQXA,  IQXA2,          
!     >                IQXA4, IQXIN, IQXOUT,IQXAW4,IQXAW2,IQXAW, IQXFAC,         
!     >                IY0,   IYL8,  IYL4,  IYL2,  IYL,   IY3L2, IY7L4,          
!     >                IY158, IY2L,  IYS2,  IY3S2, IY5S2, IXA,   IXA2,           
!     >                IXA4,  IXIN,  IXOUT, IXAW4, IXAW2, IXAW,  IXFAC,          
!     >                IQY5S, IQY15S,IY0LT                                             




  implicit none
  private

  INTEGER,public::IQY0,  IQYL8, IQYL4, IQYL2, IQYL,  IQY3L2,IQY7L4,&
                  IQY158,IQY2L, IQYS2, IQY3S2,IQY5S2,IQXA,  IQXA2,&          
                  IQXA4, IQXIN, IQXOUT,IQXAW4,IQXAW2,IQXAW, IQXFAC,&
                  IY0,   IYL8,  IYL4,  IYL2,  IYL,   IY3L2, IY7L4,&          
                  IY158, IY2L,  IYS2,  IY3S2, IY5S2, IXA,   IXA2,&           
                  IXA4,  IXIN,  IXOUT, IXAW4, IXAW2, IXAW,  IXFAC,&          
                  IQY5S, IQY15S,IY0LT                                             



  
  public :: allocate_mod_printr, deallocate_mod_printr


contains

  subroutine allocate_mod_printr
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_printr


  subroutine deallocate_mod_printr
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_printr



end module mod_printr
