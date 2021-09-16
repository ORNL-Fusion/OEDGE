module mod_coords

  use mod_params

!c     -*-Fortran-*-
!c
!      COMMON /COORDS/ X0S,X0L,Y0S,Y0L,XL1,YL1,XL2,YL2,WEDGAN,                   
!     >                TC,SC,TO,SO,TV,SV,GC,RP,
!     >                XST1,YST1,XST2,YST2,XST3,YST3,
!     >                P0S,P0L
!
!      REAL X0S,X0L,Y0S,Y0L,XL1(2),YL1(2),XL2(2),YL2(2),WEDGAN(2)                
!      REAL TC,SC,TO,SO,TV,SV,GC,RP                                              
!      REAL XST1(2),YST1(2),XST2(2),YST2(2),XST3(2),YST3(2) 
!      real p0s,p0l
  
  implicit none
  private

      REAL,PUBLIC:: X0S,X0L,Y0S,Y0L,XL1(2),YL1(2),XL2(2),YL2(2),WEDGAN(2)                
      REAL,PUBLIC:: TC,SC,TO,SO,TV,SV,GC,RP                                              
      REAL,PUBLIC:: XST1(2),YST1(2),XST2(2),YST2(2),XST3(2),YST3(2) 
      real,public:: p0s,p0l


  
  public :: allocate_mod_coords, deallocate_mod_coords


contains

  subroutine allocate_mod_coords
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_coords


  subroutine deallocate_mod_coords
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_coords



end module mod_coords
