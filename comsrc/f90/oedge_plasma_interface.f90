module oedge_plasma_interface

implicit none

private







public :: get_oedge_plasma,init_oedge_plasma




contains


subroutine init_oedge_plasma(filename,ierr)
  implicit none

  ! open input files 
  ! read array sizes from the geometry file


  ! Allocate storage


  

  ! load geometry data




  ! load plasma data




end subroutine init_oedge_plasma


subroutine load_oedge_geometry(fnum,ierr)
   implicit none

   ! Load tagged geometry data and allocate storage for geometry



end subroutine load_oedge_geometry



subroutine load_oedge_plasma(fnum,ierr)
  implicit none
  ! load oedge plasma file


end subroutine load_oedge_plasma



subroutine close_oedge_plasma
  implicit none

  ! deallocate any allocated storage




end subroutine close_oedge_plasma


subroutine get_oedge_plasma(r,z,ne,te,ti,vb,ef)
  implicit none

  real :: r,z,ne,te,ti,vb,ef

  ! Get the OEDGE plasma conditions at the specified R,Z location
  ! Two options - value in cell and interpolated - value in cell is quicker



end subroutine get_oedge_plasma




end module oedge_plasma_interface
