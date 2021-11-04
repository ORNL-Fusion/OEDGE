module hc_velocity_type

  !
  ! hc_v defines the velocity of the HC molecule
  !
  ! This type contains - the 3D velocity for neutrals (as a 3-vector)
  !                    - the parallel and perpendicular components for ion velocity 
  !                      (note that the perpendicular direction is not defined since precession and particle
  !                      collisions are assumed to randomize it).
  !                    - the "temperature" of the HC molecule
  !                    - separate parallel and perpendicular temperatures
  !                  
  ! Some of the options of the treatment will map the "temperature" values to and from the velocity of 
  ! the HC fragment. This is not quite accurate since "temperature" is a measure of the velocity magnitude and 
  ! distribution of an ensemble of particles. 
  !
  ! In the cases where energy (in eV) is mapped to a change in velocity deltaV - the formula used is:
  ! E (eV)  - kT = 1/2 m deltaV **2   or deltaV = (2 * 1.6e-19 * E  / (1.67e-27 * Mhc) ) ** (1/2)   
  !
! slmod begin
! Intel compiler does not like a module and derived type to have the same name,
! so I changed all instances in the code to type (hc_velocity_type1) from
! type (hc_velocity_type).  Use replace command to change to whatever suites you.
!
  type hc_velocity_type1
!
!   type hc_velocity_type
! slmod end
     ! real vr,vz,vt
     ! v(1)=vr, v(2)=vz, v(3)=vt
     real v(3)
     real vpara,vperp
     real vtot
     real t,tperp,tpara
     logical enable3D
! slmod begin
  end type hc_velocity_type1
!
!  end type hc_velocity_type1
! slmod end

contains 


  real function mag_v(v)
    implicit none
    real :: v(3) 
    
       mag_v = sqrt(sum(v**2))

     end function mag_v

  subroutine norm_v(v,normv)
    implicit none
    real :: v(3), normv(3) 
    real :: mag

    mag = mag_v(v) 
    normv = v / mag

  end subroutine norm_v


end module hc_velocity_type
